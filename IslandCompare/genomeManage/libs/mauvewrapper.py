from tempfile import NamedTemporaryFile
import subprocess
import os
import shutil
import csv

MAUVE_PATH = "/apps/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve"
MAUVE_OUTPUT_PATH = "/data/mauve"

def runMauve(sequencepaths, outputbackbonepath):
    # Parameters = path to 2 genbank files
    # Returns None
    # Creates an output file at path outputfile and backbone file at path backbonefile
    scriptFile = NamedTemporaryFile(delete=True)

    tmppaths = []

    # Copy the gbk files needed for mauve
    for gbk in sequencepaths:
        temppath = os.path.dirname(outputbackbonepath)+"/"+os.path.basename(gbk)
        tmppaths.append(temppath)
        shutil.copyfile(gbk,temppath)

    with open(scriptFile.name,'w') as script:
        script.write("#!/bin/bash\n")
        script.write(MAUVE_PATH+" --backbone-output="+outputbackbonepath+
                     ".backbone ")
        for sequence in tmppaths:
            script.write(sequence+" ")
        script.close()

    os.chmod(scriptFile.name, 0755)
    scriptFile.file.close()
    subprocess.check_call(scriptFile.name)
    scriptFile.close()

    # Delete the temporary gbk files used for mauve
    for tmp in tmppaths:
        os.remove(tmp)

    return None

def parseMauveBackbone(backbonePath):
    # input path to backbone file
    # mauveOutput will look like:
    # [{seq1:[start,end],seq2:[start,end]},{seq1:[start,end],seq2:[start,end]}]
    mauveOutput = []
    with open(backbonePath,'r') as backbone:
        tsvReader = csv.reader(backbone, delimiter='\t')
        header = tsvReader.next()
        for row in tsvReader:
            homologousRegion ={}
            for sequence in range(0,len(row)/2):
                homologousRegion[sequence]=row[sequence*2],row[sequence*2+1]
            mauveOutput.append(homologousRegion)
    return mauveOutput

def combineMauveBackbones(backbonePaths, outputfile, orderedIdList = None):
    # backbonePaths = list of backbone paths
    # will create a file at path outputfile of all the paths merged in order of backbonePaths
    # orderedIdList takes a list where the index is the index of the (sorted sequence by id)
    # and the value is the column of this value in the merged file
    with open(outputfile, 'w') as output:
        outputWriter = csv.writer(output, delimiter='\t')
        # Add the header row to the outputfile
        headerRow = []
        for outputCounter in range(len(backbonePaths)+1):
            headerRow.append("seq"+str(outputCounter)+"_leftend")
            headerRow.append("seq"+str(outputCounter)+"_rightend")
        outputWriter.writerow(headerRow)
        # Add the columns to the outputfile file by file
        for pathCounter in range(len(backbonePaths)):
            with open(backbonePaths[pathCounter]) as inputfile:
                inputReader = csv.reader(inputfile, delimiter='\t')
                # Skip the header row
                inputReader.next()
                for row in inputReader:
                    topStart = int(row[0])
                    topEnd = int(row[1])
                    bottomStart = int(row[2])
                    bottomEnd = int(row[3])
                    outputRow = []

                    if orderedIdList is None:
                        orderList = [i for i in range(len(backbonePaths)+1)]
                    else:
                        orderList = list(orderedIdList)
                        orderList.append(len(backbonePaths))

                    firstSequence = int(orderList[pathCounter]) # the column to put top sequence in
                    secondSequence = int(orderList[pathCounter+1]) # the column to put bot sequence in

                    # Place sequences in appropriate columns
                    for columnIndex in range(2*(len(backbonePaths)+1)):
                        if columnIndex == firstSequence*2:
                            outputRow.append(topStart)
                        elif columnIndex == firstSequence*2+1:
                            outputRow.append(topEnd)
                        elif columnIndex == secondSequence*2:
                            outputRow.append(bottomStart)
                        elif columnIndex == secondSequence*2+1:
                            outputRow.append(bottomEnd)
                        else:
                            outputRow.append(0)
                    outputWriter.writerow(outputRow)

### Tests

def testMauve(outputName):
    testfiledir = os.path.dirname(os.path.realpath(__file__))+"/testfiles/"
    seqpaths = [testfiledir+"AE009952.gbk",testfiledir+"BX936398.gbk"]
    runMauve(seqpaths, outputName)

def testParser():
    print parseMauveBackbone("/tmp/test1.backbone")

def testMergeMauve(outputList):
    combineMauveBackbones(outputList,"/tmp/testmerged.backbone",[3,2,1,0])

if __name__ == "__main__":
    testMauve("/tmp/1")
    testParser()

    testMauve("/tmp/2")
    testMauve("/tmp/3")
    testMauve("/tmp/4")
    testMergeMauve(["/tmp/1.backbone","/tmp/2.backbone",
                    "/tmp/3.backbone","/tmp/4.backbone"])