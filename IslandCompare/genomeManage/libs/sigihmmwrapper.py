import subprocess
import os
from tempfile import NamedTemporaryFile
from django.conf import settings
import logging

SIGIHMM_PATH = settings.SIGIHMM_PATH
SIGIHMM_EXE = settings.SIGIHMM_EXE

def runSigiHMM(emblinput,embloutput,gffoutput):
    # emblinput = embl input file
    # embloutput = embl output
    # gff = gff output file
    # blocking call to SigiHMM, returns None on completion
    scriptFile = NamedTemporaryFile(delete=True)

    with open(scriptFile.name,'w') as script:
        script.write("#!/bin/bash\n")
        script.write("/usr/bin/java "+SIGIHMM_EXE+" input="+emblinput+" output="+embloutput+" gff="+gffoutput)
        script.close()

    os.chmod(scriptFile.name, 0755)
    scriptFile.file.close()

    sp = subprocess.Popen(scriptFile.name, cwd=SIGIHMM_PATH)
    sp.wait()

    scriptFile.close()
    return None

def parseSigiGFF(gffoutput):
    # gffoutput = output from running SigiHMM
    # returns a list of dicts of start, stop of Genomic Islands
    # Putal genes are considered Islands
    # this parses sigihmm in the same way that Islandviewer does it
    if not os.path.isfile(gffoutput):
        logging.warn("No GFF File Found: Assuming no Genomic Islands for this file.")
        return []

    listPutalGenes = []
    with open(gffoutput,'r') as gff:
        currentStart = 0
        currentEnd = 0
        islandFlag = False # keeps track of whether currently parsing a GI or not
        for line in gff:
            # skip the lines with parameters
            if line[0] == '#':
                continue
            else:
                cleanedLine = ' '.join(line.split())
                geneDict = cleanedLine.split(' ')
                # at the start of a genomic island, set start and end of possible genomic island
                if geneDict[2] == 'PUTAL' and not islandFlag:
                    currentStart = geneDict[3]
                    currentEnd = geneDict[4]
                    islandFlag = True
                # continuation of current genomic island, change end and continue
                elif geneDict[2] == 'PUTAL' and islandFlag:
                    currentEnd = geneDict[4]
                # end of genomic island, append current start and end to list
                elif islandFlag:
                    listPutalGenes.append({'start':currentStart,'end':currentEnd})
                    islandFlag = False
                # not currently in a genomic island, continue to next line
                elif not islandFlag:
                    continue
                # condition not included in above reached, throw an error
                else:
                    raise Exception("Error occured in sigi wrapper, unexpected condition reached")
    return listPutalGenes

### Tests

def testRunSigiHMM():
    testfiledir = os.path.dirname(os.path.realpath(__file__))+"/testfiles/"
    runSigiHMM(testfiledir+"Pseudomonas_aeruginosa_PAO1_107.embl","/tmp/testsigi.embl","/tmp/testsigi.gff")

def testParser():
    testfile = "/vagrant/IslandCompare/genomeManage/libs/testfiles/AE009952.11.gff"
    output = parseSigiGFF(testfile)
    for line in output:
        print output
        print '\n'
