import subprocess
import os
from tempfile import NamedTemporaryFile

SIGIHMM_PATH = "/apps/Colombo_3.8"
SIGIHMM_EXE = "SigiHMM"

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
    # returns a list of dicts of start, stop, and strand (+ or -)  of Genomic Islands
    # Putal genes are considered Islands
    listPutalGenes = []
    with open(gffoutput,'r') as gff:
        for line in gff:
            if line[0] == '#':
                continue
            else:
                cleanedLine = ' '.join(line.split())
                geneDict = cleanedLine.split(' ')
                if geneDict[2] == 'PUTAL':
                    try:
                        details = geneDict[8]
                    except:
                        details = ''
                    listPutalGenes.append({'start':geneDict[3],'end':geneDict[4],
                                           'strand':geneDict[6],'details':details})
                else:
                    pass
    return listPutalGenes

def testRunSigiHMM():
    testfiledir = os.path.dirname(os.path.realpath(__file__))+"/testfiles/"
    runSigiHMM(testfiledir+"Pseudomonas_aeruginosa_PAO1_107.embl","/tmp/testsigi.embl","/tmp/testsigi.gff")

def testParser():
    testfile = "/data/sigi/AE009952.1.gff"
    output = parseSigiGFF(testfile)
    for line in output:
        print output
        print '\n'
