from django.conf import settings
from tempfile import NamedTemporaryFile
import os
import subprocess
import logging
import csv

MASH_PATH = "/apps/mash-Linux64-v1.1.1/mash" #settings.MASH_PATH
MASH_OUTPUT_PATH = "/data/mash"#settings.MASH_OUTPUT_PATH

def createCompoundSketch(fastaFileList, outputFileName):
    scriptFile = NamedTemporaryFile(delete=False)

    with open(scriptFile.name, 'w') as script:
        script.write("#!/bin/bash\n")
        script.write(MASH_PATH + " sketch -o " + outputFileName + " ")
        for fastaFile in fastaFileList:
            script.write(fastaFile+" ")
        script.close()

    os.chmod(scriptFile.name, 0755)
    scriptFile.file.close()

    logging.info("Running script: " + scriptFile.name)
    subprocess.check_call(scriptFile.name)
    scriptFile.close()

def calculateMashDistance(referenceFile, queryFastaFile, outputFile):
    scriptFile = NamedTemporaryFile(delete=True)

    with open(scriptFile.name, 'w') as script:
        script.write("#!/bin/bash\n")
        script.write(MASH_PATH + " dist " + referenceFile + " " + queryFastaFile + " > " + outputFile)
        script.close()

    os.chmod(scriptFile.name, 0755)
    scriptFile.file.close()
    subprocess.check_call(scriptFile.name)
    scriptFile.close()

def mashOutputFileParser(outputFile):
    outputList = []
    with open(outputFile, 'r') as output:
        reader = csv.reader(output, delimiter='\t')
        for row in reader:
            outputList.append({'referenceId':row[0], 'queryId':row[1], 'mashDistance':row[2], 'pValue':row[3], 'matchingHashes':row[4]})
    return outputList

##TEST

def testOutputFileParser():
    print(mashOutputFileParser("/data/mash/10/output"))

if __name__ == "__main__":
    testOutputFileParser()