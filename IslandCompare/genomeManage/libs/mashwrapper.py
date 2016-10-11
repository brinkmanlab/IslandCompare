from django.conf import settings
from tempfile import NamedTemporaryFile
import os
import subprocess

MASH_PATH = settings.MASH_PATH
MASH_OUTPUT_PATH = settings.MASH_OUTPUT_PATH

def createCompoundSketch(fastaFileList, outputFileName):
    scriptFile = NamedTemporaryFile(delete=True)

    with open(scriptFile.name, 'w') as script:
        script.write("#!/bin/bash\n")
        script.write(MASH_PATH + " sketch -o " + MASH_OUTPUT_PATH + "/" + outputFileName + " ")
        for fastaFile in fastaFileList:
            script.write(fastaFile+" ")
        script.close()

    os.chmod(scriptFile.name, 0755)
    scriptFile.file.close()
    subprocess.check_call(scriptFile.name)
    scriptFile.close()

def calculateMashDistance(referenceFile, queryFastaFile):
    scriptFile = NamedTemporaryFile(delete=True)

    with open(scriptFile.name, 'w') as script:
        script.write("#!/bin/bash\n")
        script.write(MASH_PATH + " dist " + referenceFile + " " + queryFastaFile)
        script.close()

    os.chmod(scriptFile.name, 0755)
    scriptFile.file.close()
    subprocess.check_call(scriptFile.name)
    scriptFile.close()