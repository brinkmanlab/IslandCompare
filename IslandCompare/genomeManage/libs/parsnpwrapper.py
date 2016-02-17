from tempfile import NamedTemporaryFile, mkdtemp
from shutil import copy, rmtree
import os
import subprocess

PARSNP_PATH = "/apps/Parsnp-Linux64-v1.2/parsnp"

def runParsnp(inputFiles,outputDir):
    scriptFile = NamedTemporaryFile(delete=True)
    tempDirPath = mkdtemp()

    for inputFile in inputFiles:
        copy(inputFile,tempDirPath)

    with open(scriptFile.name,'w') as script:
        script.write("#!/bin/bash\n")
        script.write(PARSNP_PATH+" -r ! -d "+tempDirPath+" -o "+outputDir)
        script.close()

    os.chmod(scriptFile.name, 0755)
    scriptFile.file.close()
    subprocess.check_call(scriptFile.name)
    scriptFile.close()

    rmtree(tempDirPath)
    return outputDir

def testRunParsnp():
    testfiledir = os.path.dirname(os.path.realpath(__file__))+"/testfiles/"
    seqpaths = [testfiledir+"Al-Hasa_1_2013.fna",testfiledir+"Al-Hasa_2_2013.fna"]
    runParsnp(seqpaths,"/tmp")
