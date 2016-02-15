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

def test1():
    testfiledir = os.path.dirname(os.path.realpath(__file__))+"/testfiles/"
    runSigiHMM(testfiledir+"Pseudomonas_aeruginosa_PAO1_107.embl","/tmp/testsigi.embl","/tmp/testsigi.gff")
