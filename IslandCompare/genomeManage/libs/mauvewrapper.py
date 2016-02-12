from tempfile import NamedTemporaryFile
import subprocess
import os
import shutil

MAUVE_PATH = "/apps/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve"
MAUVE_OUTPUT_PATH = "/data/mauve"

def runMauve(sequencepaths, outputbackbonepath):
    # Parameters = path to 2 genbank files
    # Returns None
    # Creates an output file at path outputfile and backbone file at path backbonefile
    scriptFile = NamedTemporaryFile(delete=True)

    tmppaths = []
    for gbk in sequencepaths:
        temppath = MAUVE_OUTPUT_PATH+"/"+os.path.basename(gbk)
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
    return None

### Tests

def test1():
    testfiledir = os.path.dirname(os.path.realpath(__file__))+"/testfiles/"
    seqpaths = [testfiledir+"AE009952.gbk",testfiledir+"BX936398.gbk"]
    runMauve(seqpaths, "/tmp/test1")
