from tempfile import NamedTemporaryFile, mkdtemp
from shutil import copy, rmtree
import os
import subprocess
from Bio import Phylo

PARSNP_PATH = "/apps/Parsnp-Linux64-v1.2/parsnp"

def runParsnp(inputFiles,outputDir):
    # wrapper for parsnp,
    # inputFiles = list of all fna files to run parsnp with
    # output dir = desired output directory

    # parsnp wants a directory with the input files as opposed to a list so I need to create a temporary directory
    # and copy the input files into that directory before running parsnp
    tempDirPath = mkdtemp()
    for inputFile in inputFiles:
        copy(inputFile,tempDirPath+"/"+(os.path.splitext(inputFile)[0]).split("/")[-1]+".fna")

    # create the script that will call parsnp
    scriptFile = NamedTemporaryFile(delete=True)
    with open(scriptFile.name,'w') as script:
        script.write("#!/bin/bash\n")
        script.write(PARSNP_PATH+" -r ! -d "+tempDirPath+" -o "+outputDir)
        script.close()

    # run parsnp
    os.chmod(scriptFile.name, 0777)
    scriptFile.file.close()
    subprocess.check_call(scriptFile.name)
    scriptFile.close()

    # delete the temporary directory before function completes
    rmtree(tempDirPath)
    return outputDir

def newickToArray(inputFile):
    # parses a newick file and returns an array of dicts of clades
    tree = Phylo.read(inputFile, 'newick')
    return parsePhyloTree(tree.root)

def parsePhyloTree(node):
    # recursively builds a dict representing a tree rooted at the first input node
    currentNode = {}
    currentNode['name'] = node.name
    currentNode['length'] = node.branch_length
    if len(node.clades) > 0:
        currentNode['children'] = []
        for childNode in node.clades:
            currentNode['children'].append(parsePhyloTree(childNode))
    return currentNode


### Tests

def testRunParsnp():
    testfiledir = os.path.dirname(os.path.realpath(__file__))+"/testfiles/"
    seqpaths = [testfiledir+"Al-Hasa_1_2013.fna",testfiledir+"Al-Hasa_2_2013.fna"]
    runParsnp(seqpaths,"/tmp")

def testNewickToArray():
    newickToArray("/data/parsnp/1/parsnp.tree")