from tempfile import NamedTemporaryFile
import os
import subprocess
import glob
import genomeparser
import logging

VSEARCH_EXECUTABLE = "vsearch"

def cluster(combinedFastaFilePath, clusterFilePrefix, identityThreshold):
    scriptFile = NamedTemporaryFile(delete=True)

    with open(scriptFile.name, 'w') as script:
        script.write("#!/bin/bash\n")
        script.write(VSEARCH_EXECUTABLE + " --cluster_fast " + combinedFastaFilePath + " --clusters " + clusterFilePrefix + "-" + " --id " + str(identityThreshold))
        script.close()

    os.chmod(scriptFile.name, 0755)
    scriptFile.file.close()
    subprocess.check_call(scriptFile.name)
    scriptFile.close()

def aggregateClusterFiles(targetDirectory, outputPrefix):
    outputDict = {}
    logging.info("Attempting to read: " + targetDirectory+"/"+outputPrefix+"*")
    targetFiles = glob.glob(targetDirectory+"/"+outputPrefix+"*")
    for file in targetFiles:
        clusterId = file.split("-")[-1]
        seqRecords = genomeparser.readFastaFile(file)
        for record in seqRecords:
            sequenceId = record.id
            splitSequenceId = sequenceId.split("-")
            sequenceName = splitSequenceId[0]
            if sequenceName not in outputDict.keys():
                outputDict[sequenceName] = {}
            islandId = splitSequenceId[1]
            outputDict[sequenceName][islandId] = clusterId
    logging.info(outputDict)
    return outputDict

def countNumberClusters(targetDirectory, outputPrefix):
    return(len(glob.glob(targetDirectory+"/"+outputPrefix+"*")))

## TESTS

def clusterTest():
    cluster("/data/fna/TEST.fna", "/data/fna/output", 0.9)

def aggregateClusterFilesTest():
    aggregateClusterFiles("/data/vsearch/6", "output")

if __name__ == "__main__":
    #clusterTest()
    aggregateClusterFilesTest()