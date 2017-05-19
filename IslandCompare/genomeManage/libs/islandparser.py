import logging
from django.conf import settings
import os
import sigihmmwrapper
import genomeparser
import errno

MASH_OUTPUT_PATH = settings.MASH_OUTPUT_PATH

def createSequenceGIFolder(sequenceName):
    outputDirectory = MASH_OUTPUT_PATH + "/" + sequenceName
    try:
        logging.info("Creating Directory: " + outputDirectory)
        os.mkdir(outputDirectory)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(outputDirectory):
            pass


def createFastaFilesForGis(genbankFile, sigiFile, sequenceName):
    genomicIslandList = sigihmmwrapper.parseSigiGFF(sigiFile)
    numberGenomicIslands = len(genomicIslandList)
    outputFileList = []

    for genomicIslandIndex in range(numberGenomicIslands):
        genomicIsland = genomicIslandList[genomicIslandIndex]

        islandRecord = genomeparser.getSubsequence(genbankFile, int(genomicIsland['start']), int(genomicIsland['end']), genomicIslandIndex)
        outputFileName = MASH_OUTPUT_PATH + "/" + sequenceName + "/" + str(genomicIslandIndex)
        genomeparser.writeFastaFile(outputFileName, [islandRecord])
        outputFileList.append(outputFileName)

    return outputFileList

# Tests

def testGIFolderCreation():
    createSequenceGIFolder("test")

def testFastaFilesForGis():
    createFastaFilesForGis("/vagrant/IslandCompare/genomeManage/libs/testfiles/AE009952.gbk",
                           "/vagrant/IslandCompare/genomeManage/libs/testfiles/AE009952.11.gff",
                           "test")

if __name__ == "__main__":
    testGIFolderCreation()
    testFastaFilesForGis()

