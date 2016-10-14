from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def getSubsequence(genbankFile, startPosition, endPosition, islandNumber, description=None):
    record_dict = SeqIO.index(genbankFile, "genbank")
    sequenceName = list(record_dict.keys())[0]
    if (description is not None):
        return SeqRecord(record_dict[sequenceName].seq[int(startPosition):int(endPosition)], id=sequenceName + "-" + str(islandNumber), description=description)
    else:
        return SeqRecord(record_dict[sequenceName].seq[int(startPosition):int(endPosition)], id=sequenceName + "-" + str(islandNumber))

def writeFastaFile(outputFileName, seqRecordList):
    with open(outputFileName, 'w') as outputFileHandle:
        SeqIO.write(seqRecordList, outputFileHandle, "fasta")

def readFastaFile(fastaFile):
    returnRecords = []
    for record in SeqIO.parse(fastaFile, "fasta"):
        returnRecords.append(record)
    return returnRecords

# Tests

def testFastaFileOpen():
    print(getSubsequence("/vagrant/IslandCompare/genomeManage/libs/testfiles/AE009952.gbk", 0, 100, 1).seq)

def testWriteFasta():
    writeFastaFile("/vagrant/IslandCompare/genomeManage/libs/testfiles/AE009952IslandFake.gbk", getSubsequence("/vagrant/IslandCompare/genomeManage/libs/testfiles/AE009952.gbk", 0, 100, "1"))

def testReadFasta():
    readFastaFile("/data/vsearch/output0")

if __name__ == "__main__":
    testFastaFileOpen()
    testWriteFasta()
    testReadFasta()