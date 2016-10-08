from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def getSubsequence(genbankFile, startPosition, endPosition, islandNumber, description=None):
    record_dict = SeqIO.index(genbankFile, "genbank")
    if (description is not None):
        return SeqRecord(record_dict[list(record_dict.keys())[0]].seq[startPosition:endPosition], id=list(record_dict.keys())[0] + "-" + str(islandNumber), description=description)
    else:
        return SeqRecord(record_dict[list(record_dict.keys())[0]].seq[startPosition:endPosition], id=list(record_dict.keys())[0] + "-" + str(islandNumber))

def writeFastaFile(outputFileName, seqRecordList):
    with open(outputFileName, 'w') as outputFileHandle:
        SeqIO.write(seqRecordList, outputFileHandle, "fasta")

# Tests

def testFastaFileOpen():
    print(getSubsequence("/vagrant/IslandCompare/genomeManage/libs/testfiles/AE009952.gbk", 0, 100, 1).format("fasta"))

def testWriteFasta():
    writeFastaFile("/vagrant/IslandCompare/genomeManage/libs/testfiles/AE009952IslandFake.gbk", getSubsequence("/vagrant/IslandCompare/genomeManage/libs/testfiles/AE009952.gbk", 0, 100, "1"))

if __name__ == "__main__":
    testFastaFileOpen()
    testWriteFasta()