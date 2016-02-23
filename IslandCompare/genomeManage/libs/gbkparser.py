from Bio import SeqIO
import os

def getGenesFromGbk(filePath):
    # Given a path to a gbk file, this will return all CDS
    geneList = []
    for record in SeqIO.parse(open(filePath),"genbank"):
        geneInfo = {}
        for feature in record.features:
            if feature.type=='gene':
                geneInfo['start']=feature.location.start
                geneInfo['end']=feature.location.end
                geneInfo['name']=feature.qualifiers['gene']
                try:
                    geneInfo['note']=feature.qualifiers['note']
                except:
                    pass
                geneList.append(geneInfo)
        break
    return geneList

def testGetGenesFromGbk():
    testfiledir = os.path.dirname(os.path.realpath(__file__))+"/testfiles/"
    filePath = testfiledir+"/"+"AE009952.gbk"
    print getGenesFromGbk(filePath)
