from Bio import SeqIO
import os

def getGenesFromGbk(filePath):
    # Given a path to a gbk file, this will return all CDS
    geneList = []
    for record in SeqIO.parse(open(filePath),"genbank"):
        for feature in record.features:
            toSendFlag = True
            geneInfo = {}
            if feature.type=='gene':
                geneInfo['start']=feature.location.start
                geneInfo['end']=feature.location.end
                try:
                    geneInfo['note']=feature.qualifiers['note']
                except:
                    print "No Notes Found For This Gene"
                try:
                    geneInfo['name']=feature.qualifiers['gene']
                except:
                    print "No Name Found For This Gene"
                    toSendFlag = False
                if toSendFlag:
                    geneList.append(geneInfo)
        break
    return geneList

def testGetGenesFromGbk():
    testfiledir = os.path.dirname(os.path.realpath(__file__))+"/testfiles/"
    filePath = testfiledir+"/"+"AE009952.gbk"
    print getGenesFromGbk(filePath)