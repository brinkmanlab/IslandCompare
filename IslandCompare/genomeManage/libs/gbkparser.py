from Bio import SeqIO
import os
import logging

def getGenesFromGbk(filePath):
    # Given a path to a gbk file, this will return all CDS
    geneList = []
    for record in SeqIO.parse(open(filePath),"genbank"):
        for feature in record.features:
            geneInfo = {}
            if feature.type=='gene' or feature.type=='CDS':
                # Bio.SeqIO returns 1 for (+) and  -1 for (-)
                geneInfo['strand']=feature.location.strand
                geneInfo['start']=feature.location.start
                geneInfo['end']=feature.location.end
                try:
                    geneInfo['note']=feature.qualifiers['note']
                except:
                    logging.info("No Notes Found For This Gene")
                try:
                    geneInfo['name']=feature.qualifiers['gene'][0]
                except:
                    logging.info("No Name Found For This Gene")
                    try:
                        geneInfo['name']=feature.qualifiers['locus_tag']
                    except:
                        logging.info("No Locus Found For This Gene")
                geneList.append(geneInfo)
        # Only gather data from the first genome in a gbk file
        break
    return geneList

def testGetGenesFromGbk():
    testfiledir = os.path.dirname(os.path.realpath(__file__))+"/testfiles/"
    filePath = testfiledir+"/"+"AE009952.gbk"
    print getGenesFromGbk(filePath)

if __name__ == "__main__":
    testGetGenesFromGbk()