import csv
import os

# Create a dict with contains keys as filename and value as list of genomic islands {start,end}
def parseGiFile(inputFile):
    genomeDict = dict()
    with open(inputFile, 'rU') as gifile:
        gireader = csv.reader(gifile, dialect=csv.excel_tab)
        for row in gireader:
            genomeName = row[0]
            giStart = row[1]
            giEnd = row[2]
            # if genome name is not in genomeDict then add it
            if genomeName not in genomeDict:
                genomeDict[genomeName] = list()
            # add start and end of current genome list in genomedict
            if len(row)>3:
                giColor = row[3]
                genomeDict[genomeName].append({'start':giStart,'end':giEnd,'color':giColor})
            else:
                genomeDict[genomeName].append({'start':giStart,'end':giEnd})
    return genomeDict

# Parse input file to determine if any entries have a color provided
def determineIfColorsProvided(inputFile):
    colorProvidedFlag = False
    with open(inputFile, 'rU') as gifile:
        gireader = csv.reader(gifile, dialect=csv.excel_tab)
        for row in gireader:
            if len(row)>3:
                colorProvidedFlag = True
                break
    return colorProvidedFlag

# Tests

def parseTest():
    testDict = parseGiFile(os.path.dirname(os.path.realpath(__file__))+"/testfiles/GenomicIslands")
    print(testDict)

if __name__ == "__main__":
    parseTest()