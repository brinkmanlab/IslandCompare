import numpy as np
import scipy.cluster.hierarchy as hac

def computeVector(sequence):
    countA = sequence.count('A')
    countT = sequence.count('T')
    countC = sequence.count('C')
    countG = sequence.count('G')
    return [countA, countT, countC, countG]

def computeClusters(data, method="single", criterion="maxclust"):
    array = np.array(data)
    linkageMatrix = hac.linkage(array, method=method)
    knee = np.diff(linkageMatrix[::-1, 2], 2)

    numberClusters = knee.argmax() + 2
    clusterGroups = hac.fcluster(linkageMatrix, numberClusters, criterion)

    return {"numberClusters": numberClusters, "clusterGroups": clusterGroups}

## TESTS

def computeVectorTest():
    print(computeVector("ACCTTTGGGG"))

def computeClusterTest():
    print(computeClusters([[0.1,   2.5, 5.,1],
                           [1.5,   .4, 3,1],
                           [0.3,   1, 6,1],
                           [1  ,   .8, 4,1],
                           [0.5,   0, 1  ,1],
                           [0  ,   0.5, 3,1],
                           [0.5,   0.5, 7,1],
                           [2.7,   2, 3 ,1],
                           [2.2,   3.1, 1,7],
                           [3  ,   2  , 1,3],
                           [3.2,   1.3, 4,4]]))

if __name__ == "__main__":
    computeVectorTest()
    computeClusterTest()

