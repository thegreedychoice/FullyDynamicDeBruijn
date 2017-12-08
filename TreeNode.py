from RabinKarpHash import *
import perfection.czech
from In_Out_Matrix import *
from collections import deque

class TreeNode:

    height = 0
    visited = False
    val = None
    level = None
    parent = None
    isRoot = False


    def __init__(self, val, level, parent, isRoot):

        # value inside the node (hash value)
        self.val = val
        # level of the node in relation to the root, i.e. root has level 0
        self.level = level
        # pointer to its parent node
        self.parent = parent
        # boolean value if its a root
        self.isRoot = isRoot
        self.visited = False





#TODO: take care of root case
#initialize

    #def deleteNode(kmer_string):
     #   rootFound = findRoot(kmer_string)
        #calls find root: if findroot is True then delete from in and out


    #def initializeTree(kmer_values):

    '''def findRoot(self,kmer_string):
        if self.level == 0:
            #get leaf
            #searchedNode = findLeaf(kmer_string)
        #if it is the root
        if self.isRoot == True and self.val != None:
            return True
        elif self.isRoot == True and self.parent == None:
            return False
        else:
            kmerList = list(self.val)
            # List containing the first potential letters in the kmer
            listAlphabet = ["A", "T", "C", "G"]
            # parent for the kmer
            parentDeque = deque(kmerList)
            parentDeque = self.parent.rotate(1)
            # convert back to a list
            parentList = list(collections.deque(parentDeque))
            # replace the first letter in the k_mer and hash it to check if it exists
            for x in range(0, 4):
                # replaces the first value in the kmer
                parentList[0] = listAlphabet[x]
                # hash the potential parent
                hashedValue = rabin_karp_hash(parentList)
                # check if the parent exists
                if self.parent.val == hashedValue:
                    return  self.findRoot(self.parent)

    #def findLeaf(kmer_string):
'''



    #def insertNode(self,val, 0, None, True):
    #    return __init__(self, val, newLevel, None, isRoot)


    # function to check if the node is a root
    def isRoot(self):
        return not self.parent
