from TreeNode import *
from In_Out_Matrix import *
import math
import collections


# case of leaf not done


class Forest:
    maxTreeHeight = 2

    # key is address of node and value is list typle in <value, visited>
    nodeDict = {}
    # list of all the tree roots for deletition and insertion purposes
    treeRoots = list()
    # list of all the leafs (may not be needed anymore)
    treeLeafs = list()
    str_to_mph = None
    inOutMatrix = None
    #visited list for checking duplicates
    visitedList= None

    # parentList = list()
    # number of kmers to compare and see if we have have added all of the kmers from the in-out matrix
    countNumberofKmers = 0

    def __init__(self, inoutMatrix, k, kmer, str_to_mph):
        if not inoutMatrix:
            # print("No in and out matrix")
            return

        self.inOutMatrix = inoutMatrix
        self.inOutMatrix.in_matrix = inoutMatrix.in_matrix
        self.inOutMatrix.out_matrix = inoutMatrix.out_matrix
        self.maxTreeHeight = 3 * k * math.log(5, 2)
        self.maxTreeHeight = 2
        self.str_to_mph = str_to_mph
        #print(self.inOutMatrix.in_matrix.rows)
        self.visitedList = [None] * self.inOutMatrix.in_matrix.rows

        # definition of a node: def __init__(val, level, parent, isRoot):

        # may not need id(initialNode)
        initialNode = TreeNode(kmer, 0, None, True)
        self.treeRoots.append(initialNode)

        curr = initialNode

        hashValue = self.str_to_mph.get(kmer)

        self.nodeDict.update({hashValue: curr})

        # start at index 0
        rootIndex = 0
        while rootIndex < len(self.treeRoots):
            # print(rootIndex)
            if rootIndex == 3:
                pass
            curr = self.treeRoots[rootIndex]

            self.insert_init(curr, curr, curr.val)

            self.countNumberofKmers += 1

            rootIndex += 1

            #print rootIndex

    # recursive function to build tree
    def insert_init(self, node, prev_node, kmer):


        # if the new height is overflown (above the max tree height we return and get that node)
        if node.level > self.maxTreeHeight:
            node.isRoot = True
            node.parent = None
            node.val = kmer
            node.level = 0
            self.treeRoots.append(node)
            return node

        #mark node as visited
        checkedIndex = node.val
        if node.isRoot == True:
            checkedIndex = self.str_to_mph.get(node.val)
        self.visitedList[checkedIndex] = 1

        # initialize current node
        curr = node

        # turn Node value to a list in order to rotate and update the value

        kmerList = list(kmer)
        # update the node level that we are going to add
        newLevel = node.level + 1

        # if the new level is the max height we have a leaf
        if newLevel > self.maxTreeHeight:
            # put into leaf list
            self.treeLeafs.append(curr)

        if curr.visited == True:
            return node

        # else we are not full and we have not visited the node
        else:
            # {
            # get all the headers which return true for the out matrix
            curr.visited = True
            # print(curr.val)
            arrayOut = self.inOutMatrix.return_column_headers_Out(kmer, self.str_to_mph)
            arrayIn = self.inOutMatrix.return_column_headers_In(kmer, self.str_to_mph)

            if arrayOut:
                for indexOut in arrayOut:
                    prevValOut = indexOut

                    curr, kmer_string = self.getOutNodes(curr, prevValOut, kmerList, newLevel)
                    if curr == None:
                        self.treeLeafs.append(node)
                        continue

                    # print("now here")
                    # print(kmer_string)

                    '''if arrayIn:
                        for indexIn in arrayIn:
                            prevValIn = indexIn
                            curr,kmer_string = self.getInNodes(curr,prevValIn,kmerList,newLevel)
                            return self.insert_init(curr, kmer_string)'''

                    self.insert_init(curr, node, kmer_string)

            return curr


    def insert_node(self, new_kmer):

        isAlreadyPresent = self.search(new_kmer, 0)

        if isAlreadyPresent:
            #print "Node already present. Can't be inserted!"
            return False
        else:
            #create a new node and insert it as a root

            #generate hash for this newley created node
            #hash_val = self.str_to_mph.get(new_kmer)

            hash_val = max(self.str_to_mph.values()) + 1

            new_node = TreeNode(new_kmer, 0, None, True)

            self.nodeDict.update({hash_val: new_node})
            self.str_to_mph[new_kmer] =  hash_val
            self.treeRoots.append(new_node)

            #also initialize in In and Out Matrix

            #print "Node Inserted!"
            return True

    #get the incoming nodes to a particular kmer from InMatrix
    def get_incoming_nodes(self, kmer):
        hash = self.str_to_mph.get(kmer)

        if hash == None:
            return dict()

        #check if hash has a valid entry in In Matrix and get all the incoming edges
        list = self.inOutMatrix.return_column_headers_Out(kmer, self.str_to_mph)


        nodes_to_be_updated_OutM = dict()

        for c in list:
            valueDeque = deque(kmer)  # 'C G A'
            # rotate the list right (kmer)
            valueDeque.rotate(-1)  # 'A C G'
            # convert back to list (kmer)
            valueList = list(collections.deque(valueDeque))
            # add column header to the end of the kmer value
            col_tobe_updated = valueList[0]
            valueList[0] = c  # 'T C G'
            # convert list back to a string
            kmer_string = ''.join(valueList)
            out_node_hash = self.str_to_mph.get(kmer_string)
            nodes_to_be_updated_OutM[out_node_hash] = col_tobe_updated

    #this method gives out all the nodes(hash) that had an incoming edge towards kmer
    #So the OutMatrix needs to be updated for these nodes
    def incoming(self, kmer):


        hash = self.str_to_mph.get(kmer)

        if hash == None:
            return

        # get address
        node = self.nodeDict.get(hash)
        nodes_to_be_updated_OutM = dict()

        if node != None:



            #We need to create the next adjacent outgoing in the tree as the new root for both Cases
            # get the correct outgoing node (connected in the tree) from out matrix and make it the new root node

            #if it gives an exception, then the qquery kmer is dynamically inserted
            try:
                in_edges = self.inOutMatrix.return_column_headers_In(kmer, self.str_to_mph)
            except:
                in_edges = []

            potentialStrings = {}
            k = len(kmer)


            for c in in_edges:  # c = 'T'
                valueDeque = deque(kmer)  # 'C G A'
                # rotate the list left (kmer)
                valueDeque.rotate(1)  # 'A C G'
                col_to_be_updated = valueDeque[0]
                # convert back to list (kmer)
                valueList = list(collections.deque(valueDeque))
                # add column header to the end of the kmer value
                valueList[0] = c # 'C G T'
                # convert list back to a string
                kmer_string = ''.join(valueList)
                out_node_hash = self.str_to_mph.get(kmer_string)
                nodes_to_be_updated_OutM[out_node_hash] = col_to_be_updated


        return nodes_to_be_updated_OutM

    # this method gives out all the nodes(hash) that had an outgoing edge from kmer
    #So the In matrix needs to be updated for these nodes
    def outgoing(self, kmer):


        hash = self.str_to_mph.get(kmer)

        if hash == None:
            return

        # get address
        node = self.nodeDict.get(hash)
        nodes_to_be_updated_InM = dict()

        if node != None:



            #We need to create the next adjacent outgoing in the tree as the new root for both Cases
            # get the correct outgoing node (connected in the tree) from out matrix and make it the new root node

            #if it gives an exception, then the qquery kmer is dynamically inserted
            try:
                in_edges = self.inOutMatrix.return_column_headers_Out(kmer, self.str_to_mph)
            except:
                in_edges = []

            potentialStrings = {}
            k = len(kmer)


            for c in in_edges:  # c = 'C'
                valueDeque = deque(kmer)  # 'G A T'
                col_to_be_updated = valueDeque[0]
                # rotate the list left (kmer)
                valueDeque.rotate(-1)  # 'A T G'
                # convert back to list (kmer)
                valueList = list(collections.deque(valueDeque))
                # add column header to the end of the kmer value
                valueList[k-1] = c # 'C G T'
                # convert list back to a string
                kmer_string = ''.join(valueList)
                out_node_hash = self.str_to_mph.get(kmer_string)
                nodes_to_be_updated_InM[out_node_hash] = col_to_be_updated


        return nodes_to_be_updated_InM

    def delete(self, kmer):

        isAlreadyPresent = self.search(kmer, 0)

        if isAlreadyPresent:
            hash = self.str_to_mph.get(kmer)

            if hash == None:
                return False

            # get address
            node = self.nodeDict.get(hash)

            if self.isNode_a_leaf(kmer) == False:

                #Case 1: If node is root
                if node.isRoot:

                    # delete the root node from rootsList
                    self.treeRoots.remove(node)

                else:

                    #Case2: If node is not a root and also not a leaf
                    pass


                #We need to create the next adjacent outgoing in the tree as the new root for both Cases
                # get the correct outgoing node (connected in the tree) from out matrix and make it the new root node

                #exception happens if node to be deleted is a dynamic node
                try:
                    out_edges = self.inOutMatrix.return_column_headers_Out(kmer, self.str_to_mph)
                except:
                    out_edges = []
                    return True

                potentialStrings = {}
                k = len(kmer)

                for c in out_edges:  # c = 'T'
                    valueDeque = deque(kmer)  # 'A C G'
                    # rotate the list left (kmer)
                    valueDeque.rotate(-1)  # 'C G A'
                    # convert back to list (kmer)
                    valueList = list(collections.deque(valueDeque))
                    # add column header to the end of the kmer value
                    valueList[k - 1] = c  # 'C G T'
                    # convert list back to a string
                    kmer_string = ''.join(valueList)
                    out_node_hash = self.str_to_mph.get(kmer_string)



                    # Now potential strings = Hashes of ['CGT', 'CGA']

                    # get the node for the hash value of this new potential root 'CGT' and check if it has a
                    # parent pointing to the root that needs to be deleted

                    out_node = self.nodeDict.get(out_node_hash)  # Gives the tree node for 'CGT'

                    if out_node.parent == node:
                        # make the existing root node None
                        node = None

                        # this is the new root
                        out_node.parent = None
                        out_node.level = 0
                        out_node.isRoot = True
                        out_node.val = kmer_string
                        self.treeRoots.append(out_node)


                        #remove from self.nodeDict as it is used to do search query
                        self.nodeDict[hash] = None

                        #################################################################

                        #Also delete all references of this node from in and out matrices

                        #dict of all nodes whose columns need to be set 0 in Out Matrix
                        #          { 4 : 'C' }  where 4 is the hash of 'CGA' and 'C' is column
                        nodes_to_update_OutM = self.incoming(kmer)

                        #dict of all nodes whose columns need to be set 0 in In Matrix
                        #          { 4 : 'T' }  where 4 is the hash of 'CGA' and 'T' is column
                        nodes_to_update_InM = self.outgoing(kmer)

                        #################################################################




                        #print "Node has been deleted!"
                        return True

                return False
            else:
                #if node is a leaf
                if node in self.treeLeafs:
                    self.treeLeafs.remove(node)
                node = None





                #################################################################

                # Also delete all references of this node from in and out matrices

                # dict of all nodes whose columns need to be set 0 in Out Matrix
                #          { 4 : 'C' }  where 4 is the hash of 'CGA' and 'C' is column
                nodes_to_update_OutM = self.incoming(kmer)

                # dict of all nodes whose columns need to be set 0 in In Matrix
                #          { 4 : 'T' }  where 4 is the hash of 'CGA' and 'T' is column
                nodes_to_update_InM = self.outgoing(kmer)

                #################################################################



                #print "Leaf Node deleted!"
                return True



        #Update In and Out Matrix

        else:
            #print "Node not present. Can't perform deletion!"
            return False


    def isNode_a_leaf(self, kmer):


        isAlreadyPresent = self.search(kmer, 0)


        if isAlreadyPresent:
            hash = self.str_to_mph.get(kmer)

            if hash == None:
                return False

            # get address
            node = self.nodeDict.get(hash)

            if node.isRoot:
                return False

            if node in self.treeLeafs:
                return True


            #We need to check if the outgoing vertices of the nodes are root nodes, if yes, then,
            #it is a leaf node
            out_edges = self.inOutMatrix.return_column_headers_Out(kmer, self.str_to_mph)

            if len(out_edges) == 0:
                return True

            potentialStrings = {}
            k = len(kmer)

            for c in out_edges:  # c = 'T'
                valueDeque = deque(kmer)  # 'A C G'
                # rotate the list left (kmer)
                valueDeque.rotate(-1)  # 'C G A'
                # convert back to list (kmer)
                valueList = list(collections.deque(valueDeque))
                # add column header to the end of the kmer value
                valueList[k - 1] = c  # 'C G T'
                # convert list back to a string
                kmer_string = ''.join(valueList)
                out_node_hash = self.str_to_mph.get(kmer_string)
                # potentialStrings.update({out_node_hash: kmer_string})


                # Now potential strings = Hashes of ['CGT', 'CGA']

                # get the node for the hash value of this new potential root 'CGT' and check if it has a
                # parent pointing to the root that needs to be deleted

                out_node = self.nodeDict.get(out_node_hash)  # Gives the tree node for 'CGT'

                if out_node.isRoot:
                    print "Node is a leaf"
                    return True

            return False

        else:
            print "Node not present. Can't be a leaf."
            return False




    def search(self, input_kmer,searchCount):
        #Input Kmer = ['GAT']
        #get hash
        hash = self.str_to_mph.get(input_kmer)

        if hash == None:
            return False

        #get address
        node = self.nodeDict.get(hash)

        if node == None:
            return False

        if node.isRoot and searchCount == 0:
            return True

        #get the concatonated parent

        #all incoming edges
        #['C', 'T'] from incoming matrix
        in_edges = self.inOutMatrix.return_column_headers_In(input_kmer, self.str_to_mph)

        potentialStrings = {}

        for c in in_edges:
            valueDeque = deque(input_kmer)
            # rotate the list left (kmer)
            valueDeque.rotate(1)
            # convert back to list (kmer)
            valueList = list(collections.deque(valueDeque))
            # add column header to the end of the kmer value
            valueList[0] = c
            # convert list back to a string
            kmer_string = ''.join(valueList)
            parentHash = self.str_to_mph.get(kmer_string)
            potentialStrings.update({parentHash: kmer_string})


            #Now potential strings = Hashes of ['CGA', 'GGA']
        true_parent_hash = node.parent.val

        if node.parent.isRoot == True:
            #print "Found"
            return True

        if true_parent_hash not in potentialStrings.keys():
            #print("NOT FOUND IN MATRIX")
            return False

        elif searchCount > self.maxTreeHeight:
            # print("EXCEEDED TREE HEIGHT - NOT FOUND")
            return False

        else:
            parent_string = potentialStrings[true_parent_hash]
            searchCount+=1
            return self.search(parent_string,searchCount)









    def getOutNodes(self, curr, prevVal, kmerList, newLevel):

        # get header columns marked as 1
        # prevalue is the value in the header of the out matrix: e.g. a,t,c,g
        # put it into a deque in order to manipulate the list (rotate the list)
        valueDeque = deque(kmerList)
        # rotate the list left (kmer)
        # print(valueDeque)
        valueDeque.rotate(-1)
        # print(valueDeque)

        # convert back to list (kmer)

        valueList = list(collections.deque(valueDeque))
        # add column header to the end of the kmer value
        valueList[-1] = prevVal
        # convert list back to a string
        kmer_string = ''.join(valueList)

#         if (kmer_string == 'CTA'):
# #            print "here"

        ##################
        hashValue = self.str_to_mph.get(kmer_string)
        ###################

        # check the hash of the value = hashValue
        # !!!!!still need to do
        #

        # create a new node (internal node)with hash value in value
        if self.visitedList[hashValue] == None:
            newNode = TreeNode(hashValue, newLevel, curr, False)
            self.nodeDict.update({hashValue: newNode})
            # increase the number of kmers for the check above (that is, cehcking if we have visited all the kmers)
            self.countNumberofKmers += 1
        else:
            newNode = None

        # mark the kmer as visited

        # add the node to the dictionary with the hashvalue


        # insert recursively
        # update curr pointer
        # visited for new node is
        curr = newNode
        #print("here")
        # print(kmer_string)
        # call recursively until we either reach max height or reach all the nodes visited
        return curr, kmer_string
        # }

    # loop through all the true values in in matrix to add to the forest
    def getInNodes(self, curr, prevVal, kmerList, newLevel):
        # get header columns marked as 1
        # prevalue is the value in the header of the out matrix: e.g. a,t,c,g
        # put it into a deque in order to manipulate the list (rotate the list)
        valueDeque = deque(kmerList)

        # rotate the list left (kmer)
        valueDeque.rotate(1)
        # print(kmerList)

        # convert back to list (kmer)
        valueList = list(collections.deque(valueDeque))
        # add column header to the end of the kmer value
        valueList[0] = prevVal
        # convert list back to a string
        kmer_string = ''.join(valueList)

        ##################
        hashValue = self.str_to_mph.get(kmer_string)

        ###################

        # check the hash of the value = hashValue
        # !!!!!still need to do
        #
        # create a new node (internal node)with hash value in value
        newNode = TreeNode(hashValue, newLevel, curr, False)

        # mark the kmer as visited
        newNode.visit = True
        # add the node to the dictionary with the hashvalue
        self.nodeDict.update({hashValue: newNode})
        # increase the number of kmers for the check above (that is, cehcking if we have visited all the kmers)
        self.countNumberofKmers += 1
        # insert recursively
        # update curr pointer
        curr = newNode
        return curr, kmer_string

