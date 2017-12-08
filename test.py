import GraphConstruction
import In_Out_Matrix
import StaticRabinKarp as rk

import MinimalPerfectHashing as mph
import time
from TreeNode import *
from Forest import *

def test():

    #Step 1: Input String

    fileData = 'ACGATCGATAGCTA'

    print "####################################"
    print "Input String : " + fileData
    print "####################################"


    k = 4
    print "\n\n"

    #Step 2: Create De-bruijn graph
    nodes, edges, string_hash_map, node_hash_values = GraphConstruction.construct_de_bruijn_graph(fileData, k)

    # Step 3: Construct and Print In and Out Matrix
    inOut = In_Out_Matrix()
    in_matrix, out_matrix = inOut.construct_in_out_matrix(nodes, edges, string_hash_map)
    print "\n\n"




    #Step 4:Initialize Forest
    initKmer = fileData[0:k - 1]
    print "####################################"
    print "Covering Forest Construction"
    print "####################################"
    print "\n"
    print "Forest Construction started"
    newForest = Forest(inOut, k, initKmer, string_hash_map)
    print "Forest Initialized"
    print "\n"
    print "Roots in the Forest :"
    list =  newForest.treeRoots
    for node in list:
        print node.val

    print "\n\n"
    #Step 5:

    print "Nodes in Debruijn Graph after Initialization:"
    print nodes

    print "\n\n"
    #Step 6: Membership Query

    #-------> Search for something that is not present
    print "####################################"
    print "Membership Query"
    print "####################################"
    print "\n"
    query_kmer = 'GGG'

    print "Search if '" + query_kmer + "' is present in the graph?"
    isPresent =  newForest.search(query_kmer, 0)
    if isPresent:
        print query_kmer + " : Found!"
    else:
        print query_kmer + " : Not Found!"

    print "\n\n"


    query_kmer = 'GCT'
    print "Search if '" + query_kmer + "' is present in the graph?"
    isPresent =  newForest.search(query_kmer, 0)
    if isPresent:
        print query_kmer + " : Found!"
    else:
        print query_kmer + " : Not Found!"

    print "\n\n"



    #Step 6: Dynamic Insert


    print "####################################"
    print "Dynamic Insertion"
    print "####################################"
    print "\n"
    insert_kmer = 'GAT'

    print "Insert '" + insert_kmer + "' if is not present in the graph"
    isInserted =  newForest.insert_node(insert_kmer)
    if isInserted:
        print insert_kmer + " : is successfully inserted!"
    else:
        print insert_kmer + " : already present. Insertion unsuccessful!"

    print "\n\n"

    insert_kmer = 'GGG'
    print "Insert '" + insert_kmer + "' if is not present in the graph"
    isInserted =  newForest.insert_node(insert_kmer)
    if isInserted:
        print insert_kmer + " : is successfully inserted!"
    else:
        print insert_kmer + " : already present. Insertion unsuccessful!"

    print "\n\n"

    insert_kmer = 'GGG'
    print "Insert '" + insert_kmer + "' if is not present in the graph"
    isInserted =  newForest.insert_node(insert_kmer)
    if isInserted:
        print insert_kmer + " : is successfully inserted!"
    else:
        print insert_kmer + " : already present. Insertion unsuccessful!"

    print "\n\n"

    insert_kmer = 'TTT'
    print "Insert '" + insert_kmer + "' if is not present in the graph"
    isInserted =  newForest.insert_node(insert_kmer)
    if isInserted:
        print insert_kmer + " : is successfully inserted!"
    else:
        print insert_kmer + " : already present. Insertion unsuccessful!"

    print "\n\n"

    query_kmer = 'GGG'
    print "Search if '" + query_kmer + "' is present in the graph?"
    isPresent =  newForest.search(query_kmer, 0)
    if isPresent:
        print query_kmer + " : Found!"
    else:
        print query_kmer + " : Not Found!"

    print "\n\n"


    #Step 6: Dynamic Deletion

    #-------> Search for something that is not present
    print "####################################"
    print "Dynamic Deletion"
    print "####################################"
    print "\n"
    deletion_kmer = 'CCC'

    print "Delete '" + deletion_kmer + "' if present in the graph"
    isDeleted =  newForest.delete(deletion_kmer)
    if isDeleted:
        print deletion_kmer + " : is successfully deleted!"
    else:
        print deletion_kmer + " : not present. Deletion unsuccessful!"

    print "\n\n"

    deletion_kmer = 'GGG'
    print "Delete '" + deletion_kmer + "' if present in the graph"
    isDeleted =  newForest.delete(deletion_kmer)
    if isDeleted:
        print deletion_kmer + " : is successfully deleted!"
    else:
        print deletion_kmer + " : not present. Deletion unsuccessful!"

    print "\n\n"


    print "Search a k-mer Before and After Deletion"
    print "####################################"

    query_kmer = 'GCT'
    print "Search if '" + query_kmer + "' is present in the graph?"
    isPresent =  newForest.search(query_kmer, 0)
    if isPresent:
        print query_kmer + " : Found!"
    else:
        print query_kmer + " : Not Found!"

    print "\n\n"

    deletion_kmer = 'GCT'
    print "Delete '" + deletion_kmer + "' if present in the graph"
    isDeleted =  newForest.delete(deletion_kmer)
    if isDeleted:
        print deletion_kmer + " : is successfully deleted!"
    else:
        print deletion_kmer + " : not present. Deletion unsuccessful!"

    print "\n\n"

    query_kmer = 'GCT'
    print "Search if '" + query_kmer + "' is present in the graph?"
    isPresent =  newForest.search(query_kmer, 0)
    if isPresent:
        print query_kmer + " : Found!"
    else:
        print query_kmer + " : Not Found!"

    print "\n\n"

    insert_kmer = 'GCT'
    print "Insert '" + insert_kmer + "' if is not present in the graph"
    isInserted =  newForest.insert_node(insert_kmer)
    if isInserted:
        print insert_kmer + " : is successfully inserted!"
    else:
        print insert_kmer + " : already present. Insertion unsuccessful!"

    print "\n\n"




def test_2():

    #Step 1: Input String

    fileData = 'ACGATCGATAGCTA'
    fileData = open("test2.txt", "r").read().replace('\n', '')

    print "####################################"
    print "Input String : " + fileData
    print "####################################"


    k = 4
    print "\n\n"

    #Step 2: Create De-bruijn graph
    nodes, edges, string_hash_map, node_hash_values = GraphConstruction.construct_de_bruijn_graph(fileData, k)

    # Step 3: Construct and Print In and Out Matrix
    inOut = In_Out_Matrix()
    in_matrix, out_matrix = inOut.construct_in_out_matrix(nodes, edges, string_hash_map)
    print "\n\n"




    #Step 4:Initialize Forest
    initKmer = fileData[0:k - 1]
    print "####################################"
    print "Covering Forest Construction"
    print "####################################"
    print "\n"
    print "Forest Construction started"
    newForest = Forest(inOut, k, initKmer, string_hash_map)
    print "Forest Initialized"
    print "\n"
    print "Roots in the Forest :"
    list =  newForest.treeRoots
    for node in list:
        print node.val

    print "\n\n"
    #Step 5:

    print "Nodes in Debruijn Graph after Initialization:"
    print nodes

    print "\n\n"
    #Step 6: Membership Query




if __name__ == "__main__":
    test_2()