import GraphConstruction
import In_Out_Matrix
import StaticRabinKarp as rk
import perfection.czech
import MinimalPerfectHashing as mph
import time
from TreeNode import *
from Forest import *
from test import *

def compute_permutation(permuted_set):
    new_permuted_set = set()
    list = ['A', 'C', 'G', 'T']
    for char in list:
        for str in permuted_set:
            new_permuted_set.add(char + str)
    return new_permuted_set


def main():
    k = input("Enter the value of K:")
    start_time  = time.time()
   # permuted_set = set(['A', 'C', 'G', 'T'])
   # for i in range(k-1):
    #    permuted_set = compute_permutation(permuted_set)
    #print(permuted_set)
    #print(len(permuted_set))
    #hash_set = StaticRabinKarp.static_rabin_karp(permuted_set)
    #print(hash_set)
    #perfect_hash = perfection.make_hash(hash_set.values())
    fileData = ""
    firstK = ""
    '''
    with open('5mb.fastq', 'r') as f:
        for count, line in enumerate(f, start=3):
            if count % 4 == 0:
                fileData += line.replace('\n', '')
                # print(line)
    '''
    #firstK is the first K characters needed to create the first root, to use just the out matrix and not the in matrix
    with open('5mb.fastq', 'r') as f:
        for count, line in enumerate(f, start=3):
            if count % 4 == 0:
                firstK += line.replace('\n', '')
                if len(firstK) >= k:
                    break

    firstK = (firstK[:k]) if len(firstK) > k else firstK

    print "FirstK"
    print firstK

    # fileData = open("test.txt", "r").read().replace('\n', '')
    # print(fileData)
    fileData = "ACGATGATCAGTAGCATGATCAGT"
    fileData = "ACGATCGATAGCTA"
    #fileData = open("test.txt", "r").read().replace('\n', '')
    print "De-bruijn Construction started"
    nodes, edges, string_hash_map, node_hash_values = GraphConstruction.construct_de_bruijn_graph(fileData, k)
    print "De-bruijn Construction completed"
    # print(nodes)
    # print(string_hash_map)
    # print("No. of nodes in de brujin graph", len(nodes))
   # print(edges)

    inOut = In_Out_Matrix()
    in_matrix, out_matrix = inOut.construct_in_out_matrix(nodes, edges, string_hash_map)
    initKmer = fileData[0:k-1]
    print "Forest Construction started"
    newForest = Forest(inOut, k, initKmer, string_hash_map)
    print "Forest Construction ended"

    print("--- %s seconds for hashing and forest construction ---" % (time.time() - start_time))
    print("Enter the nodes you want to add, with the size k")
    s = raw_input("Enter the nodes you want to add, seprated by space, with the size k")
    start_time = time.time()
    inserted_nodes = [i for i in s.split()]
    str_to_rk = rk.static_rabin_karp(set(inserted_nodes))
    new_str_to_rk = str_to_rk.copy()
    new_str_to_rk.update(node_hash_values)
    print("computing mph")
    print(mph.mph(new_str_to_rk))
    print("--- %s seconds for INSERTION ---" % (time.time() - start_time))
    searchedNode = raw_input("Enter the node to search of size k:")


    isLeaf = newForest.isNode_a_leaf(searchedNode)
    result = newForest.search(searchedNode, 0)
    print result

    newForest.insert_node('GGG')
    result = newForest.search(searchedNode, 0)
    print result

    nodes_to_update_OutM = newForest.incoming('ATA')
    nodes_to_update_InM = newForest.outgoing('ATA')

    isDeleted = newForest.delete('CGA')
    #print isDeleted



def test1():

    #Step 1: Input String

    fileData = 'ACGATCGATAGCTA'

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

    insert_kmer = 'TTT'
    print "Insert '" + insert_kmer + "' if is not present in the graph"
    isInserted =  newForest.insert_node(insert_kmer)
    if isInserted:
        print insert_kmer + " : is successfully inserted!"
    else:
        print insert_kmer + " : already present. Insertion unsuccessful!"

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







if __name__ == "__main__":
    test()
