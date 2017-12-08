import GraphConstruction
import In_Out_Matrix
import StaticRabinKarp as rk
import perfection.czech
import MinimalPerfectHashing as mph
import time
from TreeNode import *
from Forest import *

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

    firstK = (firstK[:k] + '..') if len(firstK) > k else firstK

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

    #nodes_to_update_OutM = newForest.incoming('AGC')
    nodes_to_update_InM = newForest.outgoing('ATA')

    isDeleted = newForest.delete('CGA')
    #print isDeleted



if __name__ == "__main__":
    main()
