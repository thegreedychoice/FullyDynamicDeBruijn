import GraphConstruction
import In_Out_Matrix
import StaticRabinKarp as rk

import MinimalPerfectHashing as mph
import time
from TreeNode import *
from Forest import *

def test_simulation(fileData, k):

    #Step 1: Input String

    #fileData = 'ACGATCGATAGCTA'

    print "####################################"
    print "Input String : " + fileData
    print "####################################"


    #k = 4
    print "\n\n"

    #Step 2: Create De-bruijn graph
    nodes, edges, string_hash_map, node_hash_values = GraphConstruction.construct_de_bruijn_graph(fileData, k)

    # Step 3: Construct and Print In and Out Matrix
    inOut = In_Out_Matrix()
    in_matrix, out_matrix = inOut.construct_in_out_matrix(nodes, edges, string_hash_map)
    print "\n\n"

    print("In Matrix")
    print(inOut.in_matrix)
    print("************************************")
    print("Out Matrix")
    print(inOut.out_matrix)


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


    k = 10
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
        try:
            if node.val == None:
                print "None"
                pass
        except:
            pass
        print node.val

    print "\n\n"
    #Step 5:

    print "Nodes in Debruijn Graph after Initialization:"
    print nodes

    print "\n\n"
    #Step 6: Membership Query

def dynamic_operations(newForest, k):
    print "########################################################################################"
    print ""

    print "Test Dynamic Operations!"
    print ""
    print "########################################################################################"
    print ""

    while True:
        print "Press 1 to insert a kmer"
        print "Press 2 to delete a kmer"
        print "Press 3 to search a kmer"
        print "Press 4 to exit!"

        option = int(raw_input())

        if option == 1:
            insert_kmer = raw_input("Enter the kmer to be inserted :")
            while len(insert_kmer) != (k-1):
                print "Please enter the kmer of size : ", k-1
                insert_kmer = raw_input("Enter the kmer to be inserted :")

            print "Insert '" + insert_kmer + "' if is not present in the graph"
            isInserted = newForest.insert_node(insert_kmer)

            if isInserted:
                print insert_kmer + " : is successfully inserted!"
            else:
                print insert_kmer + " : already present. Insertion unsuccessful!"

            print "\n\n"


        elif option == 2:
            deletion_kmer = raw_input("Enter the kmer to be deleted: ")
            while len(deletion_kmer) != (k-1):
                print "Please enter the kmer of size : ", k-1
                deletion_kmer = raw_input("Enter the kmer to be deleted: ")

            print "Delete '" + deletion_kmer + "' if present in the graph"
            isDeleted = newForest.delete(deletion_kmer)
            if isDeleted:
                print deletion_kmer + " : is successfully deleted!"
            else:
                print deletion_kmer + " : not present. Deletion unsuccessful!"

            print "\n\n"

        elif option == 3:
            query_kmer = raw_input("Enter the query kmer : ")
            while len(query_kmer) != (k-1):
                print "Please enter the kmer of size : ", k-1
                query_kmer = raw_input("Enter the query kmer : ")

            print "Search if '" + query_kmer + "' is present in the graph?"
            isPresent = newForest.search(query_kmer, 0)
            if isPresent:
                print query_kmer + " : Found!"
            else:
                print query_kmer + " : Not Found!"

            print "\n\n"

        else:
            break






def interactive_test():

    print "########################################################################################"
    print ""

    print "Fully Dynamic de Bruijn implementation!"
    print ""
    print "########################################################################################"
    print ""

    while True:
        print "Press 1 to test for Input file"
        print "Press 2 to test for an Input String"
        print "Press 3 to run the project on a simulation"
        print "Press 4 to exit the program!"
        option = int(raw_input())

        if option == 1:
            print "The project can be tested for any user provided input file (fastq) "
            print "The file needs to be in the current directory"
            print "Note : If k > 5, the construction may take a long time!"
            print "\n"

            print "Enter the name of the test file (with extension) : "
            file_name = raw_input()
            k = input("Enter the value of K:")
            #start_time = time.time()

            fileData = ""
            firstK = ""

            with open(file_name, 'r') as f:
                for count, line in enumerate(f, start=3):
                    if count % 4 == 0:
                        fileData += line.replace('\n', '')
                        # print(line)

            # Step 2: Create De-bruijn graph
            nodes, edges, string_hash_map, node_hash_values = GraphConstruction.construct_de_bruijn_graph(fileData, k)

            # Step 3: Construct and Print In and Out Matrix
            inOut = In_Out_Matrix()
            in_matrix, out_matrix = inOut.construct_in_out_matrix(nodes, edges, string_hash_map)
            print "\n\n"
            print "####################################"
            print "Dynamic Hashing Implemented"
            print "####################################"
            print "\n"

            print "####################################"
            print "IN/OUT Matrices Constructed"
            print "####################################"
            print "\n"

            # Step 4:Initialize Forest
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
            list = newForest.treeRoots
            for node in list:
                try:
                    print node.val
                except:
                    pass

            print "\n\n"
            # Step 5:

            print "Nodes in Debruijn Graph after Initialization:"
            print nodes

            dynamic_operations(newForest, k)


        elif option == 2:
            print "Enter the Input String: "
            fileData = raw_input()
            k = input("Enter the value of K:")
            #start_time = time.time()

            # Step 2: Create De-bruijn graph
            nodes, edges, string_hash_map, node_hash_values = GraphConstruction.construct_de_bruijn_graph(fileData, k)

            # Step 3: Construct and Print In and Out Matrix
            inOut = In_Out_Matrix()
            in_matrix, out_matrix = inOut.construct_in_out_matrix(nodes, edges, string_hash_map)
            print "\n\n"
            print "####################################"
            print "Dynamic Hashing Implemented"
            print "####################################"
            print "\n"

            print "####################################"
            print "IN/OUT Matrices Constructed"
            print "####################################"
            print "\n"

            # Step 4:Initialize Forest
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
            list = newForest.treeRoots
            for node in list:
                try:
                    print node.val
                except:
                    pass

            print "\n\n"
            # Step 5:

            print "Nodes in Debruijn Graph after Initialization:"
            print nodes

            dynamic_operations(newForest, k)


        elif option == 3:
            fileData = "ACGATGATCAGTAGCATGATCAGT"
            k = 4
            test_simulation(fileData, k)
        else:
            break





if __name__ == "__main__":
    interactive_test()