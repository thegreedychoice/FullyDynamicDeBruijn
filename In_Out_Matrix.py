import BitArray2D
from BitArray2D import godel

class In_Out_Matrix:

    in_matrix = None
    out_matrix = None
    column_char_map = {}
    char_column_map = {}

    def __init__(self):
        self.in_matrix = None
        self.out_matrix = None
        self.column_char_map.update({'A': 0})
        self.column_char_map.update({'C': 1})
        self.column_char_map.update({'G': 2})
        self.column_char_map.update({'T': 3})
        self.column_char_map.update({'N': 4})
        self.char_column_map = dict((v, k) for k, v in self.column_char_map.iteritems())

    def construct_in_out_matrix(self, nodes, edges, string_hash_map):
        self.in_matrix = BitArray2D.BitArray2D(rows=len(nodes), columns=5)
        self.out_matrix = BitArray2D.BitArray2D(rows=len(nodes), columns=5)

        print(self.in_matrix.size())

        # for str in nodes:
        count = 0;
        for v1, v2 in edges:
            self.in_matrix.__setitem__((string_hash_map.get(v2), self.column_char_map.get(v1[0])), 1)
            self.out_matrix.__setitem__((string_hash_map.get(v1), self.column_char_map.get(v2[len(v2) - 1])), 1)
            count+=1;
            # print(count)

        print("In Matrix")
        print(self.in_matrix)
        print("************************************")
        print("Out Matrix")
        print(self.out_matrix)

        #print(string_hash_map.get("TAC"))
        #in_column_header = return_column_headers_In("TAC", in_matrix, out_matrix, string_hash_map, char_column_map)
        #out_column_header = return_column_headers_Out("TAC", in_matrix, out_matrix, string_hash_map, char_column_map)
        #print(In_column_header)
        #print(Out_column_header)
        return self.in_matrix,self.out_matrix

    def return_column_headers_In(self, kmer,str_to_mph):
        kmer_to_column_header = {}
        column_header_set = set()
        mph_val = str_to_mph.get(kmer)
        for j in range(5):
            if (self.in_matrix[godel(mph_val, j)]==1):
                column_header_set.add(self.char_column_map.get(j))

        return column_header_set


    def return_column_headers_Out(self, kmer,str_to_mph):
        kmer_to_column_header = {}
        column_header_set = set()
        mph_val = str_to_mph.get(kmer)
        print("mph val")
        print(mph_val)
        for j in range(5):
            if (self.out_matrix[godel(mph_val, j)]==1):
                column_header_set.add(self.char_column_map.get(j))

        return column_header_set




