def static_rabin_karp(permutation_set):
    static_permutation_hash_set = {}
    #print("I m computing rkp4")
    for str in permutation_set:
        val = 0
        char_num = 0
        for ch in str:
            val+=(ord(ch)-ord("A")+1)*(5**(len(str)-1-char_num))
            char_num+=1
        static_permutation_hash_set.update({str:val})

    return static_permutation_hash_set
