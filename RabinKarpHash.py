previous_hash_value = 0
first_char = ''


def rabin_karp_hash(node):
    current_hash_value = 0
    char_num = 0
    global previous_hash_value, first_char

    if previous_hash_value == 0:
        for ch in node:
            current_hash_value += (ord(ch)-ord("A")+1)*(5**(len(node)-1-char_num))
            char_num += 1
    else:
        new_char_offset = (ord(node[len(node)-1])-ord("A")+1) % 48112959837082048697
        temp1 = (previous_hash_value - ((ord(first_char)-ord("A")+1)*(5**(len(node)-1-char_num))% 48112959837082048697))
        current_hash_value = (temp1 *(5) + new_char_offset) % 48112959837082048697

    previous_hash_value = current_hash_value
    first_char = node[len(node) - 1]
    return current_hash_value
