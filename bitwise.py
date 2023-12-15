import sys
import numpy as np
import preprocess as pre

def BPM(text, pattern, Pek, k, alphabet):
    '''
    Myers bitwise DP implementation 
    Outputs only the end indices of any occurrences
    '''
    pattern_length = len(pattern)
    Pv = 2**pattern_length -1 #all 1's of length m  
    Mv = 0
    score = pattern_length
    text_length = len(text)

    matches = []
    ten_m = 1 << (pattern_length-1)
    for j in range(text_length):
        Ek = Pek[alphabet.index(text[j])]
        Xv = Ek | Mv
        Xh = (((Ek & Pv) + Pv) ^Pv) | Ek

        Ph = Mv | ~ (Xh|Pv) 
        Mh = Pv & Xh

        if Ph & ten_m:
            score += 1
        elif Mh & ten_m:
            score -= 1       

        Ph <<=  1
        Mh <<=  1 
        Pv = Mh | ~(Xv | Ph)
        Mv = Ph & Xv

        if score <= k:
            matches.append(j)
    return matches


def shift_or(text, pattern, alphabet, Pek):
    '''
    EXACT patten matching using bitwise parallelism, sped  up version of KMP.
    Basis for many other bitwise algorithms such as BMP, BPD
    Inputs: quite self explanatory, Pek is the precomputed bit masks for each character in alphabet
    Output: list of indices where matches end 
    '''
    m = len(pattern)
    D = 0

    match_indices = []
    for j in range(len(text)):
        idx = alphabet.index(text[j])
        ek = Pek[idx]

        D = ((D<<1 | 1)) & ek #if the m-1 bit is set and the cuent tj matches the m'th bit, then we have a match!
        if (D >> (m-1)) & 1: #if the mth bit is set, match!
            match_indices.append(j)

    return match_indices



if __name__ == '__main__':
    num_args = 3
    assert(len(sys.argv) == num_args+1)
    path_to_sequences =  sys.argv[1]
    path_to_alphabet = sys.argv[2]
    k = int(sys.argv[3])

    text, pattern = pre.get_pattern_and_text(path_to_sequences)
    alphabet = pre.get_alphabet(path_to_alphabet)

    B = pre.precomp_alphabet_bit_masks(pattern, alphabet)

    match_indices = []
    if k == 0:
        match_indices = shift_or(text, pattern, alphabet, B)
    else:
        match_indices = BPM(text, pattern, B, k, alphabet)

    print('\n--------------------\nTotal matches: \n', len(match_indices))

    # print("----------Matches at indices:--------")
    # print
    # print(match_indices)
    # print


    