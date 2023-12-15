import numpy as np

#Returns the DNA sequence of a file path
def get_DNA(path_to_file):
    with open(path_to_file, 'r') as f:
        lines = f.readlines()
        text = lines[0].strip('\n').upper()
    return text


#Returns both the pattern and text from a file path
def get_pattern_and_text(path_to_file):
    with open(path_to_file, 'r') as f:
        lines = f.readlines()
        text = lines[0].strip('\n').upper()
        pattern = lines[1].strip('\n').upper()

    return text, pattern

#Returns the scoring matrix given a file path
def get_scoring_matrix(path_to_mat):
    mat_dict = {}
    with open(path_to_mat, 'r') as f:
        l1 = f.readline()
        l1 = l1.split()
        # print(l1)
        l1 = l1[1:]
        dim = len(l1)
        delta_mat = np.zeros((dim, dim))
        for (k, w) in enumerate(l1):
            # print(w, k)
            mat_dict[w] = k
        
        #mat should be symmetic, so could speed this up if necessay     
        for (j,l) in enumerate(f): 
            # print(l)
            l = l.split()[1:]
            for i, w in enumerate(l):
                delta_mat[i, j] = w

    return delta_mat, mat_dict


# returns the alphabet contained in an input file
def get_alphabet(path_to_alphabet):
    with open(path_to_alphabet, 'r') as f:
        alphabet = f.readline().strip('\n').upper()
        alphabet = list(set(alphabet)) #make sure only unique characters
        return alphabet


# '''
# Precomputes the m length bit masks for each character in the alphabet
# i.e. for each letter, there are m bits, where the ith bit is set if the letter occurs at that position
# in the pattern.
# p = ACGT      p = GATTACA
# A : 0001      A : 1010010
# C : 0010      C : 0100000
# G : 0100      G : 0000001
# T : 1000      T : 0001100
# Output: B, an np array of length alphabet with the appropriate int bit masks, where the ith index of B is 
# the ith index of the alphabet, alphabet
# '''   
def precomp_alphabet_bit_masks(pattern, alphabet):
    B = np.zeros((len(alphabet)), dtype=int) #initialized to zeos so no need to worry about values that dont appear in the text
    for idx, char in enumerate(pattern): 
        curr_value = B[alphabet.index(char)] #find the current bit values
        new_value= curr_value | 1 << idx #place a bit in the appropriate index
        B[alphabet.index(char)] = new_value


    return B

    



