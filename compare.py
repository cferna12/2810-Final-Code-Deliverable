import sys
from dynamic import *
from preprocess import *
from bitwise import *
import time
import random
import matplotlib.pyplot as plt
import json


#Sanity check for run-times: not used in final version
def virus_comparisions(max_k):
    print('--------------------------------------------------------------------------------')
    virus_file_paths = 'virus_data/sars-cov.txt'
    scoring_matrix_path = 'extra_data/standard.m'
    dna = get_DNA(virus_file_paths)
    print('len dna: ', len(dna))
    delta_mat, char_to_idx = get_scoring_matrix(scoring_matrix_path)
    motifs = ['CTAAC', 'TCTAAAC', 'ACGAAC', 'CTTAACAA', 'GATTACATGCAACCTTAACAA']


    alphabet_path = 'extra_data/alphabet.txt'
    alphabet = get_alphabet(alphabet_path)
    
    # gen_pat(dna, alphabet)
    u_times = {motif: [] for motif in motifs}
    c_times = {motif: [] for motif in motifs}
    bpm_times = {motif: [] for motif in motifs}
    for k in range(1,max_k+1):
        for motif in motifs:
            c_times, u_times, bpm_times = compute_times(motif, dna, char_to_idx, alphabet, delta_mat, c_times, u_times, bpm_times, k)
    print('ctimes: ', c_times)
    print('utimes: ', u_times)
    print('bpm times: ', bpm_times)



#Computes and returns the run times for each algorithm base on the given pattern, DNA, and k 
def compute_times(pattern, dna, char_to_idx, alphabet, delta_mat, k):
    Pek = precomp_alphabet_bit_masks(pattern, alphabet)
    c_start = time.time()
    C , C_back, naive_indices = dynamic_Sellers(dna, pattern, delta_mat, char_to_idx, k, 1)
    c_end = time.time()
    c_elapsed = c_end - c_start
    # print('c elapsed: ', c_elapsed)
    # c_times.append(c_elapsed)

    u_start = time.time()
    U, U_back, u_indices = dynamic_Ukkonen_cutoff(dna, pattern, delta_mat, char_to_idx, k, 1)
    u_end= time.time()
    u_elapsed = u_end - u_start

    # print('u elapsed: ', u_elapsed)
    # u_times.append(u_elapsed)

    bpm_start= time.time()
    bpm_indices = BPM(dna, pattern, Pek, k, alphabet)
    bpm_end= time.time()
    bpm_elapsed = bpm_end - bpm_start

    # bpm_times.append(bpm_elapsed)
    # print('bpm elapsed :', bpm_elapsed)

    #Sanity check to ensure theyre working
    # print("total indices sellers: ", len(naive_indices))
    # print('total ukkonene: ', len(u_indices))
    # print('bpm: ', len(bpm_indices))
    return c_elapsed, u_elapsed, bpm_elapsed


#Generates a random slice from the dna and randomly changes up to k positions in the slice
#Returns the new pattern 
def gen_pat(dna, alphabet, length, k):
    pat = ''
    # length = 30
    index = random.randint(0, len(dna))
    # print('rand index in dna: ', index)

    pat = (dna[index:index+length+1])
    # print('slice from dna: ', pat)
    
    # rand_k= random.randint(1, 5)
    # print('rand k: ', rand_k)
    indices_to_change = random.sample(range(length), k)
    # print('indices to change: ', indices_to_change)
    pat = list(pat)
    for idx in indices_to_change:
        new_char = random.randint(0, len(alphabet)-1)
        # print('new char', new_char)
        pat[idx] = alphabet[new_char]
    
    # for i in range()
    pat = ''.join(pat)
    # print('new pat: ', pat)

    return pat


#Calculates and saves the runtimes for varying m and k 
# Saved into jsons file to not have to run everytime 
def test_random_pattern(path_to_dna, path_to_alphabet):
   
    dna = get_DNA(path_to_dna)

    alphabet = get_alphabet(path_to_alphabet)
    # print('alphabet: ', alphabet)
    scoring_matrix_path = 'extra_data/standard.m'
    delta_mat, char_to_idx = get_scoring_matrix(scoring_matrix_path)
  
    k = 5

    num_pats = 3
    
    c_length_to_k = {}
    u_length_to_k = {}
    bpm_length_to_k = {}
    for m in range(10,40, 10): #loop though lengths
        m_c_times = []
        m_u_times = []
        m_bpm_times = []
        for k in range(1, int(m/2)+1,1): # loop thouhg k values fo current length
            c_times = []
            u_times = []
            bpm_times = []
            for i in range(num_pats): #calculate time with 5 different patterns and get avg
                pattern = (gen_pat(dna, alphabet, m, k))
            # pattern = 
                
                c_time, u_time,bpm_time = compute_times(pattern, dna, char_to_idx, alphabet, delta_mat,  k)
                c_times.append(c_time)
                u_times.append(u_time)
                bpm_times.append(bpm_time)
                # print(c_times, u_times, bpm_times)
            c_avg = np.average(c_times)
            # print('c avg: ', c_avg)
            #add the avevge time to the current k index in the time lists
            u_avg = np.average(u_times)
            bpm_avg = np.average(bpm_times)
            m_c_times.append(c_avg)
            m_u_times.append(u_avg)
            m_bpm_times.append(bpm_avg)
        #save the list of times (with k as indices), keyed on the length of the pattern
        c_length_to_k[m] = m_c_times
        u_length_to_k[m] = m_u_times
        bpm_length_to_k[m] = m_bpm_times

    # print(c_length_to_k)
    # print(u_length_to_k)
    # print(bpm_length_to_k)
    with open('jsons/c_times2.json', 'w') as f:
        json.dump(c_length_to_k, f)
    with open('jsons/u_times2.json', 'w') as f:
        json.dump(u_length_to_k, f)
    with open('jsons/bpm_times2.json', 'w') as f:
        json.dump(bpm_length_to_k, f)
    #gen3 seuences fo m = 8, 12, 16, 20, 24, 32
    #
    

#Plots the run xtimes of each algorithm
#Note that using sh does not seem to work in   
def plot_times():
    c_json = 'jsons/c_times.json'
    u_json = 'jsons/u_times.json'
    bpm_json = 'jsons/bpm_times.json'

    with open(c_json) as f:
        c_data = json.load(f)
    with open(u_json, 'r') as f:
        u_data = json.load(f)
    with open(bpm_json, 'r') as f:
        bpm_data = json.load(f)
    # print('c_data : ', c_data)

    plt.figure(figsize=(10, 6))
    for key in c_data.keys():
        k_values = range(1, len(c_data[str(key)])+1)
        # x1 = k_values

        check = str(key)
        y1 = c_data[check]
        # Data points of line 2
        # x2 = [1, 2, 3, 4, 5]
        y2 = u_data[check]
        # Data points of line 3
        # x3 = [1, 2, 3, 4, 5]
        y3 = bpm_data[check]
        # Plotting all lines with specifying labels
        # print('cecL : ', check)
        if check == '10':
            # print('hi')
            plt.plot(k_values, y1, 'g-', label='DP')
            plt.plot(k_values, y2, 'r-', label='Cut-Off')
            plt.plot(k_values, y3, 'b-', label='BPM')
        elif check == '20':
            plt.plot(k_values, y1, 'g--', label='DP m = 20')
            plt.plot(k_values, y2, 'r--', label='Cut-Off m = 20')
            plt.plot(k_values, y3, 'b--', label='BPM m = 20')
        else:
            plt.plot(k_values, y1, 'g:', label='DP m = 30')
            plt.plot(k_values, y2, 'r:', label='Cut-Off m = 30')
            plt.plot(k_values, y3, 'b:', label='BPM m = 30')


    # Adding legend, x and y labels, and titles for the lines
    plt.legend()
    plt.xlabel('k')
    plt.ylabel('Time in log(seconds)')
    plt.yscale('log')
    plt.title('Runtimes for m=10, 20, 30')
    # Displaying the plot
    plt.show() 
   


if __name__ == '__main__':
    path_to_dna = 'virus_data/mers.txt'
    path_to_alphabet = 'extra_data/alphabet.txt'
    test_random_pattern(path_to_dna, path_to_alphabet)
    # plot_times()
    
    # virus_comparisions(max_k)

