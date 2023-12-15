import sys
import numpy as np
import preprocess as pre

def dynamic_Sellers(text, pattern, delta_mat, mat_dict, k, gap_pen):
    '''
    Implementation of the first approximate pattern matching algorithm based on the edit distance dynamic programming matrix
    Runs in O(mn) time as other naive DP algorithms would and O(mn) space 
    '''
    x_dim = len(text)+1
    y_dim = len(pattern)+1

    table = np.zeros((y_dim, x_dim))
    back = np.zeros((y_dim,x_dim))

    for row in range(1, y_dim): #base case first col, first row already handled by np.zero
        table[row, 0] = table[row-1, 0] + gap_pen
        back[row, 0] = 2 #down

    #Computes V(i,j) fo all pais
    for row in range(1, y_dim) :
        for col in range(1, x_dim):
   
            down = table[row-1, col] + gap_pen
            right = table[row, col-1] + gap_pen
            diagonal_value = delta_mat[mat_dict[text[col-1]], mat_dict[pattern[row-1]]]

            diagonal = (table[row-1, col-1] + diagonal_value) 

            value = min(diagonal, down, right)

            back[row, col] = np.argmin([right, down, diagonal])+1
            table[row, col] = value


    max_indices = np.argwhere(table[y_dim-1] <= k) #Find entries in last row with less than k total errors
    max_indices = max_indices.flatten().tolist()
    return table, back, max_indices



def dynamic_Ukkonen_cutoff(text, pattern, delta_mat, mat_dict, k, gap_pen):
    '''
    Implementation of the cutoff heurstic proposed by Ukkonen in 1985, which uses diagonal montonicity to discard diagonals/swaths of cells 
    that would be > k anyways
    Runs in O(mn) worst case and O(kn) average case. 
    Used often in more advanced/efficeint bitwise algorithms such as BPM as BPD
    Input: delta_mat is the scoring matrix, mat_dict is a dict from char ->index in scoring matrix
    Output: The filled out DP and backtracking matries and a list of indices where matches end
    '''
    n = len(text)
    m = len(pattern)

    U = np.zeros((m+1, n+1))
    back = np.zeros((m+1,n+1))

    top = min(k+1, m) #Keeps track of current diagonal

    #could do (1, top) too probs
    for i in range(top+1): #Initialize the first column
        U[i, 0] = i


    position_matches = []
    for col in range(1, n+1):
        for row in range(1, top+1):
            down = U[row-1, col] + gap_pen
            right = U[row, col-1] + gap_pen
            diagonal_value = delta_mat[mat_dict[text[col-1]], mat_dict[pattern[row-1]]]

            diagonal = (U[row-1, col-1] + diagonal_value) 
            value = min(diagonal, down, right)

            back[row, col] = np.argmin([right, down, diagonal])+1
            U[row, col] = value

        while U[top, col]> k: #while the current diagonal value is still greater than k, move up the diagonal
            top -= 1
        if top == m: #if we've reached the final diagonal, then there was a match
            # print('match ending at position j', col)
            position_matches.append(col)
        else:
            top += 1
            U[top, col] = k+1

    # print(U)
    return U, back, position_matches



def backtrack(row, col, back, s1, s2, text, pattern):
    '''
    Backtrack through the table to find the alignment configuration and start index
    Outputs the sequence alignment in reverse order
    '''
    if(row == 0): #if at top row, exit recursion
        return s1, s2, col
    # print("CUrr: rOW: ", row)
    # print('curr col: ', col)
    new_row = row
    new_col = col
    if(back[row, col]) == 3: #Diagonal
        # print('diagonal')
        new_row = row-1
        new_col = col-1
        # print('newx: ', new_x)
        s1 += text[new_col]
        s2 += pattern[new_row] 

    # down
    elif(back[row, col]) ==2: #Up
        # print('down')
        new_row = row -1
        s1 += '-'
        s2 += pattern[new_row]  
    
    #right
    elif(back[row, col] ==1): #Left
        # print('right')
        new_col = col-1
        s1 += text[new_col]
        s2 += '-'

    return backtrack(new_row, new_col, back, s1, s2, text, pattern)

if __name__ == '__main__':
    num_args = 5
    assert(len(sys.argv) == num_args+1)
    path_to_file =  sys.argv[1]
    path_to_mat = sys.argv[2]

    k = int(sys.argv[3])
    gap_pen = int(sys.argv[4])
    ukkonen = (sys.argv[5])

    text, pattern = pre.get_pattern_and_text(path_to_file)
    delta_mat, mat_to_dict = pre.get_scoring_matrix(path_to_mat)

    # k =3
    # gap_pen = -1
    final_row = len(pattern)
    pattern_matches = []

    C = np.zeros((1))
    C_back = np.zeros((1))
    C_indices =[]
    

    if(ukkonen == 'U'):
        C, C_back, C_indices = dynamic_Ukkonen_cutoff(text, pattern, delta_mat, mat_to_dict, k, gap_pen)
        print('------------------- Ukkonen alignments: ---------')

    else:
        C, C_back, C_indices = dynamic_Sellers(text, pattern, delta_mat, mat_to_dict, k, gap_pen)
    
        print('------------------- Sellers alignments: ------------')

    for i, max_idx in enumerate(C_indices):
        s1, s2, start_idx = backtrack(final_row, max_idx, C_back, '', '', text, pattern)
        curr_k = C[final_row, max_idx]
        
        pattern_matches.append([start_idx, curr_k, s1[::-1], s2[::-1]])
        # print('back table: ', C_back)
    print('Total matches: ', len(pattern_matches))

    
        
    
    

    
    
    
    

