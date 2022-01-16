import numpy as np
import matplotlib.pyplot as plt
# ============================================================================
# function definition
# ============================================================================
def calculate_score(test_sequence):
    score = 0
    for str_index in range(1, len(test_sequence)):
        row = state_dict[test_sequence[str_index-1]]
        col = state_dict[test_sequence[str_index]]
        score += beta_matrix[row, col]
    score /= len(test_sequence)
    score = round(score, 2)
    return score
# ============================================================================
# Start
# ============================================================================
# file open
# ============================================================================
with open('CpG.txt') as cpg_file:
    cpg_sequence = cpg_file.readlines()
    cpg_sequence = cpg_sequence[0]
with open('Non_CpG.txt') as non_cpg_file:
    non_cpg_sequence = non_cpg_file.readlines()
    non_cpg_sequence = non_cpg_sequence[0]
with open('test_sequence.txt') as test_file:
    test_sequence = test_file.readlines()
    test_sequence = test_sequence[0]


# ============================================================================
# generate a+, a- transition matrix
# ============================================================================
state_list = ['A', 'T', 'C', 'G']
state_dict = {'A': 0, 'T' : 1, 'C' : 2, 'G' : 3}
positve_matrix = np.zeros((4,4))
negative_matrix = np.zeros((4,4))
for pre_char_index in range(len(state_list)): # row
    pre_char_positive_total = cpg_sequence.count(state_list[pre_char_index])
    pre_char_negative_total = non_cpg_sequence.count(state_list[pre_char_index])
    for next_char_index in range(len(state_list)): # col
        sub = state_list[pre_char_index] + state_list[next_char_index]
        positve_matrix[pre_char_index, next_char_index] = cpg_sequence.count(sub) / pre_char_positive_total
        negative_matrix[pre_char_index, next_char_index] = non_cpg_sequence.count(sub) / pre_char_negative_total
print('a+ transirion matrix = ')
print(positve_matrix)
print('a- transirion matrix = ')
print(negative_matrix)

# ============================================================================
# generate beta log likelihood ratios matrix 
# ============================================================================

beta_matrix = np.zeros((4,4))
beta_matrix = np.log(np.divide(positve_matrix, negative_matrix))
print('beta log likelihood ratios matrix = ')
print(beta_matrix)


# ============================================================================
# test
# ============================================================================
test_score = calculate_score(test_sequence)
print(test_score)
if test_score>=0:
    print('Belongs to the CpG islands sequence!')
else:
    print('Not belongs to the CpG islands sequence!')

# ============================================================================
# visualization
# cut cpg and non_cpg sequence into substring list every 70 char
# ============================================================================
n = 30
cpg_test_list = [cpg_sequence[i:i+n] for i in range(0, len(cpg_sequence), n)]
non_cpg_test_list = [non_cpg_sequence[i:i+n] for i in range(0, len(non_cpg_sequence), n)]
positve_score_list = []
negative_score_list = []
for item_index in range(len(cpg_test_list)):
    positve_score_list.append(calculate_score(cpg_test_list[item_index]))
for item_index in range(len(non_cpg_test_list)):
    negative_score_list.append(calculate_score(non_cpg_test_list[item_index]))

# print(positve_score_list)
# print(negative_score_list)
plt.hist(positve_score_list, alpha=0.5)
plt.hist(negative_score_list, alpha=0.5)
plt.title('length-normalised scores for divided training sequence')
plt.xlabel('Bits')
plt.show()
# plt.hist(data, bins=bins, alpha=0.5)
# =====================================================================
