from tkinter import W
from unittest import result
import numpy as np

import matplotlib.pyplot as plt

with open('Life_X_Query_Seq.txt') as f:
    lines = f.readlines()
life_x_name, life_x_seq = lines

with open('100_known_species_Seq.txt') as f:
    species_lines = f.readlines()

score = []
# parameters ===========================
match = 2
mismatch = -1
gap_penalty = -3
# =======================================
count = 0
best = 0
best_loc = (0, 0)
best_name = ''
update_flag = 0
str1 = ''
str2 = ''

# visualization ===========================
vis_row_length = 0
vis_col_length = 0
vis_loc = []
# =========================================

while count < len(species_lines):
    print( str(int(count/2 + 1))+ '/' +  str(int(len(species_lines)/2)))
    species_name = species_lines[count]
    species_seq = species_lines[count+1]
    # local alignment
    mt = np.zeros((len(life_x_seq)+1, len(species_seq)+1))
    for row in range(1, len(mt)):
        for col in range(1, len(mt[0])):
            if species_seq[col-1] == life_x_seq[row-1]:
                s = match
            else:
                s = mismatch
            mt[row, col] = max(0, mt[row-1, col-1]+s, mt[row,
                               col-1]+gap_penalty, mt[row-1, col]+gap_penalty)
            # new best
            if mt[row][col] >= best:
                best = mt[row][col]
                best_loc = (row, col)
                str1 = ''
                str2 = ''
                update_flag = 1
                best_name = species_name
                vis_row_length = len(mt)
                vis_col_length = len(mt[0])
                vis_loc = []
    # traceback
    if update_flag == 1:
        best_loc_row, best_loc_col = best_loc
        str1 += species_seq[best_loc_col - 1]
        str2 += life_x_seq[best_loc_row - 1]
        loc_row, loc_col = best_loc
        vis_loc.append((loc_row, loc_col))
        while mt[loc_row][loc_col] != 0:
            d = mt[loc_row-1][loc_col-1]
            h = mt[loc_row][loc_col-1]
            v = mt[loc_row-1][loc_col]
            if max(d, h, v) == d:
                str1 += life_x_seq[loc_row - 2]
                str2 += life_x_seq[loc_row - 2]
                loc_col = loc_col - 1
                loc_row = loc_row - 1
            elif max(d, h, v) == h:
                str1 += species_seq[loc_col - 2]
                str2 += '-'
                loc_col = loc_col - 1
            elif max(d, h, v) == v:
                str1 += '-'
                str2 += life_x_seq[loc_row - 2]
                loc_row = loc_row - 1
            vis_loc.append((loc_row, loc_col))
        str1 = str1[::-1]
        str2 = str2[::-1]
        update_flag = 0
    count += 2
result = 'Most similar species: '+ best_name + '\n' + 'Alignment: ' +  '\n'  + str1 + '\n' + str2
print(best_name + '\n'  + str1 + '\n' + str2)
output = './local_alignment_result.txt'
with open(output, W) as output_file:
    output_file.write(result)
# print(vis_loc)
# visualize
vis_mt = np.zeros((vis_row_length, vis_col_length))
for point in vis_loc:
    vis_mt[point[0], point[1]] = 1
plt.spy(vis_mt, markersize = 2)
plt.title('Best local alignment trackback visualizaiton')
plt.xlabel(best_name + ' sequence')
plt.ylabel('life X sequence')
plt.show()