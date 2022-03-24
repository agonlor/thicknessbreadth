#!/usr/bin/python3
# Script for plotting the thickness-breadth diagrams
# given a JSON file containing the thickness-breadth pairs

import sys
import json
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    sys.exit("Usage: python3 tb-diag.py filename")

with open(sys.argv[1]) as f:
    tb_pairs = json.load(f)

mark = ['o', '^', 's']
color = ['r', 'g', 'b']
for q in range(3):
    x = [pair['thickness'] for pair in tb_pairs['pairs'] if pair['dimension'] == q]
    y = [pair['breadth']   for pair in tb_pairs['pairs'] if pair['dimension'] == q]
    plt.scatter(x, y, marker = mark[q], facecolors='None', edgecolors=color[q])
plt.axhline(y=0, color='k', linestyle='--') # dashed line at y=0 to dinstinguish the pair (t, inf)
plt.xlabel('Thickness')
plt.ylabel('Breadth')
plt.xlim(left=0) # This is not working, probably because the following line
plt.axis('square')
plt.savefig('tb-diagram.png')
# plt.show()