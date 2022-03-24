#!/usr/bin/python3
# Script for making a OBJ mesh file with the 
# thickness-breadth pairs

import sys
import json
import math
# import matplotlib.pyplot as plt

def write_balls(dimension, measure, tb_pairs, divs=16):
    """
    Write 6 OBJ files with the thicknes and breadth balls of each dimension
    """
    suffix = f't{dimension}' if measure == 'thickness' else f'b{dimension}'
    fn = f'{tb_pairs["meta"]["filename"][:-4]}_{suffix}.obj'
    center = measure + '_center'
    offset = 0
    inc = 3.14159265359 / divs
    with open(fn, 'w') as f:
        for pair in tb_pairs['pairs']:
            if pair['dimension'] == dimension and pair[measure] > -1:
                print(pair)
                vertices = (2*divs) * (divs-1) + 2 # number of vertices in the sphere
                # write vertices
                r = pair[measure] - 0.5
                p = [x + 0.5 for x in pair[center]]
                f.write(f"v {p[0]}  {p[1]}  {p[2] + r}\n") # top vertex
                f.write(f"v {p[0]}  {p[1]}  {p[2] - r}\n") # bottom vertex
                for u in range(1, divs): # u in {pi/4, 2*pi/4, 3*pi/4}
                    r_sin_u = r * math.sin(u * inc)
                    z       = r * math.cos(u * inc)
                    for v in range(2*divs): # v in {0, pi/4, ..., 7*pi/4}
                        x = r_sin_u * math.cos(v * inc)
                        y = r_sin_u * math.sin(v * inc)
                        f.write(f"v {p[0]+x}  {p[1]+y}  {p[2]+z}\n")
                # write triangle faces incident to the top and bottom vertices
                for v in range(2*divs):
                    a = offset + 3 + v
                    b = offset + 3 + (v+1)%(2*divs)
                    c = offset + vertices - 2*divs + 1 + (v+1)%(2*divs)
                    d = offset + vertices - 2*divs + 1 + v
                    f.write(f"f {offset+1} {a} {b}\n")
                    f.write(f"f {offset+2} {c} {d}\n")
                # write square faces
                for u in range(1,divs-1):
                    for v in range(2*divs):
                        a = offset + 2 + 2*divs*(u-1) + v + 1
                        b = offset + 2 + 2*divs*u     + v + 1
                        c = offset + 2 + 2*divs*u     + (v+1)%(2*divs) + 1
                        d = offset + 2 + 2*divs*(u-1) + (v+1)%(2*divs) + 1
                        f.write(f"f {a} {b} {c} {d}\n")
                offset += vertices

    print(f'{fn} written')

if len(sys.argv) < 2:
    sys.exit("Usage: python3 tb-balls.py filename")

with open(sys.argv[1]) as f:
    tb_pairs = json.load(f)

for q in range(3):
    write_balls(q, 'thickness', tb_pairs)
    write_balls(q, 'breadth',   tb_pairs)