#!/usr/bin/env python

import argparse
import gudhi as gd

#points = gd.read_points_from_off_file(off_file = args.file)

from sys import argv
from sys import float_info

from math import sqrt

f=open(argv[1],"r")

if len(argv)>2:
    threshold=float(argv[2])
else:
    threshold=float_info.max

points=[]

for line in f:
    vals = line.split(" ")
    no_vals = len(vals)
    density=float(vals[no_vals-1])
    if density<=threshold:
        points.append([float(vals[i]) for i in range(no_vals-1)])

f.close()

print("Compute alpha persistence on ",len(points), "points")
alpha_complex = gd.AlphaComplex(points = points)

simplex_tree = alpha_complex.create_simplex_tree()

print("Number of simplices=", simplex_tree.num_simplices())

diag = simplex_tree.persistence()

ppairs=[]

for dot in diag:
    birth = dot[1][0]
    death = dot[1][1]
    if death<float_info.max:
        ppairs.append([birth,death])
    else:
        print("Found infinite point")

ppairs.sort(key=lambda x : x[0])

output = open("Gudhi_pers_pairs.txt","w")

for p in ppairs:
    print(sqrt(p[0]),sqrt(p[1]),file=output)

output.close()



