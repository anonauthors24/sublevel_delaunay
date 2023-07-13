import tadasets 

from sys import argv

from random import random

n = int(argv[1])
d = int(argv[2])
amb = int(argv[3])



#print("n=",n, " d=",d, " amb=",amb)
if(amb>d+1):
    points = tadasets.dsphere(n,d,ambient=amb,noise=0.01)
else:
    points = tadasets.dsphere(n,d,noise=0.01)

for p in points:
    for coor in p:
        print(coor,end=' ')
    print(random())
