import numpy as np
import matplotlib.pyplot as plt
from sys import exit

q = 10.
Lx = 1
Ly = 1

Lr = 2
LR = 5

c = [ 'red', 'blue', 'green', 'black' ]
def xy2R(x,y,q):
    return x/(1+q) + q*y/(1+q)
def xy2r(x, y,q):
    return y - x

def xy2Rr(xy,q):
    return [ xy2R(xy[0], xy[1],q), xy2r(xy[0],xy[1],q) ]

def Rr2x(R,r, q):
    return R - q*r/(1+q)
def Rr2y(R,r, q):
    return R + r/(1+q)
def Rr2xy(Rr,q):
    return [ Rr2x(Rr[0], Rr[1],q),Rr2y(Rr[0], Rr[1],q) ]


Rrcorners = [ [1,1], [-1,1], [-1,-1], [1,-1] ]
xycorners = [ Rr2xy(xyi,q) for xyi in Rrcorners ]

for i in range(4):
    j = i+ 1
    if j ==4: j=0
    #plt.subplot(1,2,1)
    plt.plot( [ Rrcorners[i][0], Rrcorners[j][0] ],  [ Rrcorners[i][1], Rrcorners[j][1] ], color=c[i])
for i in range(4):
    j = i+ 1
    if j ==4: j=0
    #plt.subplot(1,2,2)
    plt.plot( [ xycorners[i][0], xycorners[j][0] ],  [ xycorners[i][1], xycorners[j][1] ], color=c[i])
    #plt.plot( [ xy0corners[i][0], xy0corners[j][0] ],  [ xy0corners[i][1], xy0corners[j][1] ], color=c[i], ls=":")
#plt.subplot(1,2,1)
#plt.axhline(-2)
#plt.axhline(2)
#plt.axvline(2)
#plt.axvline(-2)
#plt.subplot(1,2,2)
#plt.axhline(-2)
#plt.axhline(2)
#plt.axvline(2)
#plt.axvline(-2)

d = 2
for i in range(len(Rrcorners)):
    Rrcorners[i][0] += d

xycorners = [ Rr2xy(xyi,q) for xyi in Rrcorners ]

for i in range(4):
    j = i+ 1
    if j ==4: j=0

    plt.plot( [ Rrcorners[i][0], Rrcorners[j][0] ],  [ Rrcorners[i][1], Rrcorners[j][1] ], color=c[i])
    plt.plot( [ xycorners[i][0], xycorners[j][0] ],  [ xycorners[i][1], xycorners[j][1] ], color=c[i])

d = 2
for i in range(len(Rrcorners)):
    Rrcorners[i][1] += d

xycorners = [ Rr2xy(xyi,q) for xyi in Rrcorners ]

for i in range(4):
    j = i+ 1
    if j ==4: j=0

    plt.plot( [ Rrcorners[i][0], Rrcorners[j][0] ],  [ Rrcorners[i][1], Rrcorners[j][1] ], color=c[i])
    plt.plot( [ xycorners[i][0], xycorners[j][0] ],  [ xycorners[i][1], xycorners[j][1] ], color=c[i])

d = 2
for i in range(len(Rrcorners)):
    Rrcorners[i][1] += d

xycorners = [ Rr2xy(xyi,q) for xyi in Rrcorners ]

for i in range(4):
    j = i+ 1
    if j ==4: j=0

    plt.plot( [ Rrcorners[i][0], Rrcorners[j][0] ],  [ Rrcorners[i][1], Rrcorners[j][1] ], color=c[i])
#    plt.plot( [ xycorners[i][0], xycorners[j][0] ],  [ xycorners[i][1], xycorners[j][1] ], color=c[i])
 
plt.show()







