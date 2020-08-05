import numpy as np
import matplotlib.pyplot as plt
from sys import exit


def norm_array(x,y,l,r):
    i = 0;
    while x[i] < l: i += 1
    j = len(y)
    while x[j-1] > r : j -= 1

    return x[i:j], y[i:j]


