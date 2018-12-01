#!/usr/bin/env sage -python

import numpy as np
from sage.all import *
from functions import a4 as func

def create_func(F,f,q,x,n,m):
    input = [0 for y in F.variable_names()];
    input[0] = q;
    z = var('z');
    input[n] = z;
    input[m] = x;
    g = f(input);
    return g;

# bound = 1;

# # for q in range(1,100):
# for i in range(-25,50):
#     x = (2.0*float(i) - 50.0)/50.0;
#     if not (x == 1.0):
#         h = create_func(func.F,func.f,10,x,1,2);
#         p = complex_plot(h/h.abs(),(-bound,bound),(-bound,bound));
#         show(p);

x = var('x');
y = var('y');
p = contour_plot(func.f(sqrt(100),x,y,0,0),(x,-2,2),(y,-2,2),contours=np.arange(-10, 10, 0.1));
show(p);