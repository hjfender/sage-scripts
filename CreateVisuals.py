#!/usr/bin/env sage -python

import numpy as np
from sage.all import *
from functions import a5 as func

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

contour_lines = np.arange(0, 2.05, 0.025)

lab = True
x = var('x');
y = var('y');
# p = contour_plot(func.f(sqrt(100),x,y,0.0,0.0,0.0), (x,-2,2), (y,-2,2), contours=contour_lines, fill=True, cmap='jet', colorbar=True);
# t1 = text("q=10", (-1.8,1.9), rgbcolor=(1,1,1))
# t2 = text("f_A5(x,y,i,i,i)", (-1.55,1.75), rgbcolor=(1,1,1))
# t3 = text("i="+'{0:.5}'.format(0.0005), (-1.65,1.6), rgbcolor=(1,1,1))
# p = p+t1+t2+t3;
# show(p);

# for i in range(1,101):
#         p = contour_plot(func.f(sqrt(i),x,y,0,0,0), (x,-2,2), (y,-2,2), contours=np.arange(0, 5.5, 0.5), cmap='jet', colorbar=True);
#         p.save('visuals/'+str(i)+'.png');

for i in range(0,1000):
    j = float(i)/2000.0;
    if not (j == 1.0):
        p = contour_plot(func.f(sqrt(100),x,y,j,j,j), (x,-2,2), (y,-2,2), contours=contour_lines, fill=True, cmap='jet', colorbar=True);
        t1 = text("q=100", (-1.75,1.9), rgbcolor=(1,1,1))
        t2 = text("f_A5(x,y,j,j,j)", (-1.55,1.75), rgbcolor=(1,1,1))
        t3 = text("j="+'{:05.4f}'.format(j), (-1.65,1.6), rgbcolor=(1,1,1))
        p = p+t1+t2+t3;
        p.save('visuals/'+str(i)+'.png');