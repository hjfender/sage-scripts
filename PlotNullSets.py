#!/usr/bin/env sage -python

import numpy as np
from sage.all import *
from functions import d3 as func

x1 = func.F.gen(1)
x2 = func.F.gen(2)
x3 = func.F.gen(3)
# x4 = func.F.gen(4)

g = func.f.numerator()

start = 50
l = 1
end = start + l

r = 4

for i in range(start,end):
    p = Graphics()
    p += implicit_plot3d(g(sqrt(i),x1,x2,x3),(x1,-r,r),(x2,-r,r),(x3,-r,r),plot_points=100)
    p += text3d("q="+str(i),(-3.5,3.5,3.5),color='black',fontsize='large')
    p.save('visuals/d3q'+str(i)+'.png');