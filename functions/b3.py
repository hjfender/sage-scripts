from weyl_group_multiple_dirichlet_series_v4.field_of_rational_functions_with_weyl_action import *

R = RootSystem(['B', 3]);

F = FieldOfRationalFunctionsWithWeylAction(R);

sqrtq = F.CF.gen(0);
x1 = F.CF.gen(1);
x2 = F.CF.gen(2);
x3 = F.CF.gen(3);
f1 = (sqrtq**36*x1**9*x2**15*x3**16 - sqrtq**36*x1**9*x2**14*x3**16 - sqrtq**36*x1**8*x2**15*x3**16 + sqrtq**36*x1**8*x2**14*x3**16 + sqrtq**34*x1**9*x2**13*x3**16 + sqrtq**34*x1**7*x2**15*x3**16 - sqrtq**34*x1**8*x2**13*x3**16 - sqrtq**34*x1**7*x2**14*x3**16 - sqrtq**34*x1**7*x2**13*x3**16 + sqrtq**34*x1**7*x2**12*x3**16 + sqrtq**34*x1**6*x2**13*x3**16 - sqrtq**32*x1**9*x2**12*x3**15 - sqrtq**34*x1**6*x2**12*x3**16 + 2*sqrtq**32*x1**7*x2**13*x3**16 + sqrtq**32*x1**9*x2**12*x3**14 + sqrtq**32*x1**8*x2**12*x3**15 - sqrtq**32*x1**7*x2**12*x3**16 - sqrtq**32*x1**6*x2**13*x3**16 - sqrtq**32*x1**8*x2**12*x3**14 - sqrtq**32*x1**7*x2**11*x3**16 - sqrtq**32*x1**5*x2**13*x3**16 + sqrtq**30*x1**9*x2**11*x3**15 + sqrtq**32*x1**6*x2**11*x3**16 + sqrtq**32*x1**5*x2**12*x3**16 - sqrtq**30*x1**9*x2**12*x3**13 - sqrtq**30*x1**9*x2**11*x3**14 - sqrtq**30*x1**8*x2**11*x3**15 - sqrtq**30*x1**7*x2**12*x3**15 + sqrtq**30*x1**7*x2**11*x3**16 + sqrtq**30*x1**5*x2**13*x3**16 + sqrtq**30*x1**9*x2**12*x3**12 + sqrtq**30*x1**8*x2**12*x3**13 + sqrtq**30*x1**8*x2**11*x3**14 + sqrtq**30*x1**7*x2**12*x3**14 + sqrtq**30*x1**6*x2**12*x3**15 - sqrtq**30*x1**8*x2**12*x3**12 - sqrtq**30*x1**6*x2**12*x3**14 + sqrtq**30*x1**7*x2**10*x3**15 - sqrtq**30*x1**6*x2**10*x3**16 - 2*sqrtq**30*x1**5*x2**11*x3**16 - sqrtq**30*x1**4*x2**12*x3**16 + sqrtq**28*x1**9*x2**11*x3**13 - sqrtq**30*x1**7*x2**10*x3**14 - sqrtq**30*x1**6*x2**10*x3**15 + sqrtq**28*x1**7*x2**11*x3**15 + sqrtq**30*x1**5*x2**10*x3**16 + sqrtq**30*x1**4*x2**11*x3**16 - sqrtq**28*x1**9*x2**11*x3**12 - sqrtq**28*x1**8*x2**11*x3**13 - sqrtq**28*x1**7*x2**12*x3**13 + sqrtq**30*x1**6*x2**10*x3**14 - sqrtq**28*x1**7*x2**11*x3**14 - sqrtq**28*x1**7*x2**10*x3**15 - sqrtq**28*x1**6*x2**11*x3**15 - sqrtq**28*x1**5*x2**12*x3**15 + sqrtq**28*x1**5*x2**11*x3**16 + sqrtq**28*x1**8*x2**11*x3**12 + sqrtq**28*x1**7*x2**12*x3**12 + sqrtq**28*x1**6*x2**12*x3**13 + sqrtq**28*x1**7*x2**10*x3**14 + sqrtq**28*x1**6*x2**11*x3**14 + sqrtq**28*x1**5*x2**12*x3**14 - sqrtq**28*x1**7*x2**9*x3**15 + sqrtq**28*x1**6*x2**10*x3**15 + sqrtq**28*x1**4*x2**12*x3**15 - sqrtq**28*x1**6*x2**12*x3**12 + sqrtq**28*x1**7*x2**10*x3**13 + sqrtq**28*x1**7*x2**9*x3**14 - sqrtq**28*x1**6*x2**10*x3**14 - sqrtq**28*x1**4*x2**12*x3**14 + sqrtq**28*x1**6*x2**9*x3**15 + sqrtq**28*x1**5*x2**10*x3**15 - sqrtq**28*x1**4*x2**10*x3**16 - sqrtq**28*x1**7*x2**10*x3**12 - sqrtq**28*x1**6*x2**10*x3**13 + sqrtq**26*x1**7*x2**11*x3**13 - sqrtq**28*x1**6*x2**9*x3**14 - sqrtq**28*x1**5*x2**10*x3**14 + sqrtq**26*x1**7*x2**9*x3**15 - sqrtq**28*x1**4*x2**10*x3**15 + sqrtq**26*x1**5*x2**11*x3**15 + sqrtq**26*x1**9*x2**10*x3**11 + sqrtq**28*x1**6*x2**10*x3**12 - sqrtq**26*x1**7*x2**11*x3**12 - sqrtq**26*x1**7*x2**10*x3**13 - sqrtq**26*x1**6*x2**11*x3**13 - sqrtq**26*x1**5*x2**12*x3**13 - sqrtq**26*x1**7*x2**9*x3**14 + sqrtq**28*x1**4*x2**10*x3**14 - sqrtq**26*x1**5*x2**11*x3**14 - sqrtq**26*x1**5*x2**10*x3**15 - sqrtq**26*x1**4*x2**11*x3**15 - sqrtq**26*x1**9*x2**10*x3**10 - sqrtq**26*x1**8*x2**10*x3**11 + sqrtq**26*x1**7*x2**10*x3**12 + sqrtq**26*x1**6*x2**11*x3**12 + sqrtq**26*x1**5*x2**12*x3**12 - sqrtq**26*x1**7*x2**9*x3**13 + sqrtq**26*x1**6*x2**10*x3**13 + sqrtq**26*x1**4*x2**12*x3**13 + sqrtq**26*x1**5*x2**10*x3**14 + sqrtq**26*x1**4*x2**11*x3**14 - sqrtq**26*x1**6*x2**8*x3**15 - 2*sqrtq**26*x1**5*x2**9*x3**15 + sqrtq**26*x1**4*x2**10*x3**15 - sqrtq**26*x1**3*x2**11*x3**15 + sqrtq**26*x1**8*x2**10*x3**10 + sqrtq**26*x1**7*x2**9*x3**12 - sqrtq**26*x1**6*x2**10*x3**12 - sqrtq**26*x1**4*x2**12*x3**12 + sqrtq**26*x1**6*x2**9*x3**13 + sqrtq**26*x1**5*x2**10*x3**13 + sqrtq**26*x1**6*x2**8*x3**14 + 2*sqrtq**26*x1**5*x2**9*x3**14 - sqrtq**26*x1**4*x2**10*x3**14 + sqrtq**26*x1**3*x2**11*x3**14 + sqrtq**26*x1**5*x2**8*x3**15 + sqrtq**26*x1**4*x2**9*x3**15 + sqrtq**26*x1**3*x2**10*x3**15 - sqrtq**26*x1**6*x2**9*x3**12 - sqrtq**26*x1**5*x2**10*x3**12 + sqrtq**24*x1**7*x2**9*x3**13 - sqrtq**26*x1**4*x2**10*x3**13 + sqrtq**24*x1**5*x2**11*x3**13 - sqrtq**26*x1**5*x2**8*x3**14 - sqrtq**26*x1**4*x2**9*x3**14 - sqrtq**26*x1**3*x2**10*x3**14 + sqrtq**24*x1**5*x2**9*x3**15 + sqrtq**24*x1**3*x2**11*x3**15 + sqrtq**24*x1**9*x2**10*x3**9 + sqrtq**24*x1**7*x2**10*x3**11 - sqrtq**24*x1**7*x2**9*x3**12 + sqrtq**26*x1**4*x2**10*x3**12 - sqrtq**24*x1**5*x2**11*x3**12 - sqrtq**24*x1**5*x2**10*x3**13 - sqrtq**24*x1**4*x2**11*x3**13 - sqrtq**24*x1**5*x2**9*x3**14 - sqrtq**24*x1**3*x2**11*x3**14 - sqrtq**24*x1**8*x2**10*x3**9 - sqrtq**24*x1**7*x2**10*x3**10 - sqrtq**24*x1**6*x2**10*x3**11 + sqrtq**24*x1**5*x2**10*x3**12 + sqrtq**24*x1**4*x2**11*x3**12 - sqrtq**24*x1**6*x2**8*x3**13 - 2*sqrtq**24*x1**5*x2**9*x3**13 + sqrtq**24*x1**4*x2**10*x3**13 - sqrtq**24*x1**3*x2**11*x3**13 - sqrtq**24*x1**4*x2**8*x3**15 - 2*sqrtq**24*x1**3*x2**9*x3**15 - sqrtq**24*x1**2*x2**10*x3**15 + sqrtq**24*x1**6*x2**10*x3**10 - sqrtq**24*x1**7*x2**8*x3**11 + sqrtq**24*x1**6*x2**8*x3**12 + 2*sqrtq**24*x1**5*x2**9*x3**12 - sqrtq**24*x1**4*x2**10*x3**12 + sqrtq**24*x1**3*x2**11*x3**12 + sqrtq**24*x1**5*x2**8*x3**13 + sqrtq**24*x1**4*x2**9*x3**13 + sqrtq**24*x1**3*x2**10*x3**13 + sqrtq**24*x1**4*x2**8*x3**14 + 2*sqrtq**24*x1**3*x2**9*x3**14 + sqrtq**24*x1**2*x2**10*x3**14 + sqrtq**24*x1**3*x2**8*x3**15 + sqrtq**24*x1**2*x2**9*x3**15 + sqrtq**24*x1**7*x2**8*x3**10 + sqrtq**24*x1**6*x2**8*x3**11 - sqrtq**24*x1**5*x2**8*x3**12 - sqrtq**24*x1**4*x2**9*x3**12 - sqrtq**24*x1**3*x2**10*x3**12 + sqrtq**22*x1**5*x2**9*x3**13 + sqrtq**22*x1**3*x2**11*x3**13 - sqrtq**24*x1**3*x2**8*x3**14 - sqrtq**24*x1**2*x2**9*x3**14 + sqrtq**22*x1**3*x2**9*x3**15 + sqrtq**22*x1**7*x2**10*x3**9 - sqrtq**24*x1**6*x2**8*x3**10 + sqrtq**22*x1**7*x2**8*x3**11 + sqrtq**22*x1**6*x2**9*x3**11 + sqrtq**22*x1**5*x2**10*x3**11 - sqrtq**22*x1**5*x2**9*x3**12 - sqrtq**22*x1**3*x2**11*x3**12 - sqrtq**22*x1**3*x2**9*x3**14 - sqrtq**22*x1**6*x2**10*x3**9 - sqrtq**22*x1**7*x2**8*x3**10 - sqrtq**22*x1**6*x2**9*x3**10 - sqrtq**22*x1**5*x2**10*x3**10 - 2*sqrtq**22*x1**6*x2**8*x3**11 - sqrtq**22*x1**5*x2**9*x3**11 - sqrtq**22*x1**4*x2**10*x3**11 - sqrtq**22*x1**4*x2**8*x3**13 - 2*sqrtq**22*x1**3*x2**9*x3**13 - sqrtq**22*x1**2*x2**10*x3**13 - sqrtq**22*x1**2*x2**8*x3**15 - sqrtq**22*x1**7*x2**8*x3**9 + 2*sqrtq**22*x1**6*x2**8*x3**10 + sqrtq**22*x1**5*x2**9*x3**10 + sqrtq**22*x1**4*x2**10*x3**10 - sqrtq**22*x1**5*x2**8*x3**11 + sqrtq**22*x1**4*x2**8*x3**12 + 2*sqrtq**22*x1**3*x2**9*x3**12 + sqrtq**22*x1**2*x2**10*x3**12 + sqrtq**22*x1**3*x2**8*x3**13 + sqrtq**22*x1**2*x2**9*x3**13 + sqrtq**22*x1**2*x2**8*x3**14 + sqrtq**22*x1**6*x2**8*x3**9 + sqrtq**22*x1**5*x2**8*x3**10 + 2*sqrtq**22*x1**4*x2**8*x3**11 - sqrtq**22*x1**3*x2**8*x3**12 - sqrtq**22*x1**2*x2**9*x3**12 + sqrtq**20*x1**3*x2**9*x3**13 + sqrtq**20*x1**7*x2**8*x3**9 + sqrtq**20*x1**6*x2**9*x3**9 + sqrtq**20*x1**5*x2**10*x3**9 - 2*sqrtq**22*x1**4*x2**8*x3**10 + 2*sqrtq**20*x1**5*x2**8*x3**11 + sqrtq**20*x1**4*x2**9*x3**11 + sqrtq**20*x1**3*x2**10*x3**11 - sqrtq**20*x1**3*x2**9*x3**12 - sqrtq**20*x1**6*x2**9*x3**8 - 2*sqrtq**20*x1**6*x2**8*x3**9 - sqrtq**20*x1**5*x2**9*x3**9 - sqrtq**20*x1**4*x2**10*x3**9 - 2*sqrtq**20*x1**5*x2**8*x3**10 - sqrtq**20*x1**4*x2**9*x3**10 - sqrtq**20*x1**3*x2**10*x3**10 + sqrtq**20*x1**5*x2**7*x3**11 - 3*sqrtq**20*x1**4*x2**8*x3**11 - sqrtq**20*x1**2*x2**10*x3**11 - sqrtq**20*x1**2*x2**8*x3**13 + sqrtq**20*x1**6*x2**8*x3**8 + sqrtq**20*x1**5*x2**9*x3**8 - sqrtq**20*x1**5*x2**8*x3**9 - sqrtq**20*x1**5*x2**7*x3**10 + 3*sqrtq**20*x1**4*x2**8*x3**10 + sqrtq**20*x1**2*x2**10*x3**10 - sqrtq**20*x1**5*x2**6*x3**11 - sqrtq**20*x1**4*x2**7*x3**11 - sqrtq**20*x1**3*x2**8*x3**11 + sqrtq**20*x1**2*x2**8*x3**12 + 2*sqrtq**20*x1**4*x2**8*x3**9 + sqrtq**20*x1**5*x2**6*x3**10 + sqrtq**20*x1**4*x2**7*x3**10 + sqrtq**20*x1**3*x2**8*x3**10 + sqrtq**20*x1**4*x2**6*x3**11 - sqrtq**18*x1**5*x2**7*x3**11 + sqrtq**20*x1**2*x2**8*x3**11 - sqrtq**20*x1**4*x2**8*x3**8 + 2*sqrtq**18*x1**5*x2**8*x3**9 + sqrtq**18*x1**4*x2**9*x3**9 + sqrtq**18*x1**3*x2**10*x3**9 - sqrtq**20*x1**4*x2**6*x3**10 + sqrtq**18*x1**5*x2**7*x3**10 - sqrtq**20*x1**2*x2**8*x3**10 + sqrtq**18*x1**5*x2**6*x3**11 + sqrtq**18*x1**4*x2**7*x3**11 + sqrtq**18*x1**3*x2**8*x3**11 - sqrtq**18*x1**5*x2**8*x3**8 - sqrtq**18*x1**4*x2**9*x3**8 + sqrtq**18*x1**5*x2**7*x3**9 - 3*sqrtq**18*x1**4*x2**8*x3**9 - sqrtq**18*x1**2*x2**10*x3**9 - sqrtq**18*x1**5*x2**6*x3**10 - sqrtq**18*x1**4*x2**7*x3**10 - sqrtq**18*x1**3*x2**8*x3**10 - sqrtq**18*x1**4*x2**6*x3**11 + 2*sqrtq**18*x1**3*x2**7*x3**11 - sqrtq**18*x1**2*x2**8*x3**11 + sqrtq**18*x1*x2**9*x3**11 - sqrtq**18*x1**5*x2**7*x3**8 + 2*sqrtq**18*x1**4*x2**8*x3**8 - sqrtq**18*x1**5*x2**6*x3**9 - sqrtq**18*x1**4*x2**7*x3**9 - sqrtq**18*x1**3*x2**8*x3**9 + sqrtq**18*x1**4*x2**6*x3**10 - 2*sqrtq**18*x1**3*x2**7*x3**10 + sqrtq**18*x1**2*x2**8*x3**10 - sqrtq**18*x1*x2**9*x3**10 - 2*sqrtq**18*x1**3*x2**6*x3**11 - sqrtq**18*x1**2*x2**7*x3**11 - sqrtq**18*x1*x2**8*x3**11 + sqrtq**18*x1**4*x2**7*x3**8 + sqrtq**18*x1**4*x2**6*x3**9 - sqrtq**16*x1**5*x2**7*x3**9 + sqrtq**18*x1**2*x2**8*x3**9 + 2*sqrtq**18*x1**3*x2**6*x3**10 + sqrtq**18*x1**2*x2**7*x3**10 + sqrtq**18*x1*x2**8*x3**10 + sqrtq**18*x1**2*x2**6*x3**11 - sqrtq**16*x1**3*x2**7*x3**11 - sqrtq**16*x1**6*x2**7*x3**7 + sqrtq**16*x1**5*x2**7*x3**8 + sqrtq**16*x1**5*x2**6*x3**9 + sqrtq**16*x1**4*x2**7*x3**9 + sqrtq**16*x1**3*x2**8*x3**9 - sqrtq**18*x1**2*x2**6*x3**10 + sqrtq**16*x1**3*x2**7*x3**10 + sqrtq**16*x1**3*x2**6*x3**11 + sqrtq**16*x1**6*x2**7*x3**6 + sqrtq**16*x1**6*x2**6*x3**7 + sqrtq**16*x1**5*x2**7*x3**7 - sqrtq**16*x1**4*x2**7*x3**8 - sqrtq**16*x1**4*x2**6*x3**9 + 2*sqrtq**16*x1**3*x2**7*x3**9 - sqrtq**16*x1**2*x2**8*x3**9 + sqrtq**16*x1*x2**9*x3**9 - sqrtq**16*x1**3*x2**6*x3**10 + 2*sqrtq**16*x1**3*x2**5*x3**11 + sqrtq**16*x1*x2**7*x3**11 - sqrtq**16*x1**6*x2**6*x3**6 - sqrtq**16*x1**5*x2**7*x3**6 - 2*sqrtq**16*x1**3*x2**7*x3**8 - 2*sqrtq**16*x1**3*x2**6*x3**9 - sqrtq**16*x1**2*x2**7*x3**9 - sqrtq**16*x1*x2**8*x3**9 - 2*sqrtq**16*x1**3*x2**5*x3**10 - sqrtq**16*x1*x2**7*x3**10 - sqrtq**16*x1**3*x2**4*x3**11 - sqrtq**16*x1**2*x2**5*x3**11 - sqrtq**16*x1*x2**6*x3**11 - sqrtq**16*x1**4*x2**6*x3**7 + sqrtq**16*x1**3*x2**6*x3**8 + sqrtq**16*x1**2*x2**7*x3**8 + sqrtq**16*x1**2*x2**6*x3**9 - sqrtq**14*x1**3*x2**7*x3**9 + sqrtq**16*x1**3*x2**4*x3**10 + sqrtq**16*x1**2*x2**5*x3**10 + sqrtq**16*x1*x2**6*x3**10 - sqrtq**14*x1**3*x2**5*x3**11 - sqrtq**14*x1**6*x2**7*x3**5 + sqrtq**16*x1**4*x2**6*x3**6 - sqrtq**14*x1**5*x2**6*x3**7 - sqrtq**14*x1**4*x2**7*x3**7 + sqrtq**14*x1**3*x2**7*x3**8 + sqrtq**14*x1**3*x2**6*x3**9 + sqrtq**14*x1**3*x2**5*x3**10 + sqrtq**14*x1**6*x2**6*x3**5 + sqrtq**14*x1**5*x2**7*x3**5 + sqrtq**14*x1**5*x2**6*x3**6 + sqrtq**14*x1**4*x2**7*x3**6 + 2*sqrtq**14*x1**4*x2**6*x3**7 + sqrtq**14*x1**3*x2**7*x3**7 + 2*sqrtq**14*x1**3*x2**5*x3**9 + sqrtq**14*x1*x2**7*x3**9 + sqrtq**14*x1**2*x2**4*x3**11 + sqrtq**14*x1*x2**5*x3**11 - 2*sqrtq**14*x1**4*x2**6*x3**6 - sqrtq**14*x1**3*x2**7*x3**6 + sqrtq**14*x1**4*x2**5*x3**7 - 2*sqrtq**14*x1**3*x2**5*x3**8 - sqrtq**14*x1**2*x2**6*x3**8 - sqrtq**14*x1**3*x2**4*x3**9 - sqrtq**14*x1**2*x2**5*x3**9 - sqrtq**14*x1*x2**6*x3**9 - sqrtq**14*x1**2*x2**4*x3**10 - sqrtq**14*x1*x2**5*x3**10 - sqrtq**14*x1*x2**4*x3**11 - sqrtq**14*x1**4*x2**6*x3**5 - sqrtq**14*x1**4*x2**5*x3**6 - sqrtq**14*x1**4*x2**4*x3**7 - sqrtq**14*x1**2*x2**6*x3**7 + sqrtq**14*x1**2*x2**5*x3**8 - sqrtq**12*x1**3*x2**5*x3**9 + sqrtq**14*x1*x2**4*x3**10 - sqrtq**12*x1**5*x2**6*x3**5 - sqrtq**12*x1**4*x2**7*x3**5 + sqrtq**14*x1**4*x2**4*x3**6 + sqrtq**14*x1**2*x2**6*x3**6 - sqrtq**12*x1**4*x2**5*x3**7 - sqrtq**12*x1**3*x2**6*x3**7 - sqrtq**12*x1**2*x2**7*x3**7 + sqrtq**14*x1**2*x2**4*x3**8 + sqrtq**12*x1**3*x2**5*x3**8 + 2*sqrtq**12*x1**4*x2**6*x3**5 + sqrtq**12*x1**3*x2**7*x3**5 + sqrtq**12*x1**4*x2**5*x3**6 + sqrtq**12*x1**3*x2**6*x3**6 + sqrtq**12*x1**2*x2**7*x3**6 + sqrtq**12*x1**4*x2**4*x3**7 + 2*sqrtq**12*x1**2*x2**6*x3**7 + sqrtq**12*x1**3*x2**4*x3**8 + sqrtq**12*x1**2*x2**4*x3**9 + sqrtq**12*x1*x2**5*x3**9 + sqrtq**12*x1**4*x2**5*x3**5 - sqrtq**12*x1**4*x2**4*x3**6 - 2*sqrtq**12*x1**2*x2**6*x3**6 + sqrtq**12*x1**3*x2**4*x3**7 + sqrtq**12*x1**2*x2**5*x3**7 - 2*sqrtq**12*x1**2*x2**4*x3**8 - sqrtq**12*x1*x2**4*x3**9 - sqrtq**12*x1**4*x2**4*x3**5 - sqrtq**12*x1**2*x2**6*x3**5 - sqrtq**12*x1**3*x2**4*x3**6 - sqrtq**12*x1**2*x2**5*x3**6 - 2*sqrtq**12*x1**2*x2**4*x3**7 - sqrtq**10*x1**4*x2**5*x3**5 - sqrtq**10*x1**3*x2**6*x3**5 - sqrtq**10*x1**2*x2**7*x3**5 + 2*sqrtq**12*x1**2*x2**4*x3**6 - sqrtq**10*x1**3*x2**4*x3**7 - sqrtq**10*x1**2*x2**5*x3**7 + sqrtq**10*x1**4*x2**4*x3**5 + 2*sqrtq**10*x1**2*x2**6*x3**5 + sqrtq**10*x1**3*x2**4*x3**6 + sqrtq**10*x1**2*x2**5*x3**6 - sqrtq**10*x1**3*x2**3*x3**7 + 2*sqrtq**10*x1**2*x2**4*x3**7 - sqrtq**10*x1*x2**5*x3**7 + sqrtq**10*x1**3*x2**4*x3**5 + sqrtq**10*x1**2*x2**5*x3**5 + sqrtq**10*x1**3*x2**3*x3**6 - 2*sqrtq**10*x1**2*x2**4*x3**6 + sqrtq**10*x1*x2**5*x3**6 + sqrtq**10*x1**2*x2**3*x3**7 + sqrtq**10*x1*x2**4*x3**7 + sqrtq**10*x2**5*x3**7 + sqrtq**10*x1*x2**3*x3**8 - 2*sqrtq**10*x1**2*x2**4*x3**5 - sqrtq**10*x1**2*x2**3*x3**6 - sqrtq**10*x1*x2**4*x3**6 - sqrtq**10*x2**5*x3**6 - sqrtq**10*x2**4*x3**7 - sqrtq**10*x1*x2**2*x3**8 - sqrtq**8*x1**3*x2**4*x3**5 - sqrtq**8*x1**2*x2**5*x3**5 + sqrtq**10*x2**4*x3**6 - sqrtq**8*x1**3*x2**3*x3**5 + 2*sqrtq**8*x1**2*x2**4*x3**5 - sqrtq**8*x1*x2**5*x3**5 - 2*sqrtq**8*x1*x2**3*x3**7 + sqrtq**8*x1**3*x2**3*x3**4 + sqrtq**8*x1**2*x2**3*x3**5 + sqrtq**8*x1*x2**4*x3**5 + sqrtq**8*x2**5*x3**5 + 2*sqrtq**8*x1*x2**3*x3**6 + sqrtq**8*x1*x2**2*x3**7 + sqrtq**8*x2**3*x3**7 - sqrtq**8*x2**4*x3**5 - sqrtq**8*x1*x2**2*x3**6 - sqrtq**8*x2**3*x3**6 - sqrtq**8*x1**2*x2**2*x3**4 - sqrtq**6*x1**3*x2**3*x3**3 - sqrtq**6*x1**2*x2**3*x3**4 - 2*sqrtq**6*x1*x2**3*x3**5 + sqrtq**6*x1**3*x2**3*x3**2 + sqrtq**6*x1**2*x2**2*x3**4 + sqrtq**6*x1*x2**3*x3**4 + sqrtq**6*x1*x2**2*x3**5 + sqrtq**6*x2**3*x3**5 + sqrtq**6*x1**2*x2**2*x3**3 - sqrtq**6*x1**2*x2**2*x3**2 + sqrtq**4*x1**2*x2**3*x3**3 - sqrtq**6*x2**2*x3**4 - sqrtq**4*x1**3*x2**3*x3 - sqrtq**4*x1**2*x2**3*x3**2 - sqrtq**4*x1**2*x2**2*x3**3 - sqrtq**4*x1*x2**3*x3**3 + sqrtq**4*x1**2*x2**2*x3**2 + sqrtq**4*x1*x2**3*x3**2 - sqrtq**4*x1*x2*x3**4 + sqrtq**4*x1**2*x2**2*x3 + sqrtq**4*x2**2*x3**3 + sqrtq**4*x2*x3**4 + sqrtq**2*x1**2*x2**3*x3 - sqrtq**4*x2**2*x3**2 - sqrtq**2*x1**2*x2**2*x3 - sqrtq**2*x1*x2**3*x3 + sqrtq**2*x1*x2*x3**3 - sqrtq**2*x1*x2*x3**2 - sqrtq**2*x2*x3**3 + sqrtq**2*x2**2*x3 + sqrtq**2*x2*x3**2 + x1*x2*x3 - x1*x2 - x2*x3 + 1)/(-sqrtq**36*x1**9*x2**15*x3**17 + sqrtq**36*x1**9*x2**15*x3**16 + sqrtq**36*x1**9*x2**14*x3**17 + sqrtq**36*x1**8*x2**15*x3**17 - sqrtq**36*x1**9*x2**14*x3**16 - sqrtq**36*x1**8*x2**15*x3**16 - sqrtq**36*x1**8*x2**14*x3**17 + sqrtq**36*x1**8*x2**14*x3**16 + sqrtq**34*x1**7*x2**13*x3**17 - sqrtq**34*x1**7*x2**13*x3**16 - sqrtq**34*x1**7*x2**12*x3**17 - sqrtq**34*x1**6*x2**13*x3**17 + sqrtq**32*x1**9*x2**13*x3**15 + sqrtq**34*x1**7*x2**12*x3**16 + sqrtq**34*x1**6*x2**13*x3**16 + sqrtq**34*x1**6*x2**12*x3**17 - sqrtq**32*x1**9*x2**13*x3**14 - sqrtq**32*x1**9*x2**12*x3**15 - sqrtq**32*x1**8*x2**13*x3**15 - sqrtq**34*x1**6*x2**12*x3**16 + sqrtq**32*x1**9*x2**12*x3**14 + sqrtq**32*x1**8*x2**13*x3**14 + sqrtq**32*x1**8*x2**12*x3**15 - sqrtq**32*x1**8*x2**12*x3**14 + sqrtq**30*x1**9*x2**13*x3**13 + sqrtq**30*x1**7*x2**13*x3**15 - sqrtq**30*x1**9*x2**13*x3**12 - sqrtq**30*x1**9*x2**12*x3**13 - sqrtq**30*x1**8*x2**13*x3**13 - sqrtq**30*x1**7*x2**13*x3**14 - sqrtq**30*x1**7*x2**12*x3**15 - sqrtq**30*x1**6*x2**13*x3**15 + sqrtq**30*x1**9*x2**12*x3**12 + sqrtq**30*x1**8*x2**13*x3**12 + sqrtq**30*x1**8*x2**12*x3**13 + sqrtq**30*x1**7*x2**12*x3**14 + sqrtq**30*x1**6*x2**13*x3**14 - sqrtq**30*x1**7*x2**11*x3**15 + sqrtq**30*x1**6*x2**12*x3**15 - sqrtq**30*x1**8*x2**12*x3**12 + sqrtq**30*x1**7*x2**11*x3**14 - sqrtq**30*x1**6*x2**12*x3**14 + sqrtq**30*x1**7*x2**10*x3**15 + sqrtq**30*x1**6*x2**11*x3**15 + sqrtq**28*x1**7*x2**13*x3**13 - sqrtq**30*x1**7*x2**10*x3**14 - sqrtq**30*x1**6*x2**11*x3**14 - sqrtq**30*x1**6*x2**10*x3**15 - sqrtq**28*x1**7*x2**13*x3**12 - sqrtq**28*x1**7*x2**12*x3**13 - sqrtq**28*x1**6*x2**13*x3**13 + sqrtq**30*x1**6*x2**10*x3**14 + sqrtq**28*x1**7*x2**12*x3**12 + sqrtq**28*x1**6*x2**13*x3**12 - sqrtq**28*x1**7*x2**11*x3**13 + sqrtq**28*x1**6*x2**12*x3**13 - sqrtq**28*x1**5*x2**11*x3**15 + sqrtq**28*x1**7*x2**11*x3**12 - sqrtq**28*x1**6*x2**12*x3**12 + sqrtq**28*x1**7*x2**10*x3**13 + sqrtq**28*x1**6*x2**11*x3**13 + sqrtq**28*x1**5*x2**11*x3**14 + sqrtq**28*x1**5*x2**10*x3**15 + sqrtq**28*x1**4*x2**11*x3**15 - sqrtq**26*x1**9*x2**11*x3**11 - sqrtq**28*x1**7*x2**10*x3**12 - sqrtq**28*x1**6*x2**11*x3**12 - sqrtq**28*x1**6*x2**10*x3**13 - sqrtq**28*x1**5*x2**10*x3**14 - sqrtq**28*x1**4*x2**11*x3**14 - sqrtq**28*x1**4*x2**10*x3**15 + sqrtq**26*x1**9*x2**11*x3**10 + sqrtq**26*x1**9*x2**10*x3**11 + sqrtq**26*x1**8*x2**11*x3**11 + sqrtq**28*x1**6*x2**10*x3**12 + sqrtq**28*x1**4*x2**10*x3**14 - sqrtq**26*x1**9*x2**10*x3**10 - sqrtq**26*x1**8*x2**11*x3**10 - sqrtq**26*x1**8*x2**10*x3**11 - sqrtq**26*x1**5*x2**11*x3**13 + sqrtq**26*x1**8*x2**10*x3**10 + sqrtq**26*x1**5*x2**11*x3**12 + sqrtq**26*x1**5*x2**10*x3**13 + sqrtq**26*x1**4*x2**11*x3**13 - 2*sqrtq**24*x1**7*x2**11*x3**11 - sqrtq**26*x1**5*x2**10*x3**12 - sqrtq**26*x1**4*x2**11*x3**12 - sqrtq**26*x1**4*x2**10*x3**13 + 2*sqrtq**24*x1**7*x2**11*x3**10 + 2*sqrtq**24*x1**7*x2**10*x3**11 + 2*sqrtq**24*x1**6*x2**11*x3**11 + sqrtq**26*x1**4*x2**10*x3**12 - 2*sqrtq**24*x1**7*x2**10*x3**10 - 2*sqrtq**24*x1**6*x2**11*x3**10 + sqrtq**24*x1**7*x2**9*x3**11 - 2*sqrtq**24*x1**6*x2**10*x3**11 - sqrtq**24*x1**7*x2**9*x3**10 + 2*sqrtq**24*x1**6*x2**10*x3**10 - sqrtq**24*x1**7*x2**8*x3**11 - sqrtq**24*x1**6*x2**9*x3**11 - sqrtq**22*x1**7*x2**11*x3**9 + sqrtq**24*x1**7*x2**8*x3**10 + sqrtq**24*x1**6*x2**9*x3**10 + sqrtq**24*x1**6*x2**8*x3**11 - sqrtq**22*x1**7*x2**9*x3**11 - sqrtq**22*x1**5*x2**11*x3**11 + sqrtq**22*x1**7*x2**11*x3**8 + sqrtq**22*x1**7*x2**10*x3**9 + sqrtq**22*x1**6*x2**11*x3**9 - sqrtq**24*x1**6*x2**8*x3**10 + sqrtq**22*x1**7*x2**9*x3**10 + sqrtq**22*x1**5*x2**11*x3**10 + sqrtq**22*x1**7*x2**8*x3**11 + sqrtq**22*x1**6*x2**9*x3**11 + sqrtq**22*x1**5*x2**10*x3**11 + sqrtq**22*x1**4*x2**11*x3**11 - sqrtq**22*x1**7*x2**10*x3**8 - sqrtq**22*x1**6*x2**11*x3**8 - sqrtq**22*x1**6*x2**10*x3**9 - sqrtq**22*x1**7*x2**8*x3**10 - sqrtq**22*x1**6*x2**9*x3**10 - sqrtq**22*x1**5*x2**10*x3**10 - sqrtq**22*x1**4*x2**11*x3**10 - sqrtq**22*x1**6*x2**8*x3**11 + 2*sqrtq**22*x1**5*x2**9*x3**11 - sqrtq**22*x1**4*x2**10*x3**11 + sqrtq**22*x1**6*x2**10*x3**8 + sqrtq**22*x1**6*x2**8*x3**10 - 2*sqrtq**22*x1**5*x2**9*x3**10 + sqrtq**22*x1**4*x2**10*x3**10 - 2*sqrtq**22*x1**5*x2**8*x3**11 - 2*sqrtq**22*x1**4*x2**9*x3**11 + 2*sqrtq**22*x1**5*x2**8*x3**10 + 2*sqrtq**22*x1**4*x2**9*x3**10 + 2*sqrtq**22*x1**4*x2**8*x3**11 - sqrtq**20*x1**5*x2**9*x3**11 - 2*sqrtq**22*x1**4*x2**8*x3**10 + sqrtq**20*x1**5*x2**9*x3**10 + sqrtq**20*x1**5*x2**8*x3**11 + sqrtq**20*x1**4*x2**9*x3**11 + sqrtq**20*x1**5*x2**9*x3**9 - sqrtq**20*x1**5*x2**8*x3**10 - sqrtq**20*x1**4*x2**9*x3**10 + sqrtq**20*x1**5*x2**7*x3**11 - sqrtq**20*x1**4*x2**8*x3**11 + sqrtq**20*x1**3*x2**9*x3**11 - sqrtq**20*x1**5*x2**9*x3**8 - sqrtq**20*x1**5*x2**8*x3**9 - sqrtq**20*x1**4*x2**9*x3**9 - sqrtq**20*x1**5*x2**7*x3**10 + sqrtq**20*x1**4*x2**8*x3**10 - sqrtq**20*x1**3*x2**9*x3**10 - sqrtq**20*x1**5*x2**6*x3**11 - sqrtq**20*x1**4*x2**7*x3**11 - sqrtq**20*x1**3*x2**8*x3**11 - sqrtq**20*x1**2*x2**9*x3**11 + sqrtq**18*x1**7*x2**9*x3**7 + sqrtq**20*x1**5*x2**8*x3**8 + sqrtq**20*x1**4*x2**9*x3**8 + sqrtq**20*x1**4*x2**8*x3**9 + sqrtq**20*x1**5*x2**6*x3**10 + sqrtq**20*x1**4*x2**7*x3**10 + sqrtq**20*x1**3*x2**8*x3**10 + sqrtq**20*x1**2*x2**9*x3**10 + sqrtq**20*x1**4*x2**6*x3**11 + sqrtq**20*x1**2*x2**8*x3**11 - sqrtq**18*x1**7*x2**9*x3**6 - sqrtq**18*x1**7*x2**8*x3**7 - sqrtq**18*x1**6*x2**9*x3**7 - sqrtq**20*x1**4*x2**8*x3**8 - sqrtq**20*x1**4*x2**6*x3**10 - sqrtq**20*x1**2*x2**8*x3**10 + sqrtq**18*x1**7*x2**8*x3**6 + sqrtq**18*x1**6*x2**9*x3**6 + sqrtq**18*x1**6*x2**8*x3**7 + sqrtq**18*x1**3*x2**7*x3**11 - sqrtq**18*x1**6*x2**8*x3**6 - sqrtq**18*x1**3*x2**7*x3**10 - sqrtq**18*x1**3*x2**6*x3**11 - sqrtq**18*x1**2*x2**7*x3**11 + sqrtq**16*x1**7*x2**7*x3**7 + sqrtq**16*x1**5*x2**9*x3**7 + sqrtq**16*x1**5*x2**7*x3**9 + sqrtq**18*x1**3*x2**6*x3**10 + sqrtq**18*x1**2*x2**7*x3**10 + sqrtq**18*x1**2*x2**6*x3**11 - sqrtq**16*x1**7*x2**7*x3**6 - sqrtq**16*x1**5*x2**9*x3**6 - sqrtq**16*x1**7*x2**6*x3**7 - sqrtq**16*x1**6*x2**7*x3**7 - sqrtq**16*x1**5*x2**8*x3**7 - sqrtq**16*x1**4*x2**9*x3**7 - sqrtq**16*x1**5*x2**7*x3**8 - sqrtq**16*x1**5*x2**6*x3**9 - sqrtq**16*x1**4*x2**7*x3**9 - sqrtq**18*x1**2*x2**6*x3**10 + sqrtq**16*x1**7*x2**6*x3**6 + sqrtq**16*x1**6*x2**7*x3**6 + sqrtq**16*x1**5*x2**8*x3**6 + sqrtq**16*x1**4*x2**9*x3**6 + sqrtq**16*x1**6*x2**6*x3**7 - sqrtq**16*x1**5*x2**7*x3**7 + sqrtq**16*x1**4*x2**8*x3**7 + sqrtq**16*x1**5*x2**6*x3**8 + sqrtq**16*x1**4*x2**7*x3**8 + sqrtq**16*x1**4*x2**6*x3**9 - sqrtq**16*x1**6*x2**6*x3**6 + sqrtq**16*x1**5*x2**7*x3**6 - sqrtq**16*x1**4*x2**8*x3**6 + sqrtq**16*x1**5*x2**6*x3**7 + sqrtq**16*x1**4*x2**7*x3**7 - sqrtq**16*x1**4*x2**6*x3**8 - sqrtq**16*x1**5*x2**6*x3**6 - sqrtq**16*x1**4*x2**7*x3**6 - sqrtq**16*x1**4*x2**6*x3**7 + 2*sqrtq**14*x1**5*x2**7*x3**7 + sqrtq**16*x1**4*x2**6*x3**6 - 2*sqrtq**14*x1**5*x2**7*x3**6 - 2*sqrtq**14*x1**5*x2**6*x3**7 - 2*sqrtq**14*x1**4*x2**7*x3**7 + 2*sqrtq**14*x1**5*x2**6*x3**6 + 2*sqrtq**14*x1**4*x2**7*x3**6 - sqrtq**14*x1**5*x2**5*x3**7 + 2*sqrtq**14*x1**4*x2**6*x3**7 - sqrtq**14*x1**3*x2**7*x3**7 - sqrtq**14*x1**3*x2**5*x3**9 + sqrtq**14*x1**5*x2**5*x3**6 - 2*sqrtq**14*x1**4*x2**6*x3**6 + sqrtq**14*x1**3*x2**7*x3**6 + sqrtq**14*x1**5*x2**4*x3**7 + sqrtq**14*x1**4*x2**5*x3**7 + sqrtq**14*x1**3*x2**6*x3**7 + sqrtq**14*x1**2*x2**7*x3**7 + sqrtq**14*x1**3*x2**5*x3**8 + sqrtq**14*x1**3*x2**4*x3**9 + sqrtq**14*x1**2*x2**5*x3**9 - sqrtq**14*x1**5*x2**4*x3**6 - sqrtq**14*x1**4*x2**5*x3**6 - sqrtq**14*x1**3*x2**6*x3**6 - sqrtq**14*x1**2*x2**7*x3**6 - sqrtq**14*x1**4*x2**4*x3**7 - sqrtq**14*x1**2*x2**6*x3**7 + sqrtq**12*x1**3*x2**7*x3**7 - sqrtq**14*x1**3*x2**4*x3**8 - sqrtq**14*x1**2*x2**5*x3**8 - sqrtq**14*x1**2*x2**4*x3**9 + sqrtq**14*x1**4*x2**4*x3**6 + sqrtq**14*x1**2*x2**6*x3**6 - sqrtq**12*x1**3*x2**7*x3**6 - sqrtq**12*x1**3*x2**6*x3**7 - sqrtq**12*x1**2*x2**7*x3**7 + sqrtq**14*x1**2*x2**4*x3**8 + sqrtq**12*x1**3*x2**6*x3**6 + sqrtq**12*x1**2*x2**7*x3**6 - 2*sqrtq**12*x1**3*x2**5*x3**7 + sqrtq**12*x1**2*x2**6*x3**7 + 2*sqrtq**12*x1**3*x2**5*x3**6 - sqrtq**12*x1**2*x2**6*x3**6 + 2*sqrtq**12*x1**3*x2**4*x3**7 + 2*sqrtq**12*x1**2*x2**5*x3**7 - sqrtq**10*x1**5*x2**5*x3**5 - 2*sqrtq**12*x1**3*x2**4*x3**6 - 2*sqrtq**12*x1**2*x2**5*x3**6 - 2*sqrtq**12*x1**2*x2**4*x3**7 + sqrtq**10*x1**5*x2**5*x3**4 + sqrtq**10*x1**5*x2**4*x3**5 + sqrtq**10*x1**4*x2**5*x3**5 + 2*sqrtq**12*x1**2*x2**4*x3**6 - sqrtq**10*x1**5*x2**4*x3**4 - sqrtq**10*x1**4*x2**5*x3**4 - sqrtq**10*x1**4*x2**4*x3**5 - sqrtq**10*x1*x2**5*x3**7 + sqrtq**10*x1**4*x2**4*x3**4 + sqrtq**10*x1*x2**5*x3**6 + sqrtq**10*x1*x2**4*x3**7 + sqrtq**10*x2**5*x3**7 - sqrtq**8*x1**5*x2**5*x3**3 - sqrtq**8*x1**3*x2**5*x3**5 - sqrtq**10*x1*x2**4*x3**6 - sqrtq**10*x2**5*x3**6 - sqrtq**10*x2**4*x3**7 + sqrtq**8*x1**5*x2**5*x3**2 + sqrtq**8*x1**5*x2**4*x3**3 + sqrtq**8*x1**4*x2**5*x3**3 + sqrtq**8*x1**3*x2**5*x3**4 + sqrtq**8*x1**3*x2**4*x3**5 + sqrtq**8*x1**2*x2**5*x3**5 + sqrtq**10*x2**4*x3**6 - sqrtq**8*x1**5*x2**4*x3**2 - sqrtq**8*x1**4*x2**5*x3**2 - sqrtq**8*x1**4*x2**4*x3**3 - sqrtq**8*x1**3*x2**4*x3**4 - sqrtq**8*x1**2*x2**5*x3**4 + sqrtq**8*x1**3*x2**3*x3**5 - sqrtq**8*x1**2*x2**4*x3**5 + sqrtq**8*x1**4*x2**4*x3**2 - sqrtq**8*x1**3*x2**3*x3**4 + sqrtq**8*x1**2*x2**4*x3**4 - sqrtq**8*x1**3*x2**2*x3**5 - sqrtq**8*x1**2*x2**3*x3**5 - sqrtq**6*x1**3*x2**5*x3**3 + sqrtq**8*x1**3*x2**2*x3**4 + sqrtq**8*x1**2*x2**3*x3**4 + sqrtq**8*x1**2*x2**2*x3**5 + sqrtq**6*x1**3*x2**5*x3**2 + sqrtq**6*x1**3*x2**4*x3**3 + sqrtq**6*x1**2*x2**5*x3**3 - sqrtq**8*x1**2*x2**2*x3**4 - sqrtq**6*x1**3*x2**4*x3**2 - sqrtq**6*x1**2*x2**5*x3**2 + sqrtq**6*x1**3*x2**3*x3**3 - sqrtq**6*x1**2*x2**4*x3**3 + sqrtq**6*x1*x2**3*x3**5 - sqrtq**6*x1**3*x2**3*x3**2 + sqrtq**6*x1**2*x2**4*x3**2 - sqrtq**6*x1**3*x2**2*x3**3 - sqrtq**6*x1**2*x2**3*x3**3 - sqrtq**6*x1*x2**3*x3**4 - sqrtq**6*x1*x2**2*x3**5 - sqrtq**6*x2**3*x3**5 + sqrtq**6*x1**3*x2**2*x3**2 + sqrtq**6*x1**2*x2**3*x3**2 + sqrtq**6*x1**2*x2**2*x3**3 + sqrtq**6*x1*x2**2*x3**4 + sqrtq**6*x2**3*x3**4 + sqrtq**6*x2**2*x3**5 - sqrtq**6*x1**2*x2**2*x3**2 - sqrtq**6*x2**2*x3**4 + sqrtq**4*x1*x2**3*x3**3 - sqrtq**4*x1*x2**3*x3**2 - sqrtq**4*x1*x2**2*x3**3 - sqrtq**4*x2**3*x3**3 + sqrtq**2*x1**3*x2**3*x3 + sqrtq**4*x1*x2**2*x3**2 + sqrtq**4*x2**3*x3**2 + sqrtq**4*x2**2*x3**3 - sqrtq**2*x1**3*x2**3 - sqrtq**2*x1**3*x2**2*x3 - sqrtq**2*x1**2*x2**3*x3 - sqrtq**4*x2**2*x3**2 + sqrtq**2*x1**3*x2**2 + sqrtq**2*x1**2*x2**3 + sqrtq**2*x1**2*x2**2*x3 - sqrtq**2*x1**2*x2**2 - x1*x2*x3 + x1*x2 + x1*x3 + x2*x3 - x1 - x2 - x3 + 1);
f2 = F(f1);