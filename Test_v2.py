#!/usr/bin/env sage

from weyl_group_multiple_dirichlet_series_v2.field_of_rational_functions_with_weyl_action import *

R = RootSystem(['A', 2])
F = FieldOfRationalFunctionsWithWeylAction(R)

W = F.weyl_group

for w in W.list():
    g = F.act(F.one(),w)
    g = F.reduce_function(g)
    # if len(w.reduced_word()) < 2: g = F.reduce_function(g)
    print ("\t" + str(w.reduced_word()) + ": " + str(g))