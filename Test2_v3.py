#!/usr/bin/env sage -python

from weyl_group_multiple_dirichlet_series_v3.field_of_rational_functions_with_weyl_action import *

#############################################################################
#Change the Cartan Type to see the effects on the field of rational functions
R = RootSystem(['A', 2])
#############################################################################

#Initialize the field of rational functions off of the given root system
F = FieldOfRationalFunctionsWithWeylAction(R)

print "==========f=========="
f = F.invariant_function()
print f
print f.parent()
F.pretty_print_CF(f)

print "==========g=========="
g = F(f)
print g
print g.parent()
F.pretty_print(g)