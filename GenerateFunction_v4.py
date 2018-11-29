#!/usr/bin/env sage -python

from weyl_group_multiple_dirichlet_series_v3.field_of_rational_functions_with_weyl_action import *

#############################################################################
#Change the Cartan Type to see the effects on the field of rational functions
cartan = 'F'
n = 4

R = RootSystem([cartan, n])
#############################################################################

#Initialize the field of rational functions off of the given root system
F = FieldOfRationalFunctionsWithWeylAction(R)

f = F.invariant_function()
filename = "output/" + cartan + str(n) + ".txt"
func_file = open(filename,"w+")
func_file.write(str(f))
func_file.close()