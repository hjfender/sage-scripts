#!/usr/bin/env sage -python

import sys
# sys.path.append('/Applications/sage-8.4/src')
# sys.path.append('/Applications/sage-8.4/src/sage')
# sys.path.append('/Applications/sage-8.4/src/sage/misc')
# sys.path.append('/Applications/sage-8.4/src/sage/env')

# import os
# os.environ['SAGE_ROOT'] = '/Applications/sage-8.4'
# os.environ['SAGE_SRC'] = '/Applications/sage-8.4/src'
# os.environ['SAGE_DOC_SRC'] = '/Applications/sage-8.4/src/doc'
# os.environ['SAGE_LOCAL'] = '/Applications/sage-8.4/local'
# os.environ['DOT_SAGE'] = '/Users/henry/.sage/'

from weyl_group_multiple_dirichlet_series_v2.field_of_rational_functions_with_weyl_action import *

R = RootSystem(['A', 2])
F = FieldOfRationalFunctionsWithWeylAction(R)

W = F.weyl_group

for w in W.list(): 
    g = F.act(F.one(),w)
    g = F.reduce_function(g)
    print ("\t" + str(w.reduced_word()) + ": " + str(g))