reset()

from weyl_group_multiple_dirichlet_series.variable_helper import *
from weyl_group_multiple_dirichlet_series.rational_function_with_action import *
from weyl_group_multiple_dirichlet_series.field_of_rational_functions_with_action import *

#################################################################
#Change the Cartan Type to see the effects on the polynomial ring
R = RootSystem(['A', 2])
#################################################################

#Initialize the polynomial ring off of the given root system
F = FieldOfRationalFunctionsWithAction(R)
print "1) Initializing " + str(F) + "\n"

L = F.lattice
W = F.weyl_group

print "===============Miscellaneous Functionality of Polynomial Ring===============\n"

#Choose a random element of the root lattice for example purposes
a = L.random_element()
print "2) Choosing random element of the root lattice 'a' =\n" + str(a) + "\n"

#Choose a random element of the root lattice for example purposes
w = W.random_element()
print "3) Choosing random element of the corresponding weyl group 'w' =\n" + str(w) + "\n(reduced word: " + str(w.reduced_word()) + ")\n"

#Construct complex monomial
m = F.x(a)
print "4) Constructing complex monomial 'm' for 'a' (described on page 8):\n" + str(m) + "\n"

#Find phi(w)
print "5) Find positive roots sent to negative roots by 'w' (page 7):\n" + str(F.phi(w)) + "\n"

#Load vector of generators of the polynomial ring
v = F._variable_helper.variables
print "6) Loading vector of generators of the polynomial ring:\n" + str(v) + "\n"

print "===============Weyl Action on Polynomial Ring===============\n"

#Act by simple reflections on the generators
print "7) Acting by simple reflections on the generators of the polynomial ring (Eq. 3.8 on page 8)..."
for i in range(0,len(v)):
    print "Reflection " + str(i) + ": " + str(F._variable_helper.act_on_variables("s"+str(i)))

print

#Act by epsilon on the generators
print "8) Acting by epsilon on the generators of the polynomial ring (Eq. 3.10 on page 9)..."
for i in range(0,len(v)):
    print "Action " + str(i) + ": " + str(F._variable_helper.act_on_variables("e"+str(i)))

print

#Demonstrate the usefulness of act_on_variables
print "9) Checking the relations (3.9) on page 9..."
sr = W.simple_reflections()
identity = W.random_element_of_length(0)
for i in range(0,len(v)):
    for j in range(i,len(v)):
        if (i != j and (sr[i+1] * sr[j+1])**3 == identity):
            print str(i) + " and " + str(j) + " are adjacent..." 
            print ("\ts" + str(i) + "s" + str(j) + "s" + str(i) + "x" \
                + " = " + str(F._variable_helper.act_on_variables("s"+str(i)+"s"+str(j)+"s"+str(i))))
            print ("\ts" + str(j) + "s" + str(i) + "s" + str(j) + "x" \
                + " = " + str(F._variable_helper.act_on_variables("s"+str(j)+"s"+str(i)+"s"+str(j))))
        elif i != j:
            print str(i) + " and " + str(j) + " are not adjacent..."
            print ("\ts" + str(i) + "s" + str(j) + "x" \
                + " = " + str(F._variable_helper.act_on_variables("s"+str(j)+"s"+str(i))))
            print ("\ts" + str(j) + "s" + str(i) + "x" \
                + " = " + str(F._variable_helper.act_on_variables("s"+str(i)+"s"+str(j))))
        else:
            print str(i) + " and " + str(j) + " are equal..."
            print ("\ts" + str(i) + "s" + str(j) + "x" \
                + " = " + str(F._variable_helper.act_on_variables("s"+str(i)+"s"+str(j))))

print

print "10) Checking the relations (3.11) on page 9..."
for i in range(0,len(v)):
    for j in range(0,len(v)):
        if (i != j and (sr[i+1] * sr[j+1])**3 == identity):
            print str(i) + " and " + str(j) + " are adjacent..." 
            print ("\ts" + str(i) + "e" + str(j) + "x" \
                + " = " + str(F._variable_helper.act_on_variables("e"+str(j)+"s"+str(i))))
            print ("\te" + str(i) + "e" + str(j) + "s" + str(i) + "x" \
                + " = " + str(F._variable_helper.act_on_variables("s"+str(i)+"e"+str(j)+"e"+str(i))))
        else:
            print str(i) + " and " + str(j) + " are not adjacent..."
            print ("\ts" + str(i) + "e" + str(j) + "x" \
                + " = " + str(F._variable_helper.act_on_variables("e"+str(j)+"s"+str(i))))
            print ("\te" + str(j) + "s" + str(i) + "x" \
                + " = " + str(F._variable_helper.act_on_variables("s"+str(i)+"e"+str(j))))

print

#Act on the monomial by simple reflections
print "10) Acting on 'm' by simple reflections (Eq. 3.13 on page 9)..."
for u in W.simple_reflections():
    print "Reflection " + str(u.reduced_word()) + ": " + str(m.act(u))

print

#Act by weyl group element on the monomial
wm = m.act(w)
print "11) Acting by 'w' on 'm':\n" + str(wm) + "\n"

print "===============Functionality for Constructing Invariant Polynomial===============\n"

#Create delta polynomial
print "12) Creating delta polynomial (page 12):\n" + str(F.delta()) + "\n"

#Create j polynomial
print "13) Creating j polynomial for 'w' (page 12):\n" + str(F.j(w)) + "\n"

#Create delta polynomial
f0 = F.f0()
print "14) Creating f0 polynomial (eq. 3.20 page 12):\n" + str(f0) + "\n"

#Create delta polynomial
f = F.invariant_function()
print "15) Creating invariant polynomial (eq. 3.21 page 12):\n" + str(f)