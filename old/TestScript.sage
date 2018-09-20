reset('R,F,a,w,m,v,f0,f')

#################################################################
#Change the Cartan Type to see the effects on the polynomial ring
R = RootSystem(['A', 2])

#Initialize the polynomial ring off of the given root system
F = ComplexFieldOfRationalFunctionsWithWeylAction(R)
print "1) Initializing " + str(F) + "\n"

print "===============Miscellaneous Functionality of Polynomial Ring===============\n"

#Choose a random element of the root lattice for example purposes
a = F.lattice.random_element()
print "2) Choosing random element of the root lattice 'a' =\n" + str(a) + "\n"

#Choose a random element of the root lattice for example purposes
w = F.weyl_group.random_element()
print "3) Choosing random element of the corresponding weyl group 'w' =\n" + str(w) + "\n(reduced word: " + str(w.reduced_word()) + ")\n"

#Construct complex monomial
m = F.x(a)
print "4) Constructing complex monomial 'm' for 'a' (described on page 8):\n" + str(m) + "\n"

#Find phi(w)
print "5) Find positive roots sent to negative roots by 'w' (page 7):\n" + str(F.phi(w)) + "\n"

#Load vector of generators of the polynomial ring
v = F.v()
print "6) Loading vector of generators of the polynomial ring:\n" + str(v) + "\n"

print "===============Weyl Action on Polynomial Ring===============\n"

#Act by simple reflections on the generators
print "7) Acting by simple reflections on the generators of the polynomial ring (Eq. 3.8 on page 8)..."
for i in range(0,len(v)):
    print "Reflection " + str(i) + ": " + str(F.sigma_action(i, v))

print

#Act by simple reflections on the generators
print "8) Acting by epsilon on the generators of the polynomial ring (Eq. 3.10 on page 9)..."
for i in range(0,len(v)):
    print "Reflection " + str(i) + ": " + str(F.epsilon_action(i, v))

print

#Act by weyl group element on the generators
print "9) Acting by 'w' on the generators of the polynomial ring:\n" + str(F.w_action(w, v)) +"\n"

#Act on the monomial by simple reflections
print "10) Acting on 'm' by simple reflections (Eq. 3.13 on page 9)..."
for i in range(0,len(v)):
    print "Reflection " + str(i) + ": " + str(F.weyl_action_by_simple_reflection(i, m))

print

#Act by weyl group element on the monomial
wm = F.weyl_action(w, m)
print "11) Acting by 'w' on 'm':\n" + "Value saved in variable wm\n"

print "===============Functionality for Constructing Invariant Polynomial===============\n"

#Create delta polynomial
print "12) Creating delta polynomial (page 12):\n" + str(F.delta()) + "\n"

#Create j polynomial
print "13) Creating j polynomial for 'w' (page 12):\n" + str(F.j(w)) + "\n"

#Create delta polynomial
f0 = F.f0()
print "14) Creating f0 polynomial (eq. 3.20 page 12):\n" + "Value saved in variable f0\n"

#Create delta polynomial
f = F.weyl_invariant_polynomial()
print "15) Creating invariant polynomial (eq. 3.21 page 12):\n" + "Value saved in variable f"