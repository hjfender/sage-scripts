#Change the Cartan Type to see the effects on the polynomial ring
R = RootSystem(['A', 3])

#Initialize the polynomial ring off of the given root system
F = ComplexPolynomialRingWithWeylAction(R)
print "Initializing " + str(F) + "\n"

#Choose a random element of the root lattice for example purposes
a = F.lattice.random_element()
print "Choosing random element of the root lattice 'alpha' = " + str(a) + "\n"

#Construct complex monomial
m = F.x(a)
print "Constructing complex monomial 'm' for 'alpha' (described on page 8): " + str(m) + "\n"

#Load vector of generators of the polynomial ring
v = F.v()
print "Loading vector of generators of the polynomial ring: " + str(v) + "\n"

#Act by simple reflections on the generators
print "Acting by simple reflections on the generators of the polynomial ring (Eq. 3.8 on page 8)..."
for i in range(0,len(v)):
    print "Reflection " + str(i) + ": " + str(F.sigma_action(i, v))

print

#Act by simple reflections on the generators
print "Acting by epsilon on the generators of the polynomial ring (Eq. 3.10 on page 9)..."
for i in range(0,len(v)):
    print "Reflection " + str(i) + ": " + str(F.epsilon_action(i, v))

print

#Choose a random element of the root lattice for example purposes
w = F.weyl_group.random_element()
print "Choosing random element of the corresponding weyl group 'w' = " + str(w) + "\n"

#Act by weyl group element on the generators
print "Acting by 'w' on the generators of the polynomial ring " + str(F.w_action(w, v))

#Act on the monomial by simple reflections
print "Acting on 'm' by simple reflections (Eq. 3.13 on page 9)..."
for i in range(0,len(v)):
    print "Reflection " + str(i) + ": " + str(F.weyl_action_by_simple_reflection(i, m))

print

#Act by weyl group element on the monomial
print "Acting by 'w' on 'm' " + str(F.weyl_action(w, m)) + "\n"

#Create delta polynomial
print "Creating delta polynomial (page 12): " + str(F.delta()) + "\n"

#Create j polynomial
print "Creating j polynomial for 'w' (page 12): " + str(F.j(w)) + "\n"

#Create delta polynomial
print "Creating f0 polynomial (eq. 3.20 page 12): " + str(F.f0()) + "\n"

#Create delta polynomial
print "Creating invariant polynomial (eq. 3.21 page 12): " + str(F.weyl_invariant_polynomial())