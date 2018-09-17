#Change the Cartan Type to see the effects on the polynomial ring
R = RootSystem(['A', 2])

#Initialize the polynomial ring off of the given root system
F = ComplexPolynomialRingWithWeylAction(R)
print "Initializing " + str(F) + "..."

#Choose a random element of the root lattice for example purposes
a = F.lattice.random_element()
print "Choosing random element of root lattice 'alpha' = " + str(a)

#Construct complex monomial
m = F.x(a)
print "Construct complex monomial for 'alpha' (described on page 8): " + str(m)

