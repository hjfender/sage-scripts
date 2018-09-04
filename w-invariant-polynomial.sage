#generate complex monomial based off of a member of the root lattice
#
#I assume that 'latticeMember' is a linear combination of the basis elements of
#the root lattice
def get_complex_monomial(latticeMember):
    n = latticeMember.parent().dimension()
    F = PolynomialRing(CC, 'x', n)
    p = 1
    for i in range(0,n):
        p = p*(F.gen(i))^(latticeMember.coefficient(i+1))
    return p

#generate complex delta polynomial based off of a given lattice
#for a root system
def get_delta_polynomial(rootLattice):
    p = 1
    q = var('q')
    for a in rootLattice.positive_roots():
        p = p * (1 - q^sum(a.coefficients()) * get_complex_monomial(a))
    return p

#generate half sum of positive roots of lattice
def half_sum_positive_roots(rootLattice):
    return sum(rootLattice.positive_roots()).map_coefficients(lambda x : (1/2)*x)

#compute j from scratch
def calc_j(weylGroupElement, rootLattice):
    return get_delta_polynomial(rootLattice)/get_delta_polynomial(rootLattice)

#return j based on shortcuts in the paper
# def j(weylGroupElement, rootLattice):
#     r = weylGroupElement.reduced_word()
#     F = PolynomialRing(CC, 'x', rootLattice.dimension())
#     q = var('q')
#     if len(r) == 0:
#         return 1
#     elif len(r) == 1:
#         return -q*F.gen(r.reduced_word()[0]-1)^2
#     else:
#         a = 
#         return (-1)^len(r) * q^

#Weyl group action on polynomial ring

#Sigma action
def sigma_action(i, f, W, F):
    n = len(F.gens())
    q = var('q') #this converts everything to a symbolic ring
    simpleReflections = W.simple_reflections()
    inputs = []
    for j in range(0,n):
        if (simpleReflections[i+1] * simpleReflections[j+1])^3 == W.random_element_of_length(0):
            inputs.append(F.gen(i)*F.gen(j)*sqrt(q))
        elif i == j:
            inputs.append(1/(q*F.gen(j)))
        else:
            inputs.append(F.gen(j))
    return f(inputs)

#Epsilon action
#i - ith action
#f - polynomial
#W - Weyl group
def epsilon_action(i, f, W, F):
    n = len(F.gens())
    simpleReflections = W.simple_reflections()
    inputs = []
    for j in range(0,n):
        if (simpleReflections[i+1] * simpleReflections[j+1])^3 == W.random_element_of_length(0):
            inputs.append(-F.gen(j))
        else:
            inputs.append(F.gen(j))
    return f(inputs)

#i - ith action
#f - polynomial
#W - Weyl group
def weyl_action_by_simple_reflection(i, f, W, F):
    fplus = (f + epsilon_action(i, f, W, F))/2
    fminus = (f - epsilon_action(i, f, W, F))/2
    q = var('q') #this converts everything to a symbolic ring
    return (-(1-q*F.gen(i))/(q*F.gen(i)*(1-F.gen(i))))*sigma_action(i, fplus, W, F) + (1/(F.gen(i)*sqrt(q)))*sigma_action(i, fminus, W, F)

#w - Weyl group element
def weyl_action(w, f, W, F):
    word = w.reduced_word()
    word.reverse()
    for i in word:
        f = weyl_action_by_simple_reflection(i-1,f,W,F)
    return f