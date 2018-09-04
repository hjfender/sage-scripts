#for testing purposes
F1 = PolynomialRing(CC, 'x', 2)
V1 = VectorSpace(CC, 2)
R1 = RootSystem(['A',2])
L1 = R1.root_lattice()
W1 = L1.weyl_group()

#generate complex monomial based off of a member of the root lattice
#
#I assume that 'm' is a linear combination of the basis elements of
#the root lattice
def get_complex_monomial(m):
    n = m.parent().dimension()
    F = PolynomialRing(CC, 'x', n)
    p = 1
    for i in range(0,n):
        k = m.coefficient(i+1)
        p = p*(F.gen(i))^(m.coefficient(i+1))
    return p

#generate complex delta polynomial based off of a given lattice
#for a root system
def get_delta_polynomial(L):
    p = 1
    q = var('q')
    for a in L.positive_roots():
        cfs = a.coefficients()
        p = p * (1 - q^sum(cfs) * get_complex_monomial(a))
    return p

#generate half sum of positive roots of lattice
def half_sum_positive_roots(L):
    return sum(L.positive_roots()).map_coefficients(lambda x : ZZ((1/2)*x))

#compute j from scratch
#w - Weyl Group element
#L - root lattice
# def calc_j(w, L, F):
#     d = get_delta_polynomial(L)
#     inputs = {}
#     for i in range(0, L.dimension()):
#         inputs[F.gen(i)] = weyl_action(w, F.gen(i), L.weyl_group(), F)
#     print(inputs)
#     return d/substitute(d, inputs)

#return j based on shortcuts in the paper
def j(w, L):
    rw = w.reduced_word()
    F = PolynomialRing(CC, 'x', L.dimension())
    q = var('q')
    if len(rw) == 1:
        return -q*F.gen(rw[0]-1)^2
    else:
        a = half_sum_positive_roots(L)
        a = a - w.inverse().action(a)
        b = 2 * a
        cfs = a.coefficients()
        return (-1)^len(rw) * q^(sum(cfs)) * get_complex_monomial(2*a)

#Weyl group action on polynomial ring

#Sigma action
def sigma_action(i, f, W, F):
    n = len(F.gens())
    q = var('q') #this converts everything to a symbolic ring
    sr = W.simple_reflections()
    inputs = {}
    for j in range(0,n):
        if i == j:
            inputs[F.gen(j)] = (1/(q*F.gen(j)))
        elif (sr[i+1] * sr[j+1])^3 == W.random_element_of_length(0):
            inputs[F.gen(j)] = F.gen(i)*F.gen(j)*sqrt(q)
        else:
            inputs[F.gen(j)] = F.gen(j)
    return substitute(f, inputs)

#Epsilon action
#i - ith action
#f - polynomial
#W - Weyl group
def epsilon_action(i, f, W, F):
    n = len(F.gens())
    sr = W.simple_reflections()
    inputs = {}
    for j in range(0,n):
        if (sr[i+1] * sr[j+1])^3 == W.random_element_of_length(0):
            inputs[F.gen(j)] = -F.gen(j)
        else:
            inputs[F.gen(j)] = F.gen(j)
    return substitute(f, inputs)

#i - ith action
#f - polynomial
#W - Weyl group
def weyl_action_by_simple_reflection(i, f, W, F):
    fplus = (f + epsilon_action(i, f, W, F))/2
    fminus = (f - epsilon_action(i, f, W, F))/2
    q = var('q')
    return (-(1-q*F.gen(i))/(q*F.gen(i)*(1-F.gen(i))))*sigma_action(i, fplus, W, F) + (1/(F.gen(i)*sqrt(q)))*sigma_action(i, fminus, W, F)

#w - Weyl group element
def weyl_action(w, f, W, F):
    word = w.reduced_word()
    word.reverse()
    for i in word:
        f = weyl_action_by_simple_reflection(i-1,f,W,F)
    return f

#THE POLYNOMIAL PROPER
def f_0(L):
    f0 = 0
    W = L.weyl_group()
    F = PolynomialRing(CC, 'x', L.dimension())
    for w in W.list():
        f0 += j(w,L) * weyl_action(w, 1, W, F)
    return f0

def W_invariant_polynomial(L):
    return f_0(L)/get_delta_polynomial(L)

#helpers - needed because f can be a member of CC[x1,x2,...] or SR
#and those rings behave differently with respect to substitution
def substitute(f, inputs):
    if f.parent() == SR:
        vrs = f.variables()
        for v in vrs:
            for k in inputs:
                if str(k) == str(v):
                    f = f.subs({v : inputs[k]})
                    break
        return f
    else:
        return f.subs(inputs)