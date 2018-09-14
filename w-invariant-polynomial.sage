#for testing purposes
F = PolynomialRing(CC, 'x', 2)
V = VectorSpace(CC, 2)
R = RootSystem(['A',2])
L = R.root_lattice()
W = L.weyl_group()

#generate complex monomial based off of a member of the root lattice
#
#I assume that 'm' is a linear combination of the basis elements of
#the root lattice
def x(a):
    n = a.parent().dimension()
    p = 1
    for i in range(0,n):
        p = p*(var('x'+str(i)))^(a.coefficient(i+1))
    return p

#generate the complex delta polynomial based off of a given lattice
#for a root system
def delta(L):
    p = 1
    q = var('q')
    for a in L.positive_roots():
        cfs = a.coefficients()
        p = p * (1 - q^sum(cfs) * x(a))
    return p

#generate half sum of positive roots of lattice
def half_sum_positive_roots(L):
    return sum(L.positive_roots()).map_coefficients(lambda x : ZZ((1/2)*x))

#generate a list of symbolic variables n long
def v(n):
    v = []
    for i in range(0,n):
        v.append(var('x'+str(i)))
    return v

#compute j from scratch
#w - Weyl Group element
#L - root lattice
def calc_j(w, L):
    d = delta(L)
    vrs = v(L.dimension())
    return d/evaluate(d, w_action(w, vrs, W))

#return j based on shortcuts in the paper
def j(w, L):
    rw = w.reduced_word()
    q = var('q')
    if len(rw) == 1:
        return -q*var('x'+str(rw[0]-1))^2
    else:
        a = half_sum_positive_roots(L)
        a = a - w.inverse().action(a)
        b = 2 * a
        cfs = a.coefficients()
        return (-1)^len(rw) * q^(sum(cfs)) * x(2*a)

#Weyl group action on polynomial ring

#Sigma action
def sigma_action(i, v, W):
    q = var('q')
    sr = W.simple_reflections()
    w = []
    for j in range(0, len(v)):
        if i == j:
            w.append(1/(q*v[j]))
        elif (sr[i+1] * sr[j+1])^3 == W.random_element_of_length(0):
            w.append(v[j]*v[i]*sqrt(q))
        else:
            w.append(v[j])
    return w

#W action
def w_action(w, v, W):
    word = w.reduced_word()
    word.reverse()
    for i in word:
        v = sigma_action(i-1, v, W)
    return v

#Epsilon action
def epsilon_action(i, v, W):
    sr = W.simple_reflections()
    w = []
    for j in range(0, len(v)):
        if (sr[i+1] * sr[j+1])^3 == W.random_element_of_length(0):
            w.append(-v[j])
        else:
            w.append(v[j])
    return w

def f_plus(i, v, f, W):
    return (evaluate(f, v) + evaluate(f, epsilon_action(i, v, W)))/2

def f_minus(i, v, f, W):
    return (evaluate(f, v) - evaluate(f, epsilon_action(i, v, W)))/2

#i - ith action
#f - polynomial
#W - Weyl group
#F - corresponding polynomial ring
def weyl_action_by_simple_reflection(i, f, W, F):
    q = var('q')
    vrs = v(len(F.gens()))
    A = (-(1-q*vrs[i])/((q*vrs[i])*(1-vrs[i])))
    B = (1/(vrs[i]*sqrt(q)))
    return A*f_plus(i, vrs, evaluate(f, sigma_action(i, vrs, W)), W) + B*f_minus(i, vrs, evaluate(f, sigma_action(i, vrs, W)), W)

#w - Weyl group element
def weyl_action(w, f, W, F):
    word = w.reduced_word()
    word.reverse()
    for i in word:
        f = weyl_action_by_simple_reflection(i-1, f, W, F)
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

#helpers

#needed because sigma and epsilon actions are defined on vectors (which here is an array)
#the symbolic ring has no way to evaluate a function on an array especially with 'q' thrown into the mix
def evaluate(f, v):
    if f.parent() == SR:
        n = len(f.variables())
        if var('q') in f.variables():
            n = n - 1
        for i in range(0, n):
            x = var('x'+str(i))
            f = f.function(x)
            f = f(v[i])
        return f
    else:
        return f(v)