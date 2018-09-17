import sage.combinat.root_system.root_system

class ComplexPolynomialRingWithWeylAction:
    """Class comment"""

    def __init__(self, root_system):
        self.root_system = root_system
        self.lattice = root_system.root_lattice()
        self.weyl_group = self.lattice.weyl_group()
        self.polynomial_ring = PolynomialRing(CC, 'x', self.lattice.dimension())

    def __str__(self):
        return str(self.polynomial_ring) + " endowed with an action by " + str(self.weyl_group)

    def __repr__(self):
        return str(self)

    #Generate monomial based off of a member of the root lattice
    #I assume that 'm' is a linear combination of the basis elements of the root lattice
    def x(self, a):
        v = self.v()
        p = 1
        for i in range(0,len(v)):
            p = p*(v[i])^(a.coefficient(i+1))
        return p

    #Generate half sum of positive roots of lattice
    def half_sum_positive_roots(self):
        return sum(self.lattice.positive_roots()).map_coefficients(lambda x : ZZ((1/2)*x))
    
    #Action on symbolic variables by simple reflections
    def sigma_action(self, i, v):
        q = var('q')
        sr = self.weyl_group.simple_reflections()
        w = []
        for j in range(0, len(v)):
            if i == j:
                w.append(1/(q*v[j]))
            elif (sr[i+1]*sr[j+1])^3 == self.weyl_group.random_element_of_length(0):
                w.append(v[j]*v[i]*sqrt(q))
            else:
                w.append(v[j])
        return w

    #Action on symbolic variables by weyl group element
    def w_action(self, w, v):
        word = w.reduced_word()
        word.reverse()
        u = v
        for i in word:
            u = self.sigma_action(i-1, u)
        return u

    #Action on symbolic variables switching sign of adjacent variables
    def epsilon_action(self, i, v):
        sr = self.weyl_group.simple_reflections()
        w = []
        for j in range(0, len(v)):
            if (sr[i+1] * sr[j+1])^3 == self.weyl_group.random_element_of_length(0):
                w.append(-v[j])
            else:
                w.append(v[j])
        return w

    #Helpful for the weyl action
    def f_plus(self, i, f):
        v = self.v()
        return (self.evaluate(f, v) + self.evaluate(f, self.epsilon_action(i, v)))/2

    def f_minus(self, i, f):
        v = self.v()
        return (self.evaluate(f, v) - self.evaluate(f, self.epsilon_action(i, v)))/2

    #Component of the weyl action
    def weyl_action_by_simple_reflection(self, i, f):
        q = var('q')
        v = self.v()
        A = (-(1-q*v[i]))/((q*v[i])*(1-v[i]))
        B = (1/(v[i]*sqrt(q)))
        return A*self.evaluate(self.f_plus(i,f),self.sigma_action(i,v)) + B*self.evaluate(self.f_minus(i,f),self.sigma_action(i,v))

    #w - Weyl group element
    def weyl_action(self, w, f):
        word = w.reduced_word()
        # word.reverse() #right action or left action???
        for i in word:
            f = self.weyl_action_by_simple_reflection(i-1, f)
        return self.evaluate(f,self.v())

    #Helpful for finding the invariant polynomial

    #Generate the complex delta polynomial based for the root system
    def delta(self):
        p = 1
        q = var('q')
        for a in self.lattice.positive_roots():
            cfs = a.coefficients()
            p = p * (1 - q^sum(cfs) * self.x(2*a))
        return p

    #Find j based on shortcuts in the paper
    def j(self, w):
        rw = w.reduced_word()
        v = self.v()
        q = var('q')
        if len(rw) == 1:
            return -q*v[rw[0]-1]^2
        else:
            a = self.half_sum_positive_roots()
            a = a - w.inverse().action(a)
            cfs = a.coefficients()
            return (-1)^len(rw) * q^(sum(cfs)) * self.x(2*a)

    def f0(self):
        f0 = []
        one = self.one()
        for w in self.weyl_group.list():
            g = self.j(w) * self.weyl_action(w, one)
            f0.append(g)
        return sum(f0).simplify_full()

    #Find the invariant polynomial
    def weyl_invariant_polynomial(self):
        return (self.f0()/self.delta()).simplify_full()

    #Generate a list of symbolic variables one for each dimension of the lattice
    def v(self):
        v = []
        for i in range(0,self.lattice.dimension()):
            v.append(var('x'+str(i)))
        return v

    #Create the constant symbolic function 1
    def one(self):
        v = self.v()
        one = SR(1)
        return one.function(*v)

    #evaluate a function in the ring agains a list of variables representing [x0,x1,...,xn]
    def evaluate(self, f, v):
        if f.parent() == SR:
            n = len(f.variables())
            vrs = self.v()
            if var('q') in f.variables():
                n = n - 1
            for i in range(0, n):
                x = vrs[i]
                f = f.function(x)
                f = f(v[i])
            return f
        else:
            try:
                return f(v)
            except:
                return f(*v)
    
    #needed because sage can't seem to do it by itself
    def simplify(self, f):
        inputs = [0 for i in range(0,len(xs))]
        for i in range(0,len(xs)):
            inputs[i] = v[i]
            if i != 0:
                inputs[i-1] = 0
            g = self.evaluate(f,inputs)
            g_simp = g.simplify_full()
            print str(g) + " and " + str(g_simp)
            f = ((f/g)*g_simp).simplify_full()
            print "Sub " + str(f)
        # return f

#FOR TESTING PURPOSES
R = RootSystem(['A', 2])
F = ComplexPolynomialRingWithWeylAction(R)