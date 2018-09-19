import sage.combinat.root_system.root_system

class ComplexPolynomialRingWithWeylAction:
    ###############################################################
    #Class related methods
    ###############################################################

    def __init__(self, root_system):
        self.root_system = root_system
        self.lattice = root_system.root_lattice()
        self.weyl_group = self.lattice.weyl_group()

    def __str__(self):
        return ("Complex polynomial Ring in " + self.lattice.dimension()
                " variables endowed with an action by " + str(self.weyl_group))

    def __repr__(self):
        return str(self)

    ###############################################################
    #Methods related to the creation and manipulation
    #of the poynomial ring generators
    ###############################################################

    def variables(self):
        v = []
        for i in range(0,self.lattice.dimension()):
            v.append(var('x'+str(i)))
        return v

    def sigma_action(self, i, v):
        if len(v) != self.lattice.dimension():
            raise ValueError('Invalid input vector sizes.')
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

    def w_action(self, w, v):
        if len(v) != self.lattice.dimension():
            raise ValueError('Invalid input vector sizes.')
        word = w.reduced_word()
        word.reverse()
        u = v
        for i in word:
            u = self.sigma_action(i-1, u)
        return u

    def epsilon_action(self, i, v):
        if len(v) != self.lattice.dimension():
            raise ValueError('Invalid input vector sizes.')
        sr = self.weyl_group.simple_reflections()
        w = []
        for j in range(0, len(v)):
            if (sr[i+1] * sr[j+1])^3 == self.weyl_group.random_element_of_length(0):
                w.append(-v[j])
            else:
                w.append(v[j])
        return w

    ###############################################################
    #Methods related to the weyl action on polynomials
    ###############################################################

    #The following two methods are components of the method immediately following them
    def f_plus(self, i, f):
        v = self.v()
        return (self.evaluate(f, v) + self.evaluate(f, self.epsilon_action(i, v)))/2

    def f_minus(self, i, f):
        v = self.v()
        return (self.evaluate(f, v) - self.evaluate(f, self.epsilon_action(i, v)))/2

    def weyl_action_by_simple_reflection(self, i, f):
        q = var('q')
        v = self.v()
        A = (-(1-q*v[i]))/((q*v[i])*(1-v[i]))
        B = (1/(v[i]*sqrt(q)))
        return A*self.evaluate(self.f_plus(i,f),self.sigma_action(i,v)) + B*self.evaluate(self.f_minus(i,f),self.sigma_action(i,v))

    def weyl_action(self, w, f):
        word = w.reduced_word()
        for i in word:
            f = self.weyl_action_by_simple_reflection(i-1, f)
        return self.evaluate(f,self.v())

    ################################################################
    #Methods related to the construction of the invariant polynomial
    ################################################################

    def delta(self):
        p = 1
        q = var('q')
        for a in self.lattice.positive_roots():
            cfs = a.coefficients()
            p = p * (1 - q^sum(cfs) * self.x(2*a))
        return p

    def j(self, w):
        rw = w.reduced_word()
        v = self.v()
        q = var('q')
        if len(rw) == 1:
            return -q*v[rw[0]-1]^2
        else:
            a = sum(self.phi(w))
            cfs = a.coefficients()
            return (-1)^len(rw) * q^(sum(cfs)) * self.x(2*a)

    def f0(self):
        f0 = []
        one = self.one()
        for w in self.weyl_group.list():
            g = self.j(w) * self.weyl_action(w, one)
            f0.append(g)
        return sum(f0)

    def weyl_invariant_polynomial(self):
        return (self.f0()/self.delta())

    ################################################################
    #Miscellaneous useful methods
    ################################################################
    def x(self, a):
        v = self.v()
        p = 1
        for i in range(0,len(v)):
            p = p*(v[i])^(a.coefficient(i+1))
        return p

    def phi(self, w):
        w_roots = []
        for a in self.lattice.positive_roots():
            wa = w.action(a)
            for b in self.lattice.negative_roots():
                if wa == b:
                    w_roots.append(a)
                    break
        if len(w_roots) == 0:
            w_roots.append(0*F.lattice.gens()[0])
        return w_roots

    def one(self):
        v = self.v()
        one = SR(1)
        return one.function(*v)

    ###############################################################
    #Rational Function Class to handle weyl action
    ###############################################################
    class RationalFunction:

        def __init__(self, rf, w, vrs):
            self.rational_function = rf
            self.weyl_group_element = w
            self.variables = vrs
            self.q = var('q')

        def evaluate(self, inputs):
            if self.polynomial.parent() == SR:
                f = self.polynomial
                n = len(self.variables)
                if var('q') in f.variables():
                    n = n - 1
                for i in range(0, n):
                    x = self.variables[i]
                    f = f.function(x)
                    f = f(inputs[i])
                return f
            else:
                try:
                    return f(inputs)
                except:
                    return f(*inputs)

        def c(self, i):
            return (1/2)*((q*self.variables[i] - 1)/(q*self.variables[i]*(1 - self.variables[i])) + 1/(sqrt(q)*self.variables[i]))

        def d(self, i):
            return (1/2)*((q*self.variables[i] - 1)/(q*self.variables[i]*(1 - self.variables[i])) - 1/(sqrt(q)*self.variables[i]))

    class Variables:

        def __init__(self, n, adjacency_map):
            self.q = var('q')
            self.vars = []
            for i in range(0,n):
                self.list.append(var('x'+str(i)))
            self.vars_after_sigma = []
            self.vars_after_epsilon = []