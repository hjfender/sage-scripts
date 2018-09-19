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
    

        def c(self, i):
            return (1/2)*((q*self.variables[i] - 1)/(q*self.variables[i]*(1 - self.variables[i])) + 1/(sqrt(q)*self.variables[i]))

        def d(self, i):
            return (1/2)*((q*self.variables[i] - 1)/(q*self.variables[i]*(1 - self.variables[i])) - 1/(sqrt(q)*self.variables[i]))