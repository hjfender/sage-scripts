import sage.combinat.root_system.root_system

class ComplexFieldOfRationalFunctionsWithWeylAction:
    """
    A class encapsulating the relevant entities used to find a special invariant
    polynomial used to construct a Multiple Dirichlet Series. Those entities are:
        -A root system (passed in through the constructor)
        -A root lattice (corresponding to the root system)
        -A weyl group (corresponding to the root lattice)
        -A complex polynomial ring in n variables (n is the dimension of the root lattice)

    Example:
        sage: R = RootSystem(['A',2])
        sage: F = ComplexPolynomialRingWithWeylAction(R)
        sage: F
        Multivariate Polynomial Ring in x0, x1 over Complex Field with 53 bits of precision
        endowed with an action by Weyl Group of type ['A', 2] (as a matrix group acting on
        the root lattice)
        sage: F.root_system
        Root system of type ['A', 2]
        sage: F.lattice
        Root lattice of the Root system of type ['A', 2]
        sage: F.weyl_group
        Weyl Group of type ['A', 2] (as a matrix group acting on the root lattice)
        sage: F.polynomial_ring
        Multivariate Polynomial Ring in x0, x1 over Complex Field with 53 bits of precision
    """

    ###############################################################
    #Class related methods
    ###############################################################

    def __init__(self, root_system):
        self.root_system = root_system
        self.lattice = root_system.root_lattice()
        self.weyl_group = self.lattice.weyl_group()
        self.polynomial_ring = PolynomialRing(CC, 'x', self.lattice.dimension())

    def __str__(self):
        return str(self.polynomial_ring) + " endowed with an action by " + str(self.weyl_group)

    def __repr__(self):
        return str(self)

    ###############################################################
    #Methods related to the creation and manipulation
    #of the poynomial ring generators
    ###############################################################

    def v(self):
        """
        Generate a list of symbolic variables one for each generator
        of the polynomial ring. Used to in the actions below.

        Example:
            sage: F.v() #F is based on RootSystem(['A',2])
            [x0, x1]
        """
        v = []
        for i in range(0,self.lattice.dimension()):
            v.append(var('x'+str(i)))
        return v

    def sigma_action(self, i, v):
        """
        Action on a list of inputs 'v', length equal to the lattice dimension,
        by the ith generator of the weyl group, i.e. the ith simple reflection,
        as defined in (3.8) on page 8.

        Example:
            sage: F.sigma_action(0,F.v()) #F is based on RootSystem(['A',2])
            [1/(q*x0), sqrt(q)*x0*x1]

            sage: F.sigma_action(0,[1,1])
            [1/q, sqrt(q)]
        """
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
        """
        Action on symbolic variables 'v' by weyl group element 'w'.

        Example:
            sage: w = F.weyl_group.random_element() #F is based on RootSystem(['A',2])
            sage: w
            [ 0 -1]
            [-1  0]
            sage: w.reduced_word()
            [1, 2, 1]
            sage: F.w_action(w,F.v())
            [1/(q*x1), 1/(q*x0)]
        """
        if len(v) != self.lattice.dimension():
            raise ValueError('Invalid input vector sizes.')
        word = w.reduced_word()
        word.reverse()
        u = v
        for i in word:
            u = self.sigma_action(i-1, u)
        return u

    def epsilon_action(self, i, v):
        """
        Action on symbolic variables 'v' switching sign of variables adjacent to 'xi',
        as defined in (3.10) on page 9. (Adjacency is defined on page 6)
         
        Example:
            sage: F.epsilon_action(0,F.v()) #F is based on RootSystem(['A',4])
            [-x0, -x1, x2, x3]

            sage: F.epsilon_action(0,[1,1,1,1])
            [-1, -1, 1, 1]
        """
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
        """
        Act on a polynomial 'f' by the ith generator of the weyl group, i.e. the ith simple reflection,
        as defined in (3.13) on page 9.

        Example:
            sage: a = F.lattice.random_element() #F is based on RootSystem(['A',2])
            sage: a
            -2*alpha[1]
            sage: m = F.x(a)
            sage: m
            x0^(-2)
            sage: F.weyl_action_by_simple_reflection(0,m)
            -(q*x0 - 1)*q*x0/(x0 - 1)
        """
        q = var('q')
        v = self.v()
        A = (-(1-q*v[i]))/((q*v[i])*(1-v[i]))
        B = (1/(v[i]*sqrt(q)))
        return A*self.evaluate(self.f_plus(i,f),self.sigma_action(i,v)) + B*self.evaluate(self.f_minus(i,f),self.sigma_action(i,v))

    def weyl_action(self, w, f):
        """
        Act on a polynomial 'f' by an element of the Weyl group 'w'

        Example:
            sage: a = F.lattice.random_element() #F is based on RootSystem(['A',2])
            sage: a
            -2*alpha[1]
            sage: m = F.x(a)
            sage: m
            x0^(-2)
            sage: w = F.weyl_group.random_element()
            sage: w
            [ 0 -1]
            [-1  0]
            sage: w.reduced_word()
            [1, 2, 1]
            sage: F.weyl_action(w,m).simplify_full()
            -((q^3*x0^3 - q^3*x0^2)*x1^4 - (q^2*x0^3 - q^2*x0^2 + q^2*x0 - q^2)*x1^2 + q*x0 - (((q^3 - q^2)*x0^3 - (q^3 - q^2)*x0^2)*x1^3 - ((q^3 - 2*q^2 + q)*x0^2 - (q^2 - 2*q + 1)*x0)*x1^2 - ((q - 1)*x0 - q + 1)*x1)*sqrt(q) - q)/(((q^2*x0^5 - q^2*x0^4)*x1^4 + q*x0^3 - q*x0^2 - (q^2*x0^5 - q^2*x0^4 + q*x0^3 - q*x0^2)*x1^2)*sqrt(q))
        """
        word = w.reduced_word()
        for i in word:
            f = self.weyl_action_by_simple_reflection(i-1, f)
        return self.evaluate(f,self.v())

    ################################################################
    #Methods related to the construction of the invariant polynomial
    ################################################################

    def delta(self):
        """
        Compute the 'Delta' polynomial defined on page 12.

        Example:
            sage: F.delta() #F is based on RootSystem(['A',2])
            -(q^2*x0^2*x1^2 - 1)*(q*x0^2 - 1)*(q*x1^2 - 1)
        """
        p = 1
        q = var('q')
        for a in self.lattice.positive_roots():
            cfs = a.coefficients()
            p = p * (1 - q^sum(cfs) * self.x(2*a))
        return p

    def j(self, w):
        """
        Compute the 'j' polynomial for a weyl group element 'w' defined
        on page 12 based on Lemma 3.8.

        Example:
            sage: w = F.weyl_group.random_element()
            sage: w
            [ 0 -1]
            [-1  0]
            sage: w.reduced_word()
            [1, 2, 1]
            sage: F.j(w)
            -q^4*x0^4*x1^4
        """
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
        """
        Compute the f0 part of the invariant polynomial as defined in (3.20) on page 12.
        """
        f0 = []
        one = self.one()
        for w in self.weyl_group.list():
            j = self.j(w)
            h = self.weyl_action(w, one)
            g = j * h
            f0.append(g)
        return sum(f0)

    def weyl_invariant_polynomial(self):
        """
        Compute the invariant polynomial as defined in (3.21) on page 12.
        """
        return (self.f0()/self.delta())

    ################################################################
    #Miscellaneous useful methods
    ################################################################

    def x(self, a):
        """
        Generate monomial based off of a member of the root lattice
        as defined on page 8.

        Example:
            sage: a = F.lattice.random_element() #F is based on RootSystem(['A',2])
            sage: a
            -2*alpha[1]
            sage: m = F.x(a)
            sage: m
            x0^(-2)
        """
        v = self.v()
        p = 1
        for i in range(0,len(v)):
            p = p*(v[i])^(a.coefficient(i+1))
        return p

    def phi(self, w):
        """
        Find positive roots sent to negative roots by a weyl group element as
        defined on page 7.

        Example:
            sage: w = F.weyl_group.random_element() #F is based on RootSystem(['A',2])
            sage: w
            [ 0 -1]
            [-1  0]
            sage: w.reduced_word()
            [1, 2, 1]
            sage: F.phi(w)
            [alpha[1], alpha[2], alpha[1] + alpha[2]]
        """
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
        """
        Create the constant symbolic function 1 for the polynomial ring. This is used
        in the method f0 to prevent 'coercion' errors.
        """
        v = self.v()
        one = SR(1)
        return one.function(*v)

    def evaluate(self, f, v):
        """
        Special method to evaluate a given polynomial against a list of inputs [x0,x1,...,xn].
        This method is needed because the symbolic variable 'q' in the paper forces us to do
        most of our computations in the symbolic ring instead of the polynomial ring encapsulated
        in this class. Evaluation of symbolic expressions against inputs is different than evaluation
        in the polynomial ring. This method serves as a common entry point for evaluation and thus
        standardizes it within the context of this class. Outside of this class this method should
        NOT be used.
        """
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