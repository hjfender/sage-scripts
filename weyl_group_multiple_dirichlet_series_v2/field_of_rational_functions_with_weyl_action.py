from sage.all import *
from sage.rings.fraction_field_element import make_element

class FieldOfRationalFunctionsWithWeylAction(sage.rings.fraction_field.FractionField_generic):
    """
    """

    ###############################################################
    #Class related methods
    ###############################################################
    def __init__(self, root_system):
        self.root_system = root_system
        self.lattice = root_system.root_lattice()
        self.weyl_group = self.lattice.weyl_group()
        #construct the function field
        variables = ['sqrtq']
        for i in range(1,self.lattice.dimension() + 1):
            variables.append('x'+str(i))
        sage.rings.fraction_field.FractionField_generic.__init__(self, PolynomialRing(CC,variables))

    def __str__(self):
        return (super(FieldOfRationalFunctionsWithWeylAction, self).__str__() \
                + " endowed with an action by " + str(self.weyl_group))

    def __repr__(self):
        return (super(FieldOfRationalFunctionsWithWeylAction, self).__repr__() \
                + " endowed with an action by " + str(self.weyl_group))

    ################################################################
    #Methods related to the construction of the invariant polynomial
    ################################################################
    def delta(self):
        """
        """
        p = self(1)
        q = self('sqrtq')**2
        for a in self.lattice.positive_roots():
            cfs = a.coefficients()
            p = p * (1 - q ** sum(cfs) * self.x(2*a))
        return p

    def j(self, w):
        """
        """
        rw = w.reduced_word()
        variables = self.variable_names()
        j = 0
        q = self('sqrtq')**2
        if len(rw) == 1:
            j = -q*self(variables[rw[0]])**2
        else:
            a = sum(self.phi(w))
            cfs = a.coefficients()
            j = (-1)**len(rw) * q**(sum(cfs)) * self.x(2*a)
        return j

    def f0(self):
        """
        """
        f0 = self(0)
        for w in self.weyl_group.list():
            g = self.j(w) * self.act(self.one(),w)
            g.reduce()
            f0 += g
        # f0.reduce()
        return f0

    def invariant_function(self):
        """
        """
        function = self.f0()/self.delta()
        # function.reduce()
        return function


    ################################################################
    #Miscellaneous Rational Functions
    ################################################################
    def x(self, a):
        """
        """
        variables = self.variable_names()
        p = 1
        for i in range(1,len(variables)):
            p = p * self(variables[i])**(a.coefficient(i))
        return p

    #################################################################
    #Miscellaneous Functions Related to Weyl Group
    ################################################################
    def phi(self, w):
        """
        """
        w_roots = []
        for a in self.lattice.positive_roots():
            wa = w.action(a)
            for b in self.lattice.negative_roots():
                if wa == b:
                    w_roots.append(a)
                    break
        if len(w_roots) == 0:
            #we have to do the following to ensure that the sum of the
            #roots, in the method j, is a root lattice element
            w_roots.append(0*self.lattice.gens()[0]) 
        return w_roots

    ###############################################################
    #Methods related to the Weyl action on variables
    ###############################################################
    def _sigma(self, i, v):
        """
        """
        w = [v[0]]
        sr = self.weyl_group.simple_reflections()
        for j in range(1, len(v)):
            if i == j:
                w.append(1/(self('sqrtq')**2*v[j]))
            elif (sr[i] * sr[j])**3 == self.weyl_group.random_element_of_length(0):
                w.append(v[j]*v[i]*self('sqrtq'))
            else:
                w.append(v[j])
        return w

    def _epsilon(self, i, v):
        """
        """
        w = [v[0]]
        sr = self.weyl_group.simple_reflections()
        for j in range(1, len(v)):
            if i != j and (sr[i] * sr[j])**3 == self.weyl_group.random_element_of_length(0):
                w.append(-v[j])
            else:
                w.append(v[j])
        return w

    def act_on_variables(self,pattern_string):
        """
        """
        v = [self(x) for x in self.variable_names()]
        for i in range(0,len(pattern_string),2):
            if pattern_string[i] == 's':
                n = int(pattern_string[i+1])
                v = self._sigma(n,v)
            elif pattern_string[i] == 'e':
                n = int(pattern_string[i+1])
                v = self._epsilon(n,v)
            else:
                raise ValueError("Incorrect syntax in pattern string!")
        for x in v:
            x.reduce()
        return v

    ###############################################################
    #Methods related to the Weyl action on functions
    ###############################################################
    def _c(self, i):
        """
        Compute the 'c' polynomial defined on page 9.
        """
        q = self('sqrtq')**2
        xi = self('x'+str(i))
        A = (q*xi - 1)/(q*xi*(1 - xi))
        B = 1/(self('sqrtq')*xi)
        c = (A + B)/2
        return c

    def _d(self, i):
        """
        Compute the 'd' polynomial defined on page 9.
        """
        q = self('sqrtq')**2
        xi = self('x'+str(i))
        A = (q*xi - 1)/(q*xi*(1 - xi))
        B = 1/(self('sqrtq')*xi)
        d = (A - B)/2
        return d

    def simple_action(self, f, i):
        """
        """
        vars1 = self.act_on_variables("s"+str(i))
        vars2 = self.act_on_variables("s"+str(i)+"e"+str(i))
        action = self._c(i) * f(vars1) + self._d(i) * f(vars2)
        return action

    def act(self, f, w):
        """
        """
        word = w.reduced_word()
        for i in word:
            f = self.simple_action(f, i)
        return f

    ###############################################################
    #Methods to handle simplification of rational functions
    ###############################################################
    def reduce_function(self, function):
        """
        """
        g = self._gcd_numerator_denominator(function.numerator(), function.denominator())
        if not g.is_unit():
            num, _ = function.numerator().quo_rem(g)
            den, _ = function.denominator().quo_rem(g)
        else:
            num = function.numerator()
            den = function.denominator()

        if not den.is_one() and den.is_unit():
            num *= den.inverse_of_unit()
            den  = den.parent().one()
        
        return make_element(self, num, den)

    def _gcd_numerator_denominator(self, num, den):
        """
        """
        x = num.parent().gens()[-1]
        uninumpoly = num.polynomial(x)
        uninumbase = uninumpoly.base_ring()
        return self(self._gcd_polynomials(uninumpoly,den.polynomial(x)))

    def _gcd_polynomials(self, f, g):
        """
        """
        if f.degree() < g.degree():
            A, B = g, f
        else:
            A, B = f, g

        if B.is_zero():
            return A

        a = b = f.parent().zero()
        for c in A.coefficients():
            a = a.gcd(c)
            if a.is_one():
                break
        for c in B.coefficients():
            b = b.gcd(c)
            if b.is_one():
                break

        d = a.gcd(b)
        print ("A; %s \n  B: %s \n a: %s \n b: %s \n d: %s" % (A,B,a,b,d))
        A = A // a
        B = B // b
        g = h = 1

        delta = A.degree() - B.degree()
        _, R = A.pseudo_quo_rem(B)

        while R.degree() > 0:
            A = B
            print ("A: %s \n R: %s \n g: %s \n h: %s \n delta: %s" % (A, R, g, h, delta))
            B = R // (g * h ** delta)
            print "B: " + str(B)
            g = A.leading_coefficient()
            h = h * g ** delta // h ** delta
            delta = A.degree() - B.degree()
            _, R = A.pseudo_quo_rem(B)

        if R.is_zero():
            b = f.parent().zero()
            for c in B.coefficients():
                b = b.gcd(c)
                if b.is_one():
                    break
            return d * B // b

        return f.parent()(d)

    # def _pseudo_quo_rem(self,polynomial,other):
    # """
    # """
    # if other.is_zero():
    #     raise ZeroDivisionError("Pseudo-division by zero is not possible")

    # # if other is a constant, then R = 0 and Q = self * other^(deg(self))
    # if other in self.base_ring():
    #     return (self * other**(A.degree()), self.zero())

    # R = polynomial
    # B = other
    # Q = self.zero()
    # e = polynomial.degree() - other.degree() + 1
    # d = B.leading_coefficient()

    # while not R.degree() < B.degree():
    #     c = R.leading_coefficient()
    #     diffdeg = R.degree() - B.degree()
    #     Q = d*Q + self(c).shift(diffdeg)
    #     R = d*R - c*B.shift(diffdeg)
    #     e -= 1

    # q = d**e
    # return (q*Q,q*R)