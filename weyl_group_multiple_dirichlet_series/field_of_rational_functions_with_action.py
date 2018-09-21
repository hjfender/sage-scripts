from sage.all import *
from weyl_group_multiple_dirichlet_series.rational_function_with_action import *
from weyl_group_multiple_dirichlet_series.variable_helper import *

class FieldOfRationalFunctionsWithAction:
    """
    A class encapsulating the relevant entities used to find a special invariant
    function used to construct a Multiple Dirichlet Series. Those entities are:
        -A root system (passed in through the constructor)
        -A root lattice (corresponding to the root system)
        -A weyl group (corresponding to the root lattice)

    Example:
        sage: RS = RootSystem(['A',2])
        sage: F = FieldOfRationalFunctionsWithAction(RS)
        sage: F
        Field of rational functions in 2 variables endowed with an action by Weyl Group
        of type ['A', 2] (as a matrix group acting on the root lattice)
        sage: F.root_system
        Root system of type ['A', 2]
        sage: F.lattice
        Root lattice of the Root system of type ['A', 2]
        sage: F.weyl_group
        Weyl Group of type ['A', 2] (as a matrix group acting on the root lattice)
    """

    ###############################################################
    #Class related methods
    ###############################################################
    def __init__(self, root_system):
        self.root_system = root_system
        self.lattice = root_system.root_lattice()
        self.weyl_group = self.lattice.weyl_group()
        self._variable_helper = VariableHelper(self.weyl_group)

    def __str__(self):
        return ("Field of rational functions in " + str(self.lattice.dimension()) \
                + " variables endowed with an action by " + str(self.weyl_group))

    def __repr__(self):
        return str(self)

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
            p = p * (1 - q**sum(cfs) * self.x(2*a).original_function)
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
        v = self._variable_helper.variables
        q = self._variable_helper.q
        j = 0
        if len(rw) == 1:
            j = -q*v[rw[0]-1]**2
        else:
            a = sum(self.phi(w))
            cfs = a.coefficients()
            j = (-1)**len(rw) * q**(sum(cfs)) * self.x(2*a).original_function
        return j

    def f0(self):
        """
        Compute the f0 part of the invariant polynomial as defined in (3.20) on page 12.
        """
        f0 = []
        for w in self.weyl_group.list():
            one = self.one()
            g = self.j(w) * one.act(w)
            f0.append(g)
        return sum(f0)

    def invariant_function(self):
        """
        Compute the invariant polynomial as defined in (3.21) on page 12.
        """
        return RationalFunctionWithAction((self.f0()/self.delta()), self.weyl_group)

    ################################################################
    #Miscellaneous Rational Functions
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
        v = self._variable_helper.variables
        p = 1
        for i in range(0,len(v)):
            p = p*(v[i])**(a.coefficient(i+1))
        return RationalFunctionWithAction(p, self.weyl_group)

    def one(self):
        """
        Create the constant symbolic function 1 for the field. This is used
        in the method f0 to prevent 'coercion' errors.
        """
        v = self._variable_helper.variables
        one = SR(1)
        return RationalFunctionWithAction(one.function(*v), self.weyl_group)

    ################################################################
    #Miscellaneous Functions Related to Weyl Group
    ################################################################
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
            #we have to do the following to ensure that the sum of the
            #roots, in the method j, is a root lattice element
            w_roots.append(0*self.lattice.gens()[0]) 
        return w_roots

    