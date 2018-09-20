from sage.all import *
from weyl_group_multiple_dirichlet_series.rational_function_with_action import *
from weyl_group_multiple_dirichlet_series.variable_factory import *

class FieldOfRationalFunctionsWithAction:
    """

    """

    def __init__(self, root_system):
        self.root_system = root_system
        self.lattice = root_system.root_lattice()
        self.weyl_group = self.lattice.weyl_group()
        self._variable_factory = VariableFactory(self.weyl_group)

    def __str__(self):
        return ("Field of rational functions in " + str(self.lattice.dimension()) \
                + " variables endowed with an action by " + str(self.weyl_group))

    def __repr__(self):
        return str(self)

    def delta(self):
        p = 1
        q = var('q')
        for a in self.lattice.positive_roots():
            cfs = a.coefficients()
            p = p * (1 - q**sum(cfs) * self.x(2*a))
        return p

    def j(self, w):
        rw = w.reduced_word()
        v = self._variable_factory.variables
        q = self._variable_factory.q
        if len(rw) == 1:
            return -q*v[rw[0]-1]**2
        else:
            a = sum(self.phi(w))
            cfs = a.coefficients()
            return (-1)**len(rw) * q**(sum(cfs)) * self.x(2*a)

    def f0(self):
        f0 = []
        one = self.one()
        for w in self.weyl_group.list():
            g = self.j(w) * one.act(w)
            f0.append(g)
        return sum(f0)

    def invariant_function(self):
        return RationalFunctionWithAction(self.f0()/self.delta(), self.weyl_group)

    def x(self, a):
        v = self._variable_factory.variables
        p = 1
        for i in range(0,len(v)):
            p = p*(v[i])**(a.coefficient(i+1))
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
        v = self._variable_factory.variables
        one = SR(1)
        return RationalFunctionWithAction(one.function(*v), self.weyl_group)