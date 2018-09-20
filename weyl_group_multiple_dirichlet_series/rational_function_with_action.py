from sage.all import *
from weyl_group_multiple_dirichlet_series.variable_factory import *

class RationalFunctionWithAction:
    """

    """

    def __init__(self, f, weyl_group):
        self.original_function = f
        self.transformed_function = f
        self.weyl_group_element = weyl_group.random_element_of_length(0)
        self._variable_factory = VariableFactory(weyl_group)
        self._weyl_action = WeylAction(self._variable_factory)
        self.is_transformed = True

    def __str__(self):
        return str(self.transformed_function)

    def __repr__(self):
        return str(self)

    def _transform(self):
        if self.is_transformed == False:
            word = self.weyl_group_element.reduced_word()
            for i in word:
                self.transformed_function = self._weyl_action.simple_action(self, i-1)
            self.is_transformed = True
        return self.transformed_function

    def act(self, w):
        self.is_transformed = False
        self.transformed_function = self.original_function
        self.weyl_group_element = self.weyl_group_element * w
        return self._transform()

    def _evaluate(self, inputs, f):
        if f.parent() == SR:
            xs = self._variable_factory.variables
            for i in range(0, len(xs)):
                x = xs[i]
                f = f.function(x)
                f = f(inputs[i])
            return f
        else:
            try:
                return f(inputs)
            except:
                return f(*inputs)

    def evaluate(self, inputs, evaluate_original_function=False):
        if evaluate_original_function:
            return self._evaluate(inputs, self.original_function)
        else:
            return self._evaluate(inputs, self.transformed_function)

class WeylAction:
    """

    """

    def __init__(self, variable_factory):
        self._variable_factory = variable_factory

    def c(self, i):
        q = self._variable_factory.q
        xi = self._variable_factory.variables[i]
        A = (q*xi - 1)/(q*xi*(1 - xi))
        B = 1/(sqrt(q)*xi)
        c = (A + B)/2
        return c

    def d(self, i):
        q = self._variable_factory.q
        xi = self._variable_factory.variables[i]
        A = (q*xi - 1)/(q*xi*(1 - xi))
        B = 1/(sqrt(q)*xi)
        d = (A - B)/2
        return d

    def simple_action(self, f, i):
        vars1 = self._variable_factory.act_on_variables("s"+str(i))
        vars2 = self._variable_factory.act_on_variables("s"+str(i)+"e"+str(i))
        action = self.c(i) * f.evaluate(vars1) + self.d(i) * f.evaluate(vars2)
        return action