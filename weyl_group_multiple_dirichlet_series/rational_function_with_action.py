from sage.all import *

class RationalFunctionWithAction:
    """

    """

    def __init__(self, f, w, vrs):
        self.function = f
        self.weyl_group_element = w

    def action(self):
        return f

    def evaluate(self, inputs):
        if self.f.parent() == SR:
            f = self.f
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