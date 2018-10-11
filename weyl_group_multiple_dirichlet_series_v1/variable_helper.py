from sage.all import *

class VariableHelper:
    """
    A class class to handle the creation and manipulation of the field of fraction variables.
    The weyl group of the root lattice plays an important role in the manipulation of variables
    and must be passed in the constructor

    Example:
        sage: RS = RootSystem(['A',2])
        sage: W = RS.root_lattice().weyl_group()
        sage: V = VariableHelper(W)
        sage: V.q
        q
        sage: V.variables
        [x0, x1]
    """
    
    ###############################################################
    #Class related methods
    ###############################################################
    def __init__(self, weyl_group):
        self._weyl_group = weyl_group

        #variable involved in the weyl action
        self.q = SR.var('q')

        #create polynomial ring variables
        self.variables = []
        for i in range(0,len(self._weyl_group.gens())):
            self.variables.append(SR.var('x'+str(i)))

    ###############################################################
    #Variable related methods
    ###############################################################
    def act_on_variables(self,pattern_string):
        """
        A helper method to ease composition of sequences of the following two methods.

        Note: The pattern string acts left to right! This is the opposite of what is
        presented in the paper. So given 'si' is the ith simple reflection, 'ej' is the
        jth epsilon action, and 'x' is the vector of variables, the expression "si ej x" as
        appearing in (3.11) on page 9 would be written in the code as "V.act_on_variables("ejsi")."

        Example:
            sage: V.act_on_variables("e0s1") #V is based on RootSystem(['A',2])
            [x0/(sqrt(q)*x1), -1/(q*x1)]

            sage: V.act_on_variables("s1e0")
            [-x0/(sqrt(q)*x1), -1/(q*x1)]

            sage: V.act_on_variables("s1e1s1")
            [x0, -x1]
        """
        v = list(self.variables)
        for i in range(0,len(pattern_string),2):
            if pattern_string[i] == 's':
                n = int(pattern_string[i+1])
                v = self._sigma(n,v)
            elif pattern_string[i] == 'e':
                n = int(pattern_string[i+1])
                v = self._epsilon(n,v)
            else:
                raise ValueError("Incorrect syntax in pattern string!")
        return v

    def _sigma(self, i, v):
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
        w = []
        sr = self._weyl_group.simple_reflections()
        for j in range(0, len(v)):
            if i == j:
                w.append(1/(self.q*v[j]))
            elif (sr[i+1] * sr[j+1])**3 == self._weyl_group.random_element_of_length(0):
                w.append(v[j]*v[i]*sqrt(self.q))
            else:
                w.append(v[j])
        return w

    def _epsilon(self, i, v):
        """
        Action on symbolic variables 'v' switching sign of variables adjacent to 'xi',
        as defined in (3.10) on page 9. (Adjacency is defined on page 6)
         
        Example:
            sage: F.epsilon_action(0,F.v()) #F is based on RootSystem(['A',4])
            [-x0, -x1, x2, x3]

            sage: F.epsilon_action(0,[1,1,1,1])
            [-1, -1, 1, 1]
        """
        w = []
        sr = self._weyl_group.simple_reflections()
        for j in range(0, len(v)):
            if i != j and (sr[i+1] * sr[j+1])**3 == self._weyl_group.random_element_of_length(0):
                w.append(-v[j])
            else:
                w.append(v[j])
        return w

    ###############################################################
    #Function related methods
    ###############################################################
    def _evaluate(self, inputs, f):
        """
        Special method to evaluate a given rational function against a list of inputs [x0,x1,...,xn].
        This method is needed because the symbolic variable 'q' in the paper forces us to do
        most of our computations in the symbolic rings. Evaluation of symbolic expressions against
        inputs is different than evaluation in the polynomial ring. This method serves as a common
        entry point for evaluation and thus standardizes it for the project.
        """
        if f.parent() == SR:
            xs = self.variables
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