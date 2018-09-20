from sage.all import *

class VariableFactory:
    """

    """
    
    def __init__(self, weyl_group):
        self._weyl_group = weyl_group

        #variable involved in the weyl action
        self.q = SR.var('q')

        #list polynomial ring variables
        self.variables = []
        for i in range(0,len(self._weyl_group.gens())):
            self.variables.append(SR.var('x'+str(i)))

    #pattern string acts left to right!!!!
    def act_on_variables(self,pattern_string):
        v = self.variables
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
        w = []
        sr = self._weyl_group.simple_reflections()
        for j in range(0, len(v)):
            if (sr[i+1] * sr[j+1])**3 == self._weyl_group.random_element_of_length(0):
                w.append(-v[j])
            else:
                w.append(v[j])
        return w