class VariableFactory:
    """

    """
    
    def __init__(self, weyl_group):
        self._weyl_group = weyl_group

        #variable involved in the weyl action
        self.q = var('q')

        #list polynomial ring variables
        self.variables = []
        for i in range(0,len(self._weyl_group.gens())):
            self.variables.append(var('x'+str(i)))

        #the following are for amortization purposes
        self._sigma = {}
        self._epsilon = {}

    def sigma(self, i):
        w = self._sigma.get(i, [])
        if len(w) == 0:
            sr = self._weyl_group.simple_reflections()
            v = self.variables
            for j in range(0, len(v)):
                if i == j:
                    w.append(1/(self.q*v[j]))
                elif (sr[i+1]*sr[j+1])^3 == self._weyl_group.random_element_of_length(0):
                    w.append(v[j]*v[i]*sqrt(self.q))
                else:
                    w.append(v[j])
            self._sigma[i] = w
        return w

    def epsilon(self, i):
        w = self._epsilon.get(i, [])
        if len(w) == 0:
            sr = self._weyl_group.simple_reflections()
            v = self.variables
            for j in range(0, len(v)):
                if (sr[i+1] * sr[j+1])^3 == self._weyl_group.random_element_of_length(0):
                    w.append(-v[j])
                else:
                    w.append(v[j])
            self._epsilon[i] = w
        return w

        
