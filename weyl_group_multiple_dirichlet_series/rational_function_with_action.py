from sage.all import *
from weyl_group_multiple_dirichlet_series.variable_helper import *

class RationalFunctionWithAction:
    """
    A class encapulating a rational function and handling the weyl action on it. A instance
    is created with a base rational function which is maintained as 'original_function' and
    a weyl group which allows for group calculations. A 'weyl_group_element' is set to the identity.
    This variable changes as the function as acted upon and the actual action on the polynomial is
    maintained as 'transformed_function'.
    
    This approach allows a drastic simplification of weyl action computations for functions that have already
    had actions performed on them. Instead of computing the action of 'w' directly on the function each time,
    we first determine the interaction of the previous actions with the current one by multiplying self.weyl_group_element
    on the right by 'w'. We are then free to act with the resulting element on the 'original_function' and update
    'transformed_function.' This amounts to travelling around the top and right sides of the following commutative
    diagram instead of the left and bottom sides...

    F x W x W ----mult----> F x W
        |                     |
       act                   act
        |                     |
        v                     v
      F x W ------act------>  F

    W - weyl group, F - field of rational functions

    Note: This class uses the supporting class below.
    """

    ###############################################################
    #Class related methods
    ###############################################################
    def __init__(self, f, weyl_group):
        self.original_function = f
        self.transformed_function = f
        self.weyl_group_element = weyl_group.random_element_of_length(0)
        self._variable_helper = VariableHelper(weyl_group)
        self._weyl_action = WeylAction(self._variable_helper)
        self.is_transformed = True

    def __str__(self):
        return str(self.transformed_function)

    def __repr__(self):
        return str(self)

    ###############################################################
    #Action related methods
    ###############################################################
    def _transform(self):
        """
        Handle the weyl action as illustrated in the class comment in conjunction with
        the method below.
        """
        if self.is_transformed == False:
            word = self.weyl_group_element.reduced_word()
            for i in word:
                self.transformed_function = self._weyl_action.simple_action(self, i-1)
            self.transformed_function = self.transformed_function.simplify_full()
            self.is_transformed = True
        return self.transformed_function

    def act(self, w):
        """
        Handle the weyl action as illustrated in the class comment in conjunction with
        the method above.

        Example:
            sage: w = W.random_element()
            
            sage: w
            [ 0 -1]
            [ 1 -1]

            sage: f
            x1^2/x0^2

            sage: f.act(w)
            -1/4*(q^3*x0^4*x1^2*(1/(q*x0*x1) + (q^(3/2)*x0*x1 + 1)/((sqrt(q)*x0*x1 + 1)*q^(3/2)*x0*x1))
            - q^3*x0^4*x1^2*(1/(q*x0*x1) - (q^(3/2)*x0*x1 + 1)/((sqrt(q)*x0*x1 + 1)*q^(3/2)*x0*x1)))*(1/(sqrt(q)*x1)
            + (q*x1 - 1)/(q*(x1 - 1)*x1)) - 1/4*(q^3*x0^4*x1^2*(1/(q*x0*x1) + (q^(3/2)*x0*x1 - 1)/((sqrt(q)*x0*x1 - 1)*q^(3/2)*x0*x1))
            - q^3*x0^4*x1^2*(1/(q*x0*x1) - (q^(3/2)*x0*x1 - 1)/((sqrt(q)*x0*x1 - 1)*q^(3/2)*x0*x1)))*(1/(sqrt(q)*x1)
            - (q*x1 - 1)/(q*(x1 - 1)*x1))
        """
        self.is_transformed = False
        self.transformed_function = self.original_function
        self.weyl_group_element = self.weyl_group_element * w
        return self._transform()

    ###############################################################
    #Function related methods
    ###############################################################
    def evaluate(self, inputs, evaluate_original_function=False):
        """
        Using the evaluate function provided in the variable helper evaluate the function
        against a set of inputs. There are two functions one can evaluate against, the 
        original and the transformed (after action has been applied) function, which one is
        specified by the optional parameter.

        Example:
            sage: f.original_function
            1/(x0^4*x1^12)

            sage: f.transformed_function
            -1/4*(q^6*x1^4*(1/(q*x0*x1) + (q^(3/2)*x0*x1 + 1)/((sqrt(q)*x0*x1 + 1)*q^(3/2)*x0*x1))/x0^8 - q^6*x1^4*(1/(q*x0*x1)
            - (q^(3/2)*x0*x1 + 1)/((sqrt(q)*x0*x1 + 1)*q^(3/2)*x0*x1))/x0^8)*(1/(sqrt(q)*x1) + (q*x1 - 1)/(q*(x1 - 1)*x1))
            - 1/4*(q^6*x1^4*(1/(q*x0*x1) + (q^(3/2)*x0*x1 - 1)/((sqrt(q)*x0*x1 - 1)*q^(3/2)*x0*x1))/x0^8
            - q^6*x1^4*(1/(q*x0*x1) - (q^(3/2)*x0*x1 - 1)/((sqrt(q)*x0*x1 - 1)*q^(3/2)*x0*x1))/x0^8)*(1/(sqrt(q)*x1) - (q*x1 - 1)/(q*(x1 - 1)*x1))

            sage: f.evaluate([2,2])
            -1/512*(q^6*(1/q + (4*q^(3/2) + 1)/(q^(3/2)*(4*sqrt(q) + 1))) - q^6*(1/q - (4*q^(3/2) + 1)/(q^(3/2)*(4*sqrt(q) + 1))))*((2*q - 1)/q + 1/sqrt(q))
            + 1/512*(q^6*(1/q + (4*q^(3/2) - 1)/(q^(3/2)*(4*sqrt(q) - 1))) - q^6*(1/q - (4*q^(3/2) - 1)/(q^(3/2)*(4*sqrt(q) - 1))))*((2*q - 1)/q - 1/sqrt(q))

            sage: f.evaluate([2,2], True)
            1/65536
        """
        if evaluate_original_function:
            return self._variable_helper._evaluate(inputs, self.original_function)
        else:
            return self._variable_helper._evaluate(inputs, self.transformed_function)

class WeylAction:
    """
    A class class encapsulating the components necessary to perform the simple action
    on a rational function as defined on page 9.

    Note: This class is used in a supporting manner in the class above. It shouldn't be used
    independently.
    """

    ###############################################################
    #Class related methods
    ###############################################################
    def __init__(self, variable_helper):
        self._variable_helper = variable_helper

    ###############################################################
    #Action related methods
    ###############################################################
    def _c(self, i):
        """
        Compute the 'c' polynomial defined on page 9.
        """
        q = self._variable_helper.q
        xi = self._variable_helper.variables[i]
        A = (q*xi - 1)/(q*xi*(1 - xi))
        B = 1/(sqrt(q)*xi)
        c = (A + B)/2
        return c

    def _d(self, i):
        """
        Compute the 'd' polynomial defined on page 9.
        """
        q = self._variable_helper.q
        xi = self._variable_helper.variables[i]
        A = (q*xi - 1)/(q*xi*(1 - xi))
        B = 1/(sqrt(q)*xi)
        d = (A - B)/2
        return d

    def simple_action(self, f, i):
        """
        Compute the simple action on a rational function as defined in (3.14) on page 9.

        Note: the 'f' argument must be an instance of RationalFunctionWithAction because this
        method makes use of the evaluate method defined in that class.

        Example:
            sage: WA = WeylAction(V) #V is a variable helper based off of ['A', 2]
            sage: f
            1/(x0^4*x1^12)
            sage: WA.simple_action(f,0)
            -1/2*(1/(sqrt(q)*x0) + (q*x0 - 1)/(q*(x0 - 1)*x0))/(q^2*x0^8*x1^12) + 1/2*(1/(sqrt(q)*x0) - (q*x0 - 1)/(q*(x0 - 1)*x0))/(q^2*x0^8*x1^12)
        """
        vars1 = self._variable_helper.act_on_variables("s"+str(i))
        vars2 = self._variable_helper.act_on_variables("s"+str(i)+"e"+str(i))
        action = self._c(i) * f.evaluate(vars1) + self._d(i) * f.evaluate(vars2)
        return action