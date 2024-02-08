class num_diff:

    def derivative(g, a, method='central', p=0.01):
        '''Compute the difference formula for f'(a) with step size h.

        Parameters
        ----------
        f : function
             Vectorized function of one variable
        a : number
             Compute derivative at x = a
        method : string
             Difference formula: 'forward', 'backward' or 'central'
        h : number
             Step size in difference formula

        Returns
        -------
        float
             Difference formula:
                  central: f(a+h) - f(a-h))/2h
                  forward: f(a+h) - f(a))/h
                  backward: f(a) - f(a-h))/h            
        '''
        if method == 'central':
            return (g(a + p) - g(a - p))/(2*p)
        elif method == 'forward':
            return (g(a + p) - g(a))/p
        elif method == 'backward':
            return (g(a) - g(a - p))/p
        else:
            raise ValueError(
                "Method must be 'central', 'forward' or 'backward'.")

    def derivative2(g, a, method='central', p=0.01):  # derivada de segunda ordem
        '''Compute the difference formula for f'(a) with step size h.

        Parameters
        ----------
        f : function
             Vectorized function of one variable
        a : number
             Compute derivative at x = a
        method : string
             Difference formula: 'forward', 'backward' or 'central'
        h : number
             Step size in difference formula

        Returns
        -------
        float
             Difference formula:
                  central: f(a+h) - f(a-h))/2h
                  forward: f(a+h) - f(a))/h
                  backward: f(a) - f(a-h))/h            
        '''
        if method == 'central':
            return (g(a + p) - 2*g(a)+g(a - p))/(p**2)

        elif method == 'forward':
            return (g(a + 2*p) - 2*g(a+p) + g(a))/p

        elif method == 'backward':
            return (g(a) - 2*g(a - p)+g(a-2*p))/p
        else:
            raise ValueError(
                "Method must be 'central', 'forward' or 'backward'.")