# from numpy import array, log, exp
# from scipy.special import gammaln

# def log_factorial(x):
#     """Returns the logarithm of x!
#     Also accepts lists and NumPy arrays in place of x."""
#     return gammaln(array(x)+1)

# def multinomial_(xs, ps):
#     n = sum(xs)
#     xs, ps = array(xs), array(ps)
#     result = log_factorial(n) - sum(log_factorial(xs)) + sum(xs * log(ps))
#     return exp(result)