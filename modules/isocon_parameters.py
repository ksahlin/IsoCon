from __future__ import print_function
class Parameters(object):
    """docstring for Parameters
        This class is uset to create one object containing
        all parameters and other arguments that will be
        used, and passed across all modules.
        It inherits from argument parser object
        (all the pre specified arguments), but adds on the
        parameters that needs to be inferred at runtime."""
    def __init__(self, **kwargs):
        super(Parameters, self).__init__()
        # kwargs is a dict of the keyword args passed to the function
        # we get then from the argparse object.
        # The reason for having this class is to add additional
        # parameters that might need to be inferred during
        # runtime.
        for key, value in kwargs.items():
            print("{0}: {1}".format(key, value))
            setattr(self, key, value)
