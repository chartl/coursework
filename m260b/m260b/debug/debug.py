"""\
Package holding the debug decorator and its config namespace

"""
from itertools import islice
from pprint import pprint

DEBUG=False

def _sub(arg, limit=10):
    """Subset the really large arguments so they don't destroy stdout"""
    if isinstance(arg, list):
        return [_sub(q, limit) for q in arg[:limit]]
    if isinstance(arg, dict):
        return dict({k: _sub(v, limit) for k, v in islice(arg.iteritems(), limit)})
    if isinstance(arg, tuple):
        return tuple([_sub(etry, limit) for etry in arg])
    if isinstance(arg, basestring):
        return arg[:2*limit]
    return arg


def debug(fxn):
    def wrap(*args, **kwargs):
        if not DEBUG:
            return fxn(*args, **kwargs)
        ok_args = [_sub(arg) for arg in args]
        print('Entered function {} with args {}'.format(fxn, ok_args))
        pprint({k: _sub(v) for k, v in kwargs.iteritems()})
        result = fxn(*args, **kwargs)
        print('Exited function {} with result {}'.format(fxn, _sub(result, limit=50)))
        return result
    return wrap

