"""Utilities"""
import numpy as np

def sum_log_probarray(logprobs):
    assert all([u <= 0 for u in logprobs]), logprobs
    max_val, max_idx = max(logprobs), np.argmax(logprobs)
    arrsum = np.sum(np.exp(np.array(logprobs) - max_val))
    result = max_val + np.log(arrsum)
    return min(0, max_val + np.log(arrsum))


def sum_log_prob(*args):
    return sum_log_probarray(args)
