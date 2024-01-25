import numpy as np

polymer_len = 7

def log_f(x, polymer_length):
    return -0.38*(np.sin(0.5*np.pi*x/polymer_length))**2

def get_log_expected_singular_values(log_vals, n, show_ans_tups=False):
    """ get the expected singular values where arbitrary numbers
        of quanta can be put into a given mode
        log_vals: logarithms of the value for 1 quanta in a given mode.
                  Should be sorted in increasing order!
        n: we'll compute the largest (least negative) n singular values """
    def H(tup): # energy function
        return sum([tup[i]*log_vals[i] for i in range(dim)])
    dim = len(log_vals)
    ans_tups = [tuple([0]*dim)]
    for i in range(n):
        # get a list of candidates for the next largest tuple
        # a candidate must not be an existing tuple
        # all candidates are increments of some existing tuple
        candidate_tups = []
        for tup in ans_tups:
            for j in range(dim):
                lst = list(tup)
                lst[j] += 1
                inc_tup = tuple(lst)
                if inc_tup not in ans_tups:
                    candidate_tups.append(inc_tup)
                    break
        # find and add to the list the largest of the candidates
        i_best = 0
        H_best = -np.inf
        for i, tup in enumerate(candidate_tups):
            if H(tup) > H_best:
                i_best = i
                H_best = H(tup)
        ans_tups.append(candidate_tups[i_best])
    if show_ans_tups:
        for tup in ans_tups: print(tup)
    return [H(tup) for tup in ans_tups]


log_lin_S = log_f(1 + np.arange(polymer_len-1), polymer_len)
get_log_expected_singular_values(log_lin_S, 18, True)


