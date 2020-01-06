import numpy as np

class CategoricalMM:
    # Final Model class... TODO implement
    def __init__(self, pi, theta):
        self.pi = np.array([])
        self.theta = np.array([])

    def p(self, x):
        """ Return probability """
        logpx = x @ np.log(self.theta.T)
        p =  np.sum(self.pi * np.exp(logpx))
        return p

    def simple_sample(self, draws, replace = True):
        k = np.random.choice(range(len(self.pi)), p = self.pi)
        theta = self.theta[k]

        return np.random.choice(range(len(theta)), p = theta, size = draws, replace = replace)

    @staticmethod
    def ML_est_pi(Z):
        K = np.max(Z)+1
        pi = np.array([np.sum(Z == k) for k in range(K)])
        pi /= pi.sum()

        return pi
    

class AlleleMM(CategoricalMM):
    pass