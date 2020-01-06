import numpy as np
from scipy.special import loggamma


def vectorized_np_choice(probabilities):
    cumulative = probabilities.cumsum(axis=1)
    uniform = np.random.rand(len(cumulative), 1)
    samples = (uniform < cumulative).argmax(axis=1)
    return samples


def log_categorical(X, theta):
    return X @ np.log(theta.T)  # Computes for all theta_k at once


class GibbsSampler:
    """Samples pi, theta, and Z for a categorical mixture model
    
    Returns:
        pi -- Mixture fractions
        theta -- categorical probabilities for all K clusters
        Z -- Cluster assignment
    """
    def __init__(self, seed=None):
        np.random.seed(seed)
        self.pi = np.array([])
        self.theta = np.array([])
        self.Z = np.array([])
        self.num_clusters = 0

    def fit(self, X, num_clusters, num_steps=100, alpha=1, gamma=1):
        """ Infers the parameters for a defined number of clusters
        
        Arguments:
            X {np.array} -- Data
            num_clusters {int} -- Number of clusters
        
        Keyword Arguments:
            num_steps {int} -- Maximal number of steps (default: {100})
            num_burn_in_steps {int} -- [description] (default: {20})
            alpha {float} -- Dirichlet Paramter (default: {1})
            gamma {fliat} -- Dirichlet Parameter (default: {1})
        """
        # Initialization
        ll_list = []
        self.num_clusters = num_clusters
        self.Z = np.random.choice(self.num_clusters, len(X))

        for _ in range(num_steps):
            # Sampling
            self._sample_pi(alpha)
            self._sample_theta(X, gamma)
            self._sample_Z(X)
            # Evaluation
            log_likelihood = sum([np.sum(X[self.Z == k] @ np.log(
                self.theta[k])) for k in range(self.num_clusters)])
            ll_list += [log_likelihood]

    def _sample_pi(self, alpha):
        """Samples a vector pi of length num_clusters using a Dirichlet
        distribution and hyperparameter alpha."""
        alpha_prime = alpha + \
            np.bincount(
                self.Z, minlength=self.num_clusters)  # Bincount counts the number of occurences for each value in Z
        self.pi = np.random.dirichlet(alpha_prime)

    def _sample_Z(self, X):
        """Samples the most likely vector Z of length num_datapoints by using
        X, pi and theta."""
        p_Z_not_normalized = np.exp(
            np.log(self.pi) + log_categorical(X, self.theta))
        p_Z_normalizing_const = np.sum(p_Z_not_normalized, axis=1)
        p_Z = p_Z_not_normalized / np.expand_dims(p_Z_normalizing_const, 1)
        self.Z = vectorized_np_choice(p_Z)

    def _sample_theta(self, X, gamma):
        count = np.array([np.sum(X[self.Z == k], axis=0)
                          for k in range(self.num_clusters)])
        gamma_prime = gamma + count
        self.theta = np.array([np.random.dirichlet(gamma_prime[k])
                               for k in range(self.num_clusters)])


class PiCollapsedNonparametricGibbsSampler:
    """Samples parameters for an infinte Categorical Mixture Model by sampling from a Dirichlet Process
    using Gibbs sampling.
    
    Arguments:
        seed {integer} -- Seed for pseudo random generator

    Returns:
        theta {np.array} -- Categorical Probabilities for K clusters.
        Z {np.array} -- Cluster assignment for training data.
        ll_list {list} -- Log liklihood progression during sampling.
    """
    def __init__(self, seed=None):
        np.random.seed(seed)

    def fit(self, X, delta=10, report_status=True, stop_at_conv=True, num_steps=100, num_burn_in_steps=20, alpha=1, gamma=1):
        """ Infers the parameters
        
        Arguments:
            X {np.array} -- Matrix of Data, where each row repesent one data point and each column the number of times
                            a catgory was observed in the data  
        
        Keyword Arguments:
            delta {int} -- If the loglikelihood does not change in a defined interval by delta the algorithm stops (default: {10})
            report_status {bool} -- print status at current iteration (default: {True})
            stop_at_conv {bool} -- stop if converged (default: {True})
            num_steps {int} -- maximal number of steps (default: {100})
            num_burn_in_steps {int} -- minimal number of steps (default: {20})
            alpha {int} -- Dirchichlet Parameter for pi, which can be seen as prior knowledge for the number of clusters (default: {1})
            gamma {int} -- Dirichlet Paramters for theta (default: {1})
        """
        I = len(X[0])
        # Initialization
        self.ll_list = []
        self.K_seen = 1
        self.Z = np.zeros(len(X))
        self.theta = np.zeros((self.K_seen, I))
        self.theta[0] = np.random.dirichlet(I*[gamma])

        # Useful private counts to have:
        self._logpx_n_unseen = np.zeros(len(X))
        for n in range(len(X)):
            first_gamma = np.sum(loggamma(gamma + X[n]))
            second_gamma = loggamma(np.sum([gamma]*I))
            third_gamma = loggamma(np.sum(gamma+X[n]))
            four_gamma = np.sum(loggamma([gamma]*I))
            self._logpx_n_unseen[n] = first_gamma + \
                second_gamma - third_gamma - four_gamma

        self._m_k = np.zeros(self.K_seen)
        for k in range(self.K_seen):
            self._m_k[k] = np.sum(self.Z == k)

        self._gamma_km = np.zeros((self.K_seen, I))
        for k in range(self.K_seen):
            for m in range(I):
                self._gamma_km[k, m] = gamma + np.sum(X[self.Z == k][:, m])

        # Start:
        i = 1
        while i <= num_burn_in_steps or (stop_at_conv and not self._isConverged(delta=delta)) and i <= num_burn_in_steps + num_steps:
            self._sample_Z(X, alpha, gamma)
            self._sample_theta()
            self._print_status(X, i, report_status=report_status)
            i += 1

    def _print_status(self, X, i, report_status=True):
        
        if not report_status:
            return
        text = "Iteration: {:d}; Current clusters: {:d}; Likelihood: {:10.3f}"
        print(text.format(i, self.K_seen, float(self.loglikelihood(X))))

    def loglikelihood(self, X):
        """ Returns the liklihood of the current model
        
        Arguments:
            X {np.array} -- Data
        
        Returns:
            {integer} -- Liklihood
        """
        log_likelihood = sum([np.sum(X[self.Z == k] @ np.log(self.theta[k])) for k in range(self.K_seen)])
        self.ll_list += [log_likelihood]
        return log_likelihood

    def _update_counts(self, new_z, z_old, X_n, gamma):
                    # If we observe a new cluster:
        if new_z >= self.K_seen:
            # Change word counts and K_seen
            self._m_k = np.append(self._m_k, [0])
            self._m_k[new_z] += 1
            self.K_seen += 1

            # Add new theta
            self._gamma_km = np.vstack(
                (self._gamma_km[:self.K_seen-1], X_n + gamma))
            theta_new = np.random.dirichlet(self._gamma_km[new_z])
            self.theta = np.vstack((self.theta[:self.K_seen-1], theta_new))
        # If old cluster change cluster counts to new.
        else:
            self._m_k[new_z] += 1
            self._gamma_km[new_z] += X_n

    def _remove_empty_clusters(self, z_old):
        if self._m_k[z_old] == 0:
            self._m_k = np.delete(self._m_k, z_old)
            self._gamma_km = np.delete(self._gamma_km, z_old, axis=0)
            self.theta = np.delete(self.theta, z_old, axis=0)
            self.K_seen -= 1
            self.Z[self.Z > z_old] -= 1

    def _sample_theta(self):
        for k in range(self.K_seen):
            self.theta[k] = np.random.dirichlet(self._gamma_km[k])

            self.theta[k] = np.random.dirichlet(self._gamma_km[k])

    def _sample_Z(self, X, alpha, gamma):
        """Samples the most likely vector Z of length num_datapoints by using X
        and the current cluster assignments."""
        for n in range(len(X)):
            z_old = int(self.Z[n])
            self._m_k[z_old] -= 1
            self._gamma_km[z_old] -= X[n]

            # For seen clusters:
            p_zn_Zn_seen = self._m_k/(len(X)-1 + alpha)
            logp_xn_theta_seen = X[n] @ np.log(self.theta.T)

            # For unseen clusters:
            p_zn_Zn_unseen = alpha/(len(X)-1+alpha)
            logpx_n_unseen = self._logpx_n_unseen[n]

            p_zn_seen = p_zn_Zn_seen * np.exp(np.array(logp_xn_theta_seen, dtype = np.longfloat))
            p_zn_unseen = p_zn_Zn_unseen * np.exp(np.array(logpx_n_unseen,dtype = np.longfloat))

            p_zn = np.append(p_zn_seen, p_zn_unseen) 
            p_zn /= np.sum(p_zn, dtype=np.longfloat)
            
            new_z = np.random.choice(range(self.K_seen + 1), p=np.array(p_zn,dtype = np.float64))

            # Update counting variables and differentiate new/old cluster
            self._update_counts(new_z, z_old, X[n], gamma)

            # Set new cluster, remove created empty ones.
            self.Z[n] = new_z
            self._remove_empty_clusters(z_old)

    def _isConverged(self, delta=10, interval=10):
        if len(self.ll_list) > 2*interval:
            start = len(self.ll_list) - 2*interval
            ll_before = self.ll_list[start:start+interval]
            ll_now = self.ll_list[start+interval:]
            return delta > np.abs(np.mean(ll_before) - np.mean(ll_now))
        else:
            return False
