import numpy as np
from scipy.stats import loggamma

def vectorized_np_choice(probabilities):
    cumulative = probabilities.cumsum(axis=1)
    uniform = np.random.rand(len(cumulative), 1)
    samples = (uniform < cumulative).argmax(axis=1)
    return samples


def log_categorical(X, theta):
    return X @ np.log(theta.T) # Computes for all theta_k at once

    
    
class GibbsSampler:
    def __init__(self, seed=None):
        np.random.seed(seed)
        self.pi = np.array([])
        self.theta = np.array([])
        self.Z = np.array([])
        self.num_clusters = 0
    
    def fit(self, X, num_clusters, num_steps, num_burn_in_steps, alpha=1, gamma=1):
        # Initialization
        log_likelihood_progression = []
        self.num_clusters = num_clusters
        self.Z = np.random.choice(self.num_clusters, len(X))

        for _ in range(num_steps):
            # Sampling
            self._sample_pi(alpha)
            self._sample_theta(X, gamma)
            self._sample_Z(X)
            # Evaluation
            log_likelihood = sum([np.sum(X[self.Z==k] @ np.log(self.theta[k])) for k in range(self.num_clusters)])
            log_likelihood_progression += [log_likelihood]
        return log_likelihood_progression

    def _sample_pi(self, alpha):
        """Samples a vector pi of length num_clusters using a Dirichlet distribution and hyperparameter alpha."""
        alpha_prime = alpha + np.bincount(self.Z, minlength=self.num_clusters) # Bincount counts the number of occurences for each value in Z
        self.pi = np.random.dirichlet(alpha_prime)

    def _sample_Z(self, X):
        """Samples the most likely vector Z of length num_datapoints by using X, pi and theta"""
        p_Z_not_normalized = np.exp(np.log(self.pi) + log_categorical(X, self.theta))
        p_Z_normalizing_const = np.sum(p_Z_not_normalized, axis=1)
        p_Z = p_Z_not_normalized / np.expand_dims(p_Z_normalizing_const, 1)
        self.Z = vectorized_np_choice(p_Z)
        
    def _sample_theta(self,X, gamma):
        count = np.array([np.sum(X[self.Z==k], axis=0) for k in range(self.num_clusters)])
        gamma_prime = gamma + count
        self.theta = np.array([np.random.dirichlet(gamma_prime[k]) for k in range(self.num_clusters)])
    

""" Sample from a dirichlet process """
class PiCollapsedNonparametricGibbsSampler:
    
    def __init__(self, seed=None):
        np.random.seed(seed)
    
    def fit(self, X, num_steps, num_burn_in_steps, alpha=1, gamma=1):
        # Initialization
        log_likelihood_progression = []
        I = len(X[0])
        # Starting with one cluster
        self.K_seen = 1
        self.Z = np.zeros(len(X))
        self.theta = np.zeros((self.K_seen,I))
        self.theta[0] = np.random.dirichlet(I*[gamma])
        self.logpx_n_unseen = np.zeros(len(X))

        for n in range(len(X)):
            first_gamma = np.sum(loggamma(gamma + X[n]))
            second_gamma = loggamma(np.sum([gamma]*I))
            third_gamma = loggamma(np.sum(gamma+X[n]))
            four_gamma = np.sum(loggamma([gamma]*I))
            self.logpx_n_unseen[n] = first_gamma + second_gamma - third_gamma - four_gamma
        
        self.m_k = np.zeros(self.K_seen)
        for k in range(self.K_seen):
            self.m_k[k] = np.sum(self.Z == k)
            
        self.gamma_km = np.zeros((self.K_seen,I))
        for k in range(self.K_seen):
            for m in range(I):
                self.gamma_km[k,m] = gamma + np.sum(X[self.Z == k][:,m])
            
        for _ in range(num_steps):
            self._sample_Z(X, alpha, gamma)
            for k in range(self.K_seen):
                self.theta[k] = np.random.dirichlet(self.gamma_km[k])
            print(max(self.Z))
            log_likelihood = sum([np.sum(X[self.Z==k] @ np.log(self.theta[k])) for k in range(self.K_seen)])
            print(log_likelihood)
            log_likelihood_progression += [log_likelihood]

        return log_likelihood_progression
    
    def _sample_Z(self, X, alpha, gamma):
        """Samples the most likely vector Z of length num_datapoints by using X and the current cluster assignments"""
        for n in range(len(X)):
            z_old = int(self.Z[n])
            
            self.m_k[z_old] -= 1
            
            # For seen clusters:
            p_zn_Zn_seen = self.m_k/(len(X)-1 + alpha)
            logp_xn_theta_seen = X[n] @ np.log(self.theta.T)
            #print(logp_xn_theta_seen)
            # For unseen clusters:
            p_zn_Zn_unseen = alpha/(len(X)-1+alpha)
            logpx_n_unseen = self.logpx_n_unseen[n]
    
            p_zn_seen = p_zn_Zn_seen * np.exp(logp_xn_theta_seen)
            p_zn_unseen = p_zn_Zn_unseen * np.exp(logpx_n_unseen)

            p_zn = np.append(p_zn_seen, p_zn_unseen)
            #print(p_zn)
            p_zn /= np.sum(p_zn, dtype=np.longfloat)

            new_z = np.random.choice(range(self.K_seen + 1), p=p_zn)

            # Also need to sample new theta...
            # Should be done
            # Also remove emtpy clusters..
            
            if new_z >= self.K_seen:
                self.m_k = np.append(self.m_k, [0])
                self.m_k[new_z] += 1
                self.K_seen += 1
                
                # New thetas
                self.gamma_km = np.vstack((self.gamma_km[:self.K_seen-1], X[n] + gamma))
                self.gamma_km[int(self.Z[n])] -= X[n]
                theta_new = np.random.dirichlet(self.gamma_km[new_z])
                self.theta = np.vstack((self.theta[:self.K_seen-1], theta_new))

            else:                
                self.m_k[new_z] += 1
                self.gamma_km[int(self.Z[n])] -= X[n]
                self.gamma_km[new_z] += X[n]
                                
            self.Z[n] = new_z
            
            
            if self.m_k[z_old] == 0:
                self.m_k = np.delete(self.m_k,z_old)
                self.gamma_km = np.delete(self.gamma_km,z_old, axis = 0)
                self.theta = np.delete(self.theta,z_old, axis = 0)
                self.K_seen -= 1
                self.Z[self.Z > z_old] -= 1
                        





