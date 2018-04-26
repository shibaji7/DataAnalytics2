import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

np.random.seed(1)


def f(x,m):
    """The function to predict."""
    #return x * np.sin(x)
    return 1-(np.exp(x/(100*m)))

def mm(y):
    m = (y - min(y))/(max(y)-min(y))
    return m
M = [5.0,5.5,6.0,6.5,7.0,7.5]
for m in M:
# ----------------------------------------------------------------------
#  First the noiseless case
    X = np.atleast_2d([10, 25, 65, 90]).T

# Observations
    y = f(X,m).ravel()

# Mesh the input space for evaluations of the real function, the prediction and
# its MSE
    x = np.atleast_2d(np.linspace(5, 100, 1000)).T

# Instanciate a Gaussian Process model
    kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
    gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)

# Fit to data using Maximum Likelihood Estimation of the parameters
    gp.fit(X, y)

# Make the prediction on the meshed x-axis (ask for MSE as well)
    y_pred, sigma = gp.predict(x, return_std=True)

# Plot the function, the prediction and the 95% confidence interval based on
# the MSE
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(x, f(x,m), 'g:', label=u'$C(M\geq %s)$'%str(m))
    plt.plot(X, y, 'r.', markersize=10, label=u'Observations')
    plt.plot(x, y_pred, 'k-', label=u'Prediction')
    plt.fill(np.concatenate([x, x[::-1]]),
                     np.concatenate([y_pred - 1.9600 * sigma,
                                                 (y_pred + 1.9600 * sigma)[::-1]]),
                              alpha=.5, fc='b', ec='None', label='95% confidence interval')
    
    ax.set_yticklabels(np.arange(10,400,20))
    plt.xlabel('$\delta in [^o]$')
    plt.ylabel('$Count(M\geq m)$')
    plt.legend(loc='lower left')
    plt.savefig("out/tgp%s.png"%str(m))
