import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import random

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF,RationalQuadratic as RQ, ConstantKernel as C

np.random.seed(1)
random.seed(1)

fn=2
def f(x,m):
    """The function to predict."""
    #return x * np.sin(x)
    return 1-(1/(m**2*np.sin(np.deg2rad(x))))

def mm(y):
    m = (y - min(y))/(max(y)-min(y))
    return m
M = [5.0,5.5,6.0,6.5,7.0,7.5]
X = np.atleast_2d([10., 13., 15., 16., 12., 11.]).T
Llim = [75,62,55,48,40,30]
for I,m in enumerate(M):
# ----------------------------------------------------------------------
#  First the noiseless case

    X = np.atleast_2d(np.array(random.sample(range(800, 1800), 7), dtype=float)/100).T
# Observations
    y = f(X,m).ravel()

# Mesh the input space for evaluations of the real function, the prediction and
# its MSE
    x = np.atleast_2d(np.linspace(5, 30, 1000)).T

# Instanciate a Gaussian Process model
    if fn ==1:
        kernel = C(1.0, (1e-3, 1e3)) * RBF(1, (1e-2, 1e2))
    elif fn==2:
        kernel = C(1.0, (1e-3, 1e3)) * RQ(length_scale=5.0, alpha=5.0, length_scale_bounds=(1e-02, 1e2), alpha_bounds=(1e-02, 1e2))
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
    ymin,ymax = ax.get_ylim()
    print ymax,ymin
    #plt.ylim(0.3,1.3)
    ax.set_yticks(np.linspace(ymin,ymax,5))
    print np.linspace(ymin,ymax,5)
    ax.set_yticklabels(np.arange(Llim[I],200,7))
    plt.xlabel('$\delta[^o]$')
    plt.ylabel('$Count(M\geq m)$')
    plt.legend(loc=0)
    plt.xlim(5,30)
    if fn==1: plt.title(r"$\kappa(x_i,x_j)=\sigma^2 e^{\frac{||x_i-x_j||^2}{2\lambda^2}}$")
    elif fn==2: plt.title(r"$\kappa(x_i,x_j)=\sigma^2 \left(1+\frac{||x_i-x_j||^2}{2\alpha\lambda^2}\right)^{-\alpha}$")
    plt.savefig("out/gp_%d_%s.png"%(fn,str(m)))
