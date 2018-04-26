
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

import gpplt

np.random.seed(1)

def dd_sort(x,y):
    df = pd.DataFrame()
    df["x"] = x
    df["y"] = y
    df = df.sort_values(by=["x"])
    return np.array(df["x"].tolist()),np.array(df["y"].tolist())

time_window = np.hstack( (1, np.arange(5, 105, step=5) ) )
tw_cols = list(time_window.astype('str'))
pv = pd.read_csv('../results/panamint_valley_gr.csv', index_col=[0])
X = pv[["M","dip","Ddot"]]
XX = pd.DataFrame()
Y = []
for I in time_window:
    X["tw"] = ([I] * 12000)
    XX = XX.append(X)
    Y = Y + pv[str(I)].tolist()
    pass
Y = np.array(Y)

N = 5
idx = np.random.randint(len(Y), size=N)
X_d = XX.as_matrix()[idx,:]
y_d = Y[idx]
print X_d.shape,len(y_d)

kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
gp.fit(X_d, y_d)

np.random.seed(300)
idt = np.random.randint(len(Y), size=100)
X_t = XX.as_matrix()[idt,:]
y_t = Y[idt]
y_pred, sigma = gp.predict(X_t, return_std=True)

print y_pred,sigma

f,axes = plt.subplots(nrows=1,ncols=1,figsize=(8,6),dpi=120)
ax = axes
x,y = dd_sort(X_t[:,1],y_t)
x,yp = dd_sort(X_t[:,1],y_pred)
ax.plot(x,y,"ro")
ax.plot(x,yp-(1.96*sigma),"b-.")
ax.plot(x,yp+(1.96*sigma),"k-.")




f.savefig("out/gp1.png")
