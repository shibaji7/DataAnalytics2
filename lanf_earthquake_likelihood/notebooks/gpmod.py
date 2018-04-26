
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
    print x.shape,y.shape
    df["x"] = x
    df["y"] = y
    df = df.sort_values(by=["x"])
    return np.array(df["x"].tolist()),np.array(df["y"].tolist())

pv = pd.read_csv('../results/panamint_valley_gr.csv', index_col=[0])
XY = pv[["M","dip","Ddot","35"]]
print XY.shape
XY = XY[XY.M==6.5]
print XY.shape
X = np.array(XY["dip"].tolist())
Y = XY["35"]
Y = np.array(Y.tolist())

N = 10
idx = np.random.randint(len(Y), size=N)
X_d = X[idx]
y_d = Y[idx]
print X_d.shape,len(y_d)

kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
gp.fit(X_d.reshape(len(X_d),1), y_d)

np.random.seed(300)
idt = np.random.randint(len(Y), size=200)
X_t = X[idt]
y_t = Y[idt]
y_pred, sigma = gp.predict(X_t.reshape(len(X_t),1), return_std=True)

#print y_pred,sigma

f,axes = plt.subplots(nrows=1,ncols=1,figsize=(8,6),dpi=120)
ax = axes
x,y = dd_sort(X_d,y_d)
ax.semilogy(x,np.abs(y),"ro")
x,yp = dd_sort(X_t,y_pred)
print len(x),len(yp)
ax.semilogy(x,np.abs(yp-(1.96*sigma)),"b-.")
ax.semilogy(x,np.abs(yp+(1.96*sigma)),"k-.")




f.savefig("out/gp2.png")
