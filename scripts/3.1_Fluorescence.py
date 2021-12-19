# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 11:10:21 2021

@author: marcu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.util import describe, make_func_code


class LeastSquares:
    """
    Generic least-squares cost function with error.
    """

    errordef = Minuit.LEAST_SQUARES # for Minuit to compute errors correctly

    def __init__(self, model, x, y, err):
        self.model = model  # model predicts y for given x
        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.err = np.asarray(err)

    def __call__(self, *par):  # we accept a variable number of model parameters
        ym = self.model(self.x, *par)
        return np.sum((self.y - ym) ** 2 / self.err ** 2)


class BetterLeastSquares(LeastSquares):
    def __init__(self, model, x, y, err):
        super().__init__(model, x, y, err)
        self.func_code = make_func_code(describe(model)[1:])
        

## Load data
big_path = "C:\\Users\\marcu\\OneDrive\\Desktop\\PraktikumIII\\e+e-_Annihilation\\"
path = big_path + "data\\3.1\\Dynode3.1\\800V\\800V_ALL0000_CurveData.txt"
df = pd.read_csv(path);

# df[!, :Time] = df[!, :Time] * 1e9; # nanoseconds
# df[!, :Time] = df[!, :Time] .- minimum(df[!, :Time]);
# df[!, :Time] = df[!, :Time] * u"ns"
# df[!, :Height] = abs.(df[!, :Height] * u"V") .|> u"mV"
# df[!, :σ] .= 0


## Cost function
def line(t, A, a, B, b, C, c):  # simple straight line model with explicit parameters
    return A*np.exp(-t/a) + B*np.exp(-t/b) - C*np.exp(-t/c)

t_data = df["Time"].to_numpy()
t_data = (t_data - t_data[0]) * 1e9

v_data = df["Height"] * -1000
v_err = np.full_like(v_data, 1)

start = 180
# start = 278

with plt.style.context("science"):
    plt.figure(dpi=350, figsize=(9,4))
    plt.xlim(0, 1000)
    plt.tight_layout()
    plt.errorbar(t_data, v_data, yerr=v_err, capsize=0.8, errorevery=5, label="Data")
    plt.xlabel("Time (ns)")
    plt.ylabel("Height (mV)")
    plt.grid()



lsq = BetterLeastSquares(line, t_data, v_data, v_err)

m = Minuit(lsq, A=120.0, a=230, B = 520.0, b = 5.0, C=1200.0, c=5.0)
for i in m.parameters:
    m.limits[i] = (1e-2, None)

# print(m.limits)

m.migrad()
# m.hesse() -> fails

plt.plot(t_data[start:], line(t_data[start:], *m.values), label="Fit")
plt.scatter([t_data[start]], [v_data[start]], color="red", s=2, zorder=3, label="Starting point of fit")
plt.annotate("$A e^{-t / \\tau_s} + B e^{-t / \\tau_f} - C e^{-t / \\tau}$", (600, 83), fontsize=12)
plt.annotate("$A=$"+f"{round(m.values[0], ndigits=2)}" + "$\\pm$" + f"{round(m.errors[0], ndigits=2)}", (600, 73))
plt.annotate("$\\tau_s=$"+f"{round(m.values[1], ndigits=2)}"+ "$\\pm$" + f"{round(m.errors[1], ndigits=2)}", (800, 73))
plt.annotate("$B=$"+f"{round(m.values[2], ndigits=2)}"+ "$\\pm$" + f"{round(m.errors[2], ndigits=2)}", (600, 63))
plt.annotate("$\\tau_f=$"+f"{round(m.values[3], ndigits=2)}"+ "$\\pm$" + f"{round(m.errors[3], ndigits=2)}", (800, 63))
plt.annotate("$C=$"+f"{round(m.values[4], ndigits=2)}"+ "$\\pm$" + f"{round(m.errors[4], ndigits=2)}", (600, 53))
plt.annotate("$\\tau=$"+f"{round(m.values[5], ndigits=2)}"+ "$\\pm$" + f"{round(m.errors[5], ndigits=2)}", (800, 53))
plt.legend() 

print(m)

from matplotlib import gridspec
fig = plt.figure(figsize=(12, 15), dpi=900)
# fig.tight_layout()
fig = gridspec.GridSpec(3, 1, height_ratios=[1,1,1])

x1 = plt.subplot(fig[0])
m.draw_profile("a")
plt.title(f"$\\tau_s = {round(m.values[1], ndigits=2)} \\pm {round(m.errors[1], ndigits=2)}$ ns")
plt.xlabel("$\\tau_s $")

x2 = plt.subplot(fig[1])
m.draw_profile("b")
plt.title(f"$\\tau_f = {round(m.values[3], ndigits=2)} \\pm {round(m.errors[3], ndigits=2)}$ ns")
plt.xlabel("$\\tau_f $")

x3 = plt.subplot(fig[2])
m.draw_profile("c")
plt.title(f"$\\tau = {round(m.values[5], ndigits=2)} \\pm {round(m.errors[5], ndigits=2)}$ ns")
plt.xlabel("$\\tau $")


