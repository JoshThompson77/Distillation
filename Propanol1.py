import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

xtrue = [0, 0.0575,0.061, 0.1455, 0.2285, 0.3095, 0.3125, 0.42, 0.4355, 0.519, 0.631, 0.7305, 0.7675,
0.8585, 0.91, 1]

ytrue = [0, 0.11, 0.111, 0.2325, 0.351, 0.4435, 0.45, 0.5545, 0.5725, 0.66, 0.748, 0.8225, 0.8495,
0.9175, 0.9525, 1]

x = np.linspace(0, 1, num=1000)

def antoine_eq(a, b, T, c):
    Psat = 10 ** (a - b/(c + T))
    return Psat

def calc_k(P, a, b, T, c):
    Psat = antoine_eq(a, b, T, c)
    K = Psat/ P
    return K

def calc_alpha(P, a, b, T, c, a2, b2, c2):
    K1 = calc_k(P, a, b, T, c)
    K2 = calc_k(P, a2, b2, T, c2)
    alpha = K1/K2
    return alpha

def eqline(P, a, b, T, c, a2, b2, c2):
    alpha = calc_alpha(P, a, b, T, c, a2, b2, c2)
    yeq = alpha * x/ (1 + x * (alpha - 1))
    return yeq

def eqline_input(xinput, P, a, b, T, c, a2, b2, c2, yhor1):
    alpha = calc_alpha(P, a, b, T, c, a2, b2, c2)
    yeq = alpha * xinput / (1 + xinput * (alpha - 1)) - yhor1
    return yeq

def yxline_input(xinput):
    y = xinput
    return y

def yxline():
    y=x
    return y

def rline(R, xd):
    xr = np.linspace(0, xd, 1000)
    yr = R/(R + 1)* xr + xd/(R + 1)
    return yr, xr

def rline_input(xinput, R, xd):
    yr = R / (R + 1) * xinput + xd / (R + 1)
    return yr



def stages(xd, xb, P, a, b, T, c, a2, b2, c2, R):
    xinput = xd
    yhor1 = xd
    plt.vlines(xinput, yhor1, 1)
    counter = 0
    while xinput > xb:
        yhor1 = rline_input(xinput, R, xd)
        xhor = fsolve(eqline_input, xinput, args=(P, a, b, T, c, a2, b2, c2, yhor1))
        plt.hlines(yhor1, xhor, xinput)


        yvert = rline_input(xhor, R, xd)
        plt.vlines(xhor, yvert, yhor1)
        xinput = xhor
        counter += 1
        if counter > 35:
            break

    print(counter)

    return



def eqplot(xd, xb, P, a, b, T, c, a2, b2, c2, R):
    y = yxline()
    yeq = eqline(P, a, b, T, c, a2, b2, c2)
    yr, xr = rline(R, xd)
    line, = plt.plot(xr, yr, label="R line")
    line, = plt.plot(x, y, label="xy line")
    line, = plt.plot(x, yeq, label= "Equilibrium")
    stages(xd, xb, P, a, b, T, c, a2, b2, c2, R)
    plt.legend()
    plt.title("Reflux ratio of 1")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()

    return


eqplot(0.706016, 0.259504, 760, 8.00308, 1505.52, 83, 211.6, 8.37895, 1788.02, 227.438, 1)





