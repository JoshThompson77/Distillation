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

def stages(xd, xb, P, a, b, T, c, a2, b2, c2):
    xinput = xd
    yhor1 = xd
    plt.vlines(xinput, yhor1, 1)
    while xinput > xb:
        yhor1 = yxline_input(xinput)
        xhor = fsolve(eqline_input, xinput, args=(P, a, b, T, c, a2, b2, c2, yhor1))
        plt.hlines(yhor1, xhor, xinput)


        yvert = xhor
        plt.vlines(xhor, yvert, yhor1)
        xinput = xhor

    return



def eqplot(xd, xb, P, a, b, T, c, a2, b2, c2):
    y = yxline()
    yeq = eqline(P, a, b, T, c, a2, b2, c2)
    line, = plt.plot(x, y, label="xy line")
    line, = plt.plot(x, yeq, label= "Equilibrium")
    stages(xd, xb, P, a, b, T, c, a2, b2, c2)
    plt.legend()
    plt.title("Infinite reflux stages")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()

    return


eqplot(.93, .26, 760, 8.00308, 1505.52, 83, 211.6, 8.37895, 1788.02, 227.438)





