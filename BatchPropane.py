import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#x = np.linspace(0, 1, num=1000)

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

def rline(R, xd, xinput):
    yr = R / (R + 1) * xinput + 1 / (R + 1) * xd
    return yr


def findxw(n0, V, xd, xw0, P, a, b, T, c, a2, b2, c2):
    t = 10/60
    step = 1
    n = n0 - t*V
    xw = (n*xd - n0*(xd - xw0))/n

    R = .01

    xinput = xd

    yr = R / (R + 1) * xinput + 1 / (R + 1) * xd

    xhor1 = fsolve(eqline_input, xinput, args=(P, a, b, T, c, a2, b2, c2, yr))

    yvert = rline(R, xd, xhor1)

    xhor2 = fsolve(eqline_input, xinput, args=(P, a, b, T, c, a2, b2, c2, yvert))

    yvert2 = rline(R, xd, xhor2)

    xhor3 = fsolve(eqline_input, xinput, args=(P, a, b, T, c, a2, b2, c2, yvert2))

    yvert3 = rline(R, xd, xhor3)

    if xhor3 > xw :
        plt.hlines(yr, xhor1, xd)
        plt.vlines(xhor1, yvert, yr)
        plt.hlines(yvert, xhor2, xhor1)
        plt.vlines(xhor2, yvert2, yvert)
        plt.hlines(yvert2, xhor3, xhor2)
        plt.vlines(xhor3, yvert3, yvert2)
    else:
        R = R + .01


    return

def eqplot(xd, xb, P, a, b, T, c, a2, b2, c2):
    y = yxline()
    yeq = eqline(P, a, b, T, c, a2, b2, c2)
    line, = plt.plot(x, y, label="xy line")
    line, = plt.plot(x, yeq, label= "Equilibrium")
    stages(xd, xb, P, a, b, T, c, a2, b2, c2)
    plt.legend()
    plt.title("Theoretical minimum stages")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()

    return

eqplot(.90, .1, 760, 8.00308, 1505.52, 83, 211.6, 8.37895, 1788.02, 227.438)