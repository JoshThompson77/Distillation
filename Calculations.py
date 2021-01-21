import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

x = np.linspace(0, 1, num=1000)


def calcK1(Psat, P, A):
    gamma1 = 10 ** (A * (1 - x) ** 2)
    K = Psat/ P
    return K

def  calk2(Psat, P, A):
    gamma2 = 10 ** (A * x ** 2)
    K2 = Psat /P
    return K2

def yxline():
    y = x
    return y

def alpha12( Psat1, P, Psat2, A):
    K1 = calcK1(Psat1, P, A)
    K2 = calk2(Psat2, P, A)
    alpha12 = K1/K2
    return alpha12

def eqline(Psat1, P, Psat2, A):
    alpha = alpha12(Psat1, P, Psat2, A)
    gamma1 = 10 ** (A * (1 - x) ** 2)
    gamma2 = 10 ** (A * x ** 2)
    yeq = (alpha * x) / ((alpha * x) + (gamma1 / gamma2) * (1 - x))
    return yeq

def qline(q, zf, Psat1, P, Psat2, A):
    xq = np.linspace(.2, zf, num= 100)
    yq = q/(q-1) * x - zf/(q-1)
    yeq = eqline(Psat1, P, Psat2, A)
    xint1, yint1 = findintersect(yq, yeq, x)
    xq = np.linspace(xint1, zf, num=1000)
    yq = yq = q/(q-1) * xq - zf/(q-1)
    return xq, yq

def rlineLV(L, V, D, xd):
    y = L/V * x + D/V* xd
    return y

def rlineR(R, q, zf, xd):
    global xint, yint
    xint = (1/(R+1)* xd + zf/(q-1))/(q/(q-1)-R/(R+1))
    yint= R / (R + 1) * xint + 1 / (R + 1) * xd
    xr = np.linspace(xint, xd, num=1000)
    yr = R / (R + 1) * xr + 1 / (R + 1) * xd
    print()
    return xr, yr, xint, yint

def slineVB(xb, q, zf, R, xd):
    xint2 = fsolve(qrintersect, .5, args=(q, zf, R, xd))
    yint2 = q/(q-1) * xint2 - zf/(q-1)
    VB = 1/((yint2 - xb)/(xint2 -xb)-1)
    xs = np.linspace(xb, xint2, num= 500)
    ys = (VB + 1)/ VB * xs - 1/VB * xb
    return xs, ys

# def slineLV(L, V, xb, B, q, zf):
#     ys = L/V * x - B/V * xb
#     yq = q / (q - 1) * x - zf / (q - 1)
#     xint, yint = findintersect(yq, ys, x)
#     return xs, ys

def findintersect(yline1, yline2, xline):
    """Finds the intersection of two lines with arrays of the same size"""
    xint1 = 0
    yint1 = 0
    for i in range(1000):
        if round(yline2[i], 3) >= round(yline1[i], 3):
            yint1 = yline1[i]
            xint1 = xline[i]
            break

    return xint1, yint1

def newqline(q, zf, xinput):
    nxq = xinput
    nyq = q / (q - 1) * nxq - zf / (q - 1)
    return nxq, nyq

def yeq(xd, Psat1, P, Psat2):
    alpha = alpha12(Psat1, P, Psat2)
    yeq = alpha * xd / (1 + xd*(alpha - 1))
    return yeq

def rline(R, xd, xinput):
    yr = R / (R + 1) * xinput + 1 / (R + 1) * xd
    return yr

def equiline(xh, Psat1, P, Psat2, A):
    global yhor1
    alpha = alpha12(Psat1, P, Psat2, A)
    gamma1 = 10 ** (A * (1 - xh) ** 2)
    gamma2 = 10 ** (A * (xh) ** 2)
    equilline = (alpha * xh) / ((alpha * xh) + (gamma1 / gamma2) * (1 - xh))-yhor1
    return equilline

def intercept(R, xd, zf, q):
    xint = (1 / (R + 1) * xd + zf / (q - 1)) / (q / (q - 1) - R / (R + 1))
    return xint

def sline(xb, xinput):
    VB = 1/((yint - xb) / (xint - xb)-1)
    ys = (VB + 1) / VB * xinput - 1 / VB * xb
    return ys


def stages(R, q, zf, xd, Psat1, P, Psat2, xb, A):

    xpoint = xd
    global yhor1

    xints = intercept(R, xd, zf, q)
    while  xpoint > xb:
        if xpoint > xints :
            yhor1 = rline(R, xd, xpoint)
            xhor = fsolve(equiline, xints, args=(Psat1, P, Psat2, A))
            plt.hlines(yhor1, xhor, xpoint)

            xpoint = xhor
            if xpoint > xints :
                yvert = rline(R, xd, xpoint)
                plt.vlines(xpoint, yvert, yhor1)
            else:
                yvert = sline(xb, xpoint)
                plt.vlines(xpoint, yvert, yhor1)
        else:
            yhor1 = yvert
            xhor = fsolve(equiline, .1, args=(Psat1, P, Psat2, A))
            plt.hlines(yhor1, xhor, xpoint)

            xpoint = xhor
            yvert = sline(xb, xpoint)
            plt.vlines(xpoint, yvert, yhor1)
    plt.show()
    return

# def f(x):
#     return 2*x + 2
#
# num = fsolve(f, 2)
#
# print(num * 4)
# num = num[0]
#
# print(num)

# stages(1, .5, .5, .7, 640, 760, 23.8, .1, 1)

# print(yeq(0.14916246, 640, 760, 23.8))

def qrintersect(xinput,q, zf, R, xd):
    yq = q / (q - 1) * xinput - zf / (q - 1)
    qrintersection = R / (R + 1) * xinput + 1 / (R + 1) * xd - yq
    return  qrintersection


def McCabeThielePlot (Psat1, P, Psat2, A, q, zf, R, xd, xb):
    y = yxline()
    yeq = eqline(Psat1, P, Psat2, A)
    xq, yq = qline(q, zf, Psat1, P, Psat2, A)
    xr, yr, xint, yint = rlineR(R, q, zf, xd)
    xs, ys = slineVB(xb, q, zf, R, xd)
    #nxq, nyq = newqline(q, zf)

    plt.plot(x, y, x, yeq, xq, yq, xr, yr, xs, ys)
    stages(R, q, zf, xd, Psat1, P, Psat2, xb, A)
    plt.show()

McCabeThielePlot(200, 243, 42.4, .7947, .5, .5, 1, .7, .3)






