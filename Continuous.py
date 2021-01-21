import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

x = np.linspace(0, 1, num=100)
# x = .7

def eqfitcurve():
    yfit = 16.2768*x**5 - 47.294*x**4 + 51.822*x**3 - 25.995*x**2 + 6.198*x + 0.0153
    return yfit

def eqfitcurveX(xinput, q, zf):
    yfit = 16.2768 * xinput ** 5 - 47.294 * xinput **4 + 51.822 * xinput ** 3 - 25.995 * xinput**2 + 6.198 * xinput + \
           0.0153 - q / (q - 1) * xinput + zf / (q - 1)
    return yfit

def eqfitcurveY(xinput, yhor1):
    yfit = 16.2768 * xinput ** 5 - 47.294 * xinput **4 + 51.822 * xinput ** 3 - 25.995 * xinput**2 + 6.198 * xinput + \
           0.0153 - yhor1
    return yfit

def yxline():
    y = x
    return y

def qline(q, zf):
    xinput = .5
    xint = fsolve(eqfitcurveX, xinput, args=(q, zf))
    xq = np.linspace(xint, zf, num=100)
    yq = q/(q-1) * xq - zf/(q-1)
    return yq, xq

qline(.5, .5)

def qlineX(xinput, q, zf, R, xd):
    yq = q / (q - 1) * xinput - zf / (q - 1)- R/(R+1) * xinput - 1/(R+1) * xd
    return yq

def rline(q, zf, xd, R):
    xinput = .5
    xint = fsolve(qlineX, xinput, args=(q, zf, R, xd))
    xr = np.linspace(xint, xd, num=1000)
    yr = R / (R + 1) * xr + 1 / (R + 1) * xd
    return yr, xr

def slineVB(xb, q, zf, R, xd):
    xint2 = fsolve(qlineX, .5, args=(q, zf, R, xd))
    yint2 = q/(q-1) * xint2 - zf/(q-1)
    VB = 1/((yint2 - xb)/(xint2 -xb)-1)
    xs = np.linspace(xb, xint2, num= 500)
    ys = (VB + 1)/ VB * xs - 1/VB * xb
    return xs, ys

def slineX(xinput, xb, q, zf, R, xd):
    xint2 = fsolve(qlineX, .5, args=(q, zf, R, xd))
    yint2 = q / (q - 1) * xint2 - zf / (q - 1)
    VB = 1 / ((yint2 - xb) / (xint2 - xb) - 1)
    ys = (VB + 1) / VB * xinput - 1 / VB * xb
    return ys

def rlineX(xinput, xd, R):
    yr = R / (R + 1) * xinput + 1 / (R + 1) * xd
    return yr

def stages(R, q, zf, xd, xb):

    xpoint = xd

    xints = fsolve(qlineX, xd, args=(q, zf, R, xd))
    while  xpoint > xb:
        if xpoint > xints :
            yhor1 = rlineX(xpoint, xd, R)
            xhor = fsolve(eqfitcurveY, xints, args=(yhor1))
            plt.hlines(yhor1, xhor, xpoint)

            xpoint = xhor
            if xpoint > xints :
                yvert = rlineX(xpoint, xd, R)
                plt.vlines(xpoint, yvert, yhor1)
            else:
                yvert = slineX(xpoint, xb, q, zf, R, xd)
                plt.vlines(xpoint, yvert, yhor1)
        else:
            yhor1 = yvert
            xhor = fsolve(eqfitcurveY, xints, args=(yhor1))
            plt.hlines(yhor1, xhor, xpoint)

            xpoint = xhor
            yvert = slineX(xpoint, xb, q, zf, R, xd)
            plt.vlines(xpoint, yvert, yhor1)
    plt.show()
    return

def McabeThielePlot(q, zf, xd, R, xb):
    yfit = eqfitcurve()
    y = yxline()
    yq, xq = qline(q, zf)
    yr, xr = rline(q, zf, xd, R)
    xs, ys = slineVB(xb, q, zf, R, xd)
    plt.plot(x, y, x, yfit, xr, yr, xq, yq, xs, ys)
    stages(R, q, zf, xd, xb)
    plt.show()
    return

McabeThielePlot(.5, .5, .8, 2, .2)