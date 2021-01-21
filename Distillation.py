#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 07:29:02 2019

@author: Hunter McBrayer
"""
#Import Modules

import numpy as np
import matplotlib.pyplot as plt

#Calculation Functions
def create_yxLine():
    """Creates the y=x line"""
    global x,y
    x = np.linspace(0,1,num=1000)
    y = x
    return y,x
def create_EQLine(alphaBA,A):
    """Calculates the y and x values for the Equilibrium line"""
    global xEQ, yEQ
    xEQ = np.linspace(0,1,num=1000)
    gamma1 = 10**(A*(1-xEQ)**2)
    gamma2 = 10**(A*(xEQ)**2)
    yEQ = (alphaBA*xEQ)/((alphaBA*xEQ) + (gamma1/gamma2)*(1-xEQ))
    return yEQ, xEQ
def create_qLine(q,zF):
    """Calculates the y and x values for the q line"""
    global xq, yq
    xq = np.linspace(0,zF,num=500)
    yq = (q/(q-1))*xq - (zF/(q-1))
    return yq, xq
def findintersect(yline1, yline2, xline):
    """Finds the intersection of two lines with arrays of the same size"""
    global xint, yint
    for i in range(500):
        if round(yline1[i], 2) == round(yline2[i],2):
            yint = yline1[i]
            xint = xline[i]
    return xint, yint
def calc_zFEQ(yq,xq,yEQ):
    """Finds the Equilibrium Feed mole fraction"""
    global xFEQ, yFEQ
    findintersect(yq,yEQ,xq)
    xFEQ = xint
    yFEQ = yint
    return yFEQ, xFEQ
def calc_Rmin(xD,yFEQ,xFEQ):
    """Calculates the Minimum Reflux Ratio"""
    global Rmin
    LVmin = (yFEQ - xD)/(xFEQ - xD)
    Rmin = LVmin/(1-LVmin)
    return Rmin
def create_RLine(RRmin,Rmin,xD):
    """Calculated the y and x values for the Reflux line"""
    global xR, yR
    R = RRmin * Rmin
    xR = np.linspace(0,1,num=1000)
    yR = ((R/(R+1))*xR) + ((1/(R+1))*xD)
    findintersect(yq, yR, xR)
    xR = np.array([xint, xD])
    yR = ((R/(R+1))*xR) + ((1/(R+1))*xD)
    return yR, xR
def calc_VBmin(xB,yFEQ,xFEQ):
    """Calculates the Minimum Reboil"""
    global VBmin
    LVmax = (yFEQ - xB)/(xFEQ - xB)
    VBmin = 1/(LVmax - 1)
    return VBmin
def calc_VVBminOp(yR,yq,xR,xB,VBmin):
    """Calculates the Ratio of VB to VBmin"""
    global VVBminOp
    m = (yint - xB)/(xint - xB)
    VB = 1/(m-1)
    VVBminOp = (VB/VBmin)
    return VVBminOp
def create_VBLine(VVBminOp,VBmin,xB):
    """Calculated the y and x values for the Reboil line"""
    global yVB, xVB
    VB = VVBminOp*VBmin
    xVB = np.linspace(.05,xint,num=500)
    yVB = ((VB + 1)/VB)*xVB - (1/VB)*xB
    return yVB, xVB
def fill_arrays(alphaBA,A,zF,q, xD, RRmin, xB):
    """Fills the arrays of every y and x value for each line"""
    create_yxLine()
    create_EQLine(alphaBA, A)
    create_qLine(q, zF)
    calc_zFEQ(yq, xq, yEQ)
    calc_Rmin(xD, yFEQ, xFEQ)
    create_RLine(RRmin, Rmin, xD)
    calc_VBmin(xB, yFEQ, xFEQ)
    calc_VVBminOp(yR, yq, xR, xB, VBmin)
    create_VBLine(VVBminOp, VBmin, xB)
    return x,y,xEQ,yEQ,xq,yq,xR,yR,xVB,yVB
def create_plot(alphaBA,A,zF, q, xD, RRmin, xB):
    """Creates a plot for the Distillation system"""
    fill_arrays(alphaBA,A,zF,q, xD, RRmin, xB)
    plt.plot(x,y,xEQ,yEQ,xq,yq,xR,yR,xVB,yVB)
    plt.xlabel('$x_B$')
    plt.ylabel('$y_B$')
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.show()
def print_values(alphaBA,A,zF, q, xD, RRmin, xB):
    """Prints Rmin, R, VB, VBmin, and VVBminOp"""
    fill_arrays(alphaBA,A,zF,q, xD, RRmin, xB)
    print('Minimum Reboil Ratio:',VBmin,'\nVB over VBmin at Operation:',VVBminOp
          ,'\nReboil Ratio:',VVBminOp*VBmin,'\nMinimum Reflux Ratio',Rmin
          ,'\nReflux Ratio:',RRmin*Rmin
          , '\nThe feed mole fraction of vapor at equilibrium:',yFEQ
          ,'\nThe feed mole fraction of liquid at equilibrium:',xFEQ)


create_plot(3, -.9, .5, .3, .735, 1.5, .05)