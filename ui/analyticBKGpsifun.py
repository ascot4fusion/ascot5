#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math as m
from numpy import log as l
from numpy import power as p

def psi(x,y,i):
    return {
        0 : 1.0,
        1 : p(x,2),
        2 : p(x,2)*l(x)-p(y,2),
        3 : p(x,4)-4*p(x,2)*p(y,2),
        4 : 3.0*p(x,4.0)*l(x)-9.0*p(x,2.0)*p(y,2.0)-12.0*p(x,2.0)*l(x)*p(y,2.0)+2.0*p(y,4.0),
        5 : p(x,6)-12*p(x,4)*p(y,2)+8*p(x,2)*p(y,4),
        6 : 8*p(y,6)-140*p(x,2)*p(y,4)-120*p(x,2)*l(x)*p(y,4)+180*p(x,4)*l(x)*p(y,2)+75*p(x,4)*p(y,2)-15*p(x,6)*l(x), 
        7 : y,
        8 : y*p(x,2),
        9 : p(y,3)-3*y*p(x,2)*l(x),
        10: 3*y*p(x,4)-4*p(y,3)*p(x,2),
        11: 8*p(y,5)-45*y*p(x,4)-80*p(y,3)*p(x,2)*l(x)+60*y*p(x,4)*l(x),
        }[i]

def psix(x,y,i):
    return {
        0 : 0.0,
        1 : 2*x,
        2 : 2*x*m.log(x)+x,
        3 : 4*x**3-8*x*y**2,
        4 : 12*x**3*m.log(x)+3*x**3-30*x*y**2-24*x*m.log(x)*y**2,
        5 : 6*x**5-48*x**3*y**2+16*x*y**4,
        6 : -400*x*y**4-240*x*m.log(x)*y**4+720*x**3*m.log(x)*y**2+480*x**3*y**2-90*x**5*m.log(x)-15*x**5,
        7 : 0,
        8 : 2*y*x,
        9 : -6*y*x*m.log(x)-3*y*x,
        10: 12*y*x**3-8*y**3*x,
        11: -120*y*x**3-160*y**3*x*m.log(x)-80*y**3*x+240*y*x**3*m.log(x),
        }[i]

def psixx(x,y,i):
    return {
        0 : 0,
        1 : 2.0,
        2 : 2*m.log(x)+3,
        3 : 12*x**2-8*y**2,
        4 : 36*x**2*m.log(x)+21*x**2-54*y**2-24*m.log(x)*y**2,
        5 : 30*x**4-144*x**2*y**2+16*y**4,
        6 : -640*y**4-240*m.log(x)*y**4+2160*x**2*m.log(x)*y**2+2160*x**2*y**2-450*x**4*m.log(x)-165*x**4,
        7 : 0,
        8 : 2*y,
        9 : -6*y*m.log(x)-9*y,
        10: 36*y*x**2-8*y**3,
        11: -120*y*x**2-160*y**3*m.log(x)-240*y**3+720*y*x**2*m.log(x),
        }[i]

def psiy(x,y,i):
    return {
        0 : 0,
        1 : 0,
        2 : -2*y,
        3 : -8*x**2*y,
        4 : -18*x**2*y-24*x**2*m.log(x)*y+8*y**3,
        5 : -24*x**4*y+32*x**2*y**3,
        6 : 48*y**5-560*x**2*y**3-480*x**2*m.log(x)*y**3+360*x**4*m.log(x)*y+150*x**4*y, 
        7 : 1,
        8 : x**2,
        9 : 3*y**2-3*x**2*m.log(x),
        10: 3*x**4-12*y**2*x**2,
        11: 40*y**4-45*x**4-240*y**2*x**2*m.log(x)+60*x**4*m.log(x) ,
        }[i]

def psiyy(x,y,i):
    return {
        0 : 0,
        1 : 0,
        2 : -2,
        3 : -8*x**2,
        4 : -18*x**2-24*x**2*m.log(x)+24*y**2, 
        5 : -24*x**4+96*x**2*y**2,
        6 : 240*y**4-1680*x**2*y**2-1440*x**2*m.log(x)*y**2+360*x**4*m.log(x)+150*x**4,
        7 : 0,
        8 : 0,
        9 : 6*y, 
        10: -24*y*x**2,   
        11: 160*y**3-480*y*x**2*m.log(x), 
        }[i]

def psipart(x,y,i):
    return {
        0 : 1.0/2*p(x,2)*l(x),
        1 : 1.0/8*p(x,4),
        }[i]

def psipartx(x,y,i):
    return {
        0 : x*m.log(x)+1.0/2*x,
        1 :  1.0/2*x**3
        }[i]

def psipartxx(x,y,i):
    return {
        0 : m.log(x)+3.0/2,
        1 : 3.0/2*x**2,
        }[i]

def psiparty(x,y,i):
    return 0

def psipartyy(x,y,i):
    return 0


def psiX(x,y,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,A):
    return c0*psix(x,y,0)+c1*psix(x,y,1)+c2*psix(x,y,2)+c3*psix(x,y,3)+c4*psix(x,y,4)+c5*psix(x,y,5)+c6*psix(x,y,6)+c7*psix(x,y,7)+c8*psix(x,y,8)+c9*psix(x,y,9)+c10*psix(x,y,10)+c11*psix(x,y,11)+A*psipartx(x,y,0)+(1-A)*psipartx(x,y,1) 

def psiXX(x,y,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,zeta):
    return 2*c1+c2*(2*m.log(x)+3)+c3*(12*x**2-8*y**2)+c4*(36*x**2*m.log(x)+21*x**2-54*y**2-24*m.log(x)*y**2)+c5*(30*x**4-144*x**2*y**2+16*y**4)+c6*(-640*y**4-240*m.log(x)*y**4+2160*x**2*m.log(x)*y**2+2160*x**2*y**2-450*x**4*m.log(x)-165*x**4)+cos(zeta)*m.log(x)+3.0/2*cos(zeta)+3.0/2*sin(zeta)*x**2  

def psiY(x,y,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,A):
    return c0*psiy(x,y,0)+c1*psiy(x,y,1)+c2*psiy(x,y,2)+c3*psiy(x,y,3)+c4*psiy(x,y,4)+c5*psiy(x,y,5)+c6*psiy(x,y,6)+c7*psiy(x,y,7)+c8*psiy(x,y,8)+c9*psiy(x,y,9)+c10*psiy(x,y,10)+c11*psiy(x,y,11)+A*psiparty(x,y,0)+(1-A)*psiparty(x,y,1)  

def psiYY(x,y,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,zeta):
    return -2*c2-8*c3*x**2+c4*(-18*x**2-24*x**2*m.log(x)+24*y**2)+c5*(-24*x**4+96*x**2*y**2)+c6*(240*y**4-1680*x**2*y**2-1440*x**2*m.log(x)*y**2+360*x**4*m.log(x)+150*x**4)  

def psi0(x,y,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,A):
    return c0*psi(x,y,0)+c1*psi(x,y,1)+c2*psi(x,y,2)+c3*psi(x,y,3)+c4*psi(x,y,4)+c5*psi(x,y,5)+c6*psi(x,y,6)+c7*psi(x,y,7)+c8*psi(x,y,8)+c9*psi(x,y,9)+c10*psi(x,y,10)+c11*psi(x,y,11)+A*psipart(x,y,0)+(1-A)*psipart(x,y,1)
