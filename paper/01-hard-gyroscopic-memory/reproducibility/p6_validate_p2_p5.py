#!/usr/bin/env python3
"""Deterministic P6 validation for P2/P5.

Dependencies: numpy, scipy.
No manuscript files are read or modified.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp, quad
from scipy.linalg import svd
from scipy.interpolate import CubicSpline

SIGMA3=np.array([[1.,0.],[0.,-1.]])
SIGMA1=np.array([[0.,1.],[1.,0.]])
I2=np.eye(2)
NMINUS=np.block([[I2,-I2],[np.zeros((2,2)),I2]])

def k0_norm(n):
    x,w=leggauss(n); t=(x+1)/2; w=w/2
    d=t[:,None]-t[None,:]
    K=.5*d*np.abs(d)
    M=np.sqrt(w)[:,None]*K*np.sqrt(w)[None,:]
    return svd(M,compute_uv=False)[0]

def exact_omega_control(eps,a,b,rtol=3e-11,atol=1e-13):
    def rhs(t,y):
        Phi=y.reshape(4,4)
        A=eps*(float(a(t))*SIGMA3+float(b(t))*SIGMA1)
        M=np.block([[np.zeros((2,2)),I2],[A,np.zeros((2,2))]])
        return (M@Phi).ravel()
    sol=solve_ivp(rhs,(0.,1.),np.eye(4).ravel(),method="DOP853",rtol=rtol,atol=atol)
    Phi=sol.y[:,-1].reshape(4,4)
    S=NMINUS@Phi
    C=S[:2,:2]+S[2:,2:]
    O=.5*(C-C.T)
    return O[1,0]

def top_pair(n=220):
    x,w=leggauss(n); t=(x+1)/2; w=w/2
    d=t[:,None]-t[None,:]
    K=.5*d*np.abs(d)
    sw=np.sqrt(w)
    M=sw[:,None]*K*sw[None,:]
    U,S,Vh=svd(M)
    a=CubicSpline(t,U[:,0]/sw,bc_type="natural")
    b=CubicSpline(t,Vh[0,:]/sw,bc_type="natural")
    na=np.sqrt(quad(lambda z:a(z)**2,0,1,limit=300)[0])
    nb=np.sqrt(quad(lambda z:b(z)**2,0,1,limit=300)[0])
    return S[0],lambda z:a(z)/na,lambda z:b(z)/nb

def circular_profile(N=3,ellipticity=1.0,dc=False):
    om=2*np.pi*N
    def env(t): return np.sin(np.pi*t)**2
    def env1(t): return np.pi*np.sin(2*np.pi*t)
    def env2(t): return 2*np.pi**2*np.cos(2*np.pi*t)
    def p(t): return env(t)*np.cos(om*t)+(0.3*(3*t*t-2*t**3) if dc else 0.)
    def q(t): return ellipticity*env(t)*np.sin(om*t)+(-0.2*(3*t*t-2*t**3) if dc else 0.)
    def p1(t): return env1(t)*np.cos(om*t)-om*env(t)*np.sin(om*t)+(0.3*(6*t-6*t*t) if dc else 0.)
    def q1(t): return ellipticity*(env1(t)*np.sin(om*t)+om*env(t)*np.cos(om*t))+(-0.2*(6*t-6*t*t) if dc else 0.)
    def p2(t): return env2(t)*np.cos(om*t)-2*om*env1(t)*np.sin(om*t)-om**2*env(t)*np.cos(om*t)+(0.3*(6-12*t) if dc else 0.)
    def q2(t): return ellipticity*(env2(t)*np.sin(om*t)+2*om*env1(t)*np.cos(om*t)-om**2*env(t)*np.sin(om*t))+(-0.2*(6-12*t) if dc else 0.)
    return p,q,p1,q1,p2,q2

def area(profile):
    p,q,p1,q1,_,_=profile
    val=quad(lambda t:.5*(p(t)*q1(t)-q(t)*p1(t)),0,1,limit=400)[0]
    return val+.5*(p(1)*q(0)-p(0)*q(1))

def exact_omega_tt(eps,profile):
    *_,p2,q2=profile
    return exact_omega_control(eps,lambda t:.5*p2(t),lambda t:.5*q2(t))

def main():
    print("K0 convergence")
    for n in (50,100,200,300,500):
        print(n, f"{k0_norm(n):.13f}")
    s,a,b=top_pair()
    print("cOmega",f"{1/s:.11f}")
    print("sharp exact-Jacobi ratios")
    for eps in (.03,.01,.003,.001):
        om=exact_omega_control(eps,a,b)
        Q=2*eps**2
        print(eps,f"{Q/abs(om):.11f}")
    print("P5 area tests")
    for name,prof in (
        ("circular_N3",circular_profile(3,1.0,False)),
        ("elliptic_N3",circular_profile(3,.4,False)),
        ("dc_circular_N3",circular_profile(3,1.0,True)),
    ):
        ar=area(prof)
        om=exact_omega_tt(1e-3,prof)/1e-6
        print(name,f"area={ar:.12f}",f"omega/eps2={om:.12f}",f"diff={om-ar:.3e}")
    p,q,p1,q1,p2,q2=circular_profile(3,1.0,False)
    fixed=(p,lambda t:.6*p(t),p1,lambda t:.6*p1(t),p2,lambda t:.6*p2(t))
    print("fixed_axis",f"area={area(fixed):.3e}",f"omega={exact_omega_tt(1e-2,fixed):.3e}")
    H=1.; N=3
    analytic_area=3*np.pi*N*H**2/8
    analytic_Q=np.pi**4*H**2*(3*N**4+6*N**2+1)/2
    print("benchmark",f"area={analytic_area:.12f}",f"Q={analytic_Q:.12f}")

if __name__=="__main__":
    main()
