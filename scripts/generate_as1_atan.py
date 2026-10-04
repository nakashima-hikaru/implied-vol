#!/usr/bin/env python3
"""Numerical Remez candidate, independently certified by verify_as1_atan.py.

Requires mpmath only for fitting; the exact coefficient/error checker is stdlib.
This script does not certify global optimality of the fitted real polynomial.
"""
import argparse, json, math, struct
from fractions import Fraction
from pathlib import Path
import mpmath as mp


def fit(degree, dps):
    mp.mp.dps = dps
    qmax_exact = Fraction(1, 120) * (1 + Fraction(16, 2**53))
    qmax = mp.mpf(qmax_exact.numerator) / qmax_exact.denominator
    def value(x):
        q = (x + 1) * qmax / 2
        return mp.fsum([(-q)**j / (2*j+1) for j in range(100)])
    def derivative(x):
        q = (x + 1) * qmax / 2
        return qmax/2 * mp.fsum([(-1)**j*j*q**(j-1)/(2*j+1) for j in range(1,100)])
    def cheb(x):
        p = [mp.mpf(1), x]
        for _ in range(2,degree+1): p.append(2*x*p[-1]-p[-2])
        return p[:degree+1]
    def dp(c,x):
        if degree == 0: return mp.mpf(0)
        u = [mp.mpf(1), 2*x]
        for _ in range(2,degree): u.append(2*x*u[-1]-u[-2])
        return mp.fsum([i*c[i]*u[i-1] for i in range(1,degree+1)])
    nodes = [-mp.cos(mp.pi*k/(degree+1)) for k in range(degree+2)]
    history=[]
    for iteration in range(20):
        mat=mp.matrix([cheb(x)+[mp.mpf((-1)**k)] for k,x in enumerate(nodes)])
        solution=mp.lu_solve(mat,mp.matrix([value(x) for x in nodes]))
        c=list(solution[:degree+1]);error=solution[degree+1]
        roots=[]
        for x in nodes[1:-1]:
            y=mp.findroot(lambda z: dp(c,z)-derivative(z), (x-mp.mpf('.005'),x+mp.mpf('.005')),tol=mp.mpf(10)**(-dps+30),maxsteps=100)
            assert -1<y<1
            roots.append(y)
        newnodes=[mp.mpf(-1)]+sorted(roots)+[mp.mpf(1)]
        change=max(abs(a-b) for a,b in zip(nodes,newnodes))
        history.append({'iteration':iteration,'node_change':mp.nstr(change,25),'alternating_error':mp.nstr(error,50)})
        nodes=newnodes
        if change < mp.mpf(10)**(-dps//2): break
    else: raise RuntimeError('Remez exchange did not converge')
    mat=mp.matrix([cheb(x)+[mp.mpf((-1)**k)] for k,x in enumerate(nodes)])
    solution=mp.lu_solve(mat,mp.matrix([value(x) for x in nodes]));c=list(solution[:degree+1])
    def add(a,b):
        r=[mp.mpf(0)]*max(len(a),len(b))
        for i,v in enumerate(a):r[i]+=v
        for i,v in enumerate(b):r[i]+=v
        return r
    def mul(a,b):
        r=[mp.mpf(0)]*(len(a)+len(b)-1)
        for i,v in enumerate(a):
            for j,w in enumerate(b):r[i+j]+=v*w
        return r
    t=[[mp.mpf(1)],[mp.mpf(-1),2/qmax]]
    for i in range(2,degree+1):t.append(add([2*v for v in mul(t[1],t[-1])],[-v for v in t[-2]]))
    powers=[mp.mpf(0)]*(degree+1)
    for j,cj in enumerate(c):
        for i,v in enumerate(t[j]):powers[i]+=cj*v
    rounded=[float(x) for x in powers]
    bits=[struct.unpack('>Q',struct.pack('>d',x))[0] for x in rounded]
    return {'method':'High-precision Remez exchange in a shifted Chebyshev basis; numerical candidate, not an exact global-minimax certificate',
        'degree':degree,'mpmath_dps':dps,'q_fit_max_exact':str(qmax_exact),
        'real_coefficients_ascending':[mp.nstr(v,dps) for v in powers],
        'rounded_binary64_coefficients_ascending':[x.hex() for x in rounded],
        'coefficient_bits_ascending':[f'0x{x:016x}' for x in bits],
        'alternating_extrema_q':[mp.nstr((x+1)*qmax/2,60) for x in nodes],
        'alternating_error':mp.nstr(solution[degree+1],60),'exchange_history':history,
        'certificate':'Run verify_as1_atan.py; sampled extrema are fitting diagnostics only.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--degree',type=int,default=5);ap.add_argument('--dps',type=int,default=160);ap.add_argument('--output',type=Path);args=ap.parse_args()
    result=fit(args.degree,args.dps)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['degree','coefficient_bits_ascending','alternating_error','q_fit_max_exact']},indent=2))
if __name__=='__main__':main()
