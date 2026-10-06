#!/usr/bin/env python3
"""Exact rational certificate for AS1's rounded degree-five atan polynomial.

Only the standard library is required. Numerical Remez fitting is optional and
is not used as proof. This checker proves approximation of atan(sqrt(q))/sqrt(q)
with an alternating integral series and Bernstein hulls of an exact polynomial.
It separately bounds the five source FMA roundings under RN-even/gradual
underflow; the full AS1 correction remains a separate conditional certificate.
"""
import argparse, hashlib, json, math, re, struct
from fractions import Fraction as F
from pathlib import Path

if not __debug__:
    raise RuntimeError('Verification requires Python assertions enabled; do not use -O')

U = F(1, 2**53)
FIT_MAX = F(1, 120) * (1 + 16*U)
PROOF_MAX = F(1, 50)
BITS = (0x3ff0000000000000, 0xbfd5555555554dde, 0x3fc999999947d463,
        0xbfc24923ee7916eb, 0x3fbc709e301f27e7, 0xbfb6c9831b2f66d6)
COEFFICIENTS = tuple(F.from_float(struct.unpack('>d',struct.pack('>Q',x))[0]) for x in BITS)
SERIES_DEGREE = 20
ARITHMETIC_HASHES = {
    'atan_rel_squared': 'a85891e0b81cb70468e6190e9a4bb8d7a255dfba1e508ca8530a5128cd6a309f',
    'slope': 'ebd3e15e7d9a967ebc10cf84599137f199e151303bb07b8d2c9a33d8a2dbab3c',
    'balanced': 'f80ca28bc60d08471384b78ae8bf928b2c36d7747b1e02a4dba69deda93975d2',
}
FMA_BODY_HASH = '0594d16475d930580c8097c79c0cb1dc819b28dd3bdb21e325d2e212ef2cb1b6'
CONSTANTS_HASH = '27803db098d55851d3c8e71357e72e8ea882e5d0c3a5a21729951edc5cd783b3'
CENTRAL_FILE_HASH = '61e04489aa1849afdbabe74273685f0b2062c436035c77c0d4bcae648f5a642f'
SMALL_FILE_HASHES = {
    'd0_cubic_inverse': 'd92684d6227205e627dadd2d0a9a0f0735dec34bdaf2c960228f9e676ad39286',
}


def represented_coefficients():
    return list(COEFFICIENTS)


def bernstein_bounds(p, lo, hi):
    """Exact power-to-Bernstein conversion on a closed rational interval."""
    n=len(p)-1; width=hi-lo
    shifted=[width**j * sum((p[i]*math.comb(i,j)*lo**(i-j)
                for i in range(j,n+1)),F(0)) for j in range(n+1)]
    bern=[sum((shifted[j]*F(math.comb(k,j),math.comb(n,j))
               for j in range(k+1)),F(0)) for k in range(n+1)]
    return min(bern),max(bern)


def source_q_bound():
    """Bound reached source q, including every slope rounding before Horner.

    For x=h-t and k=h²-t², max_z 4t²z/(k+z)²=t²/k.
    The rounded h-t>=7 guard implies x>=7/(1+u), while t<=RN(.7).
    kk=fma(h,h,-RN(t*t)); A=RN(kk+z); inv=RN(1/A);
    q=RN(RN(4*RN(t*t)*z)*RN(inv*inv)). All these values are normal in
    the reached AS1 domain (a>=RN(1e-100), h<=40, positive stored atoms).
    """
    tmax=F.from_float(.7); xmin=F(7)/(1+U)
    ratio=tmax*tmax/(xmin*xmin+2*xmin*tmax)
    computed=ratio*(1+U)**6/((1-U)**4*(1-U*ratio)**2)
    assert computed<FIT_MAX
    return {'status':'exact_rational_guard_and_operation_bound',
        'real_maximum_t_squared_over_k':str(ratio),
        'rounded_source_q_upper':str(computed),
        'certified_fit_domain_upper':str(FIT_MAX),
        'normality_assumption':'Actual AS1 inputs a>=RN(1e-100), h<=40 and unchanged positive normal quadrature atoms; RN-even normal-operation bounds.',
        'inequality':'q_source <= [tmax²/(xmin²+2*xmin*tmax)]*(1+u)^6/((1-u)^4*(1-u*ratio)^2) < (1/120)*(1+16u)'}


def approximation_table(cells=256):
    """Exact represented polynomial minus the Taylor-20 alternating series.

    A(q)=integral_0^1 1/(1+q*x²) dx. The finite geometric identity
    proves |A(q)-sum_0^20 (-q)^j/(2j+1)| <= q^21/43 for q>=0.
    Bernstein basis functions are nonnegative and sum to one, so coefficient
    hulls enclose every real point. No floating samples enter this certificate.
    """
    p=[F(0)]*(SERIES_DEGREE+1)
    for j,c in enumerate(COEFFICIENTS):p[j]+=c
    for j in range(SERIES_DEGREE+1):p[j]-=F((-1)**j,2*j+1)
    cuts=sorted(set([PROOF_MAX*F(j,cells) for j in range(cells+1)]+[FIT_MAX]))
    rows=[];prefix=F(0)
    for lo,hi in zip(cuts,cuts[1:]):
        lower,upper=bernstein_bounds(p,lo,hi)
        rem=hi**(SERIES_DEGREE+1)/(2*SERIES_DEGREE+3)
        error=max(abs(lower),abs(upper))+rem
        prefix=max(prefix,error)
        rows.append({'lo':str(lo),'hi':str(hi),'error_upper':str(error),
                     'prefix_error_upper':str(prefix)})
    return rows


def polynomial_error_bound(qmax, table):
    qmax=F(qmax);assert 0<=qmax<=PROOF_MAX
    if qmax==0:return abs(COEFFICIENTS[0]-1)
    for row in table:
        if qmax<=F(row['hi']):return F(row['prefix_error_upper'])
    raise AssertionError('No enclosing q interval')


def horner_rounding_bound(qmax=FIT_MAX):
    """Five-FMA absolute error for exact binary64 argument q in [0,qmax]."""
    qmax=F(qmax);lo=hi=COEFFICIENTS[-1];error=F(0);stages=[]
    for c in COEFFICIENTS[-2::-1]:
        products=[lo*0,lo*qmax,hi*0,hi*qmax]
        lo,hi=min(products)+c,max(products)+c
        propagated=error*qmax
        error=propagated+U*(max(abs(lo),abs(hi))+propagated)
        stages.append({'exact_stage_lo':str(lo),'exact_stage_hi':str(hi),
                       'absolute_rounding_error_upper':str(error)})
        # Every stage stays safely normal; underflow/overflow are unreachable.
        assert min(abs(lo),abs(hi))>F(1,100) and lo*hi>0
    return error,stages


def verify_source(path):
    # Fail closed on the whole arithmetic bodies, including slope's rounding
    # sequence and balanced's h/t construction and range guards. Only comments
    # and whitespace are ignored; matching coefficients alone is insufficient.
    central_bytes=path.read_bytes()
    got=hashlib.sha256(central_bytes).hexdigest()
    assert got==CENTRAL_FILE_HASH,('whole central source binding',got,CENTRAL_FILE_HASH)
    small_bytes=(path.parent.parent/'small.rs').read_bytes()
    got=hashlib.sha256(small_bytes).hexdigest()
    assert got in SMALL_FILE_HASHES.values(),('whole FMA-provider source binding',got)
    text=re.sub(r'/\*.*?\*/|//[^\n]*','',central_bytes.decode(),flags=re.S)
    for name,want in ARITHMETIC_HASHES.items():
        start=text.index('fn '+name+'(');at=text.index('{',start);depth=1;end=at+1
        while depth:
            depth+=(text[end]=='{')-(text[end]=='}');end+=1
        body=re.sub(r'\s+','',text[start:end])
        got=hashlib.sha256(body.encode()).hexdigest()
        assert got==want,(name,got,want)
    # The source's fma name must implement the explicit correctly rounded
    # mul_add modeled here; matching only its callers cannot establish that.
    helpers=re.sub(r'/\*.*?\*/|//[^\n]*','',small_bytes.decode(),flags=re.S)
    start=helpers.index('fn fma(');at=helpers.index('{',start);depth=1;end=at+1
    while depth:
        depth+=(helpers[end]=='{')-(helpers[end]=='}');end+=1
    got=hashlib.sha256(re.sub(r'\s+','',helpers[start:end]).encode()).hexdigest()
    assert got==FMA_BODY_HASH,('fma helper',got,FMA_BODY_HASH)
    constants=(path.parent/'constants.rs').read_bytes()
    got=hashlib.sha256(constants).hexdigest()
    assert got==CONSTANTS_HASH,('quadrature constants',got,CONSTANTS_HASH)
    # Close the declared atom-normality premise as well as binding its bits.
    # With source a>=RN1e-100 and h<=40, t>2^-341. These positive atoms,
    # below2^30, keep every intermediate in the q formula finite and normal;
    # the smallest q exceeds2^-800, including the elementary roundings.
    constants=constants.decode()
    for name in ['RULE1','RULE2','RULE3']:
        section=constants.split('static '+name+':',1)[1].split('=',1)[1].split('];',1)[0]
        bits=re.findall(r'f64::from_bits\((0x[0-9a-f]+)\)',section)
        assert len(bits)%4==0 and bits,('quadrature layout',name)
        atoms=[F.from_float(struct.unpack('>d',struct.pack('>Q',int(v,16)))[0]) for v in bits[::4]]
        assert all(F(1,2**50)<z<F(2**30) for z in atoms),('quadrature normality',name)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);ap.add_argument('--central',type=Path,default=Path(__file__).resolve().parents[1]/'src/experimental/small/central.rs');ap.add_argument('--cells',type=int,default=256);args=ap.parse_args()
    if not args.central.is_file():ap.error('central source not found; pass --central PATH')
    table=approximation_table(args.cells);fit_error=polynomial_error_bound(FIT_MAX,table)
    rounding,stages=horner_rounding_bound();result={
        'status':'exact_rational_uniform_polynomial_and_five_FMA_certificate',
        'scope':'atan helper only; full source AS1 rho composition is separate; no exact-minimax assertion.',
        'coefficient_bits_ascending':[f'0x{x:016x}' for x in BITS],
        'rounded_coefficients_exact':[str(c) for c in COEFFICIENTS],
        'source_q_bound':source_q_bound(),'proof_hull_domain':['0',str(PROOF_MAX)],
        'series_degree':SERIES_DEGREE,'Bernstein_cells':len(table),
        'fit_domain_rounded_coefficient_approximation_error_upper':str(fit_error),
        'fit_domain_rounded_coefficient_approximation_error_upper_float':float(fit_error),
        'fit_domain_five_FMA_rounding_error_upper':str(rounding),
        'fit_domain_five_FMA_rounding_error_upper_float':float(rounding),
        'fit_domain_total_atan_error_upper':str(fit_error+rounding),
        'fit_domain_total_atan_error_upper_float':float(fit_error+rounding),
        'rounded_Horner_stages':stages,'approximation_prefix_table':table,
        'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'arithmetic_body_sha256':ARITHMETIC_HASHES,
        'source_fma_body_sha256':FMA_BODY_HASH,
        'source_constants_sha256':CONSTANTS_HASH,
        'required_central_file_sha256':CENTRAL_FILE_HASH,
        'supported_small_file_sha256':SMALL_FILE_HASHES,
        'central_sha256':verify_source(args.central)}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['rounded_coefficients_exact','rounded_Horner_stages','approximation_prefix_table','source_q_bound']},indent=2))
if __name__=='__main__':main()
