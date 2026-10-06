#!/usr/bin/env python3
"""Verify the D0/odd-moment model and cubic common-finish precision transport.

Uses exact rational D0 comparisons, rigorous moment intervals and analytic
tails, source-order rounding, an exact cubic/rational identity, the seed cover,
and compensated reserves.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, csv, hashlib, json, os, re, struct, tomllib
from math import comb, factorial
from functools import lru_cache
import importlib.util, sys
sys.dont_write_bytecode = True

if not __debug__:
    raise RuntimeError("Verification requires Python assertions enabled; do not use -O")

HERE = Path(__file__).parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-root',type=Path,default=Path(__file__).resolve().parents[1],help='implied-vol repository root; defaults to the checkout containing scripts/')
parser.add_argument('--reference-root',type=Path,required=True,help='explicit archived implied-black-volatility proof repository root')
parser.add_argument('--simd-root',type=Path,help='fearless_simd 1.0.0 source directory; defaults to the Cargo registry cache')
parser.add_argument('--out',type=Path,help='result file; default source-root/target/proof/direct-d-proof.json')
args=parser.parse_args()
SOURCE=args.source_root.resolve()
ORIGINAL=args.reference_root.resolve()
OUT=args.out or SOURCE/'target/proof/direct-d-proof.json'
EXPECTED_BINDINGS = {
  "source": {
    "fn:fma": "81e7720498bb084c3eb478f27d9eb1a4fb485182561fd4060930ad675a28d7e6",
    "fn:add": "dd9e17d91093317c854c92ec0472dd05c0818ad941e1fd64f84424aea667530f",
    "fn:sub": "a313e54dd19b5d3f2a15d856ee92cd097810df1132f90d165ed76f8612a22c0c",
    "fn:mul": "319493b9c3680f1e531741b304b930ec06db33b930706c8ca901ddfa3491db59",
    "fn:div": "32e0892f20b1479d147743963a51d66cf59708d8995dd33179cb0c0393d71e96",
    "fn:scale": "6ac9a33f96ef818bb7034ad991ab5f3dd8bbecca0b51b7918983f5a51cc0cb6f",
    "fn:mills_index": "0b13019a2d9a4c38da9c0821478f57c7c47b4d6ef0258c5d5f1f8dbe2afe5cd9",
    "fn:single_mills": "51ed1db124d4054e33d693b8fa94ad8d1d6a1327edc5e96e2dd29f6094607984",
    "fn:dual_mills": "44baf834ccae9ce83592c5007be25177a573ff40aac1ea0e96d92fe7731b41e5",
    "fn:forward_D": "86ca29f9143353ca08b83c789fc797e58aae156e39c47553adfa9e995b0627e5",
    "fn:finish": "c7dc13828b7442d8d723b1f1ab85a29b9e0c1eb47d5d62e3c69f92c269ca9373",
    "fn:tensor_seed": "8e5206a4e8a00778790f539055966bb4d1d5ad871ce528e98337631195652a39",
    "fn:rank_seed": "8e96c7466df6958a92bb80e090d621a3d960ce0126ea7538935b881279d0d813",
    "fn:rank_seed_shape": "4417ffca2d35f683b7acd4f151f9a63d50e456f84a4a6440a312e7f12dd5fbf4",
    "fn:wing_raw": "a78bc8e3c69eb3b736e96d89235878066045aa8db73880f996f2c220b8e82057",
    "fn:small_rank": "3831255a115c88740b85b4bdf6526616c4e90f2307ff4495f5f669d8d9d5d31b",
    "fn:finite_rank": "696fa13a06f9780838cb97f544cd52b6169992ae2d66322a6967f1da36110b9e",
    "src/experimental/small/constants.rs": "27803db098d55851d3c8e71357e72e8ea882e5d0c3a5a21729951edc5cd783b3",
    "src/experimental/lanes.rs": "40502418e5f024c8a4a522b52e1a9b226e126906bf9aca14496fd56c336d6f32",
    "type:Pair": "3b00d0ba8421356c391f0340df8373b9fb7528bcf57b7395ce71ef881360e37b",
    "src/experimental/small/direct_d.rs": "fafdbd44da452d8f6e347ec4106253f2dd431e24fc27b2d261762236316e52f1"
  },
  "reference": {
    "scripts/prove_global_mills.py": "2dce1c923d5ba36d7109ffdc3a17f4c726e2996c2a50d6ee49663dfb436e7547",
    "scripts/prove_small_finish.py": "31377768cd0459857db2e0ac7d34350726defb03bd6e0ce7c7d9062bed2211b3",
    "scripts/prove_small_finish_reserves.py": "55642ef8855433a404c790b6fed3d2ecb0dc22933b774ce5a1985212a97dda5e",
    "notes/global-rho-small-common-finish.md": "131d57c8cc22ed1d78db5f99c78fa1c6eaeacc8037828c7392f7658194d46d1d",
    "src/small.rs": "5e74a96a7b1d2c76b7d528cd594c7e7abec6c7a5dc257929020de8dfbe584e44",
    "src/small/constants.rs": "f9efd22b8d0c88f080855f61ab5778dedbf077c9788348f5ab4f199d283b0a8d",
    "results/global_proof/small_seed/coverage.json": "1665f058e06d4372a8ef2cacd14dc0b073a3526aa17b00eca4f274d419601c1b",
    "results/global_proof/small_seed/finish.json": "5ad99fd68457c5a22291970d8b68a7483ac53b5bc65f9fbc897a6b2ba9bc5840",
    "results/global_proof/small_seed/finish_reserves.json": "8497571b0b83d7bbda0a4b385061dfcadd596cd5879a000ac8aa623d6d3fe773",
    "results/global_proof/small_seed/wing.csv": "16d108bbef32222fe84a128a8ed50c0089f20fcf7a12cd01182c34657f1f6e38",
    "results/global_proof/small_seed/small.csv": "fd3816e3c90a96b82ab2f58c9ce4772c430c8a174c7f3b0d1be38f7ee3730ff7",
    "results/global_proof/small_seed/core.csv": "ff6ba2203d608c04846824f574a41f1ccd736471a484d8e0fcd19d2c2e80f067",
    "reference/iv_large_a_cap_table_v50/scripts/exact.py": "cd66c85457348fabdc238e7395f6c4e10bc954f23a5af7ab709c43c08ac03559"
  }
}

def function_body(source,name):
    match=re.search(r'\bfn '+re.escape(name)+r'(?:\(|<)',source)
    assert match is not None, ('missing function',name)
    start=match.start();brace=source.index('{',start);depth=1;i=brace+1
    while depth:
        depth+=(source[i]=='{')-(source[i]=='}');i+=1
    return source[start:i]

small_bytes=(SOURCE/'src/experimental/small.rs').read_bytes()
assert hashlib.sha256(small_bytes).hexdigest()=='d92684d6227205e627dadd2d0a9a0f0735dec34bdaf2c960228f9e676ad39286', 'whole small source binding mismatch'
small_source=small_bytes.decode()
source_variant='direct_D_cubic_inverse'
assert 'let t = Pair::from(s * 0.5);' in function_body(small_source,'finish')
assert 't.hi <= 1.0 && (t.hi <= 0.5 || h.hi >= 2.0)' in function_body(small_source,'forward_D')
source_bindings=EXPECTED_BINDINGS['source']
for key,want in source_bindings.items():
    if key.startswith('fn:'):
        value=function_body(small_source,key[3:]).encode()
    elif key=='type:Pair':
        value=small_source[:small_source.index('fn fma(')].encode()
    else:
        value=(SOURCE/key).read_bytes()
    got=hashlib.sha256(value).hexdigest()
    assert got==want, ('source binding mismatch',key,got,want)
for key,want in EXPECTED_BINDINGS['reference'].items():
    got=hashlib.sha256((ORIGINAL/key).read_bytes()).hexdigest()
    assert got==want, ('inherited certificate binding mismatch',key,got,want)

# The adapter changes representation only: binary ops and negation delegate to
# lane-wise fearless_simd operations, with no reductions or reassociation.
# Crucially, mul_add_precise, rather than mul_add, preserves the single rounding
# used by compensated Horner. In the bound 1.0.0 source, Neon uses vfmaq_f64;
# SSE2, Fallback and WASM use f64::mul_add per lane.
# This binds the reviewed implementation; correct primitive rounding remains
# the existing conditional premise, not a new proof of a vendor/compiler FMA.
manifest=tomllib.loads((SOURCE/'Cargo.toml').read_text())
assert manifest['dependencies']['fearless_simd']['version']=='=1.0.0', 'SIMD dependency must match the reviewed version'
if args.simd_root:
    simd_root=args.simd_root.resolve()
else:
    cargo_cache=Path(os.environ.get('CARGO_HOME',Path.home()/'.cargo'))/'registry/src'
    candidates=sorted(cargo_cache.glob('*/fearless_simd-1.0.0'))
    assert len(candidates)==1, 'Pass --simd-root for the reviewed fearless_simd 1.0.0 sources'
    simd_root=candidates[0]
assert tomllib.loads((simd_root/'Cargo.toml').read_text())['package']['version']=='1.0.0'
simd_digest=hashlib.sha256()
simd_files=sorted((simd_root/'src').rglob('*.rs'))
assert simd_files, 'SIMD sources missing'
for file in simd_files:
    simd_digest.update(str(file.relative_to(simd_root)).encode()+b'\0')
    simd_digest.update(file.read_bytes()+b'\0')
SIMD_SOURCE_HASH='65450215ad42e632e66e97a0f2f8229e0934441b50b9fcca483c2cad80f5af75'
assert simd_digest.hexdigest()==SIMD_SOURCE_HASH, 'SIMD backend source binding mismatch'
SEEDS = ORIGINAL/'results/global_proof/small_seed'
u, nu = F(1,2**53), F(1,2**700)
R, QL, T = F(1,4)+128*u, 2*u, F(1,256)
RAD = R+QL
TL, TH = 4*u*T+nu, (1+u)*T+nu

def decode(v): return F(struct.unpack('<d',struct.pack('<Q',int(v,16)))[0])
def coeffs(path):
    source = path.read_text()
    s=source.split('static M_COEFF:',1)[1].split('=',1)[1].split('];',1)[0]
    vv=[decode(v) for v in re.findall(r'f64::from_bits\((0x[0-9a-f]+)\)',s)]
    assert len(vv)==2200
    rows=[[(vv[(k*25+j)*2],vv[(k*25+j)*2+1]) for j in range(25)] for k in range(44)]
    s=source.split('mills_degree:',1)[1].split('=',1)[1].split('];',1)[0]
    degrees=list(map(int,re.findall(r'\d+',s))); assert len(degrees)==44
    return rows,degrees
ROWS,DEG=coeffs(SOURCE/'src/experimental/small/constants.rs')
assert (ROWS,DEG)==coeffs(ORIGINAL/'src/small/constants.rs')

# Decode represented D0 coefficients. No fitting estimates enter this proof.
d0_source=(SOURCE/'src/experimental/small/direct_d.rs').read_text()
high_block=d0_source.split('static COEFF:',1)[1].split('static LOW:',1)[0].split('=',1)[1]
high=[decode(a) if a else F(0) for a,b in re.findall(r'f64::from_bits\((0x[0-9a-f]+)\)|(0\.0)',high_block)]
low=[decode(v) for v in re.findall(r'f64::from_bits\((0x[0-9a-f]+)\)',d0_source.split('static LOW:',1)[1].split('static DEGREE:',1)[0])]
degrees=list(map(int,re.findall(r'\d+',d0_source.split('static DEGREE:',1)[1].split('=',1)[1].split('];',1)[0])))
assert len(high)==17*16 and len(low)==51 and len(degrees)==17
d0_rows=[]
for k,n in enumerate(degrees):
    row=high[16*k:16*(k+1)]
    assert 3<=n<16 and all(v==0 for v in row[n+1:])
    assert all(abs(v)<16 for v in row) and all(abs(v)<16*u for v in low[3*k:3*k+3])
    d0_rows.append(row[:n+1])
exact_path=ORIGINAL/'reference/iv_large_a_cap_table_v50/scripts/exact.py'
spec=importlib.util.spec_from_file_location('odd_exact',exact_path)
ex=importlib.util.module_from_spec(spec);spec.loader.exec_module(ex);up=ex.ceilq
@lru_cache(maxsize=None)
def moment_intervals(h,n):
    m0=ex.mills(h);values=[m0,ex.isub((F(1),F(1)),ex.iscale(m0,h))]
    for j in range(1,n):values.append(ex.isub(ex.iscale(values[j-1],F(j)),ex.iscale(values[j],h)))
    assert all(0<a<=b for a,b in values)
    return values
def bernstein(poly,left,right):
    n=len(poly)-1
    shifted=[up(sum(poly[j]*comb(j,i)*left**(j-i)*(right-left)**i for j in range(i,n+1))) for i in range(n+1)]
    vals=[sum(shifted[j]*F(comb(i,j),comb(n,j)) for j in range(i+1)) for i in range(n+1)]
    return up(max(map(abs,vals))+F(n+1,2**400))
def enclosure(poly):return max(bernstein(poly,-RAD+2*RAD*F(i,8),-RAD+2*RAD*F(i+1,8)) for i in range(8))
S=F(2**30)*u*u
d0_records=[]
for k,row in enumerate(d0_rows):
    center=F(2*k+1,4);lower=F(k,2)-F(1,1000)
    exact=moment_intervals(center,31);lower_ms=moment_intervals(lower,31)
    degree=28
    truth=[ex.iscale(exact[i+1],F((-1)**i,factorial(i))) for i in range(degree+1)]
    polynomial=[v+(low[3*k+i] if i<3 else F(0)) for i,v in enumerate(row)]
    delta=[(polynomial[i] if i<len(polynomial) else 0)-(a+b)/2 for i,(a,b) in enumerate(truth)]
    model=up(enclosure(delta)+sum((b-a)/2*RAD**i for i,(a,b) in enumerate(truth))+lower_ms[degree+2][1]*RAD**(degree+1)/factorial(degree+1))
    # The identical high suffix and final three compensated stages satisfy
    # the inherited generic source reserve; the explicit primary suffix
    # error and centering low are propagated before that reserve is added.
    H=abs(row[-1]);E=F(0)
    for value in row[-2:2:-1]:
        mag=H*R+abs(value);E=up(RAD*E+H*QL+u*mag+nu);H=up((1+u)*mag+nu);assert H<32
    L=F(0)
    for i in range(2,-1,-1):
        E=up(RAD*E)
        # Absolute raw low bound for the source TwoSum/FMA compensation.
        L=up((1+u)*(R*L+u*H*R+u*((1+u)*H*R+abs(row[i]))+H*QL+abs(low[3*k+i]))+S)
        H=up((1+u)**2*(H*R+abs(row[i]))+4*nu)
        assert H<32
    assert L<256*u
    derivative=sum(i*abs(v)*RAD**(i-1) for i,v in enumerate(polynomial) if i)
    assert derivative<32
    d0_records.append({'model':model,'native':E+S,'raw_low':L,'raw_high':H})

# Represented reciprocal factorials are independently checked against exact
# integers. Source literals are binary64 values, not exact 1/n! constants.
series=function_body(d0_source,'odd_series')
reciprocals={1:F(1.0/6.0)}
for j,text in re.findall(r'cs\[(\d+)\]\s*=\s*cur\s*\*\s*([0-9.e+-]+);',series):reciprocals[int(j)]=F(float(text))
assert set(reciprocals)==set(range(1,16))
for j,value in reciprocals.items():assert abs(value-F(1,factorial(2*j+1)))<=u/factorial(2*j+1)

extended_records=[]
for regime,cells,count,t_squared in [('middle',range(17),6,F(1,256)),('quarter',range(17),9,F(1,16)),('half',range(17),12,F(1,4)),('three_quarter',range(4,17),14,F(9,16)),('one',range(4,17),16,F(1))]:
    for k in cells:
        D0=d0_records[k];ED0=D0['model']+D0['native']
        hhi=F(k+1,2)+128*u;hlo=max(F(0),F(k,2)-128*u)
        Hmax=hhi*hhi
        assert Hmax<80
        assert ED0<u/3, ('D0 model/rounding budget exceeded',regime,k,float(ED0/u))
        ms=moment_intervals(max(F(0),F(k,2)-F(1,1000)),2*count+1)
        I=[v[1] for v in ms]
        # H is the rounded h.high square. Hlow restores the exact square
        # up to O(u^2), including the input low. The ordinary H+3 correction
        # is a FastTwoSum when H>=3. Otherwise an explicit 4u bound covers
        # both operations, on low-h cells where FastTwoSum cannot be used.
        # QL bounds the RENORMALIZED centered argument in base(), not
        # the low of the raw h=a/s pair. Bound that division residual
        # separately: |h.low| <= (1+u)*u*|h.high|/(1-u) + nu.
        h_low=up((1+u)*u*hhi/(1-u)+nu)
        EH=up(u*Hmax+2*hhi*h_low+h_low*h_low+nu)
        beta_error=4*u if (1-u)*hlo*hlo-nu<3 else F(0)
        Eprev=ED0+D0['raw_low']
        Ecur=up((Hmax+3)*ED0+beta_error*I[1]+2*u*I[3]+S)
        reciprocal=reciprocals[1]
        cs_errors={1:up(reciprocal*Ecur+abs(reciprocal-F(1,6))*I[3]+u*reciprocal*(I[3]+Ecur)+nu)}
        for j in range(2,count):
            n=2*j-1;C=2*n+1;B=n*(n-1)
            EA=up(EH+u*((1+u)*Hmax+C+nu)+nu)
            A=up((1+u)*((1+u)*Hmax+C+nu)+nu)
            # Integer-times-previous is rounded before the fused update.
            Eproduct=up(u*B*(I[2*j-3]+Eprev)+nu)
            inside=up(A*Ecur+B*Eprev+EA*I[2*j-1]+Eproduct)
            Enext=up((1+u)*inside+u*I[2*j+1]+nu)
            assert I[2*j+1]+Enext<2**100
            Eprev,Ecur=Ecur,Enext
            reciprocal=reciprocals[j];true=F(1,factorial(2*j+1))
            cs_errors[j]=up(reciprocal*Ecur+abs(reciprocal-true)*I[2*j+1]+u*reciprocal*(I[2*j+1]+Ecur)+nu)
        # Tail uses positivity and C_(j+1)/C_j<=1/(2j+3) for h>=0.
        tail=up(I[2*count+1]*t_squared**count/factorial(2*count+1)/(1-t_squared/F(2*count+3)))
        # The source receives Pair::from(s/2), with exactly zero low. Every
        # direct-D call has t>0.001, hence scaling s by .5 is exact and
        # normal; the sole represented-square rounding has error <=u*T.
        local_TH=(1+u)*t_squared+nu;local_TL=u*t_squared+nu
        rest=Erest=F(0)
        for j in range(count-1,0,-1):
            coefficient=I[2*j+1]/factorial(2*j+1);ec=cs_errors[j]
            mag=rest*local_TH+coefficient+ec
            Erest=up(local_TH*Erest+rest*local_TL+ec+u*mag+nu)
            rest=up((1+u)*mag+nu)
        assert rest*local_TH<1
        total=up(ED0+tail+local_TH*Erest+rest*local_TL+u*local_TH*rest+nu+S)
        record={'regime':regime,'cell':k,'terms':count,'d0_model_in_u':float(D0['model']/u),'d0_native_in_u':float(D0['native']/u),'raw_low_in_u':float(D0['raw_low']/u),'tail_in_u':float(tail/u),'total_extra_in_u':float(total/u),'exact_total':str(total)}
        extended_records.append(record)
        assert total<u, ('active direct D budget exceeded',regime,k,float(total/u))
extended_by_cell={(row['regime'],row['cell']):F(row['exact_total']) for row in extended_records}
DELTA=max(F(row['exact_total']) for row in extended_records)

# Existing complete common-finish prerequisites, verified against all exact
# represented endpoints of the transported continuous seed-cover CSV files.
coverage=json.loads((SEEDS/'coverage.json').read_text())
baseline=json.loads((SEEDS/'finish.json').read_text())
rows_by_route={}
seed_rows_by_route={}
for route in ['wing','small','core']:
    file=SEEDS/(route+'.csv')
    assert hashlib.sha256(file.read_bytes()).hexdigest()==coverage['routes'][route]['csv_sha256']
    rows=list(csv.DictReader(file.open()));rows_by_route[route]=len(rows);seed_rows_by_route[route]=rows
    for row in rows:
        r={k:F(float(v)) for k,v in row.items() if k not in ['rect','depth','guard']}
        assert 0<=r['hlo']<9 and r['hhi']<9 and 0<=r['Tlo']<=r['Thi']<9
        assert F(1,100)<r['Dlo']<=r['Dhi']<2
        assert max(abs(r['nlo']),abs(r['nhi']))<F(1,50000)
        assert 0<r['radius']<F(1,50000)
    assert baseline['routes'][route]['failing_count']==0

# Full source residual normalization derivative with respect to its error
# bound ED. This follows the exact source ledger, retaining sv.low,D.low,
# small-residual FMA, scalar division and all underflow terms.
V=u+2**30*u*u
theta=(1+V)*(1+2*u)-1
assert theta<4*u
dback=1+V
dER=4*u*dback+32*u*u*(1+V)
dEn=(1+u)*(dback+dER)/(1-theta)
assert dEn<F(1001,1000)

# The complete old checker establishes ED<1e-10*Dlo and Nhat<=1/50000.
# It also establishes Ah<100,Bh<20000 and qh>.99. Here is an independent
# global source bound for En, enlarged for DELTA, to justify transport.
old_ED_upper=F(2,10**10)
back=2*V+(1+V)*(old_ED_upper+DELTA)
N=F(1,50000)
ER=4*u*(N+back)+32*u*u*(1+V)*(2+old_ED_upper+DELTA)+F(1,2**240)
En_upper=(back+ER+N*theta+u*(N+back+ER))/(1-theta)
assert En_upper<F(1,10**8)
assert N+dEn*DELTA<F(1,40000)

# Differentiate the entire update budget U(En)=Ln(N+En)*En+
# La(N+En)*EA+Lb(N+En)*EB+128u(N+En). These nonnegative rational
# majorants bound the finite HH partials on the original source guard.
z=F(1,40000); A=F(100); B=F(20000); q=F(99,100)
d=A+B*z/3; c=1+A*z/2
Ln=(1+A*z)/q+z*c*d/(q*q)
dLn=A/q+(1+A*z)*d/(q*q)+(c*d+z*A*d/2+z*c*B/3)/(q*q)+2*z*c*d*d/(q**3)
La=z*z*(F(1,2)+B*z*z/12)/(q*q)
dLa=(2*z*(F(1,2)+B*z*z/12)+B*z**3/6)/(q*q)+2*z*z*(F(1,2)+B*z*z/12)*d/(q**3)
Lb=z**3*c/(6*q*q)
dLb=(3*z*z*c+z**3*A/2)/(6*q*q)+z**3*c*d/(3*q**3)
EA=1000*u; EB=2**20*u
assert 8*u*90+2**30*u*u<EA
assert 64*u*91**2+2**30*u*u<EB
update_lipschitz=Ln+dLn*En_upper+dLa*EA+dLb*EB+128*u
assert update_lipschitz<2

# Compare both updates at the SAME computed n, A=hh.high-tt.high and
# Q=RN(RN(3*hh.high)+tt.high). Put B=A^2-Q in exact arithmetic. Sign
# symmetry of RN makes the old H3 exactly RN(A^2-Q). Thus the ideal
# cubic-minus-rational difference is
# n^4 * (A*(A^2+Q)/4 + B*(2*A^2+Q)*n/36) / (1+A*n+B*n^2/6).
# No price-series asymptotic remainder is dropped: this is a rational
# identity and transports the existing COMPLETE rational-update bound.
# Independently check that identity as an exact polynomial in A,Q,n.
def poly_add(*polys):
    out={}
    for poly in polys:
        for exponent,value in poly.items():out[exponent]=out.get(exponent,F(0))+value
    return {exponent:value for exponent,value in out.items() if value}
def poly_scale(poly,value):return {exponent:coefficient*value for exponent,coefficient in poly.items() if coefficient*value}
def poly_mul(left,right):
    out={}
    for a,x in left.items():
        for b,y in right.items():
            exponent=tuple(i+j for i,j in zip(a,b))
            out[exponent]=out.get(exponent,F(0))+x*y
    return {exponent:value for exponent,value in out.items() if value}
pA={(1,0,0):F(1)};pQ={(0,1,0):F(1)};pn={(0,0,1):F(1)};one={(0,0,0):F(1)}
pA2=poly_mul(pA,pA);pn2=poly_mul(pn,pn);pn3=poly_mul(pn2,pn);pn4=poly_mul(pn2,pn2)
pB=poly_add(pA2,poly_scale(pQ,-1))
pC=poly_add(poly_scale(pA2,2),pQ)
p_cubic=poly_add(pn,poly_scale(poly_mul(pA,pn2),F(-1,2)),poly_scale(poly_mul(pC,pn3),F(1,6)))
p_den=poly_add(one,poly_mul(pA,pn),poly_scale(poly_mul(pB,pn2),F(1,6)))
p_numerator=poly_add(pn,poly_scale(poly_mul(pA,pn2),F(1,2)))
p_difference=poly_mul(pn4,poly_add(poly_scale(poly_mul(pA,poly_add(pA2,pQ)),F(1,4)),poly_scale(poly_mul(poly_mul(pB,pC),pn),F(1,36))))
assert poly_add(poly_mul(p_cubic,p_den),poly_scale(p_numerator,-1))==p_difference

def cubic_transport(Nhat, Amag, Qmag):
    # Actual coefficient/evaluation ledger. Multiplication by 2 and -1/2
    # is exact in the reached normal range; nu covers underflow. Each
    # remaining scalar operation/FMA gets its own RN-even error bound.
    Bmag=Amag*Amag+Qmag
    den=1-Amag*Nhat-Bmag*Nhat*Nhat/6
    assert den>F(99,100) and Amag<100 and Qmag<300 and Nhat<F(1,40000)
    difference=Nhat**4*(Amag*Bmag/4+Bmag*(2*Amag*Amag+Qmag)*Nhat/36)/den
    square_error=u*Amag*Amag+nu
    sum_error=2*square_error+u*(2*(Amag*Amag+square_error)+Qmag)+3*nu
    C=(2*Amag*Amag+Qmag)/6
    coefficient_error=sum_error/6+u*(C+sum_error/6)+nu
    inner_error=Nhat*coefficient_error+u*(Amag/2+Nhat*(C+coefficient_error))+3*nu
    outer_error=Nhat*inner_error+u*(1+Nhat*(Amag/2+Nhat*C+inner_error))+nu
    evaluation_error=Nhat*outer_error+u*Nhat*(1+Nhat*(Amag/2+Nhat*C)+outer_error)+nu
    assert evaluation_error<=4*u*Nhat+100*nu
    # Bound the virtual old rational evaluation too. Its 128u*N reserve
    # is imported from the bound source ledger; H3 rounding is added
    # explicitly using the finite derivative with respect to B.
    H3_error=u*Bmag+nu
    derivative_B=Nhat**3*(1+Amag*Nhat/2)/(6*(den-H3_error*Nhat*Nhat/6)**2)
    round_difference=evaluation_error+128*u*Nhat+derivative_B*H3_error+100*nu
    return difference,round_difference

# Geometry, residual calculation and final FMA are shared. Transport the
# D error through the old rational update, then compare the two updates
# before the SAME final FMA, so its rounding is counted once.
# raw=(u+(1+u)*(ideal+update/(1-r))), eta_true>=2u, r<1/50000.
rho_increment=(1+u)*update_lipschitz*dEn*DELTA/((1-F(1,50000))*2*u)
# Diagnostic: no imposed uniform rho reserve.
# A lower bound on eta uses the inherited true-root separation. All
# accepted source geometry has H<65,T<31 and radius<1/50000, hence K<100;
# exp(K*r)<=1/(1-K*r), so (1+r)*exp(K*r)<1.003.
rr=F(1,50000)
assert F(65)/(1-rr)**3+31*(1+rr)<100
assert (1+rr)/(1-F(1,500))<F(1003,1000)
route_bounds={}
# Regenerate the complete old rational-finish bound ON EACH continuous
# seed leaf. The prior route-wide maximum is unnecessary for leaves that
# use the wider direct-D dispatch. The inherited scripts and every input
# certificate/represented coefficient are bound above before importing.
sys.path.insert(0,str(ORIGINAL/'scripts'))
spec=importlib.util.spec_from_file_location('bound_old_finish',ORIGINAL/'scripts/prove_small_finish.py')
old_finish=importlib.util.module_from_spec(spec);spec.loader.exec_module(old_finish)
assert old_finish.ROOT.resolve()==ORIGINAL
leaf_baseline={}
for route in ['wing','small','core']:
    values=[]
    for row in seed_rows_by_route[route]:
        parsed={k:(int(v) if k in ['rect','depth','guard'] else F(float(v))) for k,v in row.items()}
        records=old_finish.finish(parsed,1)
        values.append(max(r['rho'] for r in records) if records else F(0))
    assert max(values)==F(baseline['routes'][route]['max_rho']['exact'])
    leaf_baseline[route]=values
margin=4096*u
for route in ['wing','small','core']:
    old=F(baseline['routes'][route]['max_rho']['exact'])
    increments=[]
    leaf_totals=[]
    for row_index,row in enumerate(seed_rows_by_route[route]):
        hl,hh=F(float(row['hlo'])),F(float(row['hhi']))
        lo,hi=F(float(row['Tlo'])),F(float(row['Thi']))
        if hl>8+margin:continue
        cells=[k for k in range(17) if hl-margin<=F(k+1,2) and hh+margin>=F(k,2)]
        deltas=[F(0)]
        if lo<=F(1,256)+margin and hi>=F(float(.001))**2-margin:
            deltas.extend(extended_by_cell[('middle',k)] for k in cells)
        if lo<=F(1,16)+margin and hi>=F(1,256)-margin:
            deltas.extend(extended_by_cell[('quarter',k)] for k in cells)
        if lo<=F(1,4)+margin and hi>=F(1,16)-margin:
            deltas.extend(extended_by_cell[('half',k)] for k in cells if k>=0)
        if lo<=F(9,16)+margin and hi>=F(1,4)-margin:
            deltas.extend(extended_by_cell[('three_quarter',k)] for k in cells if k>=4)
        if lo<=1+margin and hi>=F(9,16)-margin:
            deltas.extend(extended_by_cell[('one',k)] for k in cells if k>=4)
        delta=max(deltas)
        assert delta<=DELTA and DELTA<u
        radius=F(float(row['radius']))
        Dlo=F(float(row['Dlo']))
        N=max(abs(F(float(row['nlo']))),abs(F(float(row['nhi']))))
        eta_lower=2*u*(1+(Dlo-N)/F(1003,1000))
        assert eta_lower>2*u
        increment=(1+u)*update_lipschitz*dEn*delta/((1-radius)*eta_lower)
        # h,T leaf bounds enclose exact seed geometry. The inherited
        # geometry reserve bounds computed A; 1000u also safely bounds Q.
        H=hh*hh
        Amag=max(abs(hl*hl-hi),abs(H-lo))+EA
        # mul() renormalizes its pair: hh.high differs from h^2 by
        # u*H plus the inherited compensated reserve. The multiplication
        # by 3 and addition in Q add 6u*H+u*T; tt.high adds u*T.
        # The audited h<9,T<9 hull is safely inside the 1000u reserve.
        assert u*(9*81+2*9)+8*S<EA
        Qmag=3*H+hi+EA
        ideal_difference,round_difference=cubic_transport(N+En_upper,Amag,Qmag)
        increment+=(1+u)*(ideal_difference+round_difference)/((1-radius)*eta_lower)
        increments.append(increment)
        leaf_totals.append(F(leaf_baseline[route][row_index])+increment)
    new=max(leaf_totals)
    extra=max(increments)
    assert new<1 and all(v<1 for v in leaf_totals), (route,float(new))
    route_bounds[route]={'failing_leaf_count':sum(v>=1 for v in leaf_totals),'old_complete_rho_bound':str(old),'old_display':float(old),
                         'rho_increment_upper':str(extra),'rho_increment_display':float(extra),
                         'maximum_rho_bound_change':str(new-old),'maximum_rho_bound_change_display':float(new-old),
                         'candidate_rho_bound':str(new),'candidate_display':float(new)}


paths=[SOURCE/'src/experimental/small.rs',SOURCE/'src/experimental/small/constants.rs',SOURCE/'src/experimental/small/direct_d.rs',Path(__file__)]
result={'status':'All exact rational assertions passed','candidate':'Leaf-conditioned direct-D error allocation with exact half-volatility argument and cubic inverse',
 'scope':'Common finish: cubic inverse on all accepted leaves; direct D for .001<t<=.5, and .5<t<=1 with h.hi>=2; all transported wing/small/core seed leaves. Conditional inherited seed/Mills/primitive contracts.',
 'method':'Rigorous moment intervals, Bernstein D0 model enclosures, source-order rounding bounds for compensated D0, initial I3, odd recurrence, reciprocal factorials, and t-square Horner with exact t=s/2 so |RN(t^2)-t^2|<=u*t^2+nu. Recomputed leaf-wise complete rational-update bounds plus finite-HH transport and exact rational cubic-update difference.',
 'max_extra_D_error_in_u':max(r['total_extra_in_u'] for r in extended_records),'cell_bounds':extended_records,
 'route_bounds':route_bounds,'seed_rows':rows_by_route,'D_absolute_increment':str(DELTA),'leaf_baseline_method':'Recomputed independently from source-bound archived prove_small_finish.py and prove_global_mills.py; each route maximum matches the bound archived finish.json exactly',
 'source_sha256':{str(p.relative_to(SOURCE)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths[:-1]},
 'script_sha256':hashlib.sha256(paths[-1].read_bytes()).hexdigest(),
 'primitive_contracts':['RN-even binary64 basic ops, correctly rounded FMA and sqrt, gradual underflow, preserved source order without unsafe contraction/reassociation',
                        'exp relative error <=u=2^-53; inherited log/log1p preprocessing contracts'],
 'conditional_inherited_premises':['The continuous seed cover, true-Mills model, residual ledger and final-three compensated reserve are imported certificates; this checker binds their reviewed hashes and audits all prerequisite endpoints, but does not regenerate those certificates.',
                                  'The D0 high suffix and three compensated stages satisfy the inherited coefficient, argument, high/low and centering bounds. The first-moment rounding reserve covers low cross products and recovery errors; all subsequent recurrence steps are bounded explicitly.'],
 'intended_source_bindings':source_bindings,
 'simd_backend':{'crate':'fearless_simd','version':'1.0.0','source_tree_sha256':SIMD_SOURCE_HASH,
                 'rounding':'Independent binary64 lanes; mul_add_precise with a single rounding; no lane reduction. Primitive correctness is conditional.'},
 'dependency_sha256':EXPECTED_BINDINGS['reference'],
 'limits':'Conditional continuous source-component proof, not a whole-domain solver certificate, end-to-end Lean proof, vendor-libm guarantee or speed claim.'}
OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'max_extra_D_error_in_u':result['max_extra_D_error_in_u'],
                  'common_finish_rho_upper':{route:row['candidate_display'] for route,row in route_bounds.items()},
                  'seed_rows':rows_by_route,'output':str(OUT)},indent=2))
