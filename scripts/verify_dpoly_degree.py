#!/usr/bin/env python3
"""Source-order degree-tail bound and finite common-finish rho transport.

No numerical solver samples occur here. Reuses the complete continuous seed
cover, original degree-18 true-Mills model and generic final-three Pair reserve.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, csv, hashlib, json, re, struct

if not __debug__:
    raise RuntimeError("Verification requires Python assertions enabled; do not use -O")

HERE = Path(__file__).parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-root',type=Path,default=Path(__file__).resolve().parents[1],help='implied-vol repository root; defaults to the checkout containing scripts/')
parser.add_argument('--reference-root',type=Path,required=True,help='explicit archived implied-black-volatility proof repository root')
parser.add_argument('--out',type=Path,help='result file; default source-root/target/proof/dpoly-degree-proof.json')
args=parser.parse_args()
SOURCE=args.source_root.resolve()
ORIGINAL=args.reference_root.resolve()
OUT=args.out or SOURCE/'target/proof/dpoly-degree-proof.json'
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
    "fn:dpoly": "0bb704f18f7da66de8acd974e322186831688c0890ea4dfceac8951b3f81f128",
    "fn:forward_D": "3140767e9507e8a00f9678497e01c096f1a85eae9dc9f629fdec3f8ccac6433a",
    "fn:finish": "2d014685a54bc1bc3b2369f8c274a7627d6dce8746c8a28bb51788b3537b9e82",
    "fn:tensor_seed": "8e5206a4e8a00778790f539055966bb4d1d5ad871ce528e98337631195652a39",
    "fn:rank_seed": "8e96c7466df6958a92bb80e090d621a3d960ce0126ea7538935b881279d0d813",
    "fn:rank_seed_shape": "4417ffca2d35f683b7acd4f151f9a63d50e456f84a4a6440a312e7f12dd5fbf4",
    "fn:wing_raw": "a78bc8e3c69eb3b736e96d89235878066045aa8db73880f996f2c220b8e82057",
    "fn:small_rank": "3831255a115c88740b85b4bdf6526616c4e90f2307ff4495f5f669d8d9d5d31b",
    "fn:finite_rank": "696fa13a06f9780838cb97f544cd52b6169992ae2d66322a6967f1da36110b9e",
    "src/experimental/small/constants.rs": "27803db098d55851d3c8e71357e72e8ea882e5d0c3a5a21729951edc5cd783b3",
    "src/experimental/lanes.rs": "b0cfa174bcddf7628d5f73c6abd60fa35a6ac15a1e6973b7e169a424cf4cef85",
    "type:Pair": "f678ef91ec0b58936b1700e7243ad65b90cb448eef1f10901b930ba08df6a20b"
  },
  "reference": {
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
    "results/global_proof/small_seed/core.csv": "ff6ba2203d608c04846824f574a41f1ccd736471a484d8e0fcd19d2c2e80f067"
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
SMALL_FILE_HASHES = {
    'dynamic_degree':'e20e0546690b1acb05e00cdec547c63ec39e047963fae64ca704385d58b460aa',
    'constant_degree_dispatch':'d3adc596e729f4e25a8721fb1e15a4895acb9ef4eade4f0c1e04e60b90b8bf53',
}
small_hash=hashlib.sha256(small_bytes).hexdigest()
assert small_hash in SMALL_FILE_HASHES.values(), ('whole small source binding mismatch',small_hash)
small_source=small_bytes.decode()
STATIC_DPOLY_HASH='38d12d548b6e933f5042482a4388f8dfa1d91a5659297ba965c85c9a86dfaa03'
STATIC_HELPER_HASH='5bf7063c437b580cb9ecd928d74b618f44af97aad44f0a9e08392a3a58d366c4'
source_bindings=EXPECTED_BINDINGS['source'].copy()
source_variant='dynamic_degree'
if hashlib.sha256(function_body(small_source,'dpoly').encode()).hexdigest()==STATIC_DPOLY_HASH:
    source_variant='constant_degree_dispatch'
    source_bindings['fn:dpoly']=STATIC_DPOLY_HASH
    source_bindings['fn:dpoly_degree']=STATIC_HELPER_HASH
assert small_hash==SMALL_FILE_HASHES[source_variant], ('source variant mismatch',source_variant,small_hash)
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
SEEDS = ORIGINAL/'results/global_proof/small_seed'
u, nu = F(1,2**53), F(1,2**700)
R, QL, T = F(1,4)+128*u, 2*u, F(1,256)
RAD = R+QL
TL, TH = 4*u*T+nu, (1+u)*T+nu
radius = RAD+F(1,16)

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

def arithmetic(k,n):
    # Exactly the nonnegative first-order source ledger in original dpoly,
    # with n in place of18. The tail excludes both x.low and tt.low; each
    # omission, both FMAs in ne, and the FMA in o are charged separately.
    stored=ROWS[k+12]; EH=abs(stored[n][0]); OH=EE=OE=F(0)
    high_bounds=[]
    for j in range(n-1,2,-1):
        inner=OH*TH+abs(stored[j][0])
        outer=EH*R+(1+u)*inner+nu
        en=RAD*EE+(TH+TL)*OE+EH*QL+OH*TL+u*inner+nu+u*outer+nu
        on=RAD*OE+EE+OH*QL+u*(OH*R+EH)+nu
        EH,OH=(1+u)*outer+nu,(1+u)*(OH*R+EH)+nu
        EE,OE=en,on
        assert EH<32 and OH<32
        high_bounds.append(max(EH,OH))
    # Final compensated stages propagate the incoming suffix errors.
    # Their secondary RN errors are exactly the inherited generic reserve.
    for j in range(2,-1,-1):
        EE,OE=RAD*EE+(TH+TL)*OE,EE+RAD*OE
    return OE,max(high_bounds,default=EH)

records=[]
DELTA=u/F(32)
for k in range(17):
    n=min(DEG[k+12]+1,18)
    if n==18:
        records.append({'k':k,'degree':n,'changed':False})
        continue
    # P18-Pn is an exact stored-coefficient polynomial. The divided
    # difference is its derivative somewhere between x-t and x+t;
    # radius=.25+130u+.0625 bounds that whole segment and indexing seams.
    omitted=sum(j*abs(ROWS[k+12][j][0])*radius**(j-1) for j in range(n+1,19))
    evaluation,high=arithmetic(k,n)
    # We deliberately do not subtract the old evaluator's positive budget.
    # New approximation/evaluation error <= old P18 model error + these.
    extra=omitted+evaluation
    assert extra<DELTA
    def short_upper(v,scale=2**24):
        return str(F(-((-v.numerator*scale)//v.denominator),scale))
    records.append({'k':k,'degree':n,'changed':True,'omitted_polynomial_tail_in_u_upper':short_upper(omitted/u),
                    'source_high_stage_arithmetic_in_u_upper':short_upper(evaluation/u),
                    'extra_D_error_in_u_upper':short_upper(extra/u),'extra_D_error_in_u_display':float(extra/u),
                    'high_stage_component_upper':short_upper(high,2**16)})

# Existing complete common-finish prerequisites, verified against all exact
# represented endpoints of the transported continuous seed-cover CSV files.
coverage=json.loads((SEEDS/'coverage.json').read_text())
baseline=json.loads((SEEDS/'finish.json').read_text())
rows_by_route={}
for route in ['wing','small','core']:
    file=SEEDS/(route+'.csv')
    assert hashlib.sha256(file.read_bytes()).hexdigest()==coverage['routes'][route]['csv_sha256']
    rows=list(csv.DictReader(file.open()));rows_by_route[route]=len(rows)
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

# Ideal HH term, geometry errors and final FMA are otherwise identical.
# raw=(u+(1+u)*(ideal+update/(1-r))), eta_true>=2u, r<1/50000.
rho_increment=(1+u)*update_lipschitz*dEn*DELTA/((1-F(1,50000))*2*u)
assert rho_increment<F(1,32)
route_bounds={}
for route in ['wing','small','core']:
    old=F(baseline['routes'][route]['max_rho']['exact'])
    new=old+F(1,32)
    assert new<1
    route_bounds[route]={'old_complete_rho_bound':str(old),'old_display':float(old),
                         'candidate_rho_bound':str(new),'candidate_display':float(new)}

dependencies=[ORIGINAL/'scripts/prove_small_finish.py',ORIGINAL/'scripts/prove_small_finish_reserves.py',
              ORIGINAL/'notes/global-rho-small-common-finish.md',SEEDS/'coverage.json',SEEDS/'finish.json',SEEDS/'finish_reserves.json']
paths=[SOURCE/'src/experimental/small.rs',SOURCE/'src/experimental/small/constants.rs',Path(__file__)]
result={'status':'All exact rational assertions passed',
 'candidate':'dpoly degree=min(mills_degree[cell]+1,18), original last three compensated stages',
 'scope':'Accepted common physical-price finish only: wing_raw, small_rank,finite_rank on their complete transported seed covers. Other forward_D branches, guards,seed expressions,exp and public adapters unchanged.',
 'method':'Exact stored-polynomial derivative tail plus source-order high-stage arithmetic; generic final-three compensated reserve transported unchanged; source residual normalization and finite HH error-budget sensitivity differentiated exactly; no finite accuracy sample used.',
 'primitive_contracts':['RN-even binary64 basic ops, correctly rounded FMA and sqrt, gradual underflow, no unsafe contraction/reassociation','exp relative error <=u=2^-53; inherited log/log1p preprocessing contracts'],
 'D_absolute_increment':str(DELTA),'max_extra_D_error_in_u':max(r.get('extra_D_error_in_u_display',0)for r in records),
 'residual_normalization_ED_derivative_upper':str(dEn),'residual_normalization_ED_derivative_display':float(dEn),
 'En_global_upper':str(En_upper),'finite_HH_update_lipschitz_upper':str(update_lipschitz),
 'finite_HH_update_lipschitz_display':float(update_lipschitz),
 'rho_increment_upper':str(rho_increment),'rho_increment_display':float(rho_increment),
 'coarse_rho_increment_used':str(F(1,32)),'route_bounds':route_bounds,'seed_rows':rows_by_route,'cell_bounds':records,
 'source_sha256':{str(p.relative_to(SOURCE)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths[:2]},
 'script_sha256':hashlib.sha256(paths[2].read_bytes()).hexdigest(),
 'source_variant':source_variant,
 'intended_source_bindings':source_bindings,
 'dependency_sha256':{str(p.relative_to(ORIGINAL)):hashlib.sha256(p.read_bytes()).hexdigest()for p in dependencies},
 'conditional_inherited_premises':['Imported continuous seed-cover, degree-18 true-Mills model, source residual ledger and generic final-three compensated reserve are established dependencies; this verifier binds their reviewed hashes and audits all CSV prerequisite endpoints, but does not regenerate those certificates.',
                                  'The reviewed baseline uses independent SIMD lanes with the reference high-stage scalar rounding order; all commonfinish source expressions except dpoly degree selection have fixed reviewed fingerprints.'],
 'limits':'Conditional continuous source-component proof under stated primitive contracts and the inherited complete seed/model/reserve certificates. It does not prove vendor libm accuracy, native Windows execution, or speed.'}
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items()if k not in ['cell_bounds','dependency_sha256','source_sha256']},indent=2))
