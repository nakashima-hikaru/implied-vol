"""Conditional continuous proof for reduced rank seeds and the cubic finish.

All input rows and kernels must first pass verify_direct_d.py source bindings.
The caller supplies freshly recomputed direct-D cell errors; no cached report
is an input. Reduced coefficients are independently regenerated exactly.
"""
from fractions import Fraction as F
from math import isqrt
import re,struct

if __name__ == '__main__':
 raise SystemExit('Run scripts/verify_direct_d.py --reference-root PATH for the integrated source-bound proof.')

def check_inverse_derivative_identity():
 # Exact polynomial check of the inverse fourth-derivative numerator.
 # A=H-T, Q=3H+T; s*dH/ds=-2H and s*dT/ds=2T.
 def add(*polys):
  out={}
  for poly in polys:
   for e,v in poly.items():out[e]=out.get(e,F(0))+v
  return {e:v for e,v in out.items() if v}
 def scale(poly,k):return {e:k*v for e,v in poly.items() if k*v}
 def mul(a,b):
  out={}
  for e,x in a.items():
   for f,y in b.items():
    k=(e[0]+f[0],e[1]+f[1]);out[k]=out.get(k,F(0))+x*y
  return {e:v for e,v in out.items() if v}
 def logarithmic_derivative(poly):return {e:(-2*e[0]+2*e[1])*v for e,v in poly.items() if (-2*e[0]+2*e[1])*v}
 H={(1,0):F(1)};T={(0,1):F(1)};A=add(H,scale(T,-1));Q=add(scale(H,3),T)
 A2=mul(A,A);A3=mul(A2,A)
 B3=add(A2,scale(A,-1),logarithmic_derivative(A))
 assert B3==add(A2,scale(Q,-1))
 B4=add(mul(A,B3),scale(B3,-2),logarithmic_derivative(B3))
 assert B4==add(A3,scale(mul(A,Q),-3),scale(H,12))
 inverse4=add(scale(A3,15),scale(mul(A,B3),-10),B4)
 assert inverse4==add(scale(A3,6),scale(mul(A,Q),7),scale(H,12))

check_inverse_derivative_identity()
u=F(1,2**53);nu=F(1,2**1075);BITS=160

def up(x):return F((x.numerator*(1<<BITS)+x.denominator-1)//x.denominator,1<<BITS)
def down(x):return -up(-x)
def sqrtlo(x):return F(isqrt((x.numerator<<(2*BITS))//x.denominator),1<<BITS)
def sqrthi(x):return sqrtlo(x)+F(1,1<<BITS)
def expup(x):
 assert 0<=x<F(1,2),('exp range',float(x))
 return up(1/(1-x))
def ia(a,b):return down(a[0]+b[0]),up(a[1]+b[1])
def im(a,b):
 v=[x*y for x in a for y in b];return down(min(v)),up(max(v))
def rn(a):
 e=up(u*max(abs(v) for v in a)+nu);return down(a[0]-e),up(a[1]+e)
def absmax(a):return max(abs(v) for v in a)
def horner(c,x):
 p=(c[-1],c[-1])
 for v in c[-2::-1]:p=ia(im(p,x),(v,v))
 return p

def decode(v):return F(struct.unpack('<d',struct.pack('<Q',int(v,16)))[0])
def table(source,name):
 body=source.split('static '+name+':',1)[1].split('=',1)[1].split('];',1)[0]
 return [decode(v) for v in re.findall(r'f64::from_bits\((0x[0-9a-f]+)\)',body)]
endpoint=1+32*u
T0,T1=F(1),endpoint
for i in range(2,11):T0,T1=T1,2*endpoint*T1-T0
chebguard=T1

def source_poly_error(c,xmag):
 # FMA Horner has at most ten rounded stages. Eleven covers rank2 y padding.
 g=11*u/(1-11*u)
 mag=sum(abs(v)*xmag**i for i,v in enumerate(c))
 return up(g*mag+64*nu),up(mag)

def seed_delta(row,route):
 k=int(row['rect']);name=f'lr_{route}_{k}';m=metadata[name];lanes=m['lanes'];half=lanes//2;rank=ranks[route][k];rect=rects[route][k]
 x=rn(im(rn((row['xl']-rect[0],row['xh']-rect[0])),(rect[1],rect[1])))
 y=rn(im(rn((row['yl']-rect[2],row['yh']-rect[2])),(rect[3],rect[3])))
 assert absmax(x)<=endpoint and absmax(y)<=endpoint
 old=oldtables[name];new=newtables[name];oldv=[];newv=[];olde=[];newe=[];d=[]
 for lane in range(lanes):
  arg=x if lane<half else y;oc=old[lane::lanes];nc=new[lane::lanes];field=m['fields'][lane]
  oldv.append(horner(oc,arg));newv.append(horner(nc,arg))
  olde.append(source_poly_error(oc,absmax(arg))[0]);newe.append(source_poly_error(nc,absmax(arg))[0])
  d.append(up(F(field['removed'])*chebguard+F(field['rounding'])*endpoint**10))
 difference=F(0);oldround=F(0);newround=F(0);oldmag=F(0);newmag=F(0);oldrange=(F(0),F(0))
 for k in range(rank):
  j=half+k;ox,oy=absmax(oldv[k]),absmax(oldv[j]);nx,ny=absmax(newv[k]),absmax(newv[j])
  difference+=d[k]*oy+d[j]*ox+d[k]*d[j]
  oldround+=olde[k]*oy+olde[j]*ox+olde[k]*olde[j]
  newround+=newe[k]*ny+newe[j]*nx+newe[k]*newe[j]
  oldmag+=(ox+olde[k])*(oy+olde[j]);newmag+=(nx+newe[k])*(ny+newe[j])
  oldrange=ia(oldrange,im(oldv[k],oldv[j]))
 g=rank*u/(1-rank*u)
 oldround=up(oldround+g*oldmag+16*nu);newround=up(newround+g*newmag+16*nu)
 delta=up(difference+oldround+newround)
 if route=='core':
  lower=2*sqrtlo(row['Tlo']);assert lower>0
  epsilon=up(delta/lower)
 else:
  # Archived seed proof establishes computed v=RN(q*factor)/q>1/64.
  lower=max(down(oldrange[0]-oldround),(F(1,64)-F(1,2**350))/(1+u))
  epsilon=up(((1+u)*delta/lower+2*u)/(1-u))
 return epsilon,delta,lower

def pmul(a,b):
 ma,ea=a;mb,eb=b
 return up(ma*mb),up(ma*eb+mb*ea+ea*eb+u*(ma+ea)*(mb+eb)+nu)
def padd(a,b):
 ma,ea=a;mb,eb=b
 return up(ma+mb),up(ea+eb+u*(ma+mb+ea+eb)+nu)
def pfma(a,b,c):
 ma,ea=a;mb,eb=b;mc,ec=c
 return up(ma*mb+mc),up(ma*eb+mb*ea+ea*eb+ec+u*((ma+ea)*(mb+eb)+mc+ec)+nu)
def pdiv(a,k):
 m,e=a;return up(m/k),up(e/k+u*(m+e)/k+nu)

def edbound(hl,hh,tl,th,Dhi):
 choices=[];margin=4096*u
 # The original source error routine supplies only unchanged tiny and Mills paths.
 if tl<=F(float(.001)):
  choices.extend(e for branch,e in oldfinish.forward_error(hl,hh,tl,min(th,F(float(.001))),Dhi) if branch=='moment')
 # Enlarged domain intersects present direct D regimes; source cell selection
 # uses h.hi, so all intersecting cells plus a 4096u geometry margin are included.
 cells=[k for k in range(17) if hl-margin<=F(k+1,2) and hh+margin>=F(k,2)]
 for regime,lo,hi in [('middle',F(float(.001)),F(1,16)),('quarter',F(1,16),F(1,4)),('half',F(1,4),F(1,2)),('three_quarter',F(1,2),F(3,4)),('one',F(3,4),F(1))]:
  if tl<=hi+margin and th>=lo-margin:
   cells2=[k for k in cells if hi<=F(1,2) or k>=4]
   if cells2:choices.append(max(DERROR[regime,k] for k in cells2)+2**30*u*u*(1+Dhi))
 # Fallback is active for t>1, or t>1/2 when h.hi<2.
 if th>=1 or (hl<=2+margin and th>=F(1,2)):
  tmin=max(tl,F(1,2));choices.extend(e for branch,e in oldfinish.forward_error(hl,hh,tmin,th,Dhi) if branch=='direct')
 assert choices
 return up(max(choices))

def finish_bound(row,epsilon):
 r0=row['radius'];N0=max(abs(row['nlo']),abs(row['nhi']));hl0=row['hlo'];hh0=row['hhi'];Tlo0=row['Tlo'];Thi0=row['Thi'];Dlo0=row['Dlo'];Dhi0=row['Dhi']
 assert epsilon<F(1,100)
 K0=hh0**2/(1-epsilon)**3+Thi0*(1+epsilon);vratio=expup(K0*epsilon)
 r=up((r0+epsilon)/(1-epsilon));N=up((N0*vratio+epsilon*vratio**2)/(1-epsilon))
 hl=down(hl0/(1+epsilon));hh=min(up(hh0/(1-epsilon)),8/(1-u))
 if hl>hh:return {'rho':F(0),'epsilon':epsilon,'skipped':True}
 Tlo=down(Tlo0*(1-epsilon)**2);Thi=up(Thi0*(1+epsilon)**2);tl=sqrtlo(Tlo);th=sqrthi(Thi)
 Dhi=up((Dhi0*vratio+epsilon*vratio**2)/(1-epsilon))
 Dlo=down((Dlo0/vratio-epsilon*vratio**2)/(1+epsilon))
 assert 0<=hl<=hh<9 and 0<=tl<=th<3
 assert hl-th>F(-3,2) and hh+th<9
 assert F(1,100)<Dlo<=Dhi<2
 assert 0<r<F(1,100) and 0<=N<F(1,100)
 H=up(hh*hh);T=Thi;A=max(abs(hl*hl-T),abs(H-Tlo));Q=up(3*H+T)
 Hstar=up(H/(1-r)**2);Tstar=up(T*(1+r)**2);Astar=max(abs(hl*hl/(1+r)**2-Tstar),abs(Hstar-Tlo*(1-r)**2));Qstar=up(3*Hstar+Tstar)
 K=up(H/(1-r)**3+T*(1+r))
 # For G=B^{-1}, G''''=-(6A^3+7AQ+12H)/(s^3 V^4).
 # Taylor's finite remainder at price B(s0), after division by s0,
 # is <=C4*N^4*exp(4Kr)/(1-r)^3 throughout the root segment.
 C4=up((6*Astar**3+7*Astar*Qstar+12*Hstar)/24)
 ideal=up(C4*N**4*expup(4*K*r)/(1-r)**3)
 ED=edbound(hl,hh,tl,th,Dhi)
 assert ED/Dlo<F(1,10**10)
 V=u+2**30*u*u;back=up(V*Dhi+(1+V)*ED)
 ER=up(4*u*(N+back)+32*u*u*(1+V)*(Dhi+ED)+F(1,2**240));theta=(1+V)*(1+2*u)-1
 En=up((back+ER+N*theta+u*(N+back+ER))/(1-theta))
 # Geometry reserve: compensated division/product highs, matching inherited
 # domain h<9,t<3. Extra 2^30u² bounds neglected DD cross terms.
 EH=up(4*u*H+2**30*u*u);ET=up(4*u*T+2**30*u*u)
 EA=up(EH+ET+u*(H+T+EH+ET)+nu)
 EQ=padd(pmul((F(3),F(0)),(H,EH)),(T,ET))[1]
 Ap=(A,EA);Qp=(Q,EQ);Hp=(H,EH)
 c2=(A/2,up(EA/2+nu))
 c3=pdiv(padd(pmul((F(2),F(0)),pmul(Ap,Ap)),Qp),F(6))
 np=(N,En)
 ds=pmul(np,pfma(np,pfma(np,c3,c2),(F(1),F(0))))
 assert ds[0]+ds[1]<F(1,2)
 raw=up(u+(1+u)*(ideal+ds[1])/(1-r))
 Kr0=(hh0**2/(1-r0)**3+Thi0*(1+r0))*r0
 eta=2*u*(1+(Dlo0-N0)/((1+r0)*expup(Kr0)))
 return {'rho':up(raw/eta),'epsilon':epsilon,'r':r,'n':N,'ideal_rho':up(ideal/((1-r)*eta)),'round_rho':up(ds[1]/((1-r)*eta)),'ED_u':up(ED/u),'hmax':hh,'tmax':th,'Dlo':Dlo,'Dhi':Dhi,'h_minus_t_min':hl-th,'h_plus_t_max':hh+th}


def audit(seed_rows, original_source, candidate_source, coefficient_metadata, old_finish, cell_bounds):
 global oldtables,newtables,metadata,rects,ranks,DERROR,oldfinish
 assert __debug__, 'Exact verification requires assertions'
 metadata=coefficient_metadata
 oldtables={name:table(original_source,name) for name in metadata}
 newtables={name:table(candidate_source,name) for name in metadata}
 for name,item in metadata.items():
  lanes=item['lanes'];x_degree=item['x_degree'];y_degree=item['y_degree']
  assert len(newtables[name])==(x_degree+1)*lanes
  assert 0<=y_degree<=x_degree<=10
  for lane in range(lanes//2,lanes):
   assert all(value==0 for value in newtables[name][lane::lanes][y_degree+1:])
 rects={route:[table(original_source,route+'_rects')[i:i+4] for i in range(0,len(table(original_source,route+'_rects')),4)] for route in ['small','core']}
 ranks={route:[int(v) for v in re.findall(r'\(&lr_'+route+r'_\d+, (\d+)\)',original_source)] for route in ['small','core']}
 DERROR={(r['regime'],r['cell']):F(r['exact_total']) for r in cell_bounds}
 oldfinish=old_finish
 # The archived seed positivity checks give R<9 for wing and factor>1/64
 # for small. q and a exceed 2^-333. These stronger floors survive epsilon
 # and preserve the inherited 2^-350 normal-scale premise.
 assert F(float(1e-100))>F(1,2**333)
 assert F(1,2**666)/(80*9)>F(1,2**700)
 assert F(99,100)*F(1,2**333)/64>F(1,2**350)
 reports={}
 for route in ['wing','small','core']:
  records=[]
  for i,raw in enumerate(seed_rows[route]):
   row={k:(int(v) if k in ['rect','depth','guard'] else F(float(v))) for k,v in raw.items()}
   epsilon=F(0) if route=='wing' else seed_delta(row,route)[0]
   if route=='core':assert (1-epsilon)*2*sqrtlo(row['Tlo'])>F(1,2**350)
   record=finish_bound(row,epsilon);record.update(index=i,rect=row['rect']);records.append(record)
   assert record['rho']<1, (route,i,float(record['rho']))
  worst=max(records,key=lambda r:r['rho'])
  print('Reduced seed',route,'leaves',len(records),'rho bound',float(worst['rho']),flush=True)
  reports[route]={'rows':len(records),'failing_leaf_count':0,'candidate_rho_bound':str(worst['rho']),'candidate_display':float(worst['rho']),
   'max_seed_relative_change':str(max(r['epsilon'] for r in records)),
   'max_root_radius':str(max(r.get('r',F(0)) for r in records)),
   'max_normalized_residual':str(max(r.get('n',F(0)) for r in records)),
   'worst':{k:str(v) if isinstance(v,F) else v for k,v in worst.items()},
   'rectangles':{str(k):{'x_degree':metadata[f'lr_{route}_{k}']['x_degree'],
                        'y_degree':metadata[f'lr_{route}_{k}']['y_degree'],
                        'maximum_rho_bound':str(max(r['rho'] for r in records if r['rect']==k))}
                 for k in sorted({r['rect'] for r in records})} if route!='wing' else {}}
 return {'status':'All represented-coefficient, source-order and continuous-leaf assertions passed',
  'method':'Source-bound per-table degrees; exact Chebyshev mode truncation; source-order seed perturbation; finite vega/geometry transport; fourth inverse derivative remainder; cubic source rounding ledger',
  'route_bounds':reports,'coefficient_bytes_before':sum(v['old_bytes'] for v in metadata.values()),
  'coefficient_bytes_after':sum(v['new_bytes'] for v in metadata.values()),
  'degree_metadata_bytes':2*len(metadata),
  'checked_reserve_domain':'0<=h<9, 0<=t<3, -1.5<h-t, h+t<9, .01<D<2, ED/D<1e-10 and seed>2^-350',
  'conditional_premises':'The hash-bound continuous old seed cover and primitive/Mills/compensated reserve contracts are inherited, not regenerated. The reduced approximation and enlarged source ledger are checked here.'}
