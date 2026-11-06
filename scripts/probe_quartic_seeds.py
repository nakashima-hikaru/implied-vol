#!/usr/bin/env python3
"""Bounded ideal-arithmetic probe of joint seed degree/rank and inverse order.

This does not change runtime code or prove an error budget. The grids evaluate
represented binary64 seed coefficients with high-precision arithmetic, then
apply ideal cubic/quartic inverse corrections to the exact Black residual.
Both mathematical grid inputs and rounded physical binary64 inputs are tested;
neither mode models runtime preprocessing, price, FMA or final rounding.

Small rectangle 1 tests 7/5 (current), 6/5 and 5/5. Optional --rank-two adds
bounded rank-two SVD candidates in small rectangles 0, 1 and 3. The search
and every ideal output are repeated at two precisions. This finite search
cannot authorize a runtime cutover without continuous seed/finish bounds and
actual binary64 validation under both FMA policies.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import struct
import time

import mpmath as mp
import seed_reduction as seeds


def exact(value):
    return mp.mpf(value.numerator)/value.denominator


def reduced_factors(original, rectangle, degrees):
    name = f'lr_small_{rectangle}'
    raw = seeds.represented_table(original,name)
    lanes = 4 if len(raw)==44 else 8
    factors = [seeds.reduce_factor(raw[k::lanes],degrees[k>=lanes//2])[0]
               for k in range(lanes)]
    return [factor[:degrees[0]+1] for factor in factors],lanes


def rank_two(factors,lanes,degrees,rank):
    flat = [factors[k][i] for i in range(degrees[0]+1) for k in range(lanes)]
    matrix = seeds.rank_polynomial(flat,lanes,rank)
    left,singular,right = mp.svd(mp.matrix([[exact(matrix[i][j])
                                           for j in range(degrees[1]+1)]
                                          for i in range(degrees[0]+1)]))
    result = []
    # Match the accepted rank-compression generation: freeze Chebyshev factors
    # once to binary64, exact power conversion, one final coefficient rounding.
    for side in range(2):
        for k in range(2):
            cheb = [F(float(mp.sqrt(singular[k])*(left[i,k] if side==0 else right[k,i])))
                    for i in range(degrees[side]+1)]
            result.append([F(float(v)) for v in seeds.chebyshev_to_power(cheb,degrees[0]+1)])
    return result


def black(a,s):
    h,t = a/s,s/2
    price = (mp.exp(-a/2)*mp.erfc((h-t)/mp.sqrt(2))
             -mp.exp(a/2)*mp.erfc((h+t)/mp.sqrt(2)))/2
    vega = mp.exp(-(h*h+t*t)/2)/mp.sqrt(2*mp.pi)
    return price,vega


def invert(a,b,initial):
    left,right = mp.mpf(0),mp.sqrt(2*a)+4
    while black(a,right)[0]<=b:
        right *= 2
    s = initial
    tolerance = mp.power(10,-mp.mp.dps+20)
    for _ in range(30):
        price,vega = black(a,s)
        correction = (price-b)/vega
        if abs(correction/s)<tolerance:
            assert abs((black(a,s)[0]-b)/(s*vega))<tolerance
            return s,vega
        if price<b:
            left = s
        else:
            right = s
        candidate = s-correction
        s = candidate if left<candidate<right else (left+right)/2
    raise ArithmeticError(('independent Black inversion did not converge',a,b))


def factor_seed(factors,rank,q,x,y):
    half = len(factors)//2
    values = [mp.polyval([exact(v) for v in factor[::-1]],x if k<half else y)
              for k,factor in enumerate(factors)]
    return q*sum(values[k]*values[half+k] for k in range(rank))


def candidates(original,rectangle,include_rank_two):
    current = seeds.DEGREES['small'][rectangle]
    degree_pairs = [(7,5),(6,5),(5,5)] if rectangle==1 else [current]
    result = []
    for degrees in degree_pairs:
        factors,lanes = reduced_factors(original,rectangle,degrees)
        rank = 4 if rectangle==1 else 3
        result.append({'label':f'rank{rank}-{degrees[0]}/{degrees[1]}',
                       'degrees':degrees,'rank':rank,'lanes':lanes,'factors':factors})
    if include_rank_two:
        factors,lanes = reduced_factors(original,rectangle,current)
        result.append({'label':f'rank2-{current[0]}/{current[1]}',
                       'degrees':current,'rank':2,'lanes':4,
                       'factors':rank_two(factors,lanes,current,4 if rectangle==1 else 3)})
    return result


def run(original,grid,precision,include_rank_two):
    rectangles = (0,1,3) if include_rank_two else (1,)
    result = {}
    with mp.workdps(precision):
        rects = seeds.represented_table(original,'small_rects')
        for rectangle in rectangles:
            models = candidates(original,rectangle,include_rank_two)
            center_x,scale_x,center_y,scale_y = map(exact,rects[4*rectangle:4*rectangle+4])
            for mode in ('mathematical','physical_binary64'):
                key = f'small_{rectangle}/{mode}'
                records = {}
                for model in models:
                    payload = b''.join(struct.pack('<d',float(v)) for factor in model['factors'] for v in factor)
                    records[model['label']] = {
                        'degrees':model['degrees'],'rank':model['rank'],
                        'two_lane_Horner_stages':(model['lanes']//4)*sum(model['degrees']),
                        'scalar_rank_FMAs':model['rank']-1,
                        'current_runtime_shape_supported':model['degrees'] in ([(6,5),(6,6),(7,5),(7,6)] if model['rank']>2 else [(6,5)]),
                        'coefficient_sha256':hashlib.sha256(payload).hexdigest(),
                        'maximum_seed_relative_error':mp.mpf(0),
                        'maximum_cubic_rho_J':mp.mpf(0),'maximum_quartic_rho_J':mp.mpf(0),
                        'cubic_grid_failures':0,'quartic_grid_failures':0,
                    }
                for i in range(grid):
                    # The original rectangle coordinates determine the box.
                    grid_x = -1+mp.mpf(2)*i/(grid-1)
                    grid_w = center_x+grid_x/scale_x
                    for j in range(grid):
                        x,w = grid_x,grid_w
                        y = -1+mp.mpf(2)*j/(grid-1)
                        q2 = center_y+y/scale_y
                        if q2<=0:
                            # Rectangle 0 has a singular zero-volatility corner.
                            # Sample a positive interior limit rather than q=0.
                            q2 = mp.mpf('1e-6')
                        q = mp.sqrt(q2)
                        b = q/(mp.sqrt(2*mp.pi)*mp.exp(w))
                        a = 2*b*mp.expm1(w)
                        if mode=='physical_binary64':
                            a,b = mp.mpf(float(a)),mp.mpf(float(b))
                            q = mp.sqrt(2*mp.pi)*(a/2+b)
                            w = mp.log1p(a/(2*b))
                            x,y = (w-center_x)*scale_x,(q*q-center_y)*scale_y
                        base = models[0]
                        initial = factor_seed(base['factors'],base['rank'],q,x,y)
                        root,vega = invert(a,b,initial)
                        eta = mp.power(2,-52)*(1+b/(root*vega))
                        for model in models:
                            record = records[model['label']]
                            s0 = factor_seed(model['factors'],model['rank'],q,x,y)
                            value,vega0 = black(a,s0)
                            n = (b-value)/(s0*vega0)
                            H,T = (a/s0)**2,s0*s0/4
                            A,Q = H-T,3*H+T
                            delta3 = n-A*n*n/2+(2*A*A+Q)*n**3/6
                            delta4 = delta3-(6*A**3+7*A*Q+12*H)*n**4/24
                            seed_error = abs((s0-root)/root)
                            cubic = abs((s0*(1+delta3)-root)/root)/eta
                            quartic = abs((s0*(1+delta4)-root)/root)/eta
                            record['maximum_seed_relative_error'] = max(record['maximum_seed_relative_error'],seed_error)
                            for name,rho in [('cubic',cubic),('quartic',quartic)]:
                                if rho>record[f'maximum_{name}_rho_J']:
                                    record[f'maximum_{name}_rho_J'] = rho
                                    record[f'worst_{name}_point'] = {'a':float(a),'b':float(b),'w':float(w),'q':float(q)}
                                record[f'{name}_grid_failures'] += int(rho>=1)
                result[key] = records
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--reference-root',type=Path,required=True)
    parser.add_argument('--grid',type=int,default=13)
    parser.add_argument('--dps',type=int,default=90)
    parser.add_argument('--check-dps',type=int,default=130)
    parser.add_argument('--rank-two',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    assert __debug__ and args.grid>=2 and args.check_dps>args.dps>=70
    started = time.monotonic()
    root = args.source_root.resolve()
    small,constants = (root/'src/experimental/small.rs').read_text(),(root/'src/experimental/small/constants.rs').read_text()
    _,original,_ = seeds.restore_and_validate(small,constants,(args.reference_root/'src/small/constants.rs').read_text())
    lower = run(original,args.grid,args.dps,args.rank_two)
    higher = run(original,args.grid,args.check_dps,args.rank_two)
    for key in lower:
        for label,record in lower[key].items():
            checked = higher[key][label]
            assert record['coefficient_sha256']==checked['coefficient_sha256'],('search precision changes candidate coefficients',key,label)
            for name in ('maximum_seed_relative_error','maximum_cubic_rho_J','maximum_quartic_rho_J'):
                assert abs(record[name]-checked[name]) <= mp.mpf('1e-20')*max(mp.mpf('1e-12'),abs(checked[name])),('precision disagreement',key,label,name)
    def serialize(value):
        if isinstance(value,mp.mpf):
            return float(value)
        raise TypeError(type(value))
    report = {'scope':'Finite ideal-arithmetic seed/inverse-series probe; represented binary64 seed coefficients, exact Black residual and high-precision corrections. No runtime preprocessing/FMA/price/final rounding, dispatcher qualification or continuous final-error proof.',
              'grid_per_rectangle_per_input_mode':args.grid**2,'precision_digits':[args.dps,args.check_dps],
              'source_sha256':{'src/experimental/small.rs':seeds.sha256(small),'src/experimental/small/constants.rs':seeds.sha256(constants)},
              'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'quartic_increment_cost_example':{'reuses':'A^2 already needed by cubic',
                    'coefficient':'c4=-fma(A,fma(.25,A^2,(7/24)*Q),.5*H)',
                    'additional_scalar_multiplies':2,'additional_scalar_FMAs_including_Horner':3,
                    'qualification':'Operation counts only. Coefficient constants, geometry/residual rounding and finite fifth inverse-derivative remainder need a new source-order ledger; CPU timing required.'},
              'results':higher,'elapsed_seconds':time.monotonic()-started}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,default=serialize)+'\n')
    print(json.dumps(report,indent=2,default=serialize))


if __name__=='__main__':
    main()
