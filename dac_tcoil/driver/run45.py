import numpy as np, json
from multiprocessing import Pool
from scipy.optimize import differential_evolution
from bwlimit import KEYS, B, met, cost as cost_bw
from tradeoff import cost as cost_i
def run(a):
    mode,out,gt,t=a
    if mode=='bw': r=differential_evolution(cost_bw,B,args=(out,gt,t),seed=3,maxiter=220,popsize=12,tol=1e-9,polish=False)
    else:          r=differential_evolution(cost_i, B,args=(out,gt,t),seed=11,maxiter=250,popsize=14,tol=1e-9,polish=False)
    return a,r.x.tolist(),met(r.x,out,gt)
if __name__=='__main__':
    jobs=[('bw','tcoil_padC',6.02,150),('bw','tcoil_drainC',6.02,150),('bw','series_shunt',6.02,150),('imin','tcoil_padC',6.02,85)]
    with Pool(4) as pl: res=pl.map(run,jobs)
    json.dump([(a,x) for a,x,m in res],open('run45.json','w'))
    for a,x,m in res:
        p=dict(zip(KEYS,x))
        print(f'{a[0]:4s} {a[1]:12s}: gain={10**(m[0]/20):.2f}x BW3={m[1]:5.1f}G ripple={m[2]:4.2f}dB Itot={m[3]:5.1f}mA Cin={m[4]:4.1f}fF | I1={p["n1"]*p["j1"]/1e3:.1f} I3={p["I2"]*1e3:.1f} Ib={p["Ib"]*1e3:.1f} mA/side n1={p["n1"]:.0f} n2={p["n2"]:.0f}')
