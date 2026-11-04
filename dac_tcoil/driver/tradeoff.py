import numpy as np, json
from multiprocessing import Pool
from scipy.optimize import differential_evolution
import chain
from optimize_driver import F
from bwlimit import KEYS, B, met
def cost(x,out,gt,bwt):
    try: gdb,bw,rip,it,cin=met(x,out,gt)
    except np.linalg.LinAlgError: return 1e6
    return it+50*max(0,abs(gdb-gt)-0.1)+8*max(0,bwt-bw)+40*max(0,rip-1.2)
def run(a):
    out,gt,bwt=a
    r=differential_evolution(cost,B,args=(out,gt,bwt),seed=11,maxiter=250,popsize=14,tol=1e-9,polish=False)
    return out,gt,bwt,r.x.tolist(),met(r.x,out,gt)
if __name__=='__main__':
    jobs=[('tcoil_padC',3.52,85),('tcoil_padC',4.86,85),('tcoil_padC',6.02,85),('tcoil_padC',6.02,82)]
    with Pool(4) as pl: res=pl.map(run,jobs)
    json.dump([(o,g,b,x) for o,g,b,x,m in res],open('tradeoff.json','w'))
    for o,g,b,x,m in res:
        p=dict(zip(KEYS,x))
        print(f'target gain {10**(g/20):.2f}x BW>={b}: gain={10**(m[0]/20):.2f}x BW3={m[1]:5.1f}G ripple={m[2]:4.2f}dB Itot={m[3]:5.1f}mA ({m[3]*1.8:.0f}mW) Cin={m[4]:4.1f}fF | I1={p["n1"]*p["j1"]/1e3:.1f} I3={p["I2"]*1e3:.1f} Ib={p["Ib"]*1e3:.1f} mA/side')
