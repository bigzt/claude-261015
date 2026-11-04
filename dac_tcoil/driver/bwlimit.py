import numpy as np, sys, json
from multiprocessing import Pool
from scipy.optimize import differential_evolution
import chain
from optimize_driver import F
KEYS = ['n1','j1','Lg','n2','I2','r3','Rs','Cs','La','Lb','ko','L1d','L2d','kd','Ib','Lin']
B = [(60,600),(25,70),(0,300e-12),(100,2000),(4e-3,17e-3),(0.4,1.0),(0,12),(0,500e-15),
     (20e-12,400e-12),(20e-12,400e-12),(0,0.7),(60e-12,250e-12),(60e-12,300e-12),(0,0.6),(0,16e-3),(0,300e-12)]
def met(x,out,gt,Q=9.0):
    p=dict(zip(KEYS,x)); p['out']=out; dac=dict(chain.DAC); dac.update(L1=p['L1d'],L2=p['L2d'],k=p['kd'])
    g,cin,dev=chain.driver_only(p,F[:1],Q); gdb=20*np.log10(abs(g[0]))
    pad,_,_=chain.responses(p,F,Q,dac); db=20*np.log10(abs(pad)); r=db-db[0]
    i=np.where(r<-3)[0]; bw=F[i[0]]/1e9 if len(i) else F[-1]/1e9
    rip=r[F<=70e9].max()-r[F<=70e9].min(); it=2*(p['n1']*p['j1']*1e-6+p['I2']+p['Ib'])*1e3
    return gdb,bw,rip,it,cin[0]*1e15
def cost(x,out,gt,imax):
    try: gdb,bw,rip,it,cin=met(x,out,gt)
    except np.linalg.LinAlgError: return 1e6
    return -bw+50*max(0,abs(gdb-gt)-0.1)+30*max(0,rip-1.2)+2*max(0,it-imax)
def run(a):
    out,gt,imax=a
    r=differential_evolution(cost,B,args=(out,gt,imax),seed=3,maxiter=200,popsize=12,tol=1e-8,polish=False)
    return out,gt,imax,r.x.tolist(),met(r.x,out,gt)
if __name__=='__main__':
    jobs=[('tcoil_padC',6.02,120),('series_shunt',6.02,120),('tcoil_padC',3.52,120),('tcoil_padC',6.02,80)]
    with Pool(4) as pl: res=pl.map(run,jobs)
    json.dump([(o,g,i,x) for o,g,i,x,m in res],open('bwlimit.json','w'))
    for o,g,i,x,m in res:
        p=dict(zip(KEYS,x))
        print(f'{o:12s} gain={g:4.2f}dB Imax={i}: -> gain={m[0]:5.2f} BW3={m[1]:5.1f}G ripple={m[2]:4.2f} Itot={m[3]:5.1f}mA Cin={m[4]:4.1f}fF | Lin={p["Lin"]*1e12:.0f}p Ib={p["Ib"]*1e3:.1f}m n1={p["n1"]:.0f} n2={p["n2"]:.0f}')
