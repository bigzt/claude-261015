"""Driver design optimisation (45 ohm on-chip load only, no external load).
Variables use current densities so every device stays inside the model range."""
import numpy as np, json, sys
from multiprocessing import Pool
from scipy.optimize import differential_evolution
import chain

F = np.concatenate([[0.1e9], np.arange(2e9, 122e9, 2e9)])
KEYS = ['n1','j1','Lg','j2','I3','j3','Rs','Cs','La','Lb','ko','L1d','L2d','kd','Ib','Lin']
B = [(40,600),(20,70),(0,300e-12),(15,70),(7e-3,16e-3),(15,70),(0,1e-6),(0,1e-21),
     (20e-12,400e-12),(20e-12,400e-12),(0,0.7),(60e-12,250e-12),(60e-12,300e-12),(0,0.45),(0,1e-9),(0,250e-12)]   # Rs=0, Cs=0, Ib=0, kd<=0.45

def to_p(x, out):
    v = dict(zip(KEYS, x))
    I2m = v['I3'] + v['Ib']; n2 = I2m*1e6/v['j2']; n3 = v['I3']*1e6/v['j3']
    p = dict(n1=v['n1'], j1=v['j1'], Lg=v['Lg'], n2=n2, I2=v['I3'], r3=n3/n2, Rs=v['Rs'], Cs=v['Cs'],
             La=v['La'], Lb=v['Lb'], ko=v['ko'], Ib=v['Ib'], Lin=v['Lin'], out=out)
    dac = dict(chain.DAC); dac.update(L1=v['L1d'], L2=v['L2d'], k=v['kd'])
    return p, dac

def met(x, out, Q=9.0, fix=chain.FIX, dac_override=None):
    p, dac = to_p(x, out)
    if dac_override: dac.update(dac_override)
    g, cin, dev = chain.driver_only(p, F[:1], Q, fix); gdb = 20*np.log10(abs(g[0]))
    pad, _, _ = chain.responses(p, F, Q, dac, fix); db = 20*np.log10(abs(pad)); r = db - db[0]
    i = np.where(r < -3)[0]; bw = F[i[0]]/1e9 if len(i) else F[-1]/1e9
    rip = r[F <= 70e9].max() - r[F <= 70e9].min()
    it = 2*(p['n1']*p['j1']*1e-6 + p['I2'] + p['Ib'])*1e3
    return dict(gdb=gdb, bw=bw, rip=rip, itot=it, cin=cin[0]*1e15, D=1+dev['m2']['gm']*p['Rs'], pk=r.max())

def cost(x, out, gt, bwt, Q=9.0):
    # nominal Q plus +-15 % (PVT spread of an EM-characterised coil)
    try: ms = [met(x, out, Q=q) for q in (Q, 0.85*Q, 1.15*Q)]
    except np.linalg.LinAlgError: return 1e6
    bw = min(m['bw'] for m in ms); rip = max(m['rip'] for m in ms)
    return ms[0]['itot'] + 50*max(0, abs(ms[0]['gdb']-gt)-0.1) + 8*max(0, bwt-bw) + 40*max(0, rip-1.3)

def run(a):
    out, gt, bwt, seed, Q = a
    r = differential_evolution(cost, B, args=(out, gt, bwt, Q), seed=seed, maxiter=220, popsize=12, tol=1e-10, polish=False)
    return a, r.x.tolist(), met(r.x, out, Q=Q)

if __name__ == "__main__":
    jobs = [('tcoil_padC', 6.02, 85, 11, 9.0), ('tcoil_padC', 6.02, 85, 23, 9.0), ('tcoil_padC', 6.02, 85, 11, 6.0), ('tcoil_padC', 6.02, 85, 11, 12.0)]
    with Pool(4) as pl: res = pl.map(run, jobs)
    json.dump([(a, x) for a, x, m in res], open('design4.json', 'w'))
    for a, x, m in res:
        v = dict(zip(KEYS, x))
        mq = [met(x, a[0], Q=q*a[4]) for q in (0.85, 1.15)]
        print(f"   Q-15%: BW3={mq[0]['bw']:.0f}G rip={mq[0]['rip']:.2f}dB | Q+15%: BW3={mq[1]['bw']:.0f}G rip={mq[1]['rip']:.2f}dB")
        print(f"Q={a[4]:4.1f} {a[0]:12s} seed{a[3]}: gain={10**(m['gdb']/20):.2f}x BW3={m['bw']:5.1f}G ripple={m['rip']:4.2f}dB Itot={m['itot']:5.1f}mA Cin={m['cin']:4.1f}fF D={m['D']:.2f} | "
              f"I1={v['n1']*v['j1']/1e3:.1f} I3={v['I3']*1e3:.1f} Ib={v['Ib']*1e3:.1f} mA/side j1/j2/j3={v['j1']:.0f}/{v['j2']:.0f}/{v['j3']:.0f}")
