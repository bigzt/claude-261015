import numpy as np, sys, json
from multiprocessing import Pool
from scipy.optimize import differential_evolution
import chain

F = np.concatenate([[0.1e9], np.arange(2e9, 112e9, 2e9)])
KEYS = ['n1', 'j1', 'Lg', 'n2', 'I2', 'r3', 'Rs', 'Cs', 'La', 'Lb', 'ko', 'L1d', 'L2d', 'kd', 'Ib']
BOUNDS = [(100, 600), (25, 70), (0, 300e-12), (300, 2000), (8e-3, 17e-3), (0.4, 1.0), (0, 12), (0, 500e-15),
          (20e-12, 400e-12), (20e-12, 400e-12), (0, 0.7), (60e-12, 250e-12), (60e-12, 300e-12), (0, 0.6), (0, 16e-3)]

def unpack(x, out):
    p = dict(zip(KEYS, x)); p['out'] = out
    dac = dict(chain.DAC); dac.update(L1=p['L1d'], L2=p['L2d'], k=p['kd'])
    return p, dac

def metrics(x, out, Q=9.0, fix=chain.FIX):
    p, dac = unpack(x, out)
    g, cin, dev = chain.driver_only(p, F[:1], Q, fix)
    gdb = 20*np.log10(abs(g[0]))
    pad, d, _ = chain.responses(p, F, Q, dac, fix)
    db = 20*np.log10(abs(pad)); rel = db - db[0]
    i3 = np.where(rel < -3)[0]; bw = F[i3[0]]/1e9 if len(i3) else F[-1]/1e9
    band = F <= 70e9; rip = rel[band].max() - rel[band].min()
    D = 1 + dev['m2']['gm']*p['Rs']
    itot = 2*(p['n1']*p['j1']*1e-6 + p['I2'] + p['Ib'])*1e3
    return dict(gdb=gdb, bw=bw, rip=rip, D=D, itot=itot, cin=cin[0]*1e15, z0=db[0], pk=rel.max())

def cost(x, out):
    try: m = metrics(x, out)
    except np.linalg.LinAlgError: return 1e6
    return (m['itot'] + 50*max(0, abs(m['gdb']-6.02)-0.1) + 5*max(0, 85-m['bw'])
            + 30*max(0, m['rip']-1.2) + 50*max(0, 1.4-m['D']))

def run(out):
    r = differential_evolution(cost, BOUNDS, args=(out,), seed=7, maxiter=200, popsize=12, tol=1e-8, polish=False)
    return out, r.x.tolist(), metrics(r.x, out)

if __name__ == "__main__":
    outs = sys.argv[1:] or ['tcoil_padC', 'tcoil_drainC', 'series_shunt', 'shunt']
    with Pool(4) as pool: res = pool.map(run, outs)
    json.dump([(o, x) for o, x, m in res], open('opt_results.json', 'w'))
    for o, x, m in res:
        print(f"{o:13s} Itot={m['itot']:5.1f}mA gain={m['gdb']:5.2f}dB BW3={m['bw']:5.1f}G ripple={m['rip']:4.2f}dB D={m['D']:4.2f} Cin={m['cin']:4.1f}fF")
