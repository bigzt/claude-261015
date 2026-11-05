import numpy as np, json, sys
from multiprocessing import Pool
from scipy.optimize import differential_evolution
import stage1 as S
from opt1 import KEYS, B
def P(x, pos):
    p = dict(zip(KEYS, x)); p['I'] = p['n1']*p['j1']*1e-6; p['padpos'] = pos; return p
def cost(x, pos):
    p = P(x, pos)
    try:
        g = S.metrics_gain(p, 'n'); c = 60*max(0, abs(g-2.0)-0.03)
        if c > 20: return 1e3 + c
        q = S.at_gain(p, 1.5, 'n'); ms = []
        for cp in (0.0, 10e-15, 25e-15):
            f = dict(S.FIX); f['Cpad'] = cp
            ms += [S.metrics(p, fix=f, kind='n'), S.metrics(q, fix=f, kind='n')]
        ms += [S.metrics(p, Q=7.65, kind='n'), S.metrics(p, Q=10.35, kind='n')]
    except np.linalg.LinAlgError: return 1e6
    bw = min(m['bw'] for m in ms); rip = max(m['rip'] for m in ms)
    return 2*p['I']*1e3 + c + 8*max(0, 85-bw) + 40*max(0, rip-1.3)
def run(a):
    pos, seed = a
    r = differential_evolution(cost, B, args=(pos,), seed=seed, maxiter=150, popsize=12, tol=1e-9, polish=False)
    return pos, seed, r.x.tolist()
if __name__ == "__main__":
    with Pool(4) as pl: res = pl.map(run, [('end', 11), ('end', 23), ('center', 11), ('center', 23)])
    json.dump(res, open('opt2.json', 'w'))
    for pos, seed, x in res:
        p = P(x, pos); q = S.at_gain(p, 1.5, 'n'); out = []
        for cp in (0.0, 10e-15, 25e-15):
            f = dict(S.FIX); f['Cpad'] = cp; a, b = S.metrics(p, fix=f, kind='n'), S.metrics(q, fix=f, kind='n')
            out.append(f"Cpad{cp*1e15:.0f}f: g2 {a['bw']:.0f}G/{a['rip']:.2f}dB g1.5 {b['bw']:.0f}G/{b['rip']:.2f}dB")
        print(f"{pos:6s} seed{seed}: I={p['I']*1e3:.1f}mA/side | " + ' | '.join(out) + f" | k={p['k']:.2f} L1={p['L1']*1e12:.0f}p L2={p['L2']*1e12:.0f}p Lp={p['Lp']*1e12:.0f}p")
