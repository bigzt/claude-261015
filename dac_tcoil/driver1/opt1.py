import numpy as np, json, sys
from multiprocessing import Pool
from scipy.optimize import differential_evolution
import stage1 as S
KEYS = ['n1', 'j1', 'j2', 'Lp', 'L1', 'L2', 'k', 'C0', 'L1d', 'L2d', 'kd', 'Lin']
B = [(60, 600), (20, 70), (20, 70), (0, 300e-12), (20e-12, 350e-12), (20e-12, 350e-12), (0, 0.7), (0, 30e-15),
     (60e-12, 250e-12), (60e-12, 300e-12), (0, 0.45), (0, 200e-12)]
KIND = 'n'
def P(x):
    p = dict(zip(KEYS, x)); p['I'] = p['n1']*p['j1']*1e-6; return p
def cost(x, bwt=85, ript=1.3):
    p = P(x)
    try:
        g = S.metrics_gain(p, KIND)
        c = 60*max(0, abs(g-2.0)-0.03)
        if c > 20: return 1e3 + c
        q = S.at_gain(p, 1.5, KIND)
        ms = [S.metrics(p, Q, kind=KIND) for Q in (9.0, 7.65, 10.35)] + [S.metrics(q, Q, kind=KIND) for Q in (9.0, 7.65, 10.35)]
    except np.linalg.LinAlgError: return 1e6
    bw = min(m['bw'] for m in ms); rip = max(m['rip'] for m in ms)
    return 2*p['I']*1e3 + c + 8*max(0, bwt-bw) + 40*max(0, rip-ript)
def run(seed):
    r = differential_evolution(cost, B, seed=seed, maxiter=160, popsize=12, tol=1e-9, polish=False)
    return seed, r.x.tolist()
if __name__ == "__main__":
    KIND = sys.argv[1] if len(sys.argv) > 1 else 'p'
    with Pool(4) as pl: res = pl.map(run, [11, 23, 37, 51])
    json.dump(res, open(f'opt1_{KIND}.json', 'w'))
    for seed, x in res:
        p = P(x); q = S.at_gain(p, 1.5, KIND); m2 = S.metrics(p, kind=KIND); m15 = S.metrics(q, kind=KIND)
        print(f"seed{seed}: I={p['I']*1e3:.1f}mA/side (gain2) {q['I']*1e3:.1f}mA (gain1.5) | gain2: BW={m2['bw']:.0f}G rip={m2['rip']:.2f} | gain1.5: BW={m15['bw']:.0f}G rip={m15['rip']:.2f} | Cin={m2['cin'][0]:.0f}fF")
