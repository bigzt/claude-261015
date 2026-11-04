"""Step C2: PAM4 transient through the behavioural full-differential chain."""
import subprocess, numpy as np, sys, itertools, time
sys.path.insert(0, '..')
import gen
S = '/tmp/claude-0/-home-user-claude-261015/059e4729-0ea3-5859-88c1-9d51e0756740/scratchpad'
IFS = 5.2e-3
SAVE = ['padp', 'padn', 'dp', 'dn', 's1p', 'yp', 'g2p', 'cs', 'xp', 'op']
def run(baud, nsym=480, seed=1, tr=1.5e-12):
    rng = np.random.default_rng(seed); sym = rng.integers(0, 4, nsym); ui = 1/baud
    # start balanced (op with equal currents, coupling caps charged symmetrically), data from t=ui
    tp, ip, inn = [0.0], [IFS/2], [IFS/2]; prev_p, prev_n = IFS/2, IFS/2
    for k in range(1, nsym):
        t = k*ui; tp += [t, t+tr]; ip += [prev_p, IFS*sym[k]/3]; inn += [prev_n, IFS*(3-sym[k])/3]; prev_p, prev_n = IFS*sym[k]/3, IFS*(3-sym[k])/3
    pw = lambda v: ' '.join(f'{a:.4e} {b:.4e}' for a, b in zip(tp, v))
    src = f'Idacp cp 0 PWL({pw(ip)})\nIdacn cn 0 PWL({pw(inn)})'
    T = nsym*ui; tag = f'tran{int(baud/1e9)}'
    net = gen.behav_full('', extra=src)
    ctrl = f"tran 0.1p {T:.4e} 0 0.1p\nlinearize\nwrdata {S}/{tag}.dat " + ' '.join(f'v({n})' for n in SAVE)
    fn = f'{S}/{tag}.cir'; open(fn, 'w').write(net + '\n.control\n' + ctrl + '\n.endc\n.end\n')
    t0 = time.time(); r = subprocess.run(['ngspice', '-b', fn], capture_output=True, text=True)
    if 'rror' in r.stdout + r.stderr: print((r.stdout + r.stderr)[-1500:])
    d = np.loadtxt(f'{S}/{tag}.dat'); t = d[:, 0]; V = {n: d[:, 2*i+1] for i, n in enumerate(SAVE)}
    return t, V, sym, ui, time.time()-t0

def eye_metrics(t, vd, sym, ui, skip=40):
    best = None
    for ph in np.linspace(0, 1, 101)[:-1]:
        for lag in range(0, 6):        # chain delay up to ~6 UI
            ts = (np.arange(skip, len(sym)-6)+ph)*ui + lag*ui
            v = np.interp(ts, t, vd); s = sym[skip-lag:len(sym)-6-lag] if False else sym[skip:len(sym)-6]
            lv = [v[s == k] for k in range(4)]
            eh = [lv[k+1].min()-lv[k].max() for k in range(3)]
            sc = min(eh)
            if best is None or sc > best[0]: best = (sc, ph, lag, eh, [x.mean() for x in lv])
    return best

if __name__ == "__main__":
    res = {}
    for baud in (91e9, 182e9):
        t, V, sym, ui, dt = run(baud)
        vd = V['padp']-V['padn']; vin = V['dp']-V['dn']
        sc, ph, lag, eh, mean = eye_metrics(t, vd, sym, ui)
        _, _, _, ehin, mein = eye_metrics(t, vin, sym, ui)
        res[baud] = (t, vd, vin, sym, ui, ph, lag, V)
        print(f'== PAM4 {baud/1e9:.0f} GBd (sim {dt:.0f}s) ==')
        print(f'  driver input : level spread {1e3*(mein[3]-mein[0]):.0f} mVppd (outer levels)   eye heights {[round(1e3*e) for e in ehin]} mV')
        print(f'  PAD diff out : level spread {1e3*(mean[3]-mean[0]):.0f} mVppd (outer levels)   eye heights {[round(1e3*e) for e in eh]} mV  '
              f'(gain {(mean[3]-mean[0])/(mein[3]-mein[0]):.2f}x)  peak-peak incl. ISI {1e3*(vd[t>40*ui].max()-vd[t>40*ui].min()):.0f} mV')
    np.save(f'{S}/tran_res.npy', np.array([0]), allow_pickle=True)
    import pickle; pickle.dump(res, open(f'{S}/tran_res.pkl', 'wb'))
    # transient 4-terminal check (182 GBd run, P side; N side symmetric)
    V = res[182e9][7]; t = res[182e9][0]; m = t > 0.1e-9
    g = lambda n: V[n][m]
    dev = {'M1p (D=Y,G=d,S=B=s1)': (g('yp'), g('dp'), g('s1p'), g('s1p')),
           'M2p (D=x,G=g2,S=B=cs)': (g('xp'), g('g2p'), g('cs'), g('cs')),
           'M3p (D=o,G=vb3,S=B=x)': (g('op'), 0.54+0*g('xp'), g('xp'), g('xp'))}
    print('== transient max |V_ab| (182 GBd, full-scale DAC) ==')
    for nm, (D_, G_, S_, B_) in dev.items():
        T = dict(D=D_, G=G_, S=S_, B=B_)
        mx = {u+v: np.max(np.abs(T[u]-T[v])) for u, v in itertools.combinations('DGSB', 2)}
        sat = np.min(S_-D_)
        print(f'  {nm:24s} worst {max(mx.values()):.2f} V ({max(mx, key=mx.get)})  min Vsd={sat:.2f} V')
    print(f"  tail (VDD-cs) max {np.max(1.8-g('cs')):.2f} V ; SF current source (VDD-s1, 2 stacked) max {np.max(1.8-g('s1p')):.2f} V -> {np.max(1.8-g('s1p'))/2:.2f} V each")
