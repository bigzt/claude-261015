"""Step A/B: ngspice linear half-circuit vs Python chain model."""
import subprocess, numpy as np, sys
sys.path.insert(0, '..')
import gen, chain
def run_ng(net, outvec='v(pad)', tag='lin'):
    S = '/tmp/claude-0/-home-user-claude-261015/059e4729-0ea3-5859-88c1-9d51e0756740/scratchpad'; fn = f'{S}/ng_{tag}.cir'; dat = f'{S}/ng_{tag}.dat'
    open(fn, 'w').write(net + f"\n.control\nrun\nwrdata {dat} {outvec}\n.endc\n.end\n")
    r = subprocess.run(['ngspice', '-b', fn], capture_output=True, text=True)
    if 'rror' in r.stdout + r.stderr: print((r.stdout + r.stderr)[-1500:])
    d = np.loadtxt(dat); return d[:, 0], d[:, 1] + 1j*d[:, 2]
def py(f, Qm, rdc):
    chain.RDC_FRAC = rdc
    s, inj, dev = chain.build(gen.P, Qm, gen.DAC); return s.solve(f, 'c', ['pad'])['pad']
def metr(f, h):
    db = 20*np.log10(abs(h)); r = db - db[0]; i = np.where(r < -3)[0]
    return db[0], (f[i[0]]/1e9 if len(i) else np.nan), r[f <= 70e9].max() - r[f <= 70e9].min()
for mode, Qm, rdc in (('none', 0, 0.15), ('const', 9.0, 1.0), ('ladder', 9.0, 0.15)):
    f, hn = run_ng(gen.linear_half(mode), tag=mode); hp = py(f, Qm, rdc)
    dmag = np.max(np.abs(20*np.log10(abs(hn)/abs(hp)))); dph = np.max(np.abs(np.angle(hn/hp, deg=True)))
    a, b = metr(f, hn), metr(f, hp)
    lab = {'none': 'Q=inf (lossless)', 'const': 'Q=9 const R', 'ladder': 'Q=9 skin (ngspice causal ladder vs python sqrt(f))'}[mode]
    print(f'{lab:52s} max|dMag|={dmag:.4f} dB  max|dPhase|={dph:.3f} deg | ngspice: DC={a[0]:.2f}dB BW3={a[1]:.1f}G rip={a[2]:.2f} | python: DC={b[0]:.2f}dB BW3={b[1]:.1f}G rip={b[2]:.2f}')
chain.RDC_FRAC = 0.15
