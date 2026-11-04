"""Step C: behavioural full-differential netlist: .op (bias + 4-terminal check), .ac, .tran PAM4."""
import subprocess, numpy as np, sys, json
sys.path.insert(0, '..')
import gen, chain
S = '/tmp/claude-0/-home-user-claude-261015/059e4729-0ea3-5859-88c1-9d51e0756740/scratchpad'
def ng(net, ctrl, tag):
    fn = f'{S}/ngb_{tag}.cir'; open(fn, 'w').write(net + '\n.control\n' + ctrl + '\n.endc\n.end\n')
    r = subprocess.run(['ngspice', '-b', fn], capture_output=True, text=True)
    out = r.stdout + r.stderr
    if 'rror' in out or 'singular' in out.lower(): print(out[-2000:])
    return out
NODES = ['cs', 'vb3'] + [q+s for s in 'pn' for q in ('d', 's1', 'Y', 'g2', 'x', 'o', 'pad', 't', 'c')]
# ---------- .op ----------
out = ng(gen.behav_full(''), 'op\n' + '\n'.join(f'print v({n})' for n in NODES), 'op')
V = {}
for line in out.splitlines():
    if line.strip().startswith('v(') and '=' in line:
        k, val = line.split('='); V[k.strip()[2:-1]] = float(val)
V['vdd'] = 1.8; V['0'] = 0.0; V['vcm'] = 0.37
print('== .op (P side) ==')
for k in ('dp', 's1p', 'yp', 'g2p', 'cs', 'xp', 'op', 'padp', 'tp', 'cp'):
    print(f'  V({k:4s}) = {V[k]:.3f} V')
print(f"  SF current = {V['yp']/25*1e3:.2f} mA/side (design 10.9), output current = {V['tp']/45*1e3:.2f} mA/side (design 15.5)")
json.dump(V, open(f'{S}/op.json', 'w'))
# ---------- .ac (differential) ----------
f0 = gen.behav_full('', extra='Iacp 0 cp AC 1\nIacn 0 cn AC 1 180')
ng(f0, f'ac lin 520 0.1G 130G\nwrdata {S}/ac_b.dat v(padp,padn) v(dp,dn)', 'ac')
d = np.loadtxt(f'{S}/ac_b.dat'); f = d[:, 0]; hc = (d[:, 1]+1j*d[:, 2])/2; hd = (d[:, 4]+1j*d[:, 5])/2
s, inj, dev = chain.build(gen.P, 9.0, gen.DAC); hp = s.solve(f, 'c', ['pad'])['pad']
def metr(h):
    db = 20*np.log10(abs(h)); r = db - db[0]; i = np.where(r < -3)[0]
    return db[0], (f[i[0]]/1e9 if len(i) else np.nan), r[f <= 70e9].max() - r[f <= 70e9].min()
a, b = metr(hc), metr(hp); gdrv = abs(hc[0]/hd[0])
print(f'== .ac differential ==\n  ngspice behavioural: chain DC={a[0]:.2f} dB(V/A)  BW3={a[1]:.1f} GHz  ripple={a[2]:.2f} dB   driver gain={gdrv:.3f}x')
print(f'  python model       : chain DC={b[0]:.2f} dB(V/A)  BW3={b[1]:.1f} GHz  ripple={b[2]:.2f} dB')
np.save(f'{S}/ac_cmp.npy', np.vstack([f, abs(hc), abs(hp), abs(hc/hd)]))
