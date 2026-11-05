"""Single-stage cascode pseudo-differential driver (user topology), half circuit (diff mode).
Idac -> DAC T-coil -> 1 pF -> [Lin] -> M1 gate (source at supply = AC gnd)
M1 -> x -> M2 cascode -> o ;  o -Lp- 45 ohm -> supply(AC gnd)
o -L1- esd(centre, C_ESD) -L2- pad ;  C0 bridging o--pad.  Gain trimmed by VCM (M1 current)."""
import sys, numpy as np
sys.path.insert(0, '../driver')
import chain, device7nm as D
from chain import Net
FIX = dict(Rdac=380., Cdac=200e-15, Cp=35e-15, Cc=1e-12, Rb=3e6, Croute_in=5e-15, RT=45., CRT=5e-15, Cesd=80e-15, Cpad=10e-15)
def build(p, Q=9.0, fix=FIX, with_dac=True, kind='p'):
    s = Net(); f = fix
    if with_dac:
        s.R('c', '0', f['Rdac']); s.C('c', '0', f['Cdac']); s.R('a', '0', 35.0); s.C('a', '0', f['Cp'])
        s.K('a', 'c', 'b', p['L1d'], p['L2d'], p['kd'], Q); s.C('a', 'b', 6e-15)
        if p['Lin'] > 1e-13: s.C('b', 'd0', f['Cc']); s.L('d0', 'd', p['Lin'], Q); s.C('d0', '0', 3e-15)
        else: s.C('b', 'd', f['Cc'])
    s.R('d', '0', f['Rb']); s.C('d', '0', f['Croute_in'])
    m1 = D.mos(p['n1'], p['j1'], kind); n2 = p['I']*1e6/p['j2']; m2 = D.mos(n2, p['j2'], kind)
    s.G('x', '0', 'd', m1['gm']); s.R('x', '0', 1/m1['gds'])
    s.C('d', '0', m1['cgs']); s.C('d', 'x', m1['cgd']); s.C('x', '0', m1['cdb'])
    s.G('o', 'x', '0', m2['gm']); s.R('o', 'x', 1/m2['gds'])
    s.C('x', '0', m2['cgs'] + m2['csb'] + D.cwell(n2))   # n-well (pmos) or DNW p-well (nmos) tied to source X
    s.C('o', '0', m2['cgd'] + m2['cdb'] + 3e-15)
    s.L('o', 't', max(p['Lp'], 1e-15), Q if p['Lp'] > 1e-13 else 0); s.R('t', '0', f['RT']); s.C('t', '0', f['CRT'])
    s.K('o', 'esd', 'pad', p['L1'], p['L2'], p['k'], Q); s.C('esd', '0', f['Cesd']); s.C('o', 'pad', p['C0'] + 1e-18); s.C('pad', '0', f['Cpad'])
    return s, dict(m1=m1, m2=m2, n2=n2)
F = np.concatenate([[0.1e9], np.arange(2e9, 122e9, 2e9)])
def metrics(p, Q=9.0, fix=FIX, kind='p'):
    s, dv = build(p, Q, fix, False, kind); r = s.solve(F[:1], 'd', ['d', 'pad']); g = abs(r['pad'][0]/r['d'][0])
    s, dv = build(p, Q, fix, True, kind); r = s.solve(F, 'c', ['pad', 'd'])
    db = 20*np.log10(abs(r['pad'])); rel = db - db[0]; i = np.where(rel < -3)[0]
    bw = F[i[0]]/1e9 if len(i) else F[-1]/1e9; rip = rel[F <= 70e9].max() - rel[F <= 70e9].min()
    sd, _ = build(p, Q, fix, False, kind); rr = sd.solve([1e9, 40e9], 'd', ['d']); cin = np.imag(1/rr['d'])/(2*np.pi*np.array([1e9, 40e9]))
    return dict(g=g, bw=bw, rip=rip, cin=cin*1e15, pk=rel.max())
def at_gain(p, gt, kind='p'):
    """lower M1 (and M2) current at fixed fins (=VCM tuning) until DC gain = gt"""
    q = dict(p); lo, hi = 3.0, p['j1']
    for _ in range(30):
        q['j1'] = 0.5*(lo+hi); q['I'] = q['n1']*q['j1']*1e-6; q['j2'] = q['I']*1e6/(p['I']*1e6/p['j2'])
        g = metrics_gain(q, kind)
        lo, hi = (q['j1'], hi) if g < gt else (lo, q['j1'])
    return q
def metrics_gain(p, kind='p'):
    s, dv = build(p, 9.0, FIX, False, kind); r = s.solve(F[:1], 'd', ['d', 'pad']); return abs(r['pad'][0]/r['d'][0])
