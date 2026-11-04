import json, numpy as np, chain, device7nm as D
import design; design.POL = 'p'
from design import met, KEYS
X = json.load(open('final_pmos.json')); OUT = 'tcoil_padC'
def row(nm, x=X, **kw):
    m = met(x, OUT, **kw)
    print(f"{nm:30s} gain={10**(m['gdb']/20):4.2f}x  BW3={m['bw']:5.1f}G  ripple(0-70G)={m['rip']:4.2f}dB")
row('nominal (Q=9)')
for Q in (7.65, 10.35): row(f'coil Q={Q}', Q=Q)
for key, vals in (('Cesd', (80e-15, 120e-15)), ('Cdac', (170e-15, 230e-15)), ('RT', (40.5, 49.5)), ('Cb_out', (3e-15, 10e-15))):
    for v in vals:
        f = dict(chain.FIX); f[key] = v; row(f'{key}={v*1e15 if v < 1e-9 else v:.4g}', fix=f)
cw = D.cwell; D.cwell = lambda n: 4*cw(n); chain.cwell = D.cwell; row('n-well caps x4 (~20fF)')
D.cwell = cw; chain.cwell = cw
saved = (D.CGS, D.CGD, D.CDB, D.CSB)
D.CGS, D.CGD, D.CDB, D.CSB = [c*1.2 for c in saved]; row('device caps +20%')
D.CGS, D.CGD, D.CDB, D.CSB = saved
D.PGM = 0.85*0.9; row('pFET gm -10%'); D.PGM = 0.85
print('gain trim via tail current:')
for I3 in (11e-3, 12.5e-3, 14e-3, 15.5e-3):
    x2 = list(X); x2[KEYS.index('I3')] = I3; m = met(x2, OUT)
    print(f"  I3={I3*1e3:4.1f} mA/side: gain={10**(m['gdb']/20):.2f}x BW3={m['bw']:.0f}G ripple={m['rip']:.2f}dB Vout_cm={I3*45:.2f}V")
