import json, numpy as np, chain, device7nm as D
from design import KEYS, met, to_p, F
X = json.load(open('final_q9.json')); OUT = 'tcoil_padC'
def row(nm, **kw):
    m = met(X, OUT, **kw)
    print(f"{nm:34s} gain={10**(m['gdb']/20):4.2f}x  BW3={m['bw']:5.1f}G  ripple(0-70G)={m['rip']:4.2f}dB  peak={m['pk']:+.2f}dB")
print('=== corners (chain Idac->PAD) ===')
row('nominal (Q=9)')
for Q in (7.65, 10.35): row(f'coil Q={Q} (9 -/+15%)', Q=Q)
for cb in (3e-15, 10e-15):
    f = dict(chain.FIX); f['Cb_out'] = cb; row(f'output Cb={cb*1e15:.0f}fF', fix=f)
for ce in (80e-15, 120e-15):
    f = dict(chain.FIX); f['Cesd'] = ce; row(f'ESD {ce*1e15:.0f}fF', fix=f)
for cd in (170e-15, 230e-15):
    f = dict(chain.FIX); f['Cdac'] = cd; row(f'Cdac {cd*1e15:.0f}fF', fix=f)
for rt in (40.5, 49.5):
    f = dict(chain.FIX); f['RT'] = rt; row(f'R_T {rt}ohm', fix=f)
# device parasitics +20 % / gm -10 %
saved = (D.CGS, D.CGD, D.CDB, D.CSB, D.GM_ID.copy())
D.CGS, D.CGD, D.CDB, D.CSB = [c*1.2 for c in saved[:4]]; row('device caps +20%')
D.CGS, D.CGD, D.CDB, D.CSB = saved[:4]; D.GM_ID[:] = saved[4]*0.9; row('device gm -10%')
D.GM_ID[:] = saved[4]
print('=== gain programmability via output tail current (per side I3) ===')
v = dict(zip(KEYS, X))
for I3 in (11e-3, 12.5e-3, 14e-3, 16e-3):
    x2 = list(X); x2[KEYS.index('I3')] = I3; m = met(x2, OUT)
    print(f"  I3={I3*1e3:4.1f} mA: gain={10**(m['gdb']/20):4.2f}x  BW3={m['bw']:5.1f}G  ripple={m['rip']:4.2f}dB  Vout_cm={1.8-I3*45:.2f}V")
