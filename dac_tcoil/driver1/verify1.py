import json, numpy as np, itertools
import stage1 as S, device7nm as D
from opt1 import KEYS, P
x = dict(zip(KEYS, dict(json.load(open('opt1_n.json')))[37]))
R = lambda v, s: round(v/s)*s
p = dict(n1=330, j1=20.2, j2=37, Lp=R(x['Lp'], 5e-12), L1=R(x['L1'], 5e-12), L2=R(x['L2'], 5e-12), k=round(x['k'], 2), C0=1e-15,
         L1d=R(x['L1d'], 5e-12), L2d=R(x['L2d'], 5e-12), kd=round(x['kd'], 2), Lin=R(x['Lin'], 5e-12))
p['I'] = p['n1']*p['j1']*1e-6
# trim j1 so that gain = 2.00
p = S.at_gain(dict(p, j1=40, I=p['n1']*40e-6, j2=37*40/20.2), 2.0, 'n'); q = S.at_gain(p, 1.5, 'n')
json.dump({k: float(v) for k, v in p.items()}, open('final1.json', 'w'))
n2 = p['I']*1e6/p['j2']
print(f"M1 {p['n1']:.0f} fin; M2 {n2:.0f} fin | gain2.0: j1={p['j1']:.1f} uA/fin I={p['I']*1e3:.2f} mA/side | gain1.5: j1={q['j1']:.1f} I={q['I']*1e3:.2f} mA/side")
print(f"Lp={p['Lp']*1e12:.0f}p  T-coil L1={p['L1']*1e12:.0f}p L2={p['L2']*1e12:.0f}p k={p['k']}  C0={p['C0']*1e15:.0f}f  Lin={p['Lin']*1e12:.0f}p  DAC T-coil {p['L1d']*1e12:.0f}/{p['L2d']*1e12:.0f}p k={p['kd']}")
def row(nm, pp, **kw):
    m = S.metrics(pp, kind='n', **kw); return f"{nm:24s} gain={m['g']:.2f}x BW3={m['bw']:5.1f}G ripple={m['rip']:.2f}dB"
print('=== corners ===')
for lab, pp in (('gain2.0', p), ('gain1.5', q)):
    print(row(f'{lab} nominal', pp), f"Cin={S.metrics(pp, kind='n')['cin'][0]:.0f}fF")
    for Q in (7.65, 10.35): print(row(f'{lab} Q={Q}', pp, Q=Q))
    for key, vals in (('Cesd', (64e-15, 96e-15)), ('Cdac', (170e-15, 230e-15)), ('RT', (40.5, 49.5)), ('Cpad', (0, 25e-15))):
        for v in vals:
            f = dict(S.FIX); f[key] = v; print(row(f'{lab} {key}={v*1e15 if v < 1e-9 else v:.4g}', pp, fix=f))
    sv = (D.CGS, D.CGD, D.CDB, D.CSB); D.CGS, D.CGD, D.CDB, D.CSB = [c*1.2 for c in sv]; print(row(f'{lab} device caps +20%', pp)); D.CGS, D.CGD, D.CDB, D.CSB = sv
# ---- bias / 4-terminal check (VDD 1.5, VB 1.2, body: M1=GND(source), M2=DNW p-well tied to X) ----
VDD, VB = 1.5, 1.2
print('=== bias & 4-terminal check (pointwise over the signal swing, static) ===')
for lab, pp, vin in (('gain2.0', p, 0.45), ('gain1.5', q, 0.45)):
    vcm = float(D.vgs(pp['j1'])); worst = 0; okd = True; msat = 9
    for v in np.linspace(-vin/2, vin/2, 41):
        for sg in (1, -1):
            vg = vcm + sg*v/2; i = pp['n1']*1e-6*D.id_per_fin(vg); j2 = i*1e6/n2
            X = VB - float(D.vgs(j2)); vo = VDD - i*45
            T1 = dict(D=X, G=vg, S=0, B=0); T2 = dict(D=vo, G=VB, S=X, B=X)
            for T in (T1, T2): worst = max(worst, max(abs(T[a]-T[b]) for a, b in itertools.combinations('DGSB', 2)))
            msat = min(msat, X - float(np.interp(np.log(i*1e6/pp['n1']), np.log(D.J), D.VDSAT)), vo - X - float(np.interp(np.log(j2), np.log(D.J), D.VDSAT)))
    i0 = pp['I']; X0 = VB - float(D.vgs(pp['j2']))
    print(f"{lab}: VCM={vcm:.3f} V  X={X0:.2f} V  Vout_cm={VDD-i0*45:.2f} V | worst |V_ab| over swing = {worst:.2f} V {'OK' if worst <= 0.9 else 'VIOLATION'} | min saturation margin = {msat*1e3:.0f} mV")
# ---- VCM / leakage sensitivity ----
print('=== sensitivity to gate-bias error (leakage x R_bias) ===')
for lab, pp in (('gain2.0', p), ('gain1.5', q)):
    vcm = float(D.vgs(pp['j1'])); g0 = S.metrics_gain(pp, 'n')
    for dv in (0.005, 0.01, 0.03):
        j = float(D.id_per_fin(vcm+dv)); qq = dict(pp, j1=j, I=pp['n1']*j*1e-6, j2=pp['j2']*j/pp['j1'])
        g1 = S.metrics_gain(qq, 'n')
        print(f"  {lab}: dVg=+{dv*1e3:.0f} mV -> current {100*(j/pp['j1']-1):+.1f}%  gain {100*(g1/g0-1):+.1f}%  Vout_cm {-(qq['I']-pp['I'])*45*1e3:+.0f} mV")
for R_ in (3e6, 300e3, 100e3):
    print(f"  R_bias={R_/1e3:.0f}k: offset per nA leakage = {R_*1e-9*1e3:.1f} mV ; HPF corner with 1pF+Cin = {1/(2*np.pi*R_*1.08e-12)/1e3:.0f} kHz")
