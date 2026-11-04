"""DC operating point & core-device reliability check (limit assumed 0.9 V per terminal pair)."""
import json, numpy as np, device7nm as D
from design import to_p
X = json.load(open('final_q9.json')); p, dac = to_p(X, 'tcoil_padC')
VDD, LIM = 1.8, 0.9
I1 = p['n1']*p['j1']*1e-6; I3 = p['I2']; n2 = p['n2']; n3 = p['n2']*p['r3']
j2 = I3*1e6/n2; j3 = I3*1e6/n3
vgs1 = float(D.vgs(p['j1'])); vgs2 = float(D.vgs(j2)); vgs3 = float(D.vgs(j3)); vgsc = vgs1
vd = lambda j: float(np.interp(np.log(j), np.log(D.J), D.VDSAT))
SW_IN, SW_G2, SW_O = 0.10, 0.10*0.78, 0.20         # per-side peak swings at 400/800 mVppd
VS2 = 0.30                                          # output-pair common source (tail headroom)
VG2 = VS2 + vgs2; VIN = VG2 + vgs1                  # SF output = M2 gate
VY = VDD - vgsc                                     # SF drain (cascode Mc gate = VDD)
VOUT = VDD - I3*45.0
VX = 0.62; VB3 = VX + vgs3
VT1M = 0.50                                         # SF tail: two stacked devices, mid node
rows = [
 ('Mc  (SF drain cascode, gate=VDD)', VDD-VY, vgsc, VDD-VDD, 'Vds=VDD-VY'),
 ('M1  (SF)', VY-VG2+SW_G2, vgs1, VIN+SW_IN-VY, ''),
 ('Mt1a (SF tail, upper)', VG2+SW_G2-VT1M, None, None, 'cascoded tail'),
 ('Mt1b (SF tail, lower)', VT1M, None, None, ''),
 ('M2  (output pair)', VX-VS2, vgs2, VG2+SW_G2-VX, 'min Vds=%.2f vs Vdsat %.2f' % (VX-VS2-0.03, vd(j2))),
 ('M3  (output cascode)', VOUT+SW_O-VX, vgs3, VB3-(VOUT-SW_O), 'min Vds=%.2f vs Vdsat %.2f' % (VOUT-SW_O-VX, vd(j3))),
 ('Mt2 (output tail)', VS2, None, None, 'Vdsat(60uA/fin)=%.2f' % vd(60)),
]
print(f"VDD={VDD} V | VCM_in={VIN:.2f} V | SF out / M2 gate={VG2:.2f} V | VY={VY:.2f} V | VS2={VS2:.2f} V | VX={VX:.2f} V | VB3={VB3:.2f} V | Vout_cm={VOUT:.2f} V (swing {VOUT-SW_O:.2f}..{VOUT+SW_O:.2f})")
print(f"currents: SF {I1*1e3:.2f} mA/side, output {I3*1e3:.1f} mA/side (tail {2*I3*1e3:.0f} mA); total {2*(I1+I3)*1e3:.1f} mA -> {VDD*2*(I1+I3)*1e3:.0f} mW")
print(f"sizes: M1/Mc {p['n1']:.0f} fins @{p['j1']:.0f}uA, M2 {n2:.0f} fins @{j2:.0f}uA, M3 {n3:.0f} fins @{j3:.0f}uA, Mt2 {2*I3*1e6/60:.0f} fins, Mt1 {I1*1e6/60:.0f} fins (x2 stacked)")
print(f"{'device':34s} {'max Vds':>8s} {'Vgs':>6s} {'max Vgd':>8s}  ok(<= {LIM} V)")
for nm, vds, vgs, vgd, note in rows:
    ok = vds <= LIM and (vgs is None or vgs <= LIM) and (vgd is None or abs(vgd) <= LIM)
    print(f"{nm:34s} {vds:8.2f} {'' if vgs is None else f'{vgs:6.2f}':>6s} {'' if vgd is None else f'{vgd:8.2f}':>8s}  {'OK' if ok else 'VIOLATION'}  {note}")
