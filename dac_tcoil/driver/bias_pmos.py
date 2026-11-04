"""4-terminal (D,G,S,B) reliability check of the PMOS-version driver. Limit: every pair <= 0.9 V.
PMOS wells tied to own source; NMOS body = GND (no DNW)."""
import json, itertools, numpy as np, device7nm as D
import design; design.POL = 'p'
from design import to_p
X = json.load(open('final_pmos.json')); p, dac = to_p(X, 'tcoil_padC')
VDD, LIM, RT = 1.8, 0.9, 45.0
I1 = p['n1']*p['j1']*1e-6; I3 = p['I2']; n2 = p['n2']; n3 = n2*p['r3']; j2 = I3*1e6/n2; j3 = I3*1e6/n3
vsg = lambda j: float(D.vgs(j)) + D.VGS_P_EXTRA
vdsat = lambda j: float(np.interp(np.log(j), np.log(D.J), D.VDSAT))
A_SF, SW_IN, SW_O, SW_X = 0.85, 0.10, 0.20, 0.06
CS = VDD - 0.30                      # output-pair common source
VG2 = CS - vsg(j2)                   # = SF output DC
VIN = VG2 - vsg(p['j1'])             # driver input VCM
RD = 25.0; VDSF = RD*I1              # SF drain lifted by Rd (bypassed)
VOUT = I3*RT                         # R_T to GND
XN = 1.20; VB3 = XN - vsg(j3)
VM1 = 1.35                           # mid node of stacked SF current source
VGT2 = VDD - vsg(60); VGU = VDD - vsg(60); VGL = VM1 - vsg(60)
g2 = (VG2 - A_SF*SW_IN, VG2 + A_SF*SW_IN); vin = (VIN - SW_IN, VIN + SW_IN)
vo = (VOUT - SW_O, VOUT + SW_O); xx = (XN - SW_X, XN + SW_X)
# device: name, {terminal: (min,max)}
dev = [
 ('M1p  SF (well=S)',            dict(D=(VDSF, VDSF), G=vin, S=g2, B=g2)),
 ('Ms1u SF current src upper',   dict(D=(VM1, VM1), G=(VGU, VGU), S=(VDD, VDD), B=(VDD, VDD))),
 ('Ms1l SF current src lower',   dict(D=g2, G=(VGL, VGL), S=(VM1, VM1), B=(VM1, VM1))),
 ('Mt2p output tail',            dict(D=(CS, CS), G=(VGT2, VGT2), S=(VDD, VDD), B=(VDD, VDD))),
 ('M2p  output pair (well=CS)',  dict(D=xx, G=g2, S=(CS, CS), B=(CS, CS))),
 ('M3p  output cascode (well=X)',dict(D=vo, G=(VB3, VB3), S=xx, B=xx)),
]
print(f"VDD={VDD}  VCM_in={VIN:.2f}  SF out/M2p gate={VG2:.2f}  CS={CS:.2f}  X={XN:.2f}  VB3={VB3:.2f}  Vout_cm={VOUT:.2f} (swing {vo[0]:.2f}..{vo[1]:.2f})  SF drain={VDSF:.2f}")
print(f"currents/side: SF {I1*1e3:.1f} mA, output {I3*1e3:.1f} mA; total {2*(I1+I3)*1e3:.1f} mA -> {VDD*2*(I1+I3)*1e3:.0f} mW")
print(f"sizes: M1p {p['n1']:.0f} fin @{p['j1']:.0f}uA, M2p {n2:.0f} fin @{j2:.0f}uA, M3p {n3:.0f} fin @{j3:.0f}uA, tail {2*I3*1e6/60:.0f} fin, SF src {I1*1e6/60:.0f} fin x2")
worst_all = 0
for nm, t in dev:
    worst = max(max(abs(a-b) for a in t[u] for b in t[v]) for u, v in itertools.combinations('DGSB', 2))
    vsd = min(t['S'][0]-t['D'][1], t['S'][1]-t['D'][0]) if 'M' in nm else 0
    worst_all = max(worst_all, worst)
    pairs = ', '.join(f"{u}{v}={max(abs(a-b) for a in t[u] for b in t[v]):.2f}" for u, v in itertools.combinations('DGSB', 2))
    print(f"  {nm:30s} worst pair {worst:.2f} V {'OK' if worst <= LIM else 'VIOLATION'} | {pairs}")
print(f"saturation margins: M2p Vsd_min={CS-xx[1]:.2f} vs {vdsat(j2):.2f}; M3p Vsd_min={xx[0]-vo[1]:.2f}... (X tracks current: at Vout max X is high) "
      f"-> M3p Vsd at Vout max={xx[1]-vo[1]:.2f} vs {vdsat(j3):.2f}; SF Vsd_min={g2[0]-VDSF:.2f} vs {vdsat(p['j1']):.2f}")
print(f"NMOS in signal path: none.  n-well/p-sub junctions: up to {VDD} V (well diode, not a MOS terminal pair).")
