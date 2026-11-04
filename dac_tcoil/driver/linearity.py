"""Large-signal linearity of the output differential pair (shared tail, Rs=0) with the
estimated FinFET I-V (device7nm.id_per_fin).  SF treated as linear gain A_sf
(taken from the small-signal model); M3 cascode passes the current unchanged."""
import json, numpy as np, chain, device7nm as D
from scipy.optimize import brentq
from design import to_p
X = json.load(open('final_q9.json')); p, dac = to_p(X, 'tcoil_padC')
n2 = p['n2']; I3 = p['I2']; ITAIL = 2*I3; RT = 45.0
# SF small-signal low-frequency gain
s, inj, dev = chain.build(p, 9.0, dac, chain.FIX, False, 'vin')
r = s.solve([0.1e9], 'd', ['d', 'g2']); A_sf = abs(r['g2'][0]/r['d'][0])
def pair(vd):                          # vd = differential gate voltage of M2 pair
    f = lambda vs: n2*1e-6*(D.id_per_fin(1.0+vd/2-vs) + D.id_per_fin(1.0-vd/2-vs)) - ITAIL
    vs = brentq(f, -1, 1.0)
    ia = n2*1e-6*D.id_per_fin(1.0+vd/2-vs); return ia, ITAIL-ia
def out_diff(vin_diff):
    ia, ib = pair(A_sf*vin_diff); return -(ia-ib)*RT, ia, ib
g0 = (out_diff(1e-4)[0]-out_diff(-1e-4)[0])/2e-4
print(f"A_sf={A_sf:.3f}  small-signal driver gain (DC, static model) = {abs(g0):.2f}x")
for vppd in (0.35, 0.40, 0.45):
    t = np.arange(256)/256; vin = vppd/2*np.sin(2*np.pi*t)
    vo = np.array([out_diff(v)[0] for v in vin]); ia = np.array([out_diff(v)[1] for v in vin])
    spec = np.abs(np.fft.rfft(vo))/128
    hd3 = 20*np.log10(spec[3]/spec[1]); comp = 20*np.log10(spec[1]/(abs(g0)*vppd/2))
    print(f"in {vppd*1e3:.0f} mVppd -> out {2*spec[1]*1e3:.0f} mVppd (fund.)  HD3={hd3:6.1f} dBc  gain compression={comp:5.2f} dB  "
          f"M2/M3 current per side {ia.min()*1e3:.1f}..{ia.max()*1e3:.1f} mA")
