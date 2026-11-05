import json, numpy as np, sys
import stage1 as S, device7nm as D
from opt1 import KEYS, P
res = json.load(open('opt1_n.json'))
def pseudo_hd3(n1, j1, gain_ss, vppd):
    """pseudo-differential pair (sources grounded): each side Id(Vcm +- v/2); static I-V; returns HD3, compression"""
    vcm = float(D.vgs(j1)); t = np.arange(512)/512; v = vppd/2*np.sin(2*np.pi*t)
    ia = n1*1e-6*D.id_per_fin(vcm + v/2); ib = n1*1e-6*D.id_per_fin(vcm - v/2)
    vo = (ia-ib)*45; sp = np.abs(np.fft.rfft(vo))/256
    g0 = n1*1e-6*(D.id_per_fin(vcm+1e-4)-D.id_per_fin(vcm-1e-4))/2e-4*45
    return 20*np.log10(sp[3]/sp[1]), 20*np.log10(sp[1]/(g0*vppd/2)), ia.min(), ia.max()
for seed, x in res:
    p = P(x); q = S.at_gain(p, 1.5, 'n')
    print(f"seed{seed}: n1={p['n1']:.0f} j1={p['j1']:.1f}uA/fin (gain1.5: {q['j1']:.1f})  n2={p['I']*1e6/p['j2']:.0f}@{p['j2']:.0f}  "
          f"Lp={p['Lp']*1e12:.0f} L1={p['L1']*1e12:.0f} L2={p['L2']*1e12:.0f} k={p['k']:.2f} C0={p['C0']*1e15:.1f}f Lin={p['Lin']*1e12:.0f} DAC={p['L1d']*1e12:.0f}/{p['L2d']*1e12:.0f}/{p['kd']:.2f}")
    for lab, pp in (('gain2.0', p), ('gain1.5', q)):
        for vin in (0.40, 0.45):
            h3, cp, imin, imax = pseudo_hd3(pp['n1'], pp['j1'], 0, vin)
            print(f"   {lab} in {vin*1e3:.0f}mVppd: HD3={h3:6.1f} dBc  compression={cp:5.2f} dB  I/side {imin*1e3:.1f}..{imax*1e3:.1f} mA")
