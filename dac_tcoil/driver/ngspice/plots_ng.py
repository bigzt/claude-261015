import pickle, numpy as np, sys, matplotlib, subprocess
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '..')
from check_tran import eye_metrics
import gen, chain
from check_linear import run_ng
S = '/tmp/claude-0/-home-user-claude-261015/059e4729-0ea3-5859-88c1-9d51e0756740/scratchpad'
res = pickle.load(open(f'{S}/tran_res.pkl', 'rb'))
for baud, (t, vd, vin, sym, ui, ph, lag, V) in res.items():
    sc, ph2, lag2, eh, mein = eye_metrics(t, -vin, sym, ui)
    print(f'{baud/1e9:.0f} GBd: driver input outer-level spread {1e3*(mein[3]-mein[0]):.0f} mVppd, eye heights {[round(1e3*e) for e in eh]} mV')
# ---- AC comparison figure ----
ac = np.load(f'{S}/ac_cmp.npy'); f = ac[0]
fl, hl = run_ng(gen.linear_half('ladder'), tag='ladder_plot')
fc, hcst = run_ng(gen.linear_half('const'), tag='const_plot')
chain.RDC_FRAC = 1.0; s, _, _ = chain.build(gen.P, 9.0, gen.DAC); hpc = s.solve(fc, 'c', ['pad'])['pad']; chain.RDC_FRAC = 0.15
n = lambda h: 20*np.log10(abs(h)/abs(h[0]))
fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
ax[0].plot(fc/1e9, n(hpc), 'k-', lw=3, alpha=.35, label='python, const-R coil')
ax[0].plot(fc/1e9, n(hcst), 'm--', lw=1.4, label='ngspice linear, const-R (match)')
ax[0].plot(f/1e9, n(ac[2]), 'b-', lw=3, alpha=.35, label='python, sqrt(f) coil (design model)')
ax[0].plot(fl/1e9, n(hl), 'g--', lw=1.6, label='ngspice linear, causal skin ladder')
ax[0].plot(f/1e9, n(ac[1]), 'r-', lw=1.8, label='ngspice behavioural full-diff')
ax[0].axhline(-3, color='k', ls=':', lw=.8); ax[0].axvline(80, color='k', ls=':', lw=.8); ax[0].set_xlim(0, 130); ax[0].set_ylim(-12, 4)
ax[0].grid(alpha=.3); ax[0].legend(fontsize=7.5, loc='lower left'); ax[0].set_xlabel('Frequency (GHz)'); ax[0].set_ylabel('chain |V_pad/I_dac| (dB, norm.)')
ax[0].set_title('AC: Python vs ngspice')
for k, (baud, (t, vd, vin, sym, ui, ph, lag, V)) in enumerate(sorted(res.items())):
    a = ax[k+1]; t0 = 40*ui; m = t > t0
    tt = ((t[m] - (ph+lag)*ui) % (2*ui)) - 0.5*ui
    a.plot(tt*1e12, vd[m]*1e3, ',', color='#dc2626', alpha=.35)
    a.set_xlim(-0.5*ui*1e12, 1.5*ui*1e12); a.set_ylim(-600, 600); a.grid(alpha=.3)
    a.set_xlabel('time (ps), 2 UI'); a.set_ylabel('V(padp)-V(padn) (mV)')
    a.set_title(f'PAM4 {baud/1e9:.0f} GBd, DAC full scale 0–5.2 mA/side (no FFE)')
plt.tight_layout(); plt.savefig('ngspice_verification.png', dpi=140); print('ok')
