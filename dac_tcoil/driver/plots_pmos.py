import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
import chain, design; design.POL = 'p'
from design import to_p
X = json.load(open('final_pmos.json')); p, dac = to_p(X, 'tcoil_padC')
f = np.linspace(0.1e9, 130e9, 651)
fig, ax = plt.subplots(1, 2, figsize=(14, 4.8))
for Q, ls, lab in ((9.0, '-', 'Q=9'), (7.65, '--', 'Q=7.65'), (10.35, ':', 'Q=10.35')):
    s, inj, dev = chain.build(p, Q, dac); r = s.solve(f, 'c', ['pad', 'd'])
    ch = 20*np.log10(abs(r['pad'])); dn = 20*np.log10(abs(r['d'])); dr = 20*np.log10(abs(r['pad']/r['d']))
    ax[0].plot(f/1e9, ch-ch[0], 'r'+ls, lw=2.2, label=f'chain Idac→PAD ({lab})')
    if Q == 9.0:
        ax[0].plot(f/1e9, dn-dn[0], color='#2563eb', lw=1.6, label='DAC node: Idac→driver input')
        ax[0].plot(f/1e9, dr-dr[0], color='#059669', lw=1.6, label=f'driver: input→PAD (DC {10**(dr[0]/20):.2f}x)')
    ph = np.unwrap(np.angle(r['pad'])); gd = -np.gradient(ph, 2*np.pi*f)*1e12
    ax[1].plot(f[f <= 100e9]/1e9, gd[f <= 100e9], 'r'+ls, lw=2, label=lab)
ax[0].axhline(-3, color='k', ls=':', lw=.8); ax[0].axvline(80, color='k', ls=':', lw=.8)
ax[0].set_xlim(0, 130); ax[0].set_ylim(-12, 4); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8, loc='lower left')
ax[0].set_xlabel('Frequency (GHz)'); ax[0].set_ylabel('dB (normalised to DC)'); ax[0].set_title('PMOS driver — full chain (45Ω to GND, ESD 100 fF)')
ax[1].axvline(80, color='k', ls=':', lw=.8); ax[1].grid(alpha=.3); ax[1].legend(fontsize=8); ax[1].set_ylim(0, 25)
ax[1].set_xlabel('Frequency (GHz)'); ax[1].set_ylabel('Group delay (ps)'); ax[1].set_title('Chain group delay')
plt.tight_layout(); plt.savefig('driver_pmos_response.png', dpi=150); plt.close()

fig, a = plt.subplots(figsize=(15.5, 8.4)); a.set_xlim(-1, 30); a.set_ylim(-7.4, 7.6); a.axis('off'); a.set_aspect('equal')
c = '#1f2937'; lw = 1.7; red = '#b91c1c'
def w(xs, ys): a.plot(xs, ys, color=c, lw=lw, solid_capstyle='round')
def gnd(x, y): w([x, x], [y, y-0.3]); [w([x-wd, x+wd], [y-0.3-0.15*i]*2) for i, wd in enumerate((0.4, 0.27, 0.13))]
def rail(x1, x2, y, t): w([x1, x2], [y, y]); a.text(x1, y+0.2, t, fontsize=9.5)
def ind_h(x1, x2, y, t):
    n = 3; ww = (x2-x1)/n; tt = np.linspace(0, np.pi, 30)
    for i in range(n): a.plot(x1+i*ww+ww/2*(1-np.cos(tt)), y+0.35*np.sin(tt), color=red, lw=lw+0.5)
    a.text((x1+x2)/2, y+0.6, t, ha='center', fontsize=10, color=red, fontweight='bold')
def ind_v(x, y1, y2, t):
    n = 3; hh = (y2-y1)/n; tt = np.linspace(0, np.pi, 30)
    for i in range(n): a.plot(x+0.35*np.sin(tt), y1+i*hh+hh/2*(1-np.cos(tt)), color=red, lw=lw+0.5)
    a.text(x+0.55, (y1+y2)/2, t, va='center', fontsize=10, color=red, fontweight='bold')
def res_v(x, y1, y2, t):
    ys = np.linspace(y1, y2, 11); pts = [(x, y1)] + [(x+(0.22 if i % 2 else -0.22), ys[i]) for i in range(1, 10)] + [(x, y2)]
    a.plot(*zip(*pts), color=c, lw=lw); a.text(x+0.35, (y1+y2)/2, t, va='center', fontsize=9.5)
def cap_v(x, y, t):
    w([x, x], [y, y-0.55]); w([x-0.4, x+0.4], [y-0.55]*2); w([x-0.4, x+0.4], [y-0.8]*2); w([x, x], [y-0.8, y-1.2]); gnd(x, y-1.2)
    a.text(x+0.5, y-0.7, t, va='center', fontsize=9.5)
def pmos(x, y, lab):   # source top (x,y+0.9), drain bottom (x,y-0.9), gate left (x-1.1,y)
    w([x, x], [y+0.9, y+0.45]); w([x, x], [y-0.9, y-0.45]); w([x-0.25, x-0.25], [y-0.55, y+0.55]); w([x-0.25, x], [y+0.45, y+0.45]); w([x-0.25, x], [y-0.45, y-0.45])
    w([x-0.45, x-0.45], [y-0.45, y+0.45]); a.add_patch(Circle((x-0.56, y), 0.11, fill=False, ec=c, lw=1.3)); w([x-1.1, x-0.67], [y, y])
    a.add_patch(FancyArrowPatch((x, y+0.45), (x-0.22, y+0.45), arrowstyle='-|>', mutation_scale=8, color=c))
    a.text(x+0.15, y, lab, va='center', fontsize=9.3, fontweight='bold')
def box(x, y, wd, ht, t, col='#6b7280'):
    a.add_patch(Rectangle((x, y), wd, ht, fill=False, ec=col, lw=1.4, ls='--')); a.text(x+wd/2, y+ht/2, t, ha='center', va='center', fontsize=8.8, color=col)
rail(6.5, 22.5, 6.6, 'VDD 1.8V')
box(-0.8, -0.9, 4.2, 1.8, 'DAC node + T-coil\nR_L=35Ω, L1=110p/L2=75p\nk=0.43, Cb≈6f, Cdac 200f', '#2563eb')
w([3.4, 4.2], [0, 0]); w([4.2, 4.2], [-0.4, 0.4]); w([4.45, 4.45], [-0.4, 0.4]); a.text(4.3, 0.6, '1pF', ha='center', fontsize=9.5); w([4.45, 5.0], [0, 0])
ind_h(5.0, 6.8, 0, 'Lin 55p'); w([6.8, 8.0], [0, 0]); a.plot(7.4, 0, 'o', color=c, ms=4)
w([7.4, 7.4], [0, -0.4]); res_v(7.4, -0.4, -2.2, '3MΩ'); a.text(7.0, -2.75, 'VCM\n0.37V', ha='right', fontsize=9)
pmos(9.1, 0, 'M1p SF\n390 fin, 28µA/fin\nwell→source')
w([9.1, 9.1], [0.9, 1.6]); a.plot(9.1, 1.6, 'o', color=c, ms=4); w([9.1, 9.1], [1.6, 2.6]); box(8.4, 2.6, 1.4, 2.0, 'I1 10.9mA\n2 stacked\nPMOS'); w([9.1, 9.1], [4.6, 6.6])
w([9.1, 9.1], [-0.9, -1.6]); w([8.5, 9.7], [-1.6, -1.6]); res_v(8.6, -1.6, -3.4, ''); a.text(8.85, -3.0, '25Ω', fontsize=8.5); gnd(8.6, -3.4)
w([9.6, 9.6], [-1.6, -2.1]); w([9.2, 10.0], [-2.1, -2.1]); w([9.2, 10.0], [-2.35, -2.35]); w([9.6, 9.6], [-2.35, -3.4]); gnd(9.6, -3.4); a.text(10.1, -2.2, 'bypass', fontsize=8.5)
w([9.1, 10.4], [1.6, 1.6]); ind_h(10.4, 12.2, 1.6, 'Lg 55p'); w([12.2, 13.1], [1.6, 1.6])
pmos(14.2, 1.6, 'M2p\n388 fin, 40µA/fin\nwell→CS'); w([14.2, 14.2], [2.5, 3.4]); a.plot(14.2, 3.4, 'o', color=c, ms=4); a.text(13.9, 3.65, 'CS 1.50V', fontsize=8.5, ha='right')
w([14.2, 17.0], [3.4, 3.4]); a.text(17.1, 3.4, '→ M2p of N side', va='center', fontsize=8.5)
w([14.2, 14.2], [3.4, 4.0]); box(13.5, 4.0, 1.4, 1.5, 'Itail\n31mA'); w([14.2, 14.2], [5.5, 6.6])
w([14.2, 14.2], [0.7, 0.0]); a.plot(14.2, 0.0, 'o', color=c, ms=4); a.text(14.35, 0.15, 'X 1.20V', fontsize=8.5)
w([14.2, 14.2], [0.0, -0.6]); pmos(14.2, -1.5, 'M3p\n282 fin, 55µA/fin\nwell→X'); w([13.1, 13.1], [-1.5, -1.5]); a.text(13.0, -1.25, 'VB3=0.54V', fontsize=8.5, ha='right')
w([14.2, 14.2], [-2.4, -3.2]); a.plot(14.2, -3.2, 'o', color=c, ms=4); a.text(13.8, -3.05, 'o', fontsize=9)
w([14.2, 15.4], [-3.2, -3.2]); ind_h(15.4, 17.4, -3.2, 'La 120p'); w([17.4, 18.8], [-3.2, -3.2]); a.plot(18.8, -3.2, 'o', color=c, ms=5)
w([18.8, 22.4], [-3.2, -3.2]); a.add_patch(Rectangle((22.4, -3.55), 0.7, 0.7, fill=False, ec=c, lw=1.5)); a.plot([22.4, 23.1], [-3.55, -2.85], color=c, lw=1); a.plot([22.4, 23.1], [-2.85, -3.55], color=c, lw=1)
a.text(22.75, -2.5, 'PAD (OUTP)', ha='center', fontsize=9.5, fontweight='bold')
a.plot(20.9, -3.2, 'o', color=c, ms=4); cap_v(20.9, -3.2, 'ESD 100fF')
w([18.8, 18.8], [-3.2, -3.9]); ind_v(18.8, -5.2, -3.9, 'Lb 155p'); res_v(18.8, -5.2, -6.4, 'R_T 45Ω'); gnd(18.8, -6.4)
a.text(22.5, 5.6, 'Notes', fontsize=10, fontweight='bold')
notes = ['• all-PMOS signal path, each n-well tied to its source', '• every D/G/S/B pair ≤ 0.76 V (limit 0.9 V)', '• no NMOS above 0.9 V (no DNW needed)',
         '• R_T 45Ω to GND; Vout_cm 0.70 V, swing 0.50–0.90 V', '• La, Lb uncoupled; Q=9 @75 GHz assumed', '• N side mirrored; tails per side 15.5 mA',
         '• total 52.8 mA @1.8 V ≈ 95 mW', '• gain 2.0x; trim via tail (needs re-tune)']
for i, t in enumerate(notes): a.text(22.5, 5.0-0.55*i, t, fontsize=9)
a.set_title('Driver, PMOS version (one side): PMOS SF → PMOS cascode CS → series/shunt peaking into 45 Ω to GND + ESD', fontsize=12)
plt.tight_layout(); plt.savefig('driver_pmos_topology.png', dpi=160); plt.close(); print('ok')
