import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch
import chain
from design import to_p
X = json.load(open('final_q9.json')); p, dac = to_p(X, 'tcoil_padC')
f = np.linspace(0.1e9, 130e9, 651)
# ---------------- responses ----------------
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
ax[0].axhline(-3, color='k', ls=':', lw=.8); ax[0].axhline(1.5, color='gray', ls=':', lw=.8); ax[0].axvline(80, color='k', ls=':', lw=.8)
ax[0].set_xlim(0, 130); ax[0].set_ylim(-12, 4); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8, loc='lower left')
ax[0].set_xlabel('Frequency (GHz)'); ax[0].set_ylabel('dB (normalised to DC)'); ax[0].set_title('Full chain (R_L=35Ω DAC node + driver, 45Ω load, ESD 100 fF)')
ax[1].axvline(80, color='k', ls=':', lw=.8); ax[1].grid(alpha=.3); ax[1].legend(fontsize=8); ax[1].set_ylim(0, 25)
ax[1].set_xlabel('Frequency (GHz)'); ax[1].set_ylabel('Group delay (ps)'); ax[1].set_title('Chain group delay')
plt.tight_layout(); plt.savefig('driver_response.png', dpi=150); plt.close()

# ---------------- schematic (one side; other side mirrored) ----------------
fig, a = plt.subplots(figsize=(15, 7.6)); a.set_xlim(-1, 29); a.set_ylim(-6.6, 7.4); a.axis('off'); a.set_aspect('equal')
c = '#1f2937'; lw = 1.7; red = '#b91c1c'
def w(xs, ys, col=c): a.plot(xs, ys, color=col, lw=lw, solid_capstyle='round')
def gnd(x, y):
    w([x, x], [y, y-0.3]); [w([x-wd, x+wd], [y-0.3-0.15*i]*2) for i, wd in enumerate((0.4, 0.27, 0.13))]
def vdd(x, y, t='VDD 1.8V'): w([x-0.5, x+0.5], [y, y]); a.text(x, y+0.2, t, ha='center', fontsize=9)
def ind_h(x1, x2, y, t, col=red):
    n = 3; ww = (x2-x1)/n; tt = np.linspace(0, np.pi, 30)
    for i in range(n): a.plot(x1+i*ww+ww/2*(1-np.cos(tt)), y+0.35*np.sin(tt), color=col, lw=lw+0.5)
    a.text((x1+x2)/2, y+0.6, t, ha='center', fontsize=10, color=col, fontweight='bold')
def ind_v(x, y1, y2, t, col=red):
    n = 3; hh = (y2-y1)/n; tt = np.linspace(0, np.pi, 30)
    for i in range(n): a.plot(x+0.35*np.sin(tt), y1+i*hh+hh/2*(1-np.cos(tt)), color=col, lw=lw+0.5)
    a.text(x+0.55, (y1+y2)/2, t, va='center', fontsize=10, color=col, fontweight='bold')
def res_v(x, y1, y2, t):
    ys = np.linspace(y1, y2, 11); pts = [(x, y1)] + [(x+(0.22 if i % 2 else -0.22), ys[i]) for i in range(1, 10)] + [(x, y2)]
    a.plot(*zip(*pts), color=c, lw=lw); a.text(x+0.35, (y1+y2)/2, t, va='center', fontsize=9.5)
def cap_v(x, y, t):
    w([x, x], [y, y-0.55]); w([x-0.4, x+0.4], [y-0.55]*2); w([x-0.4, x+0.4], [y-0.8]*2); w([x, x], [y-0.8, y-1.2]); gnd(x, y-1.2)
    a.text(x+0.5, y-0.7, t, va='center', fontsize=9.5)
def nmos(x, y, lab, gate='left'):      # drain top (x,y+0.9), source bottom (x,y-0.9), gate at (x-1, y)
    w([x, x], [y+0.9, y+0.45]); w([x, x], [y-0.9, y-0.45]); w([x-0.25, x-0.25], [y-0.55, y+0.55]); w([x-0.25, x], [y+0.45, y+0.45]); w([x-0.25, x], [y-0.45, y-0.45])
    w([x-0.45, x-0.45], [y-0.45, y+0.45]); w([x-1.0, x-0.45], [y, y])
    a.add_patch(FancyArrowPatch((x-0.2, y-0.45), (x, y-0.45), arrowstyle='-|>', mutation_scale=8, color=c))
    a.text(x+0.15, y, lab, va='center', fontsize=9.5, fontweight='bold')
def box(x, y, wd, ht, t, col='#2563eb'):
    a.add_patch(Rectangle((x, y), wd, ht, fill=False, ec=col, lw=1.4, ls='--')); a.text(x+wd/2, y+ht/2, t, ha='center', va='center', fontsize=9, color=col)
# DAC block
box(-0.8, -0.9, 4.2, 1.8, 'DAC node + T-coil\nR_L=35Ω, L1=110p/L2=105p\nk=0.45, Cb≈6f, Cdac 200f')
w([3.4, 4.2], [0, 0]); w([4.2, 4.2], [-0.4, 0.4]); w([4.45, 4.45], [-0.4, 0.4]); a.text(4.3, 0.6, '1pF', ha='center', fontsize=9.5); w([4.45, 5.0], [0, 0])
ind_h(5.0, 6.8, 0, 'Lin 65p'); w([6.8, 7.8], [0, 0]); a.plot(7.4, 0, 'o', color=c, ms=4)
w([7.4, 7.4], [0, -0.4]); res_v(7.4, -0.4, -2.2, '3MΩ'); a.text(7.4, -2.6, 'VCM=1.47V', ha='center', fontsize=9)
# SF M1 + Mc
nmos(8.8, 0, 'M1 SF\n230 fin\n28µA/fin'); w([8.8, 8.8], [0.9, 1.6]); nmos(8.8, 2.5, 'Mc 230 fin'); w([7.8, 7.8], [2.5, 2.5]); w([7.8, 7.8], [2.5, 3.6]); w([7.8, 8.8], [3.6, 3.6])
w([8.8, 8.8], [3.4, 4.2]); vdd(8.8, 4.2); a.text(7.6, 3.1, 'gate=VDD', fontsize=8, ha='right')
w([8.8, 8.8], [-0.9, -1.5]); a.plot(8.8, -1.5, 'o', color=c, ms=4)
w([8.8, 8.8], [-1.5, -2.2]); box(8.2, -3.6, 1.2, 1.4, 'I1\n6.4mA\ncasc.', '#6b7280'); gnd(8.8, -3.6)
# Lg
w([8.8, 10.0], [-1.5, -1.5]); ind_h(10.0, 11.8, -1.5, 'Lg 70p'); w([11.8, 12.6], [-1.5, -1.5]); w([12.6, 12.6], [-1.5, -1.5])
# M2, M3
nmos(13.6, -1.5, 'M2\n267 fin\n60µA/fin'); w([13.6, 13.6], [-0.6, 0.0]); a.plot(13.6, 0.0, 'o', color=c, ms=4); a.text(13.75, -0.15, 'X 0.62V', fontsize=8)
w([13.6, 13.6], [0.0, 0.6]); nmos(13.6, 1.5, 'M3\n320 fin\n50µA/fin'); w([12.6, 12.6], [1.5, 1.5]); a.text(12.5, 1.75, 'VB3=1.24V', fontsize=8.5, ha='right')
w([13.6, 13.6], [-2.4, -3.2]); w([13.6, 16.2], [-3.2, -3.2]); a.text(16.3, -3.2, '→ M2 of N side', va='center', fontsize=8.5)
w([13.6, 13.6], [-3.2, -3.8]); box(13.0, -5.2, 1.2, 1.4, 'Itail\n32mA', '#6b7280'); gnd(13.6, -5.2)
# output network
w([13.6, 13.6], [2.4, 3.2]); a.plot(13.6, 3.2, 'o', color=c, ms=4); a.text(13.2, 3.45, 'o', fontsize=9)
w([13.6, 15.0], [3.2, 3.2]); ind_h(15.0, 17.0, 3.2, 'La 105p'); w([17.0, 18.4], [3.2, 3.2]); a.plot(18.4, 3.2, 'o', color=c, ms=5)
w([18.4, 18.4], [3.2, 2.6]); cap_v(18.4, 2.6, 'ESD 100fF')
w([18.4, 20.0], [3.2, 3.2]); a.add_patch(Rectangle((20.0, 2.85), 0.7, 0.7, fill=False, ec=c, lw=1.5)); a.plot([20.0, 20.7], [2.85, 3.55], color=c, lw=1); a.plot([20.0, 20.7], [3.55, 2.85], color=c, lw=1)
a.text(20.35, 2.4, 'PAD (OUTP)', ha='center', fontsize=9.5, fontweight='bold')
w([18.4, 18.4], [3.2, 3.9]); ind_v(18.4, 3.9, 5.2, 'Lb 160p'); w([18.4, 18.4], [5.2, 5.4]); res_v(18.4, 5.4, 6.6, 'R_T 45Ω'); vdd(18.4, 6.8)
a.text(22.0, 5.8, 'Notes', fontsize=10, fontweight='bold')
notes = ['• differential; N side identical (mirror)', '• La, Lb: separate (uncoupled) coils', '• Q=9 @75 GHz assumed (retune if EM differs)',
         '• no external load: R_T is the only load', '• Vout_cm = 1.08 V, swing 0.88–1.28 V', '• gain trim: I_tail 25–32 mA → 1.6–2.0×',
         '• all core-device |Vxy| ≤ 0.66 V', '• total 44.9 mA @1.8 V ≈ 81 mW']
for i, t in enumerate(notes): a.text(22.0, 5.2-0.55*i, t, fontsize=9)
a.set_title('Driver (one side shown): SF buffer → cascode CS output stage → series/shunt peaking into 45 Ω + ESD', fontsize=12)
plt.tight_layout(); plt.savefig('driver_topology.png', dpi=160); plt.close()
print('ok')
