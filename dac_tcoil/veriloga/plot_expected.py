#!/usr/bin/env python3
"""Expected response of the 35-ohm T-coil design, same axes as the Virtuoso plot.
Reads ngspice wrdata files (freq, dB) produced from tb_tcoil_ngspice.cir."""
import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
q11, qinf, out = sys.argv[1], sys.argv[2], sys.argv[3]
a = np.loadtxt(q11); b = np.loadtxt(qinf)
fig, ax = plt.subplots(figsize=(13, 6.2), facecolor="#1b1d21"); ax.set_facecolor("#14161a")
ax.semilogx(b[:,0], b[:,1], color="#9ca3af", lw=1.4, ls="--", label="lossless coils (Q=∞): DC 29.8 dB, peaks +2.6 dB @80 GHz")
ax.semilogx(a[:,0], a[:,1], color="#ef4444", lw=2.6, label="EXPECTED  (k=+0.29, coil Q≈11 @75 GHz) — ngspice")
def mk(f, tag, dy, dx=0):
    y = np.interp(f, a[:,0], a[:,1]); ax.plot(f, y, "o", ms=9, mfc="none", mec="#fde047", mew=2)
    ax.annotate(f"{tag}: {f/1e9:.2f} GHz\n{y:.2f} dB", (f, y), xytext=(dx, dy), textcoords="offset points", ha="center",
                color="#fde047", fontsize=10.5, bbox=dict(boxstyle="round", fc="#262a31", ec="#6b7280"))
mk(1.590724e9, "A", 28); mk(36.30781e9, "B", -95, -70); mk(60.25596e9, "M1", 30); mk(80.00129e9, "M3", -95, 20)
# user's measured markers
um = [(1.590724e9, 30.90933), (36.30781e9, 27.54675), (60.25596e9, 30.41991), (80.00129e9, 21.77486)]
ax.plot(*zip(*um), "X", ms=11, color="#60a5fa", mec="white", mew=1, label="your simulation (marker readings)", zorder=5)
ax.set_xlim(1e9, 150e9); ax.set_ylim(14, 34)
ax.set_xlabel("freq (Hz)", color="#e5e7eb"); ax.set_ylabel("dB20(V(driver input)/I_dac)  [V (dB)]", color="#e5e7eb")
ax.grid(True, which="both", color="#374151", lw=0.5, alpha=.7); ax.tick_params(colors="#e5e7eb")
for s in ax.spines.values(): s.set_color("#6b7280")
ax.legend(loc="lower left", facecolor="#262a31", edgecolor="#6b7280", labelcolor="#e5e7eb", fontsize=10)
ax.set_title("DAC output T-coil — expected curve, R_L=35 Ω, L1=120 pH, L2=150 pH, k=+0.29, Cb=6 fF, Cdac=200 fF", color="#e5e7eb", fontsize=12)
plt.tight_layout(); plt.savefig(out, dpi=150, facecolor=fig.get_facecolor())
for f in (1.59e9, 10e9, 36.3e9, 60.26e9, 80e9, 100e9):
    print(f"{f/1e9:7.2f} GHz  expected {np.interp(f,a[:,0],a[:,1]):6.2f} dB   lossless {np.interp(f,b[:,0],b[:,1]):6.2f} dB")
i = np.where(a[:,1] < a[0,1]-3)[0]; print("expected -3 dB point: %.1f GHz" % (a[i[0],0]/1e9))
