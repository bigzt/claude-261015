"""5-section Foster ladder  Z = R0 + jwL0 + sum Rk||jwLk  fitted to the coil model of chain.py
Z = r75*(0.15+0.85*sqrt(f/75G)) + jwL.  Corners fixed (log spaced) -> Re(Z) linear in Rk (NNLS);
L0 chosen so the effective inductance is exact at 60 GHz."""
import numpy as np
from scipy.optimize import nnls
W0 = 2*np.pi*75e9
FC = np.array([3e9, 10e9, 30e9, 100e9, 300e9]); WC = 2*np.pi*FC
f = np.logspace(np.log10(1e9), np.log10(130e9), 120); w = 2*np.pi*f
def coeffs():
    A = np.column_stack([np.ones_like(w)] + [(w/wc)**2/(1+(w/wc)**2) for wc in WC])
    b = 0.15 + 0.85*np.sqrt(f/75e9)                           # normalised to r75
    x, _ = nnls(A, b); return x                                # [R0, R1..R5]/r75
X = coeffs()
def elements(L, Q):
    r75 = W0*L/Q; R0 = X[0]*r75; Rk = X[1:]*r75; Lk = Rk/WC
    wm = 2*np.pi*60e9; L0 = L - np.sum(Lk/(1+(wm/WC)**2))
    return R0, L0, Rk, Lk
def z(L, Q, ww):
    R0, L0, Rk, Lk = elements(L, Q)
    return R0 + 1j*ww*L0 + sum(1j*ww*Lk[i]*Rk[i]/(Rk[i]+1j*ww*Lk[i]) for i in range(len(Rk)))
if __name__ == "__main__":
    L, Q = 100e-12, 9.0; t = (W0*L/Q)*(0.15+0.85*np.sqrt(f/75e9)) + 1j*w*L; zz = z(L, Q, w)
    sel = f >= 5e9
    print('R/r75:', np.round(X, 4))
    print('Re err %.1f%% (>=5GHz), L_eff err %.1f%% (>=5GHz), |Z| err %.2f%% (all)' % (
        100*np.max(abs(zz.real[sel]/t.real[sel]-1)), 100*np.max(abs(zz.imag[sel]/t.imag[sel]-1)), 100*np.max(abs(abs(zz)/abs(t)-1))))
