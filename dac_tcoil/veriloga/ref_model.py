#!/usr/bin/env python3
"""ref_model.py -- independent numpy reference for the DAC-output T-coil.

Modified nodal analysis with the two coil currents as explicit unknowns, written
exactly like the Verilog-A branch equations in tcoil_cc.va:

    V(a)-V(c) = (r1 + s L1) i1 + s M i2        i1 : a -> c
    V(c)-V(b) =  s M i1 + (r2 + s L2) i2       i2 : c -> b

Usage:  python3 ref_model.py [--Q 9] [--k 0.29] [--L1 120 --L2 150 (pH)] [--Cb 6 (fF)]
                             [--RL 35] [--baseline]
"""
import argparse, numpy as np

def tcoil_z(f, L1=120e-12, L2=150e-12, k=0.29, Cb=6e-15, RL=35., Q=0., f0=75e9,
            Rs=380., Cdac=200e-15, Cp=35e-15, Cc=1e-12, Cin=35e-15, Rb=3e6):
    s = 2j*np.pi*np.asarray(f, float)
    M = k*np.sqrt(L1*L2)
    r1 = 2*np.pi*f0*L1/Q if Q > 0 else 0.0
    r2 = 2*np.pi*f0*L2/Q if Q > 0 else 0.0
    out = np.zeros(len(s), complex)
    for n, sn in enumerate(s):
        # unknowns: Va Vc Vb Vd i1 i2   ; 1 A injected into node c
        A = np.zeros((6, 6), complex); rhs = np.zeros(6, complex)
        A[0] = [1/RL + sn*(Cp+Cb), 0, -sn*Cb, 0, 1, 0]                 # node a
        A[1] = [0, 1/Rs + sn*Cdac, 0, 0, -1, 1];  rhs[1] = 1            # node c
        A[2] = [-sn*Cb, 0, sn*(Cb+Cc), -sn*Cc, 0, -1]                   # node b
        A[3] = [0, 0, -sn*Cc, sn*(Cc+Cin)+1/Rb, 0, 0]                   # node d
        A[4] = [1, -1, 0, 0, -(r1+sn*L1), -sn*M]                        # a-c branch
        A[5] = [0, 1, -1, 0, -sn*M, -(r2+sn*L2)]                        # c-b branch
        out[n] = np.linalg.solve(A, rhs)[3]
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--L1", type=float, default=120); ap.add_argument("--L2", type=float, default=150)
    ap.add_argument("--k", type=float, default=0.29);  ap.add_argument("--Cb", type=float, default=6)
    ap.add_argument("--RL", type=float, default=35);   ap.add_argument("--Q", type=float, default=0)
    ap.add_argument("--baseline", action="store_true", help="no inductors (a, c, b shorted)")
    a = ap.parse_args()
    L1, L2, k = (1e-18, 1e-18, 0.0) if a.baseline else (a.L1*1e-12, a.L2*1e-12, a.k)
    fp = np.array([0.25, 10, 20, 30, 36.25, 45, 60.25, 70, 80, 100])*1e9
    z = tcoil_z(fp, L1, L2, k, a.Cb*1e-15, a.RL, a.Q)
    db = 20*np.log10(np.abs(z))
    print(f"RL={a.RL} L1={L1*1e12:.0f}pH L2={L2*1e12:.0f}pH k={k} Cb={a.Cb}fF Q={a.Q or 'inf'}")
    print(" f(GHz)   dB(V/I)   rel. to 0.25G")
    for f_, d in zip(fp/1e9, db):
        print(f"{f_:7.2f}  {d:8.2f}  {d-db[0]:8.2f}")
    fine = np.linspace(0.25e9, 120e9, 480); dfine = 20*np.log10(np.abs(tcoil_z(fine, L1, L2, k, a.Cb*1e-15, a.RL, a.Q)))
    rel = dfine-dfine[0]; i = np.where(rel < -3)[0]
    print(f"-3 dB point: {fine[i[0]]/1e9:.1f} GHz" if len(i) else "-3 dB point: > 120 GHz")
