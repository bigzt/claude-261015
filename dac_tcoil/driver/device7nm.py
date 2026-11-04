"""Estimated 7nm FinFET core-NMOS per-fin model (layout-inclusive, NOT PDK data).

Numbers are order-of-magnitude estimates from public literature for an N7-class
core nFET at minimum L with local M0-M2 parasitics.  Use only for first-pass
sizing; re-size with the PDK gm/Id curves (the gm/Id + current targets transfer
directly).
"""
import numpy as np

# current density (uA/fin) -> gm/Id (1/V), intrinsic gain gm/gds, Vdsat (V)
J      = np.array([5, 10, 15, 25, 40, 60, 80], float)
GM_ID  = np.array([12.5, 10.5, 9.0, 7.0, 5.3, 4.2, 3.6])
A_INT  = np.array([10.0, 9.0, 8.0, 7.0, 6.0, 5.5, 5.0])
VDSAT  = np.array([0.06, 0.08, 0.10, 0.13, 0.17, 0.21, 0.25])
VTH, VGS_5UA = 0.30, 0.33
# capacitances per fin (fF), layout-inclusive
CGS, CGD, CDB, CSB = 0.11, 0.05, 0.045, 0.045

def _interp(j, tab):
    return np.interp(np.log(j), np.log(J), tab)

# Vgs(j) consistent with gm/Id:  dVgs = d(ln Id) / (gm/Id)
_lj = np.linspace(np.log(2.0), np.log(120.0), 400)
_vgs = VGS_5UA + np.concatenate([[0], np.cumsum(np.diff(_lj) / np.interp(0.5*(_lj[1:]+_lj[:-1]), np.log(J), GM_ID))])
_vgs -= np.interp(np.log(5.0), _lj, _vgs) - VGS_5UA

def vgs(j_uA):   return np.interp(np.log(j_uA), _lj, _vgs)
def id_per_fin(vgs_v):   # inverse of vgs(), uA/fin (large-signal I-V, saturation)
    return np.exp(np.interp(vgs_v, _vgs, _lj))

PGM = 0.85          # pFET gm per fin relative to nFET (assumption)
VGS_P_EXTRA = 0.02  # pFET |Vgs| slightly higher at same density (assumption)
def cwell(nfin):    # n-well-to-psub cap of a PMOS whose well is tied to a signal node (estimate)
    return (3.0 + 0.01*nfin)*1e-15

def mos(nfin, j_uA, kind='n'):
    """small-signal params for nfin fins biased at j_uA uA/fin"""
    I = nfin * j_uA * 1e-6
    gm = I * _interp(j_uA, GM_ID) * (PGM if kind == 'p' else 1.0)
    return dict(I=I, gm=gm, gds=gm/_interp(j_uA, A_INT), vgs=float(vgs(j_uA)) + (VGS_P_EXTRA if kind == 'p' else 0.0),
                vdsat=float(_interp(j_uA, VDSAT)),
                cgs=nfin*CGS*1e-15, cgd=nfin*CGD*1e-15, cdb=nfin*CDB*1e-15, csb=nfin*CSB*1e-15,
                ft=gm/(2*np.pi*nfin*(CGS+CGD)*1e-15))

if __name__ == "__main__":
    for j in (5, 15, 25, 40, 60, 80):
        m = mos(1, j)
        print(f"{j:3d} uA/fin: gm={m['gm']*1e6:6.1f} uS  gm/Id={m['gm']/m['I']:5.2f}  Vgs={m['vgs']:.3f} V  fT={m['ft']/1e9:5.0f} GHz")
