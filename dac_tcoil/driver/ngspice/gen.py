"""ngspice netlist generator for the PMOS-version driver chain (final_pmos.json).
  linear_half(): half-circuit, element-for-element equal to chain.build() (VCCS devices)
  behav_full():  full differential, PMOS as B-source table I-V (sat + triode), real bias
Coil loss modes: 'none' (Q=inf), 'const' (R=w0L/Q), 'ladder' (causal skin ladder)."""
import json, sys, numpy as np
sys.path.insert(0, '..')
import chain, device7nm as D, design
design.POL = 'p'
from design import to_p
import ladder_fit as LF
HERE = __file__.rsplit('/', 1)[0] if '/' in __file__ else '.'
X = json.load(open(f'{HERE}/../final_pmos.json')); P, DAC = to_p(X, 'tcoil_padC')
FX = chain.FIX; W0 = chain.W0; Q = 9.0

class NL:
    def __init__(s): s.l = []; s.k = 0
    def a(s, t): s.l.append(t)
    def nid(s, p): s.k += 1; return f'{p}{s.k}'
    def coil(s, n1, n2, L, mode, tag):
        """returns the name of the main inductor (for K coupling)"""
        if mode == 'none':
            s.a(f'L{tag} {n1} {n2} {L:.6e}'); return f'L{tag}'
        if mode == 'const':
            m = s.nid('m'); s.a(f'R{tag} {n1} {m} {W0*L/Q:.6e}'); s.a(f'L{tag} {m} {n2} {L:.6e}'); return f'L{tag}'
        R0, L0, Rk, Lk = LF.elements(L, Q)
        m = s.nid('m'); s.a(f'R{tag}_0 {n1} {m} {R0:.6e}'); prev = s.nid('m'); s.a(f'L{tag} {m} {prev} {L0:.6e}')
        for i, (r, l) in enumerate(zip(Rk, Lk)):
            nxt = n2 if i == len(Rk)-1 else s.nid('m')
            s.a(f'R{tag}_{i+1} {prev} {nxt} {r:.6e}'); s.a(f'L{tag}_{i+1} {prev} {nxt} {l:.6e}'); prev = nxt
        return f'L{tag}'
    def kpair(s, la, lb, L1, L2, k, mode, tag):
        if mode == 'ladder':   # keep mutual M = k sqrt(L1 L2) on the main sections
            _, L01, _, _ = LF.elements(L1, Q); _, L02, _, _ = LF.elements(L2, Q); k = k*np.sqrt(L1*L2)/np.sqrt(L01*L02)
        if abs(k) > 1e-6: s.a(f'K{tag} {la} {lb} {min(k, 0.999):.6f}')

def dac_tcoil(n, side, mode, sfx):
    """DAC node + T-coil + coupling + Lin. nodes: c,a,b,d0,d (suffix)"""
    c, a, b, d0, d = (x+sfx for x in 'c a b d0 d'.split())
    n.a(f'Rdac{sfx} {c} vdac 380'); n.a(f'Cdac{sfx} {c} 0 200f')
    n.a(f'RL{sfx} {a} vdac {DAC["RL"]}'); n.a(f'Cp{sfx} {a} 0 35f')
    l1 = n.coil(a, c, DAC['L1'], mode, f'1{sfx}'); l2 = n.coil(c, b, DAC['L2'], mode, f'2{sfx}')
    n.kpair(l1, l2, DAC['L1'], DAC['L2'], DAC['k'], mode, f'd{sfx}'); n.a(f'Cbd{sfx} {a} {b} {DAC["Cb"]:.4e}')
    n.a(f'Cc{sfx} {b} {d0} 1p'); n.coil(d0, d, P['Lin'], mode, f'in{sfx}'); n.a(f'Cd0{sfx} {d0} 0 3f')
    n.a(f'Rbias{sfx} {d} vcm 3Meg'); n.a(f'Croute{sfx} {d} 0 5f')

def out_net(n, o, pad, mode, sfx, rt_to='0'):
    t = 't'+sfx
    la = n.coil(o, pad, P['La'], mode, f'a{sfx}'); lb = n.coil(pad, t, P['Lb'], mode, f'b{sfx}')
    n.kpair(la, lb, P['La'], P['Lb'], P['ko'], mode, f'o{sfx}')
    n.a(f'Cbo{sfx} {o} {t} {FX["Cb_out"]:.4e}'); n.a(f'Cesd{sfx} {pad} 0 100f'); n.a(f'RT{sfx} {t} {rt_to} 45'); n.a(f'CRT{sfx} {t} 0 5f')

def devs():
    I2m = P['I2']; n2 = P['n2']; n3 = n2*P['r3']
    return dict(m1=D.mos(P['n1'], P['j1'], 'p'), mt1=D.mos(P['n1']*P['j1']/60., 60., 'p'),
                m2=D.mos(n2, I2m*1e6/n2, 'p'), m3=D.mos(n3, P['I2']*1e6/n3, 'p'), mt2=D.mos(I2m*1e6/60., 60., 'p'), n3=n3)

def linear_half(mode, fstart=0.1e9, fstop=130e9, npts=520):
    n = NL(); v = devs(); m1, mt1, m2, m3, mt2 = (v[k] for k in ('m1', 'mt1', 'm2', 'm3', 'mt2'))
    n.a('* linear half-circuit, element-for-element = chain.build(pol=p, out=tcoil_padC)')
    n.a('Vdac vdac 0 0'); n.a('Vcm vcm 0 0'); n.a('Iac 0 c AC 1')
    dac_tcoil(n, 0, mode, '')
    n.a(f'G1 Y s1 d s1 {m1["gm"]:.6e}'); n.a(f'Rds1 Y s1 {1/m1["gds"]:.6e}')
    n.a(f'Cgs1 d s1 {m1["cgs"]:.4e}'); n.a(f'Cgd1 d Y {m1["cgd"]:.4e}'); n.a(f'Cdb1 Y 0 {m1["cdb"]:.4e}'); n.a(f'Cw1 s1 0 {D.cwell(P["n1"]):.4e}')
    n.a('Rd Y 0 25'); n.a('Cbyp Y 0 2p'); n.a(f'Rt1 s1 0 {1/mt1["gds"]:.6e}'); n.a(f'Ct1 s1 0 {mt1["cdb"]+3e-15:.4e}')
    n.coil('s1', 'g2', P['Lg'], mode, 'g')
    n.a(f'G2 x s2 g2 s2 {m2["gm"]:.6e}'); n.a(f'Rds2 x s2 {1/m2["gds"]:.6e}')
    n.a(f'Cgs2 g2 s2 {m2["cgs"]:.4e}'); n.a(f'Cgd2 g2 x {m2["cgd"]:.4e}'); n.a(f'Cdb2 x 0 {m2["cdb"]:.4e}')
    n.a(f'Cs2 s2 0 {m2["csb"]+D.cwell(P["n2"]):.4e}'); n.a(f'Rs2 s2 0 {1/(1/1e-4+mt2["gds"]):.6e}'); n.a(f'Ct2 s2 0 {mt2["cdb"]+1e-18:.4e}')
    n.a(f'G3 o x 0 x {m3["gm"]:.6e}'); n.a(f'Rds3 o x {1/m3["gds"]:.6e}')
    n.a(f'Cx x 0 {m3["cgs"]+m3["csb"]+D.cwell(v["n3"]):.4e}'); n.a(f'Co o 0 {m3["cgd"]+m3["cdb"]+3e-15:.4e}')
    out_net(n, 'o', 'pad', mode, '')
    n.a(f'.ac lin {npts} {fstart:.4e} {fstop:.4e}')
    return '\n'.join(n.l)

# ---------------- behavioural full-differential netlist ----------------
def pmos_tab(j_bias):
    """ln(Id/fin) vs Vsg with slope = PGM*(gm/Id)_n  (same definition as device7nm.mos(kind='p')),
    anchored so that Vsg(j_bias) = vgs_n(j_bias)+VGS_P_EXTRA."""
    lj = D._lj; gmid = np.interp(lj, np.log(D.J), D.GM_ID)*D.PGM
    v = np.concatenate([[0], np.cumsum(np.diff(lj)/(0.5*(gmid[1:]+gmid[:-1])))])
    v += (float(D.vgs(j_bias)) + D.VGS_P_EXTRA) - np.interp(np.log(j_bias), lj, v)
    idx = np.linspace(0, len(lj)-1, 90).astype(int)
    pts = [(v[0]-0.3, lj[0]-0.3*gmid[0])] + [(v[i], lj[i]) for i in idx] + [(v[-1]+0.4, lj[-1]+0.4*gmid[-1])]
    return ','.join(f'{a:.5f},{b:.5f}' for a, b in pts)

VSD0 = {'M1': 0.64, 'M2': 0.30, 'M3': 0.50}      # design Vsd (bias_pmos.py) -> lambda factor = 1 at bias

def bpmos(n, name, d, g, s, nfin, j_bias, caps):
    """PMOS: I(s->d) = nfin*j(Vsg) * (1+lam*(Vsd-Vsd0)) * (1-exp(-Vsd/vk)); caps placed physically"""
    m = D.mos(nfin, j_bias, 'p'); lam = m['gds']/m['I']
    vk = float(np.interp(np.log(j_bias), np.log(D.J), D.VDSAT))/3; v0 = VSD0[name[:2]]
    n.a(f'B{name} {s} {d} I = {nfin*1e-6:.6e}*exp(pwl(V({s},{g}),{pmos_tab(j_bias)}))'
        f'*(1+{lam:.4f}*(max(V({s},{d}),0)-{v0}))*(1-exp(-max(V({s},{d}),0)/{vk:.4f}))')
    for (a, b, cval) in caps: n.a(f'C{name}_{a}{b} {a} {b} {cval:.4e}')

def behav_full(analysis, mode='ladder', extra=''):
    n = NL(); v = devs(); m1, mt1, m2, m3, mt2 = (v[k] for k in ('m1', 'mt1', 'm2', 'm3', 'mt2'))
    I1 = P['n1']*P['j1']*1e-6; I3 = P['I2']; n2 = P['n2']; n3 = v['n3']
    n.a('* behavioural full-differential PMOS driver chain'); n.a('Vdd vdd 0 1.8'); n.a('Vdac vdac 0 0.9'); n.a('Vcm vcm 0 0.37'); n.a('Vb3 vb3 0 0.54'); n.a('Vrs1 vrs1 0 0.91'); n.a('Vrcs vrcs 0 1.50')   # gds of ideal current sources referenced to nominal node voltage (no DC error)
    for sfx, sg in (('p', 1), ('n', -1)):
        dac_tcoil(n, 0, mode, sfx)
        d, s1, y, g2, x, o, pad = (q+sfx for q in ('d', 's1', 'Y', 'g2', 'x', 'o', 'pad'))
        bpmos(n, f'M1{sfx}', y, d, s1, P['n1'], P['j1'], [(d, s1, m1['cgs']), (d, y, m1['cgd']), (y, '0', m1['cdb']), (s1, '0', D.cwell(P['n1']))])
        n.a(f'Rd{sfx} {y} 0 25'); n.a(f'Cbyp{sfx} {y} 0 2p')
        n.a(f'Is1{sfx} vdd {s1} DC {I1:.6e}'); n.a(f'Rt1{sfx} vrs1 {s1} {1/mt1["gds"]:.6e}'); n.a(f'Ct1{sfx} {s1} 0 {mt1["cdb"]+3e-15:.4e}')
        n.coil(s1, g2, P['Lg'], mode, f'g{sfx}')
        bpmos(n, f'M2{sfx}', x, g2, 'cs', n2, I3*1e6/n2, [(g2, 'cs', m2['cgs']), (g2, x, m2['cgd']), (x, '0', m2['cdb'])])
        bpmos(n, f'M3{sfx}', o, 'vb3', x, n3, I3*1e6/n3, [(x, '0', m3['cgs']+D.cwell(n3)), (o, 'vb3', m3['cgd']), (o, '0', m3['cdb']+3e-15)])
        out_net(n, o, pad, mode, sfx)
    n.a(f'Itail vdd cs DC {2*I3:.6e}'); n.a(f'Rtail vrcs cs {1/mt2["gds"]*0.5:.6e}'); n.a(f'Ccs cs 0 {2*m2["csb"]+D.cwell(2*n2)+mt2["cdb"]:.4e}')
    n.a(extra); n.a(analysis)
    return '\n'.join(n.l)

if __name__ == "__main__":
    print(linear_half('const')[:400])
