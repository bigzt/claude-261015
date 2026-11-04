"""Half-circuit (differential-mode) small-signal model of DAC node + T-coil + driver.

Signal path (per side):
  Idac -> [DAC T-coil: R_L=35 ohm, L1/L2/k, Cb] -> 1 pF -> d (driver input, 3 Mohm to VCM)
  d -> M1 source follower (drain cascode Mc, gate at VDD) -> s1
  s1 -> Lg (inter-stage series L) -> g2
  g2 -> M2 (RC degeneration Rs || Cs at s2, split tail) -> x -> M3 cascode -> o
  o -> output network (45 ohm to VDD, T-coil / shunt L) -> pad (ESD 100 fF, external 50 ohm)
"""
import numpy as np
from device7nm import mos

W0 = 2*np.pi*75e9          # frequency at which inductor Q is specified
RDC_FRAC = 0.15            # coil loss model: R(f) = R75*(RDC_FRAC + (1-RDC_FRAC)*sqrt(f/75G))
def _rf(sw, r75):
    return r75*(RDC_FRAC + (1-RDC_FRAC)*np.sqrt(np.abs(sw)/W0))

class Net:
    def __init__(s): s.nodes = {'0': -1}; s.el = []
    def n(s, name):
        if name not in s.nodes: s.nodes[name] = len(s.nodes) - 1
        return s.nodes[name]
    def R(s, a, b, r):        s.el.append(('Y', s.n(a), s.n(b), lambda sw, r=r: np.full_like(sw, 1/r)))
    def C(s, a, b, c):        s.el.append(('Y', s.n(a), s.n(b), lambda sw, c=c: sw*c))
    def G(s, d, src, g, gm):  s.el.append(('G', s.n(d), s.n(src), s.n(g), s.n(src), gm))   # I = gm*(Vg-Vs) from d to src
    def L(s, a, b, L, Q=0):   # inductor with constant series R = w0 L / Q
        r = W0*L/Q if Q else 0.0
        s.el.append(('Y', s.n(a), s.n(b), lambda sw, L=L, r=r: 1/(_rf(sw, r) + sw*L)))
    def K(s, a, c, b, L1, L2, k, Q=0):   # centre-tapped coil a-L1-c-L2-b, k>0 aiding
        r1 = W0*L1/Q if Q else 0.0; r2 = W0*L2/Q if Q else 0.0
        s.el.append(('K', s.n(a), s.n(c), s.n(b), L1, L2, k*np.sqrt(L1*L2), r1, r2))
    def solve(s, f, inj, outs):
        sw = 2j*np.pi*np.asarray(f, float); N = len(s.nodes) - 1; nf = len(sw)
        Y = np.zeros((nf, N, N), complex)
        def st(a, b, y):
            if a >= 0: Y[:, a, a] += y
            if b >= 0: Y[:, b, b] += y
            if a >= 0 and b >= 0: Y[:, a, b] -= y; Y[:, b, a] -= y
        for e in s.el:
            if e[0] == 'Y': st(e[1], e[2], e[3](sw))
            elif e[0] == 'G':
                _, d, sr, g, sr2, gm = e
                for (row, sg) in ((d, 1), (sr, -1)):
                    if row < 0: continue
                    if g >= 0: Y[:, row, g] += sg*gm
                    if sr2 >= 0: Y[:, row, sr2] -= sg*gm
            elif e[0] == 'K':
                _, a, c, b, L1, L2, M, r1, r2 = e
                Z = np.zeros((nf, 2, 2), complex)
                Z[:, 0, 0] = _rf(sw, r1) + sw*L1; Z[:, 1, 1] = _rf(sw, r2) + sw*L2; Z[:, 0, 1] = Z[:, 1, 0] = sw*M
                Gm = np.linalg.inv(Z)
                br = [(a, c), (c, b)]
                for p in range(2):
                    for q in range(2):
                        for (n1, s1) in ((br[p][0], 1), (br[p][1], -1)):
                            for (n2, s2) in ((br[q][0], 1), (br[q][1], -1)):
                                if n1 >= 0 and n2 >= 0: Y[:, n1, n2] += s1*s2*Gm[:, p, q]
        I = np.zeros((nf, N, 1), complex); I[:, s.n(inj), 0] = 1
        V = np.linalg.solve(Y, I)[:, :, 0]
        return {o: V[:, s.n(o)] for o in outs}

# ---------------------------------------------------------------------------
DAC = dict(RL=35., L1=120e-12, L2=150e-12, k=0.29, Cb=6e-15)
FIX = dict(Rdac=380., Cdac=200e-15, Cp=35e-15, Cc=1e-12, Rb=3e6, Croute_in=5e-15,
           RT=45., CRT=5e-15, Cesd=100e-15, Rext=1e9, Cb_out=6e-15)   # no external load: only on-chip 45 ohm

def build(p, Q=9.0, dac=DAC, fix=FIX, with_dac=True, drive='dac'):
    """p: driver parameter dict. drive='dac' -> inject Idac; 'vin' -> inject at d (for Zin/gain)."""
    s = Net(); f_ = fix
    if with_dac:
        s.R('c', '0', f_['Rdac']); s.C('c', '0', f_['Cdac'])
        s.R('a', '0', dac['RL']); s.C('a', '0', f_['Cp'])
        s.K('a', 'c', 'b', dac['L1'], dac['L2'], dac['k'], Q); s.C('a', 'b', dac['Cb'])
        if p.get('Lin', 0) > 1e-13:                       # series L between coupling cap and SF gate
            s.C('b', 'd0', f_['Cc']); s.L('d0', 'd', p['Lin'], Q); s.C('d0', '0', 3e-15)
        else:
            s.C('b', 'd', f_['Cc'])
    s.R('d', '0', f_['Rb']); s.C('d', '0', f_['Croute_in'])
    # --- M1 source follower, drain cascode Mc (gate AC-ground) -------------
    m1 = mos(p['n1'], p['j1']); mc = mos(p['n1'], p['j1'])
    s.G('Y', 's1', 'd', m1['gm']); s.R('Y', 's1', 1/m1['gds'])
    s.C('d', 's1', m1['cgs']); s.C('d', 'Y', m1['cgd']); s.C('Y', '0', m1['cdb']); s.C('s1', '0', m1['csb'])
    s.G('0', 'Y', '0', mc['gm']); s.R('0', 'Y', 1/mc['gds']); s.C('Y', '0', mc['cgs'] + mc['csb'])
    mt1 = mos(p['n1']*p['j1']/60., 60.)                    # SF tail at 60 uA/fin
    s.R('s1', '0', 1/mt1['gds']); s.C('s1', '0', mt1['cdb'] + 3e-15)
    # --- inter-stage series inductor --------------------------------------
    s.L('s1', 'g2', max(p['Lg'], 1e-15), Q if p['Lg'] > 1e-13 else 0)
    # --- M2 degenerated, M3 cascode -----------------------------------------
    Ib = p.get('Ib', 0.0)                                  # PMOS bleed current into x
    I2m = p['I2'] + Ib                                     # M2 current = M3 current + bleed
    j2 = I2m*1e6/p['n2']
    m2 = mos(p['n2'], j2); n3 = p['n2']*p['r3']; m3 = mos(n3, p['I2']*1e6/n3)
    if Ib > 0:                                             # cascoded PMOS bleeder at 40 uA/fin
        mb = mos(Ib*1e6/40., 40.)
        s.R('x', '0', 5/mb['gds']); s.C('x', '0', mb['cdb'] + mb['cgd'] + 2e-15)
    s.G('x', 's2', 'g2', m2['gm']); s.R('x', 's2', 1/m2['gds'])
    s.C('g2', 's2', m2['cgs']); s.C('g2', 'x', m2['cgd']); s.C('x', '0', m2['cdb']); s.C('s2', '0', m2['csb'])
    mt2 = mos(I2m*1e6/60., 60.)
    s.R('s2', '0', 1/(1/max(p['Rs'], 1e-4) + mt2['gds']))   # Rs = half-circuit degeneration (0 = common source node, virtual ground); s.C('s2', '0', mt2['cdb'] + p['Cs'] + 1e-18)
    s.G('o', 'x', '0', m3['gm']); s.R('o', 'x', 1/m3['gds'])
    s.C('x', '0', m3['cgs'] + m3['csb']); s.C('o', '0', m3['cgd'] + m3['cdb'] + 3e-15)
    # --- output network ------------------------------------------------------
    f_o = p['out']
    s.C('pad', '0', f_['Cesd']); s.R('pad', '0', f_['Rext'])
    if f_o == 'tcoil_padC':            # o -La- pad(centre) -Lb- t ; RT at t
        s.K('o', 'pad', 't', p['La'], p['Lb'], p['ko'], Q); s.C('o', 't', f_['Cb_out'])
        s.R('t', '0', f_['RT']); s.C('t', '0', f_['CRT'])
    elif f_o == 'tcoil_drainC':        # t -La- o(centre) -Lb- pad ; RT at t
        s.K('t', 'o', 'pad', p['La'], p['Lb'], p['ko'], Q); s.C('t', 'pad', f_['Cb_out'])
        s.R('t', '0', f_['RT']); s.C('t', '0', f_['CRT'])
    elif f_o == 'series_shunt':        # RT(+La shunt) at o, o -Lb- pad
        s.R('o', 't', f_['RT']); s.L('t', '0', max(p['La'], 1e-15), Q); s.C('o', '0', f_['CRT'])
        s.L('o', 'pad', max(p['Lb'], 1e-15), Q)
    elif f_o == 'shunt':               # o = pad, RT + La to VDD
        s.R('o', 't', f_['RT']); s.L('t', '0', max(p['La'], 1e-15), Q); s.C('o', '0', f_['CRT'])
        s.L('o', 'pad', 1e-15)
    inj = 'c' if drive == 'dac' else 'd'
    return s, inj, dict(m1=m1, mc=mc, mt1=mt1, m2=m2, m3=m3, mt2=mt2, n3=n3, j2=j2)

def responses(p, f, Q=9.0, dac=DAC, fix=FIX):
    s, inj, dev = build(p, Q, dac, fix, True, 'dac')
    r = s.solve(f, inj, ['pad', 'd'])
    return r['pad'], r['d'], dev

def driver_only(p, f, Q=9.0, fix=FIX):
    """driver gain V(pad)/V(d) and input capacitance (inject current at d, no DAC)."""
    s, inj, dev = build(p, Q, DAC, fix, False, 'vin')
    r = s.solve(f, inj, ['pad', 'd'])
    zin = r['d']; gain = r['pad']/r['d']
    cin = np.imag(1/zin)/(2*np.pi*np.asarray(f))
    return gain, cin, dev
