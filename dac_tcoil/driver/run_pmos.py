import json
from multiprocessing import Pool
import design
design.POL = 'p'                      # set before fork -> workers inherit
from design import run, met, KEYS
if __name__ == "__main__":
    jobs = [('tcoil_padC', 6.02, 85, 11, 9.0), ('tcoil_padC', 6.02, 85, 23, 9.0), ('tcoil_padC', 6.02, 85, 37, 9.0), ('series_shunt', 6.02, 85, 11, 9.0)]
    with Pool(4) as pl: res = pl.map(run, jobs)
    json.dump([(a, x) for a, x, m in res], open('design_p.json', 'w'))
    for a, x, m in res:
        v = dict(zip(KEYS, x)); mq = [met(x, a[0], Q=q*9.0) for q in (0.85, 1.15)]
        print(f"{a[0]:12s} seed{a[3]}: gain={10**(m['gdb']/20):.2f}x BW3={m['bw']:.0f}G (Q-15%:{mq[0]['bw']:.0f} Q+15%:{mq[1]['bw']:.0f}) "
              f"ripple={m['rip']:.2f}dB (worst {max(m['rip'],mq[0]['rip'],mq[1]['rip']):.2f}) Itot={m['itot']:.1f}mA Cin={m['cin']:.1f}fF | "
              f"I1={v['n1']*v['j1']/1e3:.1f} I3={v['I3']*1e3:.1f} mA/side j1/j2/j3={v['j1']:.0f}/{v['j2']:.0f}/{v['j3']:.0f}")
