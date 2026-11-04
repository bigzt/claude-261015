# ngspice 验证（PMOS 版 driver，`../final_pmos.json`）

器件模型为估算的 7nm 参数（非 PDK）；ngspice 验证的是"同一器件假设下"的网络、偏置与大信号行为。

| 脚本 | 内容 |
|---|---|
| `gen.py` | 网表生成：`linear_half()` 与 `chain.build()` 逐元件对应；`behav_full()` 全差分，PMOS 为 B 源查表 I-V（含饱和过渡），真实偏置 |
| `ladder_fit.py` | 线圈趋肤损耗的 5 段因果 R-L 梯形 |
| `check_linear.py` | A/B：线性网表 vs Python |
| `check_behav.py` | C：`.op` 偏置 + 差分 `.ac` |
| `check_tran.py` | C：PAM4 瞬态（DAC 满量程 0–5.2 mA/边），眼图与瞬态四端电压 |
| `plots_ng.py` | 生成 `ngspice_verification.png` |

## 结果
| 检查 | 结果 |
|---|---|
| A 线性等效（Q=∞ / 恒定 R） | 与 Python 幅度差 ≤ 0.0001 dB，相位差 ≤ 0.001° |
| B 趋肤损耗（因果梯形） | −3 dB 83.2 GHz（Python 83.7），纹波 1.36 dB（Python 1.05） |
| C `.op` | SF 10.96 mA、输出级 15.50 mA/边；s1 0.905、CS 1.500、X 1.207、o 0.751 V |
| C `.ac` 差分 | driver 增益 1.95 倍；全链路 −3 dB 84.2 GHz；纹波 1.29 dB |
| C PAM4 91 GBd | 输入 345 mVppd → 输出 677 mVppd（2.03 倍）；三只眼高 196–201 mV |
| C PAM4 182 GBd（1 sps，无 FFE） | 输出眼闭合（眼高 −61…−77 mV），需 DSP 预加重/FFE |
| C 瞬态四端电压 | 最差 0.77 V（M3p GS）≤ 0.9 V；M2p 最小 Vsd 0.19 V（Vdsat 0.17，余量很小） |

## 发现
1. 线圈若为恒定 R（Q∝f）损耗，同一参数只有 78 GHz —— 必须用 EM 提取的 Q(f)。
2. pFET 建模方式不同（gm/Id 是否下降 15%）会使增益在 2.0–2.8 倍之间变化，纹波达 2.4 dB —— 对 pFET 参数高度敏感。
3. DAC 满量程下输出只有约 680 mVppd（输入约 340 mVppd），达不到 800 mVppd。
4. M2p 在 ISI 峰值处 Vsd 降到 Vdsat 附近，输出级电压余量是瓶颈（可把尾电流管做大以降低其 Vdsat，换约 50 mV）。
