# DAC 输出 Driver 设计（第一版，基于估算器件参数）

## 设计条件
- 工艺：7nm FinFET，core 管；VDD = 1.8 V；可靠性按任意两端 ≤ 0.9 V（**假设值，请按 PDK 核对**）
- 差分；输入 350–400 mVppd，输出 800 mVppd；driver DC 增益（输入→PAD）2.0 倍，可调到约 1.5 倍
- 负载：仅片上 45 Ω 匹配电阻（PAD 外无负载），ESD+PAD 100 fF/边
- DAC 侧：R_L = 35 Ω，Idac 单端最大 5.2 mA，Cdac 200 fF，380 Ω，1 pF 耦合
- 全链路（Idac→PAD）−3 dB ≥ 80 GHz，DC–70 GHz 纹波 ≤ 1.5 dB；线圈 Q = 9 @75 GHz，损耗按 R∝√f（DC 电阻为 75 GHz 处的 15%）

## 拓扑（每边，另一边镜像）
DAC T-coil → 1 pF → **Lin** → **M1 源跟随器**（漏端共源共栅 Mc，栅接 VDD；尾电流共源共栅）→ **Lg** → **M2 共源**（差分对共用尾电流）→ **M3 共源共栅** → **La**（漏端→PAD）→ PAD/ESD → **Lb** → **R_T 45 Ω** → VDD

## 参数
| 器件 | 值 |
|---|---|
| DAC T-coil（重调） | L1 = 110 pH，L2 = 105 pH，k = 0.45，Cb ≈ 6 fF，R_L = 35 Ω |
| Lin / Lg | 65 pH / 70 pH |
| La / Lb（不耦合） | 105 pH / 160 pH |
| M1、Mc | 230 fin，28 µA/fin，6.4 mA/边 |
| M2 | 267 fin，60 µA/fin |
| M3 | 320 fin，50 µA/fin |
| 输出尾电流 | 32 mA（16 mA/边） |
| 偏置 | VCM_in = 1.47 V，VB3 = 1.24 V，X = 0.62 V，Vout_cm = 1.08 V |
| 功耗 | 44.9 mA × 1.8 V ≈ 81 mW（不含偏置电路） |

## 仿真结果（`design.py` / `verify.py` / `linearity.py` / `bias.py`）
- 增益 2.00 倍；全链路 −3 dB 88 GHz（Q ±15%：86–90 GHz）；纹波 0.86 dB（Q ±15% 最差 1.18 dB）
- driver 输入电容约 23 fF
- 增益调节：尾电流 22–32 mA → 1.4–2.0 倍
- 估算 HD3 ≈ −53 dBc @400 mVppd（理想化 I-V，偏乐观）
- 所有 core 管任意两端电压 ≤ 0.66 V，全部饱和

## 风险与待确认
- 器件寄生 +20%：−3 dB 降到 78 GHz，纹波 2.4 dB —— 最大风险，需 PDK + 寄生提取验证
- ESD/Cdac/R_T ±15–20%：纹波约 1.4–1.6 dB，贴近指标
- 群延迟 10–20 ps，80 GHz 附近上升约 8 ps
- 线圈 Q 若为 6：需另一套参数，带宽 80–86 GHz，纹波偏紧（≈1.7 dB）
- **器件参数是估算值**，不是 PDK 数据；按 gm/Id + 电流目标用 PDK 重新定尺寸
