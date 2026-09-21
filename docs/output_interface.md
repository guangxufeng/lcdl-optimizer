# v3 输出与全过程接口

`run_two_stage` 返回 `result_schema_version=3`。统一准线、节点准线、完整输入参数、每轮申报和不确定集、全域证明界、情景调度、同一联合解的 H 诊断及反馈调度均保留。没有可用解时不会用零数组伪造结果。

## 顶层与第一阶段

| 路径 | 含义 |
| --- | --- |
| status / reason | 状态及原因 |
| feasibility_certified / robust_certified | 当前完整申报偏差集合已证明零虚拟量可行；v3 两者同义 |
| cost_optimality_certified | 可选的式（31）最坏成本上下界已闭合，独立于可行性 |
| declaration_feasible | 兼容字段，在 v3 表示整个申报集合鲁棒可行，而非只检查中心 |
| executable | 只有确认反馈已重新求解并复核通过时才为 true |
| input_case / solver_options | 完整数据及参数快照 |
| axes | 节点/用户/储能 ID、定向支路、T 个功率时刻和 T+1 个电量时刻 |
| parameters | 用户电量、用户节点、分配/路径/下游矩阵、电压灵敏度、固定无功潮流、假设 |
| stage1.U / UG / Eg_mwh | 统一准线、外部电量时段分配、外部电量 |
| stage1.L / alpha / target_dr_mw | K×T 节点准线、修正量、目标功率 |
| stage1.grid_plan_mw | T 个固定交换功率 |
| stage1.network | 节点净负荷、电压与平方、支路 P/Q、线路和电压裕度 |
| stage1.*residual* | 电量、逐时加权准线修正、功率平衡残差 |

第一阶段失败时保留已计算的统一准线参数。

## 每轮完整记录

`stage2.rounds` 中每项包含轮次、档位、申报来源、申报快照、`operational` 和适用时的升档掩码、下一轮档位。申报快照保存 rho、s、Pdec、功率上下界、电量残差、相对目标偏差、电价。预设表与实际回调分别标记来源。

`operational.feasibility` 是式（28）的全网 eta 检验。每个全域检验对象包含：

- `method`：完整多面体顶点枚举或全局对偶检验。
- `lower_bound / upper_bound / gap / bounds_consistent`：合法情景提供下界，全域求解提供上界。
- `zero_certified / positive_certified / optimality_certified`：零判定、已证明正缺口、上下界精度。
- `uncertainty`：Pdec 中心、delta 的 MW 上下界、每用户零和关系和集合类型。
- `witness`：最强已验证情景、delta、内层 y、虚拟量和调度结果。
- `scenario_results`：所有实际求解过的情景及其完整结果。
- `solver`：全局状态、运行时间、求解界，或预期与完成的全部顶点数量。

只有通过全域零判定后，才可能出现 `operational.cost`，表示式（31）的 J 检验。

eta 已证明为正时，`operational.joint_diagnostic` 保存修订式（38）的联合诊断：

- `delta_mw`：K×T，所有用户共同对应的同一个偏差情景。
- `primary_objective`：固定场景的最小总虚拟量。
- `selection_rule / secondary_objective`：保持最小总量的条件下，最小化所有节点时段缺口幅值平方和。
- `y / dispatch`：这一套联合解的完整变量、储能调度、全节点虚拟量与网络结果。
- `H_mw`：K 个用户节点跨时段缺口，同节点用户共享指标；无可用联合解时不伪造该字段。
- `selection_certified / max_constraint_violation`：二次优化是否完成及约束残差。
- `scenario_is_certified_worst`：eta 全局界是否闭合；false 时仅为已验证反例。
- `upgrade_certified`：最坏情景和联合缺口分配都验证通过，才允许自动升档。
- `model_export`：启用矩阵导出后，关联该轮 eta 矩阵、固定 delta、总量等式及二次目标，足以重建联合诊断 QP。

不再输出旧草稿的 `user_diagnostics`、`H_lower_bound / H_upper_bound`。H 是选定联合解的节点分量，不是独立节点 max-min 的上下界。可升档轮次另保存 `round.H_mw`、`needs_upgrade`、`upgrade_mask`。同节点用户的 H 不应重复累加为全网缺口。

等待新申报时，`stage2.requested_levels` 是需要的新档位；`stage2.declaration` 仍可能是最后一次旧申报，其自身 levels 是唯一有效归属。每轮集合重新验证，不假设嵌套。

## 单情景调度 dispatch

| 字段 | 维度 / 内容 |
| --- | --- |
| charge_mw / discharge_mw | S×T 正的充/放电功率 |
| mode_u / mode_z | S×T 同义模式字段；1 充电，0 放电，待机规范为 0 |
| storage_net_injection_mw | S×T，放电减充电 |
| energy_mwh / energy_with_initial_mwh | S×T 时段末电量 / S×(T+1) 含初值电量 |
| soc_fraction / terminal_energy_residual_mwh | SOC / S 个末端电量残差 |
| response_dr_mw | K×T，此情景的 Pdec+delta 或已确认反馈 |
| virtual_increase_mw / virtual_decrease_mw | N×T，v_plus、v_minus；根节点为零 |
| virtual_adjustment_mw / epsilon_mw | N×T，正负差与其绝对值 |
| node_virtual_sum_mw | N，每节点跨时段虚拟量总和，仅为此情景结果 |
| diagnostic_indicator_mw | K，此情景用户节点总量；仅 joint_diagnostic 中的值用于修订 H |
| physical_net_p_mw / physical_network | 不含虚拟调节的实际净负荷与网络状态 |
| net_p_mw / flow_mw / flow_mvar | 含虚拟调节的节点净负荷 / 支路有功与无功 |
| voltage_pu / voltage_squared_pu / *margin* | 电压及网络裕度 |
| power_balance_residual_mw / physical_power_balance_residual_mw | 含/不含虚拟量的平衡残差 |
| virtual_cycle_residual_mwh | 每节点虚拟电量，v3 不要求逐节点为零，作为过程参数记录 |
| throughput_mwh_per_storage / storage_cost_per_device | S 个吞吐电量和设备成本 |
| objective_kind / cost / cost_components | 本次目标、目标值、储能成本及虚拟量；诊断的 cost 是虚拟量，不是货币成本 |

虚拟量不属于可执行设备指令。内层 `executable` 只表示这一具体曲线满足物理校验，不能代替整个集合的鲁棒证书；流程顶层会额外要求全域证明和确认反馈。

## 状态

| 状态 | 含义 |
| --- | --- |
| robust_optimal | 鲁棒可行，且 J 成本界通过 |
| robust_feasible | 鲁棒可行，主动跳过可选成本优化 |
| robust_feasible_cost_unverified | 鲁棒可行已证明，但 J 尚未闭合 |
| ready_to_execute | 已确认反馈通过集合与物理复核 |
| stalled | eta>0，但联合 H 在用户节点均未超过容差；检查无用户节点和容差，不执行 |
| no_feasible_scheme | 当前非鲁棒且所选用户已在最高档，本规则无法继续；不证明其他档位组合无解 |
| unverified | 全域可行性或升档判断所需界未闭合；不推断安全或升档 |
| awaiting_declaration / invalid_declaration | 尚未收到新申报 / 申报格式或模型关系不合法 |
| stage1_failed | 节点准线失败，保留已算中间参数 |
| out_of_set / energy_changed | 反馈越出申报集合 / 周期电量变化，分别重建申报范围或两个阶段 |

`proven_not_robust` 是轮次 operational 的状态。单情景解的 `optimal`、`requires_upgrade`、`infeasible` 只描述该情景；自动升档必须使用经过全域最坏情景验证和联合二次优化验证的 H。

## 反馈与矩阵

`execution.feedback_check` 包含 delta_mw、集合归属、上下界裕度、反馈电量和电量残差。确认反馈的储能求解强制所有虚拟量为零。

`export_result(result, path, include_matrices=True)` 输出严格 JSON 和每轮 eta、成本及已算反馈的稀疏矩阵目录。每个目录含 A、B、C、F 的 NPZ 文件，vectors.npz 含 b、e、c。

约束为 `A@y<=b+B@delta`、`C@y=e+F@delta`，delta 为相对于该轮申报中心的 MW 偏差；反馈目录以反馈为基准，使用 delta=0。y 按 charge、discharge、virtual_increase、virtual_decrease 展平，切片与约束行分组在 JSON 中给出。

显式二进制部分单独导出：u∈{0,1}，charge≤Pch_max*u，discharge≤Pdis_max*(1−u)。完整 MIP 需要这些关系；A/B/C/F 是等价连续部分。矩阵路径相对结果 JSON 所在目录，便于迁移。`DeclarationSystem` 可直接访问矩阵和目标向量。联合诊断复用 eta 的 A/B/C/F，再按 joint_diagnostic.model_export 加入固定总缺口等式和二次目标；所有诊断变量及偏差在结果中保留。
