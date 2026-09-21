# v3 申报集合鲁棒模型输入接口

本接口依据模型 v3，使用 `schema_version=3`。`Case(dict)` 或 `Case.from_json(path)` 读取输入；JSON 为 UTF-8。N 为节点数、K 为用户数、S 为储能数、T 为等长时段数、M 为激励档位数。二维数组始终按“实体、时段”排列。最小完整输入见 `data/minimal_case.json`。

## 网络与时序

| 字段 | 格式与限制 |
| --- | --- |
| schema_version / name | 3 / 算例名称字符串 |
| bus_ids / root_bus | N 个不重复节点 ID，root_bus 必须在其中；N≥2，可非连续或乱序 |
| base_kv / dt_hours | 单电压等级 kV / 正的等长时段长度 h；T≥2，不限定 24 时段 |
| branches | N−1 条连通树支路；每条有 from_bus、to_bus、r_ohm、x_ohm |
| 支路有功界 | p_limit_mw 表示对称界；或提供 p_min_mw、p_max_mw。按输入方向解释，自动根定向时同步转换 |
| rigid_p_mw | N×T 非负刚性负荷 MW，必须已扣除 DR 部分 |
| renewable_p_mw | N×T 非负可用新能源 MW，全额消纳 |
| fixed_q_mvar | N×T 已知净无功 Mvar，包括负荷减去给定注入，可为负 |
| dr_pre_mw | K×T 未响应预测负荷 MW，非负，每用户周期电量必须为正 |
| dr_allocation | N×K，每列只有一个 1，其余为 0；每用户固定接入一个节点 |
| root_voltage_pu | 根节点电压幅值 p.u. |
| voltage_min_pu / voltage_max_pu | 标量或长度 N 的幅值界，代码内部平方 |
| grid_min_mw / grid_max_mw | 标量或长度 T 的交换功率界，正值为进口 |
| user_ids / storage_ids | 可选，长度 K / S 的不重复标识 |

根节点是纯外部交换节点，刚性负荷、无功、新能源和 DR 分配的对应行必须为零，储能也不能放在根节点。schema 3 的节点诊断要求单节点用户；跨节点聚合体应拆为多个固定节点用户。

## 储能

`storage` 为 S 个字典的列表，允许 `[]`。每项必须包含：`bus`、`p_charge_mw`、`p_discharge_mw`、`e_min_mwh`、`e_max_mwh`、`e_initial_mwh`。功率 MW，电量 MWh，均非负，初始电量在容量界内。

`eta_charge`、`eta_discharge` 默认 1，当前模型只接受 1。终端电量固定为初始电量。设备级旧字段 `cost_per_mwh` 只为兼容保留，在 schema 3 中不参与目标；每台设备必须提供非负有限 `throughput_cost_coefficient`，对应式（31）的 c_e。

## 激励、申报和成本

| 字段 | 格式与限制 |
| --- | --- |
| incentive_prices_per_mwh | 长度 M 的非负严格递增电价，M≥1，不再固定为四档 |
| response_degrees | K×M，rho∈[0,1]，每行随档位非递减 |
| deviation_factors | K×M，s∈[0,1)，不强制随档位单调 |
| storage[e].throughput_cost_coefficient | 每设备非负系数，式（31）的 c_e |

鲁棒可行性目标为最坏情况下最小全网虚拟量之和，成本单独求解为 `sum_e c_e*sum_t(charge_mw+discharge_mw)`，没有额外 dt。若成本按货币/MWh 给定，应先乘 dt。v2 草稿的 slack_penalty、全局 throughput_coefficient 必须删除，代码会明确拒绝。两张响应表是算例假设或事先收集的用户参数，不是算法拟合得到的真实响应。激励电价用于发布和记录，不另外加到式（31）的储能目标。

## 调用与真实申报接口

```python
case = Case.from_json('data/minimal_case.json')
result = run_two_stage(case, SolverOptions(time_limit=60))
```

未传回调时，每一轮根据表中 rho、s 生成式（18）的申报曲线，输出明确标记 `preset_response_table`。真实系统可使用 `declaration_provider(context)`：

```python
def provider(context):
    # 从业务系统读取此轮已经取得的申报；算法本身不会给用户发送消息。
    # context: round_index, levels, target_dr_mw, forecast_dr_mw,
    #          prices_per_mwh, previous_round
    return {
        'response_degree': [0.8],       # K 个值
        'deviation_factor': [0.1],      # K 个值
        'declared_p_mw': [[0.46, 0.34]] # 可省略，提供时须符合式（18）
    }

result = run_two_stage(case, declaration_provider=provider)
```

回调返回 `None` 表示尚未收到申报：返回 `awaiting_declaration` 和 `stage2.requested_levels`，不假装已收到反馈。业务层可以保存结果，在收到数据后用 `initial_levels` 启动下一次调用；每次调用的轮次历史独立，应由业务层关联跨调用记录。单次循环内会检查重新申报的 rho 非递减。回调异常的业务处理由调用方负责，格式/模型数据错误返回 `invalid_declaration`。

反馈必须提供 K×T 的已确认整周期负荷 MW。schema 3 的 `validate_feedback`、`dispatch_feedback` 第三个参数是完整的 `stage2.declaration`，不能只传档位，以免误用预设表替代真实申报。

```python
execution = dispatch_feedback(case, result['stage1'],
                              result['stage2']['declaration'], [[0.47, 0.33]])
# 或一次性完成：run_two_stage(case, feedback_p_mw=[[0.47, 0.33]])
```

## 显式迁移其他算例

```python
from lcdl_algorithm import with_declaration_model
new_data = with_declaration_model(old_data,
    response_degrees=rho_table, deviation_factors=s_table,
    storage_cost_coefficients=cost_per_mw_step_per_device)
case = Case(new_data)
```

这会复制原拓扑、负荷和设备参数，移除旧鲁棒集合字段，但不会推算 rho、s。换成其他拓扑时同步更新所有 N 维数组、支路端点和设备节点；换成其他时段数时更新所有 T 维数组。网状网络、多电压等级、有损储能、未知无功或实时滚动反馈需要扩展模型，不能仅靠换数据支持。

## 全域求解选项

`SolverOptions(declaration_vertex_limit=256)` 控制完整顶点枚举的安全规模上限，0 强制使用全局对偶检验。超限不会截取部分顶点后声称鲁棒，而是切换方法。`optimize_robust_cost=False` 只验证鲁棒可行性，跳过可选的式（31）。

`feasibility_tol` 是全域虚拟量零判定容差，`upgrade_tol` 是同一联合解中 H 指标的升档容差，二者均针对按时段求和的 MW 数值。`cost_abs_tol`、`cost_rel_tol` 用于成本上下界。`time_limit` 对每次内层或全局求解生效，不是整个流程总时限。`max_iterations`、`two_period_vertex_limit` 属于旧 schema 1/2 路径，不控制 v3 的完整枚举或全局对偶。

不确定量 delta 使用 MW，不再使用旧模型的无量纲 xi 或二次预算。每轮新申报中心及 s 都会重新构建完整集合。
