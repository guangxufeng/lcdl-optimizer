# 两阶段输出接口 v2

入口：`result=run_two_stage(case,opts)`；内部数组为NumPy数组。`export_result(result,path,include_matrices=True)`输出严格JSON及可选稀疏NPZ。无限/不可用界导出为null，需结合状态解释。功率MW、无功MVAr、电量MWh、时间小时、成本与输入货币一致；数组默认为实体×时间。

## 顶层字段

| 键 | 内容 |
| --- | --- |
| result_schema_version | 2，显式区别旧模型输出 |
| model / solution_method | 模型语义及精确连续等价方法 |
| status / reason | 总状态与退出原因 |
| feasibility_certified | 全域零松弛可行性证书 |
| cost_optimality_certified | 最坏成本全局界精度证书 |
| robust_certified | 上述两证书均取得 |
| input_case / solver_options | 数据与求解选项快照 |
| axes | 节点/用户/储能ID、功率时间T、状态时间T+1、已定向支路及反向标记 |
| parameters | Ek、用户分配、下游/路径矩阵、电压灵敏度、固定支路Q、建模假设 |
| stage1 / stage2 | 两阶段完整结果 |
| source_sha256 | 导出时算法文件哈希 |

## 第一阶段 stage1

- `U,UG,Eg_mwh,grid_plan_mw`：统一准线、外部电量及固定交换。
- `user_energy_mwh,total_dr_energy_mwh,unified_denominator_mwh,rigid_total_mw,renewable_total_mw`：闭式中间量。
- `alpha,L,target_dr_mw,alpha_squared`：节点修正、准线、目标功率、QP目标。
- `network`：净P/Q、支路P/Q、电压幅值/平方、线路/电压各侧裕度和越限数。
- `energy_residual_mwh,aggregate_alpha_residual_mwh,power_balance_residual_mw`：守恒残差。
- `solver_status,runtime_seconds`：QP状态与时间。

节点QP失败时返回stage1_failed，仍保留此前统一准线中间量。统一式本身未定义时stage1为空并给出原因，未求出的字段不伪造。

## 第二阶段 stage2

`levels`为停止时档位，`rounds`逐轮保存levels、operational、必要时diagnostic，以及发生升档时的upgrade_mask和next_levels。每个运行/诊断检查包含：

| 键 | 内容 |
| --- | --- |
| uncertainty | 当前价格/阈值/β/预算、申报与有效盒界、最小范数可行点和范数 |
| scenarios | 保留的ξ矩阵列表 |
| scenario_results | 逐情景ξ、状态、y、dispatch、成本、残差 |
| history | 各次迭代情景结果、主问题下界、全局可行性/成本搜索输出及原始状态 |
| worst_cost_witness | 找到成本候选偏差后的重新调度；未达到此阶段则无此字段 |
| diagnostic_indicator | 完成诊断后保留情景中各用户最大绝对d |
| lower_bound / upper_bound / gap | 已取得的成本全局界和间隙 |

history中的Phase-I松弛与y是数值可行性检验量，不是可执行物理调度。成本候选偏差未必已被证明为全域最坏情景，必须查看证书和界。全局未完成时，已计算的场景结果与参数仍保留。

## 逐情景 dispatch

访问路径示例：`result['stage2']['rounds'][0]['operational']['scenario_results'][0]['dispatch']`。

- 储能：`charge_mw,discharge_mw,storage_net_injection_mw,mode_z`，S×T。
- 电量：`energy_mwh,energy_with_initial_mwh,soc_fraction`；含初值为S×(T+1)。SOC分母为Emax，零容量设备输出0。
- 成本：`throughput_mwh_per_storage,storage_cost_per_device,cost_components,cost`。
- 响应：`xi,response_degree,squared_response_deviation,response_dr_mw`。
- 松弛：`d,auxiliary_power_mw,adjusted_dr_mw`；运行/执行d=0，诊断值不可执行。
- 网络：`grid_mw,net_p_mw,net_q_mvar,flow_mw,flow_mvar,voltage_pu,voltage_squared_pu`，以及线路/电压上下裕度和越限数。
- 残差：`power_balance_residual_mw,terminal_energy_residual_mwh`。

待机约定z=0，模式不是跨情景计划。输出不包含补损耗购电资源。

## 单独调用与反馈

```python
unified = compute_unified_directrix(case)
stage1 = compute_directrix(case, opts)
scenario = solve_scenario(case, stage1, xi, opts, diagnostic=False)
check = validate_feedback(case, stage1, levels, feedback_p_mw)
dispatch = dispatch_feedback(case, stage1, levels, feedback_p_mw, opts)
```

反馈检查输出ξ、响应度、预算/盒界余量、非负余量及零和/电量残差。周期电量变化返回energy_changed并要求重建第一阶段；仅偏差越集合返回out_of_set。正常反馈重新决定其功率和模式。

## 稀疏约束矩阵

可选导出在 `<文件名>_matrices/operational` 和 `diagnostic`：A/B/C/F.npz，以及vectors.npz中的b、e、c、quad、行缩放尺度。JSON的matrix_exports含绝对路径、y切片及ξ形状。

```text
y = [charge.ravel(), discharge.ravel(), d.ravel()]
A@y <= b+B@xi.ravel()
C@y == e+F@xi.ravel()
J = c@y + sum(quad*y*y)
```

没有h，也没有共用z。此矩阵仅在当前理想效率等前提下与原情景MIQP等价。使用 `explicit_mip=True` 可独立验证情景二进制模型。
