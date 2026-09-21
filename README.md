# LCDL Optimizer

Reusable Python two-stage load-directrix and ideal-storage robust correction for radial networks.

0.3.0 对应《模型v3 联合补偿修订版》：**节点准线 → 用户申报曲线与偏差集合 → 全域储能鲁棒校正 → 同一联合解的节点缺口诊断升档 → 确认反馈复核**。算法不依赖 IEEE33，拓扑、负荷、储能与用户参数通过数据接口传入。

新输入使用 `schema_version=3`。历史 schema 1/2 的鲁棒接口仍兼容；未发布的 v2 单曲线草稿已被替换。

## 安装运行

Python 3.10+，需要适用的 Gurobi 许可；本仓库的 MIT 许可不覆盖 Gurobi。

```bash
git clone https://github.com/guangxufeng/lcdl-optimizer.git
cd lcdl-optimizer
python -m venv .venv
# Windows PowerShell: .venv\Scripts\Activate.ps1
# Linux/macOS: source .venv/bin/activate
python -m pip install -e ".[dev]"
python examples/run_minimal.py
python -m pytest -q
```

运行自定义算例，并保存完整过程与约束矩阵：

```bash
python run_declaration.py --case data/minimal_case.json --output results/run/result.json
# feedback.json 为 K×T 已确认整周期负荷 MW 数组
python run_declaration.py --case data/minimal_case.json --feedback feedback.json
# 只检验鲁棒可行性，跳过可选最坏成本优化
python run_declaration.py --case data/minimal_case.json --skip-cost
```

## Python 接口

```python
from lcdl_algorithm import Case, SolverOptions, run_two_stage, export_result

case = Case.from_json('data/minimal_case.json')
result = run_two_stage(case, SolverOptions(time_limit=60),
                       feedback_p_mw=[[0.47, 0.33]])  # 仅用于随附最小例子
print(result['status'], result['feasibility_certified'], result['executable'])
print(result['stage1']['U'], result['stage1']['L'])
for row in result['stage2']['rounds']:
    eta = row['operational']['feasibility']
    print(row['levels'], eta['lower_bound'], eta['upper_bound'])
export_result(result, 'results/run/result.json', include_matrices=True)
```

省略 feedback 时，即使已证明鲁棒可行，顶层 executable 仍为 false，等待真实反馈。`declaration_provider(context)` 可对接实际申报；未传回调时使用数据表模拟申报并标记来源，算法不会代替用户实际申报或发送通知。

- [输入格式、单位、申报回调、迁移和求解选项](docs/data_interface.md)
- [完整输出、上下界、情景调度和矩阵导出](docs/output_interface.md)
- [公式对应、全域验证方法和修订式（38）的含义](docs/model.md)

分阶段接口包括 `compute_unified_directrix`、`compute_directrix`、`make_declaration`、`DeclaredUncertaintySet`、`DeclarationSystem`、`solve_declared_dispatch`、`declared_global_oracle`、`solve_declared_robust`、`joint_gap_diagnostic`、`run_incentive_loop`、`validate_feedback` 和 `dispatch_feedback`。

## 最小例子与判定

[合成最小算例](data/minimal_case.json)为两个节点、两个时段、一个用户、一台储能。预测负荷 [0.7,0.1] MW，目标 [0.4,0.4] MW；四档 rho=[0,0.5,0.8,1]，s=[0.1,0.1,0.1,0.05]。

| 档位 | 申报中心 MW | 全域最坏最小虚拟量 eta | 处理 |
| --- | --- | --- | --- |
| 1 | [0.7,0.1] | 0.42 | H 超过容差，升档 |
| 2 | [0.55,0.25] | 0.15 | H 超过容差，升档 |
| 3 | [0.46,0.34] | 0 | 鲁棒可行 |

eta 是全节点、全时段虚拟功率之和。第三档最坏储能成本 J=0.188；示例反馈 [0.47,0.33] MW 的实际成本为 0.14，储能先放电 0.07 MW、再充电 0.07 MW，电量 [1,0.93,1] MWh。这是合成教程，不是实测数据或性能基准。

`feasibility_certified`/`robust_certified` 表示当前整个申报集合零虚拟量可行；`cost_optimality_certified` 独立表示可选 J 成本界闭合；`executable` 要求确认反馈物理复核通过。成本尚未收敛但可行性已证明时，明确返回 `robust_feasible_cost_unverified`。

## 重要模型含义

申报中心为 Pdec=Ppre+rho(Ptar−Ppre)，不确定量 delta 使用 MW，满足 ±sPdec 盒界、每用户周期零和及实际负荷非负。新一轮中心移动时集合未必嵌套，必须重新验证。

鲁棒可行性先单独最小化虚拟量，再检查其全域最坏值；成本在零虚拟量条件下另解，不使用有限惩罚权衡安全与成本。每设备 `throughput_cost_coefficient` 对应式（31）的 c_e，按原式不额外乘 dt；从货币/MWh 报价换算时需先乘 dt。

**修订后的 H 直接读取同一套最坏情景联合补偿解。** 全网总缺口最小后，保持该最小值不变，再联合最小化节点时段缺口平方和，统一处理多解。所有 H 来自同一个 delta、同一套储能和网络解；同节点用户共享指标。H 用于升档，不代表责任归因。全过程见 `operational.joint_diagnostic`，保留情景、储能轨迹、节点缺口、H、求解状态与矩阵重建接口。

只有最坏情景及联合缺口分配均验证通过才自动升档。全局界未闭合时保留反例与诊断，返回 `unverified`。若缺口只在无用户节点或用户指标均低于容差，返回 `stalled`；若所选用户已达最高档，返回 `no_feasible_scheme`，不据此证明其他档位组合无解。

小规模集合完整枚举多面体所有顶点；较大集合用全局非凸对偶/强对偶检验。抽样、单条申报可行、没找到反例都不等于全域证明。高维检验可能超时，返回 `unverified` 和现有上下界。

## 换数据和旧版迁移

用 `with_declaration_model(old_data, rho_table, s_table, storage_cost_coefficients)` 保留物理数据并显式更换响应参数；旧 beta、阈值和 xi 不能直接解释为 rho、s。schema 3 的反馈接口第三参数是完整申报字典，不再只是档位。

支持单电压等级、固定辐射拓扑、无损 LinDistFlow、理想储能、给定净无功、单节点用户和执行前整周期反馈；支持任意有限激励档位数、非连续节点 ID、反向支路、独立线路界、不同等长时段与零台储能。网状网络、多电压等级、AC/三相、有损储能和实时非预知反馈需要扩展模型。

仓库使用 [MIT License](LICENSE)，仅包含通用算法、接口文档、合成示例和测试；不包含研究文档附件、原项目实验数据或计算结果。
