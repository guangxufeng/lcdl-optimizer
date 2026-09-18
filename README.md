# LCDL Optimizer

A reusable Python implementation of two-stage load-directrix and ideal-storage robust optimization for radial distribution networks.

可复用的“负荷准线优化—储能响应补偿”两阶段算法。输入自定义网架、负荷、新能源、储能及响应参数，输出统一准线、节点准线和逐情景储能调度，以及完整的全局验证状态。

仓库只包含通用算法、接口文档、合成最小数据、运行示例及通用测试；不依赖研究论文附件或特定IEEE网架数据。

## 安装

Python 3.10及以上，使用Gurobi求解。Gurobi及其许可不包含在本项目中，请自行取得适用于你的用途和模型规模的许可。MIT仅覆盖本仓库代码与文档，不改变第三方依赖的许可条件。

在克隆后的仓库根目录执行：

```bash
python -m venv .venv
# 激活虚拟环境：Windows PowerShell用 .venv\Scripts\Activate.ps1
# Linux/macOS用 source .venv/bin/activate
python -m pip install -e ".[dev]"
python examples/run_minimal.py
python -m pytest -q
```

## 接口示例

```python
import numpy as np
from lcdl_algorithm import (
    Case, SolverOptions, run_two_stage, export_result, dispatch_feedback,
)

case = Case.from_json("data/minimal_case.json")
options = SolverOptions(time_limit=60, max_iterations=30)
result = run_two_stage(case, options)
print(result["status"], result["reason"])

if result["feasibility_certified"]:
    # 此反馈仅用于随仓库附带的两时段、单用户示例。
    # 自定义算例时替换为用户确认的K×T负荷功率矩阵，单位MW。
    feedback = np.array([[0.44, 0.36]])
    execution = dispatch_feedback(
        case, result["stage1"], result["stage2"]["levels"], feedback, options,
    )
    result["execution"] = execution
    if execution["status"] == "optimal":
        print(execution["dispatch"]["mode_z"])

export_result(result, "results/run/complete.json", include_matrices=True)
```

最小例子的完整脚本在[examples/run_minimal.py](examples/run_minimal.py)，输入在[data/minimal_case.json](data/minimal_case.json)。预计得到U=L=[0.5,0.5]，反馈调度先放电0.04 MW再充电0.04 MW，初末电量都为1 MWh，吞吐成本0.08。该例是合成教程，不是实测数据或大规模性能基准。

## 更换数据

详见[输入接口](docs/data_interface.md)。`Case(data)`接受Python字典，`Case.from_json(path)`接受UTF-8 JSON。数据数组按“实体×时间”排列：N节点、K用户、S储能、T时段。

| 数据 | 主要字段 |
| --- | --- |
| 网络 | bus_ids、root_bus、base_kv、branches；支路阻抗Ω、有功限额MW |
| 时序 | dt_hours，N×T的rigid_p_mw、fixed_q_mvar、renewable_p_mw |
| DR用户 | N×K的dr_allocation，K×T的dr_pre_mw，周期电量须为正 |
| 储能 | 接入节点、功率MW、电量MWh、初始电量、非负吞吐成本；效率必须为1 |
| 网络界 | 电压p.u.、根电压、交换功率上下界 |
| 响应集合 | 四档价格、K×4响应阈值、beta、K×T偏差上下界、诊断惩罚 |

刚性负荷必须已经扣除DR部分，避免重复计量；净无功包括负荷减去给定无功注入。一个情景ξ是所有用户整个周期的偏差矩阵，不是一个时段的单个随机数。

## 输出和分阶段调用

完整字段和维度见[输出接口](docs/output_interface.md)。

- `result['stage1']`：U、UG、Eg、L、alpha、目标功率、固定交换、网络状态及守恒残差。
- `result['stage2']['rounds']`：逐轮档位、响应预算、运行/诊断情景、模式、功率、SOC、成本及全局搜索记录。
- `axes`和`parameters`：实体顺序、时刻、已定向支路、网络矩阵和数据快照。
- `export_result(..., include_matrices=True)`：输出JSON和A/B/C/F稀疏矩阵、成本系数、缩放尺度及变量切片。

单独调用入口：`compute_unified_directrix`、`compute_directrix`、`UncertaintySet`、`solve_scenario`、`solve_robust`、`validate_feedback`、`dispatch_feedback`。单情景接口不自动证明全域鲁棒性；模拟ξ时先用`UncertaintySet.contains`检查归属。正式升档推荐使用总入口`run_two_stage`。

## 适用范围与结果判定

当前范围：单电压等级、固定辐射拓扑、无损LinDistFlow、理想储能、已知净无功，以及执行前确认完整周期反馈。支持不同规模、非连续节点编号、乱序/反向支路、独立线路上下界及零台储能。网状网络、有损储能、AC/三相潮流和实时非预知策略需要另行扩展模型。

储能模式按情景自适应。默认利用理想效率和非负成本的严格等价性求连续模型，再恢复互斥充放电模式；可用`solve_scenario(..., explicit_mip=True)`交叉检查显式二进制模型。推导见[模型说明](docs/model.md)。

**抽样全部通过不等于鲁棒证书。**

| 字段/状态 | 含义 |
| --- | --- |
| feasibility_certified | 完整申报集合的零松弛可行性已通过 |
| cost_optimality_certified | 最坏成本的全局界满足精度 |
| robust_certified / robust_optimal | 上述两项均通过 |
| unverified | 尚未证实，例如全局搜索超时；不等于不可行，不因此升档 |
| stage1_failed | 第一阶段失败，已算中间参数保留 |
| invalid_uncertainty | 申报集合为空或不相容 |
| proven_infeasible | 当前被检验模型有不可行证明，需结合所在轮次和诊断/运行作用域解释 |
| no_feasible_scheme / stalled | 无可升级用户或诊断指标没有触发升档 |
| 反馈out_of_set / energy_changed | 重建响应范围，或重建两个阶段 |

`time_limit`是每次求解器调用的限时，不是整个流程总限时。非凸全局检查可能较慢，代码不保证任何给定规模在固定时间内收敛。两时段小问题存在精确区间顶点加速分支，一般多时段问题仍使用全局非凸搜索。

## 许可

本项目使用[MIT License](LICENSE)。论文级使用应自行核对建模假设、数据来源、求解状态和数值残差，不应将合成示例表现外推为真实系统保证。
