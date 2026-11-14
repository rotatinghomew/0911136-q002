# 不完整观测年代下的暴雨重现期与联合超越风险

本项目组织领域资料并提供当前摘要，通过本地 HTTP 供数服务连接应用。业务代码在 `stormcast`，字段与单位约定见 `docs/data-contract.md`。示例均为本地合成记录，不代表真实县域记录或业务结论。

## 环境

Python 3.11 及以上，业务运行仅使用标准库。构建工具及其依赖完整固定在 `requirements-build.txt`；先用 `python3 -m pip install -r requirements-build.txt` 安装。源码和业务测试无需访问外部服务，首次安装构建依赖需要访问 Python 包仓库。

## 当前摘要

在项目根目录运行，会临时启动真实本地 HTTP 服务、核验资料代际、计算摘要并完整打印 JSON：

```sh
python3 -m stormcast.runner --output /tmp/stormcast-summary.json
```

默认 `--mode all` 同时输出**现行经验摘要**与**极端降雨频率分析**。`--mode legacy` 只保留原有经验摘要，`--mode frequency` 只输出频率分析；`--bootstrap-replicates` 控制尾部重现期自助区间的重采样次数（默认 1000）。

## 极端降雨频率分析（面向防汛物资储备讨论）

`stormcast/frequency.py`（编排）与 `stormcast/extreme.py`（统计内核，纯标准库）在经验摘要之上提供：

- **可解释的高分位与超越概率范围**：经验 0.90–0.99 分位、历史事件 Weibull 经验重现期，以及 GPD/POT 尾部模型的 5/10/20/50/100 年重现期水平，每个水平附 90% 自助区间和年超越概率。
- **有效观测时长入账**：只用完整观测日，有效年数按完整日数折算；部分日不补零，其已知超阈只作风险提示。缺测年份和短记录乡镇不会被等长比较。
- **年代合并与敏感性诊断**：不同仪器年代在平稳性假设下合并，同时给分年代超阈速率；对多个经验分位阈值比较尾部形状，形状对阈值敏感时明确提示优先看经验分位与区间。
- **相邻乡镇联合超越风险**：基于同日共同超阈的经验联合概率、年化发生率与条件概率，并与独立假设对照，识别暴雨过程联动；不硬套多元分布。
- **显式边界**：重现期是长期平均间隔而非固定年数的准时事件；超阈样本不足只给经验结果，超过资料支撑能力的重现期（如 19 年记录下的百年值）标记 `unsupported` 并留空，不硬填数值。

方法选择、公式口径与边界的完整说明见 `docs/data-contract.md`。

## 分开运行供数与应用

```sh
python3 -m stormcast.runner --serve --port 8765
```

另一个终端可运行：

```sh
python3 -m stormcast.runner --source-url http://127.0.0.1:8765 --output /tmp/stormcast-summary.json
```

## 测试与检查

```sh
python3 -m unittest discover -s tests -t . -v
python3 -m compileall -q stormcast tests
python3 -m build --wheel --no-isolation
```

测试覆盖现有领域摘要、非法引用、真实 HTTP 换代拒绝、表清单内容核验、结果保存，以及 GPD 估计器恢复性、频率分析边界与联合超越风险。展示现有产品可直接运行上述摘要命令，不需要回放历史日志。扩展分析时应同步说明新增能力的实际运行入口，保留现行摘要调用。

## 容器运行

Dockerfile 固定官方 Python 3.11.14 slim-bookworm 镜像摘要，不依赖本机私有镜像；首次构建需下载基础镜像和固定构建工具。构建与运行实际摘要：

```sh
docker build -t local/stormcast .
docker run --rm local/stormcast
```

容器内测试：

```sh
docker run --rm --entrypoint python3 local/stormcast -m unittest discover -s tests -t . -v
```
