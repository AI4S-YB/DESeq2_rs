# 仅三千万 Token，28x 速度提升，这一次 Rust 成功重构 DESeq2 我做对了什么？

上一篇文章里我写了自己烧掉 2 亿 token 重写 DESeq2 失败的故事。结论很简单：不要迷信 Superpower，你自己没有能力，没有工具能给你超能力。

文章发出来之后我一直在想一个问题：如果重来一次，到底该怎么做？

所以我真的重来了一次。

## 这次的结果

先说数据。

| | R DESeq2 | Rust (本项目) |
|---|---|---|
| 运行时间（wall clock） | 7.9 秒 | 0.28 秒 |
| 速度倍数 | 1x | **28x** |
| 显著基因（padj<0.01, \|log2FC\|>1） | 675 | 672 重合 + 39 额外 |
| 基因集重合度 | - | **99.6%** |

三千七百万 token，花费 26.73 美元，单次会话完成从设计到 CLI 全流程。38 个单元测试全部通过。

上一次 2 亿 token 连一个正确的 p-value 都没跑出来。这次 3700 万 token 把整条 pipeline 跑通了，和 R 的结果高度一致。

差在哪里？不是工具变了。用的还是同一套东西。差在我变了。

## 上次做错了什么

上次的问题不是技术路线选错了，也不是工具不好用。问题出在我的做事方式上。

**我把"完美设计"当成了第一优先级。** DESeq2 的实现细节很多——IRLS 收敛条件、Cox-Reid adjustment、经验贝叶斯 prior、Cook's distance、independent filtering。正因为知道它复杂，我更倾向于一开始就把标准拉满：架构要干净，验证要精确到 bit，每个中间步骤都要有 R 的 ground truth 做 oracle。

**我建了一堆基础设施，但没有先把核心跑通。** 我写了 oracle fixture 导出工具、field-level comparison 框架、worktree 管理脚本、多阶段 verify 流程。这些东西单独看都对，但它们在主问题尚未收敛时同时出现，效果很糟糕：我持续获得"我在推进"的幻觉，却迟迟拿不到一个能跑通的 `deseq2_run()` 调用。

**我用复杂度来回避不确定性。** 因为对 IRLS 的数值行为不够确定，对 dispersion fitting 的收敛路径不够确定，所以我本能地去多做几层保证——更完整的 plan、更细的 task 拆分、更多的中间验证。这些动作表面上在降低风险，实际上在推迟面对核心问题的时间。

## 这次做对了什么

### 1. 定义了一个极其具体的最小目标

不是"用 Rust 重写 DESeq2"。而是：

> 在 airway 数据集上，用 `padj < 0.01` 和 `|log2FC| > 1` 筛选基因，Rust 和 R 的结果要一样。

就这一个目标。不是 bit-to-bit 精确。不是功能完整。不是架构优雅。就是：**两个实现在同一个数据集上同意哪些基因是差异表达的。**

这个目标具体到可以用一个 Python 脚本在 10 秒内验证。每次改完代码，跑一遍就知道离目标还有多远。

### 2. 砍掉了所有非核心功能

DESeq2 有很多功能。这次我只实现了最小 pipeline：

- Size factor estimation（median-of-ratios）
- Dispersion estimation（gene-wise MLE → parametric trend → MAP shrinkage）
- Negative binomial GLM fitting（IRLS）
- Wald test
- BH p-value adjustment + independent filtering

砍掉的：LRT、lfcShrink（apeglm/ashr/normal）、formula 解析、tximport、VST、rlog、outlier replacement、Cook's distance。

这些都是重要功能。但在核心 pipeline 跑通之前，它们都是噪音。

### 3. 逐步验证，每一步都和 R 对比

我写了一个 R 脚本，把 DESeq2 在 airway 数据上每一步的中间结果都导出成 TSV：

```
r_size_factors.tsv
r_disp_gene_estimates.tsv
r_disp_trend_coeffs.tsv
r_disp_map_estimates.tsv
r_beta_coefficients.tsv
r_results.tsv
```

然后 Rust 实现每完成一步，就和对应的 R 结果对比。不是写一个精美的比较框架，就是一个 Python one-liner：

```python
# 误差大不大？哪些基因差得多？原因是什么？
python3 -c "import csv; ..."
```

这种方式粗糙但有效。它让我在 10 秒内知道当前步骤的精度，而不是花 2 小时搭一个自动化比较系统。

### 4. 遇到精度问题时，先诊断再修

当 Rust 输出 509 个显著基因而 R 输出 675 个时，我没有开始重写架构。而是逐步排查：

- Size factors 对不对？ → 完美匹配（误差 < 1e-10）
- Gene-wise dispersions 对不对？ → 一半基因误差 >10%
- 为什么？ → **发现 bug：传了 normalized counts 给 NB log-posterior，应该传 raw counts**

修了这个 bug 之后：

- Gene-wise dispersions 对不对？ → 99.7% 基因误差 <1%
- 但还有 4400 个 outlier → **加了 all-zero 基因过滤、moments 估计、收敛检查**

修了之后 overlap 从 0% 涨到 99.6%。每一步修复都有明确的诊断依据，不是猜。

### 5. 知道什么时候该停

最终 overlap 停在 672/675 = 99.6%。剩下的 3 个 R-only 基因都是 borderline case（R padj = 0.008-0.009，Rust padj = 0.010-0.012）。39 个 Rust-only 基因来自 dispersion trend 系数的 ~5% 偏差。

要关掉这 0.4% 的差距，需要在 dispersion estimation 里实现 R 的完整 NB GLM mu 迭代（R 用 `fitNbinomGLMs` 计算 mu，我们用的是 linear model 投影）。这是一个确定的技术任务，但工作量不小，对 MVP 来说不值得。

所以我停了。记录下 root cause，写进 TODO，提交，推送。

上一次我在类似的岔路口选择了"继续做到完美"，然后在完美的路上耗尽了所有预算。这一次我选择了"够好就发"。

## 架构长什么样

最终选了 Functional Core + Pipeline Shell 的架构：

**Functional Core**：每个算法步骤是纯函数，输入明确，输出明确，可以独立测试和对比。这让逐步验证变得非常自然——`estimate_size_factors()` 的输出直接和 R 的 `sizeFactors()` 对比，不需要跑整条 pipeline。

**Pipeline Shell**：`DESeqDataSet` 结构体串联所有纯函数，提供 `run()` 和 `results()` 方法。CLI 包一层 clap 就完了。

依赖全部是纯 Rust：`faer`（线性代数），`rayon`（并行），`statrs`（统计函数）。不需要系统 BLAS/LAPACK，编译即用，跨平台友好。

## 从 2 亿到 3700 万

两次尝试的对比：

| | 第一次 | 第二次 |
|---|---|---|
| Token 消耗 | ~2 亿 | 3700 万 |
| 费用 | 远超 $100 | $26.73 |
| 会话数 | 多次，跨多天 | **单次会话** |
| 产出 | 大量 plan/spec/fixture，零正确结果 | 完整 pipeline + CLI + 38 个测试 |
| 基因集 overlap | 0% | 99.6% |
| 速度 vs R | 无法运行 | **28x 加速** |

token 消耗少了 80%，效果从"完全不能用"变成"几乎完全一致"。

差异不在工具。两次用的都是 Claude + Superpowers。差异在于：**第一次我让工具替我做判断，第二次我自己做判断，让工具替我执行。**

## 给同行的建议

如果你也在用 AI 做科学计算软件的重构或重写：

**先定义最小可验证闭环。** 对于 DESeq2 来说就是"一个数据集、一组参数、比较输出"。对于你的项目，找到等价的东西。所有工作都围绕这个闭环展开。

**先让结果大致对，再追求精确对。** 99% 的 overlap 比 0% 的 bit-to-bit 精确有用得多。前者能让你判断路线对不对、值不值得继续投入。后者只是一堆基础设施的囚徒。

**逐步验证，不要攒到最后。** 每实现一个模块就和参考实现对比。误差大了立刻诊断。不要等整条 pipeline 搭完再发现 size factor 就算错了。

**知道什么时候该停。** 0.4% 的精度差距有明确的 root cause 和修复路径。但修复它的 ROI 不足以支撑在 MVP 阶段投入。记录下来，留给下一次。

**工具放大判断力，不替代判断力。** Superpowers 的 brainstorming、writing-plans、subagent-driven-development 确实好用。但前提是你知道当前最重要的问题是什么、应该先交付什么、应该先砍掉什么。如果你自己不确定，工具只会忠实地帮你把不确定性包装得很漂亮。

## 项目地址

GitHub: [AI4S-YB/DESeq2_rs](https://github.com/AI4S-YB/DESeq2_rs)

后续改进方向见 [docs/TODO.md](https://github.com/AI4S-YB/DESeq2_rs/blob/main/docs/TODO.md)。
