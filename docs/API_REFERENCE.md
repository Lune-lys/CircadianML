# 📘 API参考文档

**Circadian Transcriptomics Analysis Pipeline**

本文档提供所有公共API的详细参考。

---

## 目录

- [rhythm_detection](#rhyth

m_detection)
  - [jtk_cycle](#jtk_cycle)
  - [cosinor](#cosinor)
  - [lomb_scargle](#lomb_scargle)
  - [wavelet_analysis](#wavelet_analysis)
- [phase_prediction](#phase_prediction)
  - [ml_models](#ml_models)
  - [deep_learning](#deep_learning)
  - [ensemble](#ensemble)
- [network_analysis](#network_analysis)
  - [coexpression](#coexpression)
  - [phase_coupling](#phase_coupling)
  - [community](#community)
- [visualization](#visualization)
  - [phase_plots](#phase_plots)
  - [network_viz](#network_viz)
  - [dashboard](#dashboard)
- [utils](#utils)
  - [helpers](#helpers)
  - [metrics](#metrics)

---

## rhythm_detection

### jtk_cycle

#### `detect_rhythms()`

检测基因表达的昼夜节律（JTK_CYCLE算法）。

**函数签名**:
```python
def detect_rhythms(
    data: pd.DataFrame,
    period: float = 24.0,
    lag_range: Optional[Tuple[int, int]] = None,
    alpha: float = 0.05
) -> pd.DataFrame
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| data | pd.DataFrame | - | 表达矩阵（基因×时间点） |
| period | float | 24.0 | 预期周期长度（小时） |
| lag_range | Tuple[int, int] | None | 相位搜索范围，None表示测试所有相位 |
| alpha | float | 0.05 | FDR显著性阈值 |

**返回值**:
```python
pd.DataFrame: 包含以下列的结果表
    - GeneID (str): 基因标识符
    - P (float): 原始p值
    - TAU (float): Kendall's tau系数 [-1, 1]
    - PER (float): 周期（小时）
    - LAG (float): 相位/峰时间（小时）[0, period)
    - AMP (float): 振幅（峰-谷差的一半）
    - BH.Q (float): Benjamini-Hochberg校正后的q值
    - ADJ.P (float): 调整后的p值（与BH.Q相同）
```

**异常**:
- `ValueError`: 如果data不是DataFrame或为空
- `TypeError`: 如果参数类型不正确

**示例**:
```python
import pandas as pd
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms

# 加载数据
data = pd.read_csv('expression_data.csv', index_col=0)

# 运行JTK_CYCLE
results = detect_rhythms(
    data,
    period=24.0,
    lag_range=(0, 23),
    alpha=0.05
)

# 筛选显著节律基因
rhythmic_genes = results[results['BH.Q'] < 0.05]
print(f"检测到 {len(rhythmic_genes)} 个节律基因")

# 查看top10基因
top10 = rhythmic_genes.head(10)
print(top10[['GeneID', 'BH.Q', 'LAG', 'AMP']])
```

**性能**:
- 时间复杂度: O(n × m × p)，其中n=基因数，m=时间点数，p=相位数
- 内存占用: O(n × m)
- 基准: 1000基因×24时间点 ≈ 30秒（单核）

**注意事项**:
1. 输入数据应该已经过标准化处理
2. 列名应该表示时间点（如ZT00, ZT02...）
3. 缺失值会被自动跳过，但可能影响结果
4. 建议数据至少覆盖1.5个周期（36小时）

---

#### `jtk_test()`

对单个基因运行JTK检验。

**函数签名**:
```python
def jtk_test(
    expression: np.ndarray,
    period: float,
    lag_range: Tuple[int, int]
) -> Tuple[float, float, float]
```

**参数**:
| 参数 | 类型 | 说明 |
|------|------|------|
| expression | np.ndarray | 单个基因的表达值序列 |
| period | float | 周期长度 |
| lag_range | Tuple[int, int] | 相位搜索范围（起始，结束） |

**返回值**:
```python
Tuple[float, float, float]:
    - p-value: 显著性p值
    - best_lag: 最佳相位（小时）
    - best_tau: 最大Kendall's tau值
```

**示例**:
```python
import numpy as np
from circadian_analysis.rhythm_detection.jtk_cycle import jtk_test

# 单个基因的表达数据
expression = np.array([5.2, 6.1, 7.8, 8.9, 8.5, 7.2, 5.5, 4.8, ...])

# 运行JTK检验
pval, lag, tau = jtk_test(expression, period=24.0, lag_range=(0, 23))

print(f"p-value: {pval:.4f}")
print(f"相位: {lag:.1f} 小时")
print(f"Tau: {tau:.3f}")
```

---

#### `filter_rhythmic_genes()`

筛选显著节律基因。

**函数签名**:
```python
def filter_rhythmic_genes(
    results: pd.DataFrame,
    qvalue_cutoff: float = 0.05
) -> pd.DataFrame
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| results | pd.DataFrame | - | JTK_CYCLE结果 |
| qvalue_cutoff | float | 0.05 | Q值阈值 |

**返回值**:
```python
pd.DataFrame: 筛选后的结果（副本）
```

**示例**:
```python
from circadian_analysis.rhythm_detection.jtk_cycle import (
    detect_rhythms, filter_rhythmic_genes
)

results = detect_rhythms(data)
rhythmic = filter_rhythmic_genes(results, qvalue_cutoff=0.01)
```

---

### cosinor

#### `cosinor_regression()`

对单个时间序列进行Cosinor回归拟合。

**函数签名**:
```python
def cosinor_regression(
    expression: np.ndarray,
    time: np.ndarray,
    period: float = 24.0
) -> Dict[str, float]
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| expression | np.ndarray | - | 表达值序列 |
| time | np.ndarray | - | 对应的时间点 |
| period | float | 24.0 | 周期长度（小时） |

**返回值**:
```python
Dict[str, float]: 包含以下键的字典
    - mesor: MESOR（平均表达水平）
    - amplitude: 振幅
    - acrophase: 峰相位（小时）
    - beta: cos项系数
    - gamma: sin项系数
    - rsquared: R²（拟合优度）
    - pvalue: F检验p值
    - fitted: np.ndarray, 拟合值
```

**数学模型**:
```
y(t) = M + A·cos(ωt - φ) + ε
     = M + β·cos(ωt) + γ·sin(ωt) + ε

其中:
    M = MESOR
    A = √(β² + γ²) = 振幅
    φ = arctan2(γ, β) = 峰相位
    ω = 2π/T = 角频率
```

**示例**:
```python
import numpy as np
from circadian_analysis.rhythm_detection.cosinor import cosinor_regression

# 时间和表达数据
time = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22])
expression = np.array([5.2, 6.5, 8.1, 9.2, 8.8, 7.5, 5.8, 4.5, 4.2, 4.8, 5.0, 5.1])

# Cosinor拟合
result = cosinor_regression(expression, time, period=24.0)

print(f"MESOR: {result['mesor']:.2f}")
print(f"振幅: {result['amplitude']:.2f}")
print(f"峰相位: {result['acrophase']:.2f} 小时")
print(f"R²: {result['rsquared']:.3f}")
print(f"p-value: {result['pvalue']:.4f}")

# 可视化拟合
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(time, expression, 'o', label='观测值', markersize=8)
plt.plot(time, result['fitted'], '-', label='拟合曲线', linewidth=2)
plt.xlabel('时间（小时）')
plt.ylabel('表达量')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

---

#### `detect_rhythms_cosinor()`

使用Cosinor回归检测多个基因的节律。

**函数签名**:
```python
def detect_rhythms_cosinor(
    data: pd.DataFrame,
    time: Optional[np.ndarray] = None,
    period: float = 24.0,
    alpha: float = 0.05
) -> pd.DataFrame
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| data | pd.DataFrame | - | 表达矩阵（基因×时间点） |
| time | np.ndarray | None | 时间点数组，None时从列名提取 |
| period | float | 24.0 | 周期长度 |
| alpha | float | 0.05 | 显著性阈值 |

**返回值**:
```python
pd.DataFrame: 包含以下列
    - GeneID: 基因ID
    - MESOR: 平均水平
    - Amplitude: 振幅
    - Acrophase: 峰相位（小时）
    - Rsquared: R²
    - Pvalue: p值
    - Qvalue: FDR校正后的q值
```

**示例**:
```python
from circadian_analysis.rhythm_detection.cosinor import detect_rhythms_cosinor

results = detect_rhythms_cosinor(data, period=24.0)
significant = results[results['Qvalue'] < 0.05]

# 相位分布
import matplotlib.pyplot as plt
plt.hist(significant['Acrophase'], bins=24, edgecolor='black')
plt.xlabel('峰相位（小时）')
plt.ylabel('基因数')
plt.title(f'节律基因相位分布 (n={len(significant)})')
plt.show()
```

---

#### `multicomponent_cosinor()`

多谐波Cosinor拟合（同时拟合多个周期成分）。

**函数签名**:
```python
def multicomponent_cosinor(
    expression: np.ndarray,
    time: np.ndarray,
    periods: list = [24.0, 12.0]
) -> Dict[str, float]
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| expression | np.ndarray | - | 表达值 |
| time | np.ndarray | - | 时间点 |
| periods | list | [24.0, 12.0] | 要拟合的周期列表 |

**返回值**:
```python
Dict[str, float]:
    - mesor: 平均水平
    - amplitude_24: 24h成分的振幅
    - acrophase_24: 24h成分的峰相位
    - amplitude_12: 12h成分的振幅
    - acrophase_12: 12h成分的峰相位
    - rsquared: 整体R²
```

**应用场景**:
- 检测次谐波（如12小时节律）
- 复杂节律模式分析
- 提高拟合精度

**示例**:
```python
from circadian_analysis.rhythm_detection.cosinor import multicomponent_cosinor

# 拟合24h和12h成分
result = multicomponent_cosinor(
    expression,
    time,
    periods=[24.0, 12.0, 8.0]  # 可以添加更多周期
)

print(f"24h振幅: {result['amplitude_24.0']:.2f}")
print(f"12h振幅: {result['amplitude_12.0']:.2f}")
print(f"R²: {result['rsquared']:.3f}")
```

---

### lomb_scargle

#### `detect_rhythms_ls()`

使用Lomb-Scargle周期图检测节律。

**函数签名**:
```python
def detect_rhythms_ls(
    data: pd.DataFrame,
    time: Optional[np.ndarray] = None,
    period_range: Tuple[float, float] = (18.0, 30.0),
    n_periods: int = 100
) -> pd.DataFrame
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| data | pd.DataFrame | - | 表达矩阵 |
| time | np.ndarray | None | 时间点 |
| period_range | Tuple | (18.0, 30.0) | 搜索的周期范围（小时） |
| n_periods | int | 100 | 测试的周期数 |

**返回值**:
```python
pd.DataFrame:
    - GeneID: 基因ID
    - Period: 主导周期（小时）
    - Power: Lomb-Scargle功率
    - Amplitude: 表达标准差（作为振幅估计）
```

**优势**:
- 适用于非等间隔采样数据
- 不需要预先指定周期
- 可以检测多个周期成分

**示例**:
```python
from circadian_analysis.rhythm_detection.lomb_scargle import detect_rhythms_ls

# 检测18-30小时范围内的周期
results = detect_rhythms_ls(
    data,
    period_range=(18.0, 30.0),
    n_periods=200
)

# 统计主导周期分布
print(results['Period'].describe())

# 高功率基因
high_power = results[results['Power'] > results['Power'].quantile(0.9)]
print(f"高功率基因数: {len(high_power)}")
```

---

### wavelet_analysis

#### `detect_rhythms_wavelet()`

使用小波分析检测节律。

**函数签名**:
```python
def detect_rhythms_wavelet(
    data: pd.DataFrame,
    time: Optional[np.ndarray] = None,
    period_range: Tuple[float, float] = (18.0, 30.0)
) -> pd.DataFrame
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| data | pd.DataFrame | - | 表达矩阵 |
| time | np.ndarray | None | 时间点 |
| period_range | Tuple | (18.0, 30.0) | 周期搜索范围 |

**返回值**:
```python
pd.DataFrame:
    - GeneID: 基因ID
    - Period: 主导周期
    - Power: 小波功率
```

**特点**:
- 时频联合分析
- 可检测周期变化
- 适用于非平稳信号

**示例**:
```python
from circadian_analysis.rhythm_detection.wavelet_analysis import detect_rhythms_wavelet

results = detect_rhythms_wavelet(data, period_range=(20.0, 28.0))
```

---

## phase_prediction

### ml_models

#### `train_phase_predictor()`

训练相位预测模型。

**函数签名**:
```python
def train_phase_predictor(
    X: np.ndarray,
    y: np.ndarray,
    model_type: str = 'random_forest',
    **model_params
) -> Tuple[object, object]
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| X | np.ndarray | - | 特征矩阵 (n_samples × n_features) |
| y | np.ndarray | - | 相位标签（小时）|
| model_type | str | 'random_forest' | 模型类型 |
| **model_params | dict | - | 模型超参数 |

**支持的模型类型**:
- `'random_forest'`: 随机森林回归器
- `'gradient_boosting'`: 梯度提升回归器

**返回值**:
```python
Tuple[object, object]:
    - model: 训练好的模型对象
    - scaler: StandardScaler对象（用于特征标准化）
```

**示例**:
```python
from circadian_analysis.phase_prediction.ml_models import (
    train_phase_predictor, predict_phase, extract_gene_features
)

# 准备数据
X, feature_names = extract_gene_features(expression_data, rhythmic_results)
y = rhythmic_results['LAG'].values

# 训练模型
model, scaler = train_phase_predictor(
    X, y,
    model_type='random_forest',
    n_estimators=200,
    max_depth=10,
    random_state=42
)

# 预测新数据
predictions = predict_phase(model, scaler, X_new)
```

---

#### `predict_phase()`

使用训练好的模型预测相位。

**函数签名**:
```python
def predict_phase(
    model: object,
    scaler: object,
    X: np.ndarray
) -> np.ndarray
```

**参数**:
| 参数 | 类型 | 说明 |
|------|------|------|
| model | object | 训练好的模型 |
| scaler | object | StandardScaler对象 |
| X | np.ndarray | 特征矩阵 |

**返回值**:
```python
np.ndarray: 预测的相位（小时），范围[0, 24)
```

---

#### `extract_gene_features()`

从基因表达数据提取特征。

**函数签名**:
```python
def extract_gene_features(
    data: pd.DataFrame,
    rhythmic_results: Optional[pd.DataFrame] = None
) -> Tuple[np.ndarray, list]
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| data | pd.DataFrame | - | 表达矩阵 |
| rhythmic_results | pd.DataFrame | None | 节律检测结果（可选） |

**返回值**:
```python
Tuple[np.ndarray, list]:
    - X: 特征矩阵
    - feature_names: 特征名称列表
```

**提取的特征**:
1. 原始表达值（每个时间点）
2. 统计特征：mean, std, max, min
3. 节律特征（如果提供）：amplitude, tau

---

#### `evaluate_phase_predictor()`

使用交叉验证评估相位预测模型。

**函数签名**:
```python
def evaluate_phase_predictor(
    X: np.ndarray,
    y: np.ndarray,
    model_type: str = 'random_forest',
    cv: int = 5,
    **model_params
) -> Dict[str, float]
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| X | np.ndarray | - | 特征矩阵 |
| y | np.ndarray | - | 真实相位 |
| model_type | str | 'random_forest' | 模型类型 |
| cv | int | 5 | 交叉验证折数 |
| **model_params | dict | - | 模型参数 |

**返回值**:
```python
Dict[str, float]:
    - mean_mae: 平均绝对误差的均值
    - std_mae: MAE的标准差
    - mean_phase_error: 平均相位误差（考虑周期性）
    - std_phase_error: 相位误差标准差
```

**示例**:
```python
from circadian_analysis.phase_prediction.ml_models import evaluate_phase_predictor

# 评估模型性能
metrics = evaluate_phase_predictor(
    X, y,
    model_type='random_forest',
    cv=10,
    n_estimators=100
)

print(f"MAE: {metrics['mean_mae']:.2f} ± {metrics['std_mae']:.2f} 小时")
print(f"相位误差: {metrics['mean_phase_error']:.2f} 小时")
```

---

### ensemble

#### `PhaseEnsemble`

相位预测的集成模型类。

**类定义**:
```python
class PhaseEnsemble:
    def __init__(self, model_types: List[str] = None)
    def fit(self, X: np.ndarray, y: np.ndarray, **model_params)
    def predict(self, X: np.ndarray, method: str = 'weighted_mean') -> np.ndarray
    def optimize_weights(self, X_val: np.ndarray, y_val: np.ndarray)
```

**方法详解**:

##### `__init__()`
```python
def __init__(self, model_types: List[str] = None)
```

**参数**:
- `model_types`: 要包含的模型类型列表，默认`['random_forest', 'gradient_boosting']`

##### `fit()`
训练集成中的所有模型。

```python
def fit(self, X: np.ndarray, y: np.ndarray, **model_params)
```

##### `predict()`
进行集成预测。

```python
def predict(self, X: np.ndarray, method: str = 'weighted_mean') -> np.ndarray
```

**参数**:
- `method`: 集成方法
  - `'weighted_mean'`: 加权平均
  - `'median'`: 中位数
  - `'circular_mean'`: 循环均值（考虑相位周期性）

##### `optimize_weights()`
在验证集上优化集成权重。

```python
def optimize_weights(self, X_val: np.ndarray, y_val: np.ndarray)
```

**完整示例**:
```python
from circadian_analysis.phase_prediction.ensemble import PhaseEnsemble

# 创建集成模型
ensemble = PhaseEnsemble(model_types=['random_forest', 'gradient_boosting'])

# 训练
ensemble.fit(X_train, y_train, n_estimators=100)

# 在验证集上优化权重
ensemble.optimize_weights(X_val, y_val)

# 预测
predictions = ensemble.predict(X_test, method='circular_mean')
```

---

## network_analysis

### coexpression

#### `build_coexpression_network()`

构建基因共表达网络。

**函数签名**:
```python
def build_coexpression_network(
    data: pd.DataFrame,
    threshold: float = 0.7,
    method: str = 'pearson'
) -> Tuple[pd.DataFrame, pd.DataFrame]
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| data | pd.DataFrame | - | 表达矩阵 |
| threshold | float | 0.7 | 相关系数阈值 |
| method | str | 'pearson' | 相关方法 |

**相关方法**:
- `'pearson'`: Pearson相关系数
- `'spearman'`: Spearman秩相关
- `'kendall'`: Kendall tau

**返回值**:
```python
Tuple[pd.DataFrame, pd.DataFrame]:
    - adjacency_matrix: 邻接矩阵（基因×基因）
    - edge_list: 边列表，包含以下列：
        - source: 源基因
        - target: 目标基因
        - weight: 相关系数
```

**示例**:
```python
from circadian_analysis.network_analysis.coexpression import (
    build_coexpression_network, compute_network_metrics
)

# 构建网络
adj_matrix, edge_list = build_coexpression_network(
    data,
    threshold=0.8,
    method='pearson'
)

print(f"网络规模: {len(adj_matrix)} 个节点, {len(edge_list)} 条边")

# 计算网络指标
metrics = compute_network_metrics(adj_matrix)
print(f"平均度: {metrics['Degree'].mean():.2f}")
```

---

#### `compute_network_metrics()`

计算网络拓扑指标。

**函数签名**:
```python
def compute_network_metrics(adj_matrix: pd.DataFrame) -> pd.DataFrame
```

**参数**:
- `adj_matrix`: 邻接矩阵

**返回值**:
```python
pd.DataFrame:
    - Gene: 基因名
    - Degree: 度（连接数）
    - Clustering: 聚类系数
```

---

### phase_coupling

#### `build_phase_coupling_network()`

构建相位耦合网络。

**函数签名**:
```python
def build_phase_coupling_network(
    rhythmic_results: pd.DataFrame,
    coupling_threshold: float = 3.0,
    period: float = 24.0
) -> pd.DataFrame
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| rhythmic_results | pd.DataFrame | - | 节律检测结果 |
| coupling_threshold | float | 3.0 | 最大相位差（小时） |
| period | float | 24.0 | 周期 |

**返回值**:
```python
pd.DataFrame: 边列表
    - source: 源基因
    - target: 目标基因
    - phase_diff: 相位差（小时）
    - phase_source: 源基因相位
    - phase_target: 目标基因相位
```

**示例**:
```python
from circadian_analysis.network_analysis.phase_coupling import (
    build_phase_coupling_network
)

# 构建相位耦合网络（相位差<2小时的基因对）
coupling_edges = build_phase_coupling_network(
    rhythmic_results,
    coupling_threshold=2.0
)

print(f"相位耦合的基因对: {len(coupling_edges)}")
```

---

### community

#### `greedy_modularity_communities()`

使用贪婪模块度优化检测社区。

**函数签名**:
```python
def greedy_modularity_communities(
    adj_matrix: pd.DataFrame,
    max_iterations: int = 100
) -> Dict[str, int]
```

**参数**:
- `adj_matrix`: 邻接矩阵
- `max_iterations`: 最大迭代次数

**返回值**:
```python
Dict[str, int]: 基因到社区ID的映射
```

**示例**:
```python
from circadian_analysis.network_analysis.community import (
    greedy_modularity_communities, annotate_communities
)

# 检测社区
communities = greedy_modularity_communities(adj_matrix)

# 创建注释表
comm_df = annotate_communities(communities)
print(comm_df.groupby('Community').size())
```

---

## visualization

### phase_plots

#### `plot_rhythmic_genes()`

绘制节律基因的表达曲线。

**函数签名**:
```python
def plot_rhythmic_genes(
    data: pd.DataFrame,
    rhythmic_results: pd.DataFrame,
    top_n: int = 6,
    save_path: Optional[str] = None
)
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| data | pd.DataFrame | - | 表达矩阵 |
| rhythmic_results | pd.DataFrame | - | 节律检测结果 |
| top_n | int | 6 | 绘制前N个基因 |
| save_path | str | None | 保存路径 |

**示例**:
```python
from circadian_analysis.visualization.phase_plots import plot_rhythmic_genes

plot_rhythmic_genes(
    data,
    rhythmic_results,
    top_n=9,
    save_path='top9_genes.png'
)
```

---

#### `plot_phase_distribution()`

绘制相位分布图。

**函数签名**:
```python
def plot_phase_distribution(
    rhythmic_results: pd.DataFrame,
    phase_col: str = 'LAG',
    bins: int = 24,
    save_path: Optional[str] = None
)
```

**创建**:
- 线性直方图
- 环形（polar）直方图

---

#### `plot_phase_wheel()`

创建相位轮盘图。

**函数签名**:
```python
def plot_phase_wheel(
    rhythmic_results: pd.DataFrame,
    phase_col: str = 'LAG',
    amplitude_col: Optional[str] = 'AMP',
    top_n: int = 50,
    save_path: Optional[str] = None
)
```

---

### dashboard

#### `create_analysis_dashboard()`

创建综合分析仪表板。

**函数签名**:
```python
def create_analysis_dashboard(
    data: pd.DataFrame,
    rhythmic_results: pd.DataFrame,
    save_path: Optional[str] = None
)
```

**包含**:
- Top 6节律基因表达曲线
- 相位分布（线性+环形）
- 振幅分布
- Q值分布
- 周期分布
- 统计摘要

---

## utils

### helpers

#### `generate_synthetic_data()`

生成合成的昼夜节律表达数据。

**函数签名**:
```python
def generate_synthetic_data(
    n_genes: int = 100,
    n_timepoints: int = 24,
    n_rhythmic: int = 30,
    period: float = 24.0,
    noise_level: float = 0.3,
    random_state: Optional[int] = 42
) -> pd.DataFrame
```

**参数**:
| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| n_genes | int | 100 | 总基因数 |
| n_timepoints | int | 24 | 时间点数 |
| n_rhythmic | int | 30 | 节律基因数 |
| period | float | 24.0 | 周期（小时） |
| noise_level | float | 0.3 | 噪声水平（高斯噪声的标准差） |
| random_state | int | 42 | 随机种子 |

**返回值**:
```python
pd.DataFrame: 表达矩阵（基因×时间点）
    - 索引: Gene_001, Gene_002, ...
    - 列: ZT00, ZT01, ..., ZT23
```

**生成机制**:
- 节律基因: $y(t) = baseline + amplitude \cdot \cos(2\pi t/T - \phi) + \epsilon$
- 非节律基因: $y(t) = baseline + \epsilon$

**示例**:
```python
from circadian_analysis.utils.helpers import generate_synthetic_data

# 生成测试数据
data = generate_synthetic_data(
    n_genes=200,
    n_timepoints=24,
    n_rhythmic=60,
    period=24.0,
    noise_level=0.2,
    random_state=123
)

print(data.head())
```

---

#### `normalize_expression()`

标准化表达数据。

**函数签名**:
```python
def normalize_expression(
    data: pd.DataFrame,
    method: str = "zscore"
) -> pd.DataFrame
```

**参数**:
- `method`: 标准化方法
  - `'zscore'`: Z-score标准化（每个基因）
  - `'minmax'`: Min-max标准化到[0, 1]
  - `'log2'`: Log2转换

**返回值**:
- 标准化后的DataFrame

---

#### `adjust_pvalues()`

P值多重检验校正。

**函数签名**:
```python
def adjust_pvalues(
    pvalues: np.ndarray,
    method: str = "fdr_bh"
) -> np.ndarray
```

**参数**:
- `pvalues`: P值数组
- `method`: 校正方法
  - `'fdr_bh'`: Benjamini-Hochberg FDR
  - `'bonferroni'`: Bonferroni校正

**返回值**:
- 校正后的P值

---

### metrics

#### `evaluate_rhythm_detection()`

评估节律检测性能。

**函数签名**:
```python
def evaluate_rhythm_detection(
    true_labels: np.ndarray,
    pred_labels: np.ndarray,
    pred_scores: Optional[np.ndarray] = None
) -> Dict[str, float]
```

**参数**:
- `true_labels`: 真实标签（1=节律，0=非节律）
- `pred_labels`: 预测标签
- `pred_scores`: 预测分数（可选，用于计算AUC）

**返回值**:
```python
Dict[str, float]:
    - accuracy: 准确率
    - precision: 精确率
    - recall: 召回率
    - f1_score: F1分数
    - specificity: 特异性
    - auc_roc: ROC曲线下面积（如果提供scores）
    - tp, tn, fp, fn: 混淆矩阵各值
```

---

#### `compute_phase_error()`

计算考虑周期性的相位误差。

**函数签名**:
```python
def compute_phase_error(
    true_phase: np.ndarray,
    pred_phase: np.ndarray,
    period: float = 24.0
) -> float
```

**参数**:
- `true_phase`: 真实相位
- `pred_phase`: 预测相位
- `period`: 周期

**返回值**:
- 平均绝对相位误差（小时）

**特点**:
- 考虑相位的循环性质
- 例如：23h和1h的误差是2h（不是22h）

---

## 完整工作流示例

### 示例1: 基础节律检测

```python
import pandas as pd
from circadian_analysis.utils.helpers import generate_synthetic_data
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.visualization.phase_plots import (
    plot_rhythmic_genes, plot_phase_distribution
)

# 1. 生成或加载数据
data = generate_synthetic_data(n_genes=200, n_rhythmic=50)

# 2. 检测节律
results = detect_rhythms(data, period=24.0)
rhythmic = results[results['BH.Q'] < 0.05]

print(f"检测到 {len(rhythmic)} 个节律基因")

# 3. 可视化
plot_rhythmic_genes(data, rhythmic, save_path='rhythmic_genes.png')
plot_phase_distribution(rhythmic, save_path='phase_dist.png')

# 4. 保存结果
results.to_csv('jtk_results.csv', index=False)
```

### 示例2: 相位预测

```python
from circadian_analysis.rhythm_detection.jtk_cycle import detect_rhythms
from circadian_analysis.phase_prediction.ml_models import (
    extract_gene_features, train_phase_predictor, predict_phase
)
from circadian_analysis.utils.metrics import compute_phase_error
from sklearn.model_selection import train_test_split

# 1. 检测节律并获取相位
results = detect_rhythms(data)
rhythmic = results[results['BH.Q'] < 0.05]

# 2. 提取特征
rhythmic_data = data.loc[rhythmic['GeneID'].values]
X, feature_names = extract_gene_features(rhythmic_data, rhythmic)
y = rhythmic['LAG'].values

# 3. 划分数据集
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# 4. 训练模型
model, scaler = train_phase_predictor(
    X_train, y_train,
    model_type='random_forest',
    n_estimators=200
)

# 5. 预测和评估
y_pred = predict_phase(model, scaler, X_test)
error = compute_phase_error(y_test, y_pred)

print(f"平均相位误差: {error:.2f} 小时")
```

### 示例3: 网络分析

```python
from circadian_analysis.network_analysis.coexpression import (
    build_coexpression_network, compute_network_metrics, identify_hub_genes
)
from circadian_analysis.network_analysis.community import (
    greedy_modularity_communities, annotate_communities
)

# 1. 构建网络
adj_matrix, edge_list = build_coexpression_network(
    rhythmic_data,
    threshold=0.8,
    method='pearson'
)

# 2. 网络指标
metrics = compute_network_metrics(adj_matrix)
hubs = identify_hub_genes(adj_matrix, top_n=10)

print("Top 10 Hub Genes:")
print(hubs)

# 3. 社区检测
communities = greedy_modularity_communities(adj_matrix)
comm_df = annotate_communities(communities)

print(f"\n检测到 {comm_df['Community'].nunique()} 个社区")
print(comm_df.groupby('Community').size())
```

---

## 版本信息

**当前版本**: 0.1.0
**发布日期**: 2025-11-17
**Python要求**: ≥3.8

---

## 变更日志

### v0.1.0 (2025-11-17)
- 初始发布
- 实现所有核心API
- 完整的文档覆盖

---

**文档编译日期**: 2025-11-17
**维护者**: Zhihan Luo (github.com/Lune-lys)
