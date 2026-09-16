# dem2phase 当前代码技术说明

版本日期：2026-09-13  
适用目录：仓库根目录

## 1. 当前系统定位

当前代码从真实 DEM/DSM 随机裁取地理 patch，构造低空多无人机、多垂直基线
InSAR 仿真数据。名义配置为 P 波段、500 m 飞行高度、五架 UAV、四条星形
干涉边，并在复 SLC 节点层加入热噪声、散射失相干、同步误差、LOS 轨迹误差、
姿态和亚像素配准误差。

代码现在适合以下任务：

- 多基线相位解缠和高程恢复数据生成；
- 可变有效边数、UAV 掉线和长短基线联合训练；
- 地形/地表感知相干建模；
- `sync-only / trajectory-only / coreg-only / combined` 配对消融；
- nominal、challenging、failure 三级压力测试；
- 组合故障识别和置信度研究。

它仍是可控的训练仿真器，不是原始回波级 SAR 成像器或经过实飞标定的数字孪生。
所有名称含 `unvalidated` 的参数必须通过真实飞行数据重新拟合。

## 2. 名义雷达与编队参数

| 量 | 当前值 | 说明 |
|---|---:|---|
| 中心频率 $f_c$ | 500 MHz | P 波段 |
| 波长 $\lambda$ | 0.599585 m | 由 $c/f_c$ 计算 |
| 带宽 $B_w$ | 200 MHz | 用于距离分辨率和临界基线近似 |
| 距离分辨率 $\delta_r$ | 0.74948 m | $c/(2B_w)$ |
| 飞行高度 $H$ | 500 m | 相对场景参考面 |
| 入射角 $\theta$ | 45° | 名义局部平面入射角 |
| 斜距 $R$ | 707.1068 m | 当前按 $H/\cos\theta$ 设置 |
| 航速 | 20 m/s | 已配置，当前 patch 误差模型尚未显式使用 |
| PRF | 1000 Hz | 已配置，当前行坐标尚未严格转换为慢时间 |
| UAV 数量 $N$ | 5 | UAV 1 为星形图参考节点 |
| 有效垂直基线 | 0.75/1.2/1.8/3.0 m | 允许配置范围 0.5–3.0 m |
| 模糊高程 | 199.86/124.91/83.28/49.97 m | 单站配对、双程相位模型 |

## 3. 总体数据生成链

1. 递归读取 DEM/DSM，并由 `dem_split_manifest.csv` 将每个地理瓦片固定分配到
   train、test 或 validation；同地理瓦片的跨传感器版本不得跨 split，存在逐
   DEM 对齐 WorldCover 时同时读取。
2. 在 2–4 倍 bicubic 插值 DEM 中随机选择尺度与位置，裁取 $256\times256$ patch。
3. 根据各基线计算干净连续相位和缠绕相位。
4. 计算坡度、坡向、局部入射角、粗糙度、曲率、TPI、叠掩/阴影代理。
5. 构造地形、基线、地表覆盖和随机残差共同决定的场景相干图。
6. 在 UAV 节点层生成共享散射、节点接收噪声和逐从机失相干散射。
7. 加入节点同步、LOS 轨迹、姿态和配准误差，再形成干涉边。
8. 仅从受扰节点 SLC 对估计 `coherence_observed`。
9. 随机选择有效边；组合故障模式可进一步使从 UAV/边掉线。
10. 保存一个包含 `[K,H,W]` 栈、图结构、物理量、mask 和标签的 MAT。

## 4. 高程—相位数学模型

### 4.1 基本雷达量

光速取

$$
c=299,792,458\ \mathrm{m/s},\qquad
\lambda=\frac{c}{f_c},\qquad
\delta_r=\frac{c}{2B_w}.
$$

代码使用的临界垂直基线近似为

$$
B_{\mathrm{crit}}
=\frac{\lambda R\tan\theta}{2\delta_r}.
$$

当前配置得到 $B_{\mathrm{crit}}\approx282.84$ m，因此 0.75–3 m 基线本身只产生
轻微的临界基线去相关。长基线明显变差更多来自条纹密度、地形、散射、配准和
Phase-4 误差，而不是接近临界基线。

### 4.2 高程相位比例

对第 $k$ 条边，当前远场近似的高程—相位系数为

$$
\kappa_k
=m\frac{2\pi}{\lambda}
\frac{B_{\perp,k}}{R\sin\theta},
$$

其中 $m$ 是路径贡献次数：单站图像配对取 $m=2$；一发多收的单程基线贡献
近似取 $m=1$。代码的名义 profile 为 `monostatic_pair`，因此 $m=2$。

对去除 patch 最低高程后的 DEM $h(x,y)$，

$$
\phi_{k,\mathrm{unw}}(x,y)=\kappa_k h(x,y),
$$

$$
\phi_{k,\mathrm{wrap}}(x,y)
=\arg\{\exp[j\phi_{k,\mathrm{unw}}(x,y)]\}\in[-\pi,\pi].
$$

模糊高程为

$$
h_{\mathrm{amb},k}=\frac{2\pi}{\kappa_k}
=\frac{\lambda R\sin\theta}{mB_{\perp,k}}.
$$

## 5. 地形感知相干模型

令 DEM 的地面像素间隔为 $d$，代码先计算

$$
g_x=\frac{\partial h}{\partial x},\qquad
g_y=\frac{\partial h}{\partial y},\qquad
s=\arctan\sqrt{g_x^2+g_y^2}.
$$

坡向为 $\alpha=\operatorname{atan2}(g_y,g_x)$。表面单位法向量为

$$
\mathbf n=\frac{[-g_x,-g_y,1]^T}{\sqrt{1+g_x^2+g_y^2}}.
$$

根据入射角和视线方位角构造单位视线 $\mathbf l$，局部入射角满足

$$
\theta_{\mathrm{local}}=\arccos(\mathbf n^T\mathbf l).
$$

7×7 局部粗糙度使用高程标准差：

$$
r=\sqrt{\langle h^2\rangle_7-\langle h\rangle_7^2}.
$$

曲率由 MATLAB `del2` 近似，TPI 为

$$
\mathrm{TPI}=h-\langle h\rangle_{15}.
$$

绝对尺度模式下，粗糙度和曲率归一化为

$$
\bar r=\operatorname{clip}(r/r_0,0,1),\qquad
\bar c=\operatorname{clip}(|c|/c_0,0,1).
$$

确定性地形质量为

$$
q_{\mathrm{terrain}}
=\exp\left[-w_s(s/s_0)^2-w_r\bar r-w_c\bar c\right]
\sqrt{\max(\cos\theta_{\mathrm{local}},0)}.
$$

一阶几何代理把以下区域置为最低相干：

$$
\text{layover}:\quad s_{\mathrm{range}}>\theta,
$$

$$
\text{shadow}:\quad
s_{\mathrm{range}}<-(\pi/2-\theta)
\ \text{或}\ \cos\theta_{\mathrm{local}}\le0.
$$

这只是局部坡面代理，不等价于对完整 DEM 做射线追踪。

### 5.1 基线衰减

名义模型使用

$$
d_k=\max\left(0,1-\frac{|B_{\perp,k}|}{B_{\mathrm{crit}}}\right).
$$

兼容压力模型也支持

$$
d_k=\exp\left[-\alpha
\frac{\max(B_{\perp,k}-B_{\mathrm{ref}},0)}{B_{\mathrm{ref}}}\right].
$$

### 5.2 跨边相关随机残差

令 $z_s(x,y)$ 为所有边共享的平滑标准高斯场，$z_k(x,y)$ 为第 $k$ 条边独立
平滑高斯场，则

$$
q_k=\operatorname{clip}\left[
q_{\mathrm{terrain}}
\exp(\sigma_s z_s+\sigma_e z_k),0,1\right].
$$

第 $k$ 条边的相干上界和地形相干图为

$$
\gamma_{\max,k}=\gamma_{\min}
+(\gamma_{\max}-\gamma_{\min})d_k,
$$

$$
\gamma_{\mathrm{terrain},k}=\gamma_{\min}
+q_k(\gamma_{\max,k}-\gamma_{\min}).
$$

## 6. 地表覆盖调制

WorldCover 类别通过最近邻方式与 DEM 像元中心对齐。不得对类别代码做 bicubic
或线性插值。每类配置一个因子 $f_c$，最终场景相干为

$$
\gamma_{\mathrm{scene},k}
=\max\left(\gamma_{\min},
\gamma_{\mathrm{terrain},k}f_{c(x,y)}\right).
$$

当前因子：树木 0.85、灌木 0.90、草地 0.97、农田 0.90、建设用地 0.85、
裸地 1.00、雪冰 0.75、水体 0.08、湿地 0.65、红树林 0.65、苔藓 0.90。
这些是 `initial_pband_singlepass_v1_unvalidated` 假设，不是 WorldCover 提供的
相干测量值。

## 7. 复 SLC、散射与热噪声

代码当前要求复 SLC 路径使用星形图。设参考 UAV 的共享复散射场

$$
z_0\sim\mathcal{CN}(0,1),
$$

节点接收噪声 $n_i\sim\mathcal{CN}(0,1)$，线性 SNR 为

$$
\rho_i=10^{\mathrm{SNR}_{i,\mathrm{dB}}/10}.
$$

参考节点 SLC 为

$$
S_0=z_0+\frac{n_0}{\sqrt{\rho_0}}.
$$

第 $k$ 条边对应的从节点 $i$ 为

$$
S_i
=\gamma_k z_0\exp(-j\phi_k)
+\sqrt{1-\gamma_k^2}\,z_i
+\frac{n_i}{\sqrt{\rho_i}},
$$

其中 $z_i\sim\mathcal{CN}(0,1)$ 独立，$\gamma_k$ 为场景相干。干涉相位为

$$
\hat\phi_k=\arg\{S_0S_i^*\}.
$$

考虑两个节点热噪声后的解析相干标签为

$$
\gamma_{\mathrm{true},k}
=\frac{\gamma_k}
{\sqrt{(1+1/\rho_0)(1+1/\rho_i)}}.
$$

由于多条边复用 $S_0$，它们自然共享参考节点的散射与接收噪声，不再是逐边
完全独立噪声。

### 7.1 复数多视

对窗口 $\mathcal W$，多视相位为

$$
\phi_{k,\mathrm{ML}}
=\arg\left\{\frac{1}{|\mathcal W|}
\sum_{(u,v)\in\mathcal W}S_0(u,v)S_i^*(u,v)\right\}.
$$

当前使用 7×7 boxcar。

### 7.2 仅由观测数据估计相干度

复 SLC 路径现在使用标准局部估计器：

$$
\hat\gamma_k
=\frac{|\langle S_0S_i^*\rangle_{\mathcal W}|}
{\sqrt{\langle|S_0|^2\rangle_{\mathcal W}
\langle|S_i|^2\rangle_{\mathcal W}}+\epsilon}.
$$

该量保存为 `coherence_observed`，推理时可以获得。旧版利用 noisy/clean
phase 残差反演的 `coherence_estimated` 会接触干净真值；严格训练加载器会拒绝
这种旧样本，防止输入泄漏。

### 7.3 兼容的相位 PDF 噪声

`legacy_phase_pdf` 使用单视 InSAR 相位噪声密度

$$
p(\varphi;\gamma)=
\frac{1-\gamma^2}{2\pi}
\frac{1}{1-\gamma^2\cos^2\varphi}
\left[1+
\frac{\gamma\cos\varphi\arccos(-\gamma\cos\varphi)}
{\sqrt{1-\gamma^2\cos^2\varphi}}
\right].
$$

代码离散 PDF、累积成 CDF，再以逆 CDF 逐像素采样。它仅用于兼容和消融，
不能表达共享 UAV 节点造成的跨边协方差。

## 8. Phase-4 节点误差

令归一化方位坐标 $u\in[-1/2,1/2]$，行索引为 $t$。

### 8.1 同步相位

第 $i$ 个节点的同步误差为

$$
\phi_{\mathrm{sync},i}(t)
=b_i+d_i u(t)+w_i(t)+J_i\mathbf1[t\ge t_{J,i}],
$$

其中 $b_i$ 为常相位偏置，$d_i$ 为整幅线性漂移系数，随机游走满足

$$
w_i(t)=\sum_{\tau\le t}\epsilon_i(\tau)
-\operatorname{mean}_t\sum_{\tau\le t}\epsilon_i(\tau),
\qquad \epsilon_i\sim\mathcal N(0,\sigma_w^2),
$$

$J_i$ 是按配置概率出现的稀疏相位跳变。

### 8.2 LOS 轨迹残差

节点 LOS 距离误差为

$$
\Delta R_i(t)=r_{b,i}+r_{d,i}u(t)
+a_i\sin(2\pi\nu_i\bar t+\psi_i),
$$

对应相位误差

$$
\phi_{\mathrm{traj},i}(t)
=m\frac{2\pi}{\lambda}\Delta R_i(t).
$$

当前 $\bar t\in[0,1]$ 是 patch 内归一化行坐标；尚未使用 PRF、航速和真实
合成孔径时间将其转换为严格物理慢时间。

### 8.3 姿态与亚像素配准

roll/pitch 的一阶地面位移近似为

$$
\Delta x_{\mathrm{roll}}
=\frac{H\tan(\delta\mathrm{roll})}{d_g},\qquad
\Delta y_{\mathrm{pitch}}
=\frac{H\tan(\delta\mathrm{pitch})}{d_g},
$$

其中 $d_g$ 是当前插值尺度下的地面像素间隔。再加入随机常位移与线性漂移：

$$
\Delta x_i(t)=x_{0,i}+x_{d,i}u(t)+\Delta x_{\mathrm{roll},i},
$$

$$
\Delta y_i(t)=y_{0,i}+y_{d,i}u(t)+\Delta y_{\mathrm{pitch},i}.
$$

节点先乘相位误差，再做双线性重采样：

$$
\tilde S_i(x,y)=
\mathcal I\left{S_i
\exp[j(\phi_{\mathrm{sync},i}+\phi_{\mathrm{traj},i})]
\right\}(x-\Delta x_i,y-\Delta y_i).
$$

超出源图范围的像素写零，并由 `coregistration_valid_mask` 排除。边误差由节点差
自然产生，例如同步相位为

$$
\phi_{\mathrm{sync},(i,j)}
=\phi_{\mathrm{sync},i}-\phi_{\mathrm{sync},j}.
$$

## 9. 组合故障与 mask

名义模式先随机保留 2–4 条边，并强制保留最短基线。组合故障模式额外：

- 强制保留最短和最长基线；
- 从非保护的从 UAV 中选择一个掉线；
- 将与掉线 UAV 相连的边从最终 mask 中移除；
- 按活动边统计 $\gamma<\gamma_{\mathrm{low}}$ 的像素比例；
- 根据节点跳变幅度阈值生成同步异常标签。

完整组合故障标签为

$$
y_{\mathrm{compound}}
=y_{\mathrm{dropout}}\land y_{\mathrm{sync}}
\land y_{\mathrm{lowcoh}}\land y_{\mathrm{longest}}.
$$

训练有效像素必须使用

$$
M_{kxy}=M_{\mathrm{edge},k}\land M_{\mathrm{coreg},kxy}.
$$

不能把零填充的失效边解释为零相位观测。

## 10. 地形、地表与故障均衡采样

生成阶段的接受概率为

$$
p_{\mathrm{accept}}
=\min\left(1,
p_{\mathrm{terrain\ class}}m_{\mathrm{landcover}}\right).
$$

当前地形层接受概率为 flat 0.35、rolling 0.70、steep 1.0、
geometry-hazard 1.0。达到配置面积阈值的水体、建设用地、湿地或农田可获得
土地覆盖倍率；只改变真实 crop 的接受概率，不生成虚假类别。

训练 manifest 进一步按

$$
s=(y_{\mathrm{dropout}},y_{\mathrm{sync}},y_{\mathrm{lowcoh}},K_{\mathrm{active}},
\mathrm{terrain\ class})
$$

建立 stratum。若第 $s$ 层样本数为 $n_s$、层数为 $S$、总样本为 $N$，原始
逆频率权重为

$$
w_i^{(0)}=\frac{N}{Sn_{s(i)}}.
$$

代码将权重裁剪到 $[1/r_{\max},r_{\max}]$，再归一化到均值 1。默认
$r_{\max}=10$，防止极少数层获得不稳定的巨大权重。

地理划分字段 `split_group` 使用标准化 1° 地理瓦片（如 `N27E085`），而不是源
DEM 文件名。随机裁块只在每张 DEM 内发生；同一地理瓦片的全部尺度、重叠 crop、
基线、噪声 realization、故障版本以及跨传感器版本继承同一个 split。manifest
加载和训练清单构建都会检查任一地理瓦片是否出现在多个集合，一旦发现立即报错。
逆频率权重只用 train 样本统计；test 和 validation 权重固定为 1。

当前固定划分为：

| split | DEM 数 | patch 数 | 地区 |
|---|---:|---:|---|
| train | 12 | 162 | 喜马拉雅、横断山/四川；含多源稳健性样本 |
| test | 5 | 56 | 西喜马拉雅与东喜马拉雅 |
| validation | 2 | 16 | 阿尔卑斯、安第斯 |

train 和 test 都属于喜马拉雅大区，但使用互不相同的 1° tile；validation 故意采用
少量外域山地，用于检查跨地区泛化，而不是调参时反复拟合的第二训练集。

高程源按用途分级：Copernicus GLO-30 作为主数据（120 patch），ALOS AW3D30
作为地理增广（96 patch），ASTER GDEM V3 仅作为带质量差异的稳健性数据
（18 patch）。ASTER 的 `_num.tif` 被读取为质量元数据；`num_available`、
`num_void_fraction`、`num_mean` 和 `num_min` 随 patch 保存。原目录中与当前
Copernicus tile 重叠的 ALOS N27E088、N28E086，以及完全重复的 N27E085 副本，
没有重复引入，从源头避免跨传感器泄漏和重复计数。

## 11. 随机数和可复现性

批次主种子控制裁块尺度、位置、SNR 和接受判定。组件子种子使用

$$
s_j=\operatorname{mod}(s_{\mathrm{batch}}+104729j,2147483646)+1.
$$

当前每个已保存 patch 记录五类子种子：

1. terrain/coherence residual；
2. complex SLC；
3. Phase-4 node errors；
4. random edge mask；
5. explicit failure selection。

所有带 seed 的组件在返回时恢复调用方 RNG 状态。因此相同组件 seed 可独立重放，
但不会把同批次主随机流不断重置为相同状态。批次初始/最终 RNG 状态、实际 seed
和组件计数器保存于 `rng_seed_*.mat`。

## 12. 完整可配置参数列表

主配置文件为 `configs/uav_p_500m_monostatic.json`。

### 12.1 dataset.generation

| 参数 | 当前值 | 单位/作用 |
|---|---:|---|
| `patches_per_dem` | 32 | 无 manifest 覆盖时的默认 patch 数 |
| `patch_size` | 256 | 像素 |
| `n_precomp_scales` | 3 | 预计算缩放层数 |
| `interp_scale_range` | [2,4] | DEM 放大倍数范围 |
| `rng_seed` | 42 | 批次主种子 |
| `max_wrap_count` | 40 | 任一边最大缠绕周数过滤阈值 |
| `void_fraction_threshold` | 0.05 | NUM 空洞比例阈值 |

### 12.2 dataset.storage 与 sampling

| 参数 | 当前值 | 作用 |
|---|---:|---|
| `storage.mode` | grouped_mat | 一个 patch 一个 MAT |
| `write_legacy_per_edge` | false | 是否重复保存逐边旧格式 |
| `numeric_type` | single | 大型数值栈类型 |
| `terrain_sampling.mode` | weighted | 地形分层接受 |
| `max_attempts_factor` | 12 | 最大候选数/目标数 |
| flat/rolling/steep/hazard | .35/.70/1/1 | 地形接受概率 |
| `dataset.dem_directory` | data/dem | 递归 DEM 根目录 |
| `dataset.split_manifest` | configs/dem_split_manifest.csv | 整 DEM 划分清单 |
| `require_split_assignment` | true | 禁止未分配 DEM 进入生成流程 |

`dem_split_manifest.csv` 可逐 DEM 配置 `geographic_tile`、`split`、`region`、
`quality_role`、`patches_per_dem`、经纬度和数据来源。当前每 DEM 数量根据用途为
3、6、8、12、16 或 24，不再假定每个 split 内完全相同。

### 12.3 landcover

| 参数 | 当前值 |
|---|---|
| `enabled` | true |
| `aligned_file` | 对齐 WorldCover MAT |
| `aligned_directory` | data/landcover |
| `aligned_suffix` | `_worldcover2021_aligned.mat` |
| `missing_policy` | error；缺图即停止生成 |
| `unknown_factor` | 1.0 |
| `class_codes` | 10,20,30,40,50,60,70,80,90,95,100 |
| `coherence_factors` | 0.85,0.90,0.97,0.90,0.85,1,0.75,0.08,0.65,0.65,0.90 |
| `sampling.priority_codes` | 80,50,90,40 |
| `sampling.minimum_fraction` | .01,.05,.01,.30 |
| `sampling.acceptance_multiplier` | 3,1.7,3,1.2 |

### 12.4 radar、geometry、interferometry

| 参数 | 当前值 | 单位 |
|---|---:|---|
| `center_frequency_hz` | 5e8 | Hz |
| `bandwidth_hz` | 2e8 | Hz |
| `prf_hz` | 1000 | Hz |
| `speed_mps` | 20 | m/s |
| `altitude_m` | 500 | m |
| `incidence_angle_deg` | 45 | degree |
| `look_azimuth_deg` | 90 | degree |
| `slant_range_m` | 707.1068 | m |
| `measurement_mode` | monostatic_pair | 模式 |
| `phase_path_multiplicity` | 2 | 路径贡献次数 |
| `num_uavs` | 5 | 节点数 |
| `edge_mode` | star | 图结构 |
| `reference_uav` | 1 | 参考节点 |
| `baseline_perp_m` | [.75,1.2,1.8,3] | m |
| `min_active_edges`/`max_active_edges` | 2/4 | 随机 K 范围 |

### 12.5 coherence

| 参数 | 当前值 | 作用 |
|---|---:|---|
| `min`/`max` | .01/.95 | 相干范围 |
| `baseline_decay_model` | critical_baseline | 基线衰减 |
| `terrain_normalization_model` | absolute_scales | 跨 patch 统一尺度 |
| `decay_alpha` | .35 | 指数兼容模型参数 |
| `slope_weight` | .9 | 坡度惩罚 |
| `roughness_weight` | .8 | 粗糙度惩罚 |
| `curvature_weight` | .35 | 曲率惩罚 |
| `slope_scale_deg` | 45 | degree |
| `roughness_scale_m` | 71.823859 | m |
| `curvature_scale_per_m` | .006326427 | 1/m |
| `shared_residual_std` | .12 | 跨边共享 log-residual 标准差 |
| `edge_residual_std` | .08 | 逐边 log-residual 标准差 |
| `residual_scale_px` | 24 | 平滑相关尺度 |

### 12.6 noise

| 参数 | 当前值 | 作用 |
|---|---:|---|
| `model` | complex_slc_nodes | 节点复 SLC 模型 |
| `node_snr_db_range` | [15,30] | 每节点均匀抽样 dB |
| `multilook_window` | 7 | 奇数 boxcar |
| `save_node_slc` | false | 是否保存 `[N,H,W]` 复节点数据 |

### 12.7 phase4_errors

| 参数 | 当前值 | 单位 |
|---|---:|---|
| `phase_bias_std_rad` | .02 | rad |
| `linear_drift_std_rad` | .03 | rad/patch |
| `random_walk_std_rad_per_row` | .002 | rad/√row increment |
| `jump_probability_per_node` | .1 | probability |
| `jump_std_rad` | .15 | rad |
| `los_bias_std_m` | .01 | m |
| `los_drift_std_m` | .005 | m/patch |
| `vibration_std_m` | .002 | m |
| `vibration_cycles_range` | [1,5] | cycles/patch |
| `roll_std_deg`/`pitch_std_deg` | .02/.02 | degree |
| `range_shift_std_px` | .05 | pixel |
| `azimuth_shift_std_px` | .05 | pixel |
| `linear_drift_std_px` | .02 | pixel/patch |

### 12.8 phase4_failures

| 参数 | 当前值 | 作用 |
|---|---:|---|
| `enabled` | false | 名义 profile 默认关闭 |
| `force_secondary_uav_dropout` | true | 组合 profile 强制掉一从节点 |
| `always_include_shortest` | true | 保护最短边 |
| `always_include_longest` | true | 保留压力测试长边 |
| `low_coherence_threshold` | .25 | 低相干像素阈值 |
| `low_coherence_min_fraction` | .05 | 低相干事件面积阈值 |
| `sync_jump_min_abs_rad` | .1 | 同步异常标签阈值 |

压力扫描对 Phase-4 标准差使用 1×、3×、10× 倍率；频率范围和周期范围不缩放。

## 13. MAT 数据字段

| 字段 | 形状 | 用途 |
|---|---:|---|
| `wrappedphase_withoutnoise` | `[K,H,W]` | 干净缠绕相位标签 |
| `wrappedphase_node_noise_only` | `[K,H,W]` | Phase-4 前消融观测 |
| `wrappedphase_withnoise` | `[K,H,W]` | 最终单视观测 |
| `wrappedphase_multilook` | `[K,H,W]` | 最终多视观测 |
| `unwrapped_phase` | `[K,H,W]` | 连续相位标签 |
| `coherence_terrain_only` | `[K,H,W]` | 地表调制前相干 |
| `coherence_scene` | `[K,H,W]` | 地表调制后、热噪声前相干 |
| `coherence_true` | `[K,H,W]` | 含解析热噪声因子的监督真值 |
| `coherence_observed` | `[K,H,W]` | 仅由受扰节点 SLC 估计的输入 |
| `valid_edge_mask` | `[K]` | 最终有效边 |
| `coregistration_valid_mask` | `[K,H,W]` | 重采样有效像素 |
| `edge_index` | `[2,K]` | UAV 图连接 |
| `baseline_perp_m` | `[K]` | 垂直基线 |
| `ambiguity_height_m` | `[K]` | 模糊高程 |
| `dem2phase_ratio_rad_per_m` | `[K]` | $\kappa_k$ |
| `terrain_features` | struct | 地形特征和几何代理 |
| `phase4_errors` | struct | `[N,H]` 节点误差真值和参数 |
| `failure_labels` | struct | 掉线、低相干、同步和组合标签 |
| `metadata` | struct | 来源、裁块、种子、profile 和方法 |

## 14. 主要代码入口

| 文件 | 功能 |
|---|---|
| `scripts/gen_dataset_from_dem_v13.m` | 主生成器 |
| `configs/validate_sim_config.m` | 配置与物理一致性校验 |
| `coherence/gen_coherence_map_multibaseline_terrainaware.m` | 地形相干 |
| `coherence/apply_landcover_coherence.m` | 地表覆盖调制 |
| `noise/simulate_distributed_slc_star.m` | 节点复 SLC |
| `noise/apply_distributed_uav_errors.m` | Phase-4 节点误差 |
| `dataset/apply_compound_failure_mask.m` | 掉线与组合标签 |
| `dataset/load_patch_group_for_training.m` | 严格训练加载与 mask 校验 |
| `scripts/build_failure_balanced_manifest.m` | 均衡采样 manifest |
| `scripts/run_phase4_sweep.m` | 4×3 配对压力扫描 |

## 15. 当前已知限制

1. 复 SLC 模型当前只支持共享参考节点的星形图。
2. DEM 相位采用局部远场垂直基线公式，未逐像素计算完整双基地 Tx/Rx 路径。
3. 叠掩和阴影是局部坡度代理，不是完整射线追踪。
4. Phase-4 轨迹和同步场主要沿方位行变化，尚无严格慢时间和二维空间 PSD。
5. roll/pitch 到像素位移是小角度地面投影近似，不是成像算子重聚焦。
6. 暂未加入大气、传播延迟、RFI、ADC 量化、通道增益/群时延和天线方向图。
7. WorldCover 相干因子、同步、导航和配准分布尚未由本系统实飞数据标定。
8. 当前有 19 个 DEM tile，但 validation 仍仅覆盖阿尔卑斯与安第斯各一个 tile；
   这足够做首轮外域检查，但不足以声称广泛跨地区泛化。ASTER 质量较低，只应
   用于鲁棒性训练/消融，不应作为精确高程真值。
9. `coherence_true` 是模拟器内部标签；训练输入必须使用 `coherence_observed`。
10. 现阶段生成的是聚焦后 patch 级近似，不替代原始回波仿真和完整 SAR 处理链。

## 16. 当前验证状态

- MATLAB 单元测试：41/41 通过。
- 核心新增文件 `checkcode`：0 条问题。
- 固定地理瓦片划分的数据已生成：train 162、test 56、validation 16，共 234 个
  grouped MAT；整个输出目录约 2056.8 MiB。
- 数据已统一迁移到 `coherence_observed` schema；严格加载器已在实际样本上
  对 train/test/validation 各抽样验证通过，并会拒绝旧泄漏格式。
- `failure_balanced_manifest.csv` 包含 234 行和 51 个 split×故障×地形统计层；
  train 内有 26 层，归一化训练权重范围为 0.2596--6.2308，test/validation
  权重均为 1。
- 8 张 ESA WorldCover 3° COG 已按 DEM 像元中心最近邻对齐为 19 张逐 DEM
  类别图；234 个样本均使用本 DEM 对应的 WorldCover，类别代码全集为
  10/20/30/40/50/60/70/80/90/100，未出现未知代码 0。
- 234 个 patch 的五类组件种子共 1170 个，全部唯一；批次主种子为 42，可完整
  重放同一批次，同时不会让同批次内不同 patch/组件获得相同随机实现。
- 18 个 ASTER patch 均保存有效 NUM 质量统计，其余 216 个样本明确标记为
  `num_available=false`，不会把“无 NUM”误当作高质量。
- 缺失策略已设为 `error`；新增 DEM 未准备对应的对齐 land-cover 时生成器会停止，
  不会中性回退，也不会跨 DEM 复用类别图。

## 17. Python 云端生成器迁移状态（2026-09-13）

当前 `python/dem2phase` 已覆盖生产路径
`terrain_aware + complex_slc_nodes + grouped_mat + Phase-4`，支持确定性分层随机种子、
单/多进程一致生成、原子 MAT 写入、断点续跑、严格验证和 MATLAB/Python 统计比较。
首版只支持星形共享主节点图；遇到旧 `legacy_phase_pdf` 等未移植模式会明确报错。

- Python 回归测试 9/9 通过，覆盖 1/2 worker 哈希一致、缺失文件恢复、完整任务
  快速 resume、跨平台配置哈希、Phase-4 MATLAB struct 序列化和跨语言 golden。
- 标准 Python 数据集 `data/pilot_dataset_uav_p_500m_phase4_python` 共 234 个样本：
  train 162、test 56、validation 16；19 个地理组零泄漏，1404 个组件种子唯一。
- MATLAB 严格训练加载器成功读取 234/234 个 Python 样本；每个样本的
  `phase4_errors.node_parameters` 均为 1×5 struct，不再是 cell。
- MATLAB/Python 汇总统计通过门限：平均相干度绝对差 0.01395，P05/P50/P95
  差为 0.00293/0.01636/0.01511，平均有效边数差 0.00427。
- Python WorldCover 已覆盖全部 19 个 DEM。相对现有 MATLAB 对齐图，总像素
  差异率 0.0264%，单幅最大 0.061%；差异位于浮点半像元类别边界，4 幅 ASTER
  完全一致。
- 完整数据集 resume 会先检查 generation plan 和全部接受文件；本机空操作
  resume 从约 3 分钟降至 0.762 秒，且 234 个数组哈希保持不变。
- 最新 Docker 3.11 镜像已在 Ubuntu 24.04 / Docker 29.4.0 上用真实 ALOS DEM
  验证。Windows/Linux 配置哈希与源码指纹分别完全一致为 `3bdb…13d` 和
  `a590…a9f`；容器 resume 为 0.95 秒，bind mount 文件归宿主 UID/GID 1000:1000。
- 固定 crop deterministic golden 已建立：DEM 最大差 1.1 mm，相位最大差
  1.38e-4 rad，相干度最大差 2.22e-4；卷积和线性配准达到机器精度，地表索引、
  叠掩/阴影和配准 mask 完全一致。该测试发现并修复了 Python 曲率边界采用零填充
  而非 MATLAB `del2` 延拓的问题，正式 Python 数据集已据此重新生成。

尚未完成的发布项依次为：自动化 CI；小样本 DEM（尤其安第斯）的分层统计差异
诊断。大气、RFI、ADC、通道和天线方向图等额外误差仍按计划暂缓。
