# dem2phase

面向分布式多无人机 InSAR 的 DEM→干涉相位数据集生成器，提供 MATLAB 参考实现与
可在 Linux/Docker 云端独立运行的 Python 实现。

当前生产配置采用 P 波段（500 MHz）、500 m 飞行高度、45° 入射角、共享主节点的
5-UAV 星形图，以及 0.75/1.2/1.8/3.0 m 四条有效垂直基线。对应模糊高程约为
50/80/125/200 m。生成器输出 grouped MATLAB v7 文件，面向多基线相位解缠、
故障鲁棒训练和跨基线融合。

## 当前能力

- 固定地理 tile 的 train/test/validation 划分，同一 DEM 内可随机裁 patch，但不会
  跨划分泄漏。
- 地形感知相干度：坡度、坡向、局部入射角、粗糙度、曲率、TPI、叠掩与阴影代理。
- ESA WorldCover 地表类型调制，以及共享场、边独立场和复 SLC 节点热噪声。
- 同步、LOS 轨迹、姿态、亚像素配准和多节点/多边故障标签。
- BLAKE2b + PCG64DXSM 分层随机种子；结果不依赖 worker 数量或执行顺序。
- 云端 10 倍 profile 生成 2,340 个样本；可按 DEM 在本地确定性恢复云端子集，
  无需下载对应 MAT 文件。
- MATLAB/Python grouped-MAT schema 兼容、严格验证、统计比较和固定 crop golden test。
- Docker/CPU 多进程生成和可验证的断点续跑。

## 仓库结构

```text
coherence/       MATLAB 地形与相干度模型
configs/         雷达、生成、划分与 Phase-4 配置
dataset/         grouped-MAT、划分、mask 和训练加载器
noise/           节点 SLC、热噪声、同步/轨迹/姿态/配准误差
scripts/         当前 MATLAB 生成、校准、验证和 golden 入口
tests/           MATLAB 单元测试
python/          Python 3.11+ 云端生成器、测试与 Dockerfile
docs/            数学公式、字段、参数和开发状态
archive/         不再参与当前生产路径的历史实现
```

DEM、WorldCover、ASTER NUM、生成数据集和预览输出不进入 Git；参见
`.gitignore`。数据目录应通过本地磁盘、挂载卷或对象存储单独提供。

## Python 快速开始

```bash
python -m pip install -e "python[test]"
python -m pytest python/tests

python -m dem2phase generate \
  --config configs/uav_p_500m_monostatic.json \
  --seed 42 --workers 8 \
  --output data/pilot_dataset_uav_p_500m_phase4_python

python -m dem2phase validate \
  --dataset data/pilot_dataset_uav_p_500m_phase4_python
```

完整 CLI、Docker、WorldCover 与跨语言比较说明见
[`python/README.md`](python/README.md)。

## MATLAB 快速开始

```matlab
root = pwd;
addpath(genpath(root));
results = runtests(fullfile(root, 'tests'));
assertSuccess(results);

run(fullfile(root, 'scripts', 'gen_dataset_from_dem_v13.m'));
```

MATLAB 主生成器和 Python 生成器共享
[`configs/uav_p_500m_monostatic.json`](configs/uav_p_500m_monostatic.json) 与
[`configs/dem_split_manifest.csv`](configs/dem_split_manifest.csv)。

## 跨语言 golden

```matlab
export_cross_language_golden( ...
    'python/tests/fixtures/matlab_cross_language_golden.mat')
```

```bash
python python/tools/compare_cross_language_golden.py \
  --fixture python/tests/fixtures/matlab_cross_language_golden.mat \
  --config configs/uav_p_500m_monostatic.json
```

当前 fixture 中，MATLAB/Python 的 DEM 插值最大差约 1.1 mm，相位最大差
`1.38e-4 rad`，确定性相干度最大差 `2.22e-4`；卷积与线性配准达到机器精度，
离散地表和 mask 字段完全一致。

## 当前验证状态

- Python：9 项测试通过。
- MATLAB：41 项测试通过。
- Python 全量集：234 个 patch（162/56/16），19 个地理组，零 split 泄漏。
- MATLAB 严格加载 Python grouped-MAT：234/234 通过。
- Windows 与 Linux Docker 的配置哈希和源码指纹一致。
- 实际噪声和地表参数仍属于“物理合理但未由本系统实飞数据标定”的仿真配置。

详细公式、字段、配置参数和限制见
[`docs/current_code_technical_reference_zh.md`](docs/current_code_technical_reference_zh.md)。

## 历史代码

旧的逐字段目录数据格式和早期 ALOS 批处理代码已移至
[`archive/matlab_v1/`](archive/matlab_v1/README.md)。归档代码只用于追溯，不参与
测试、Docker 构建或当前数据生成。
