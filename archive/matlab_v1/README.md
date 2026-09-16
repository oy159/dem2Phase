# MATLAB v1 archive

这里保存 2026-09 仓库整理时退出生产路径的早期 MATLAB 实现。

归档原因：

- `scripts/gen_dataset_from_dem.m` 输出按字段拆分的旧目录结构，已由 grouped-MAT
  多基线生成器 `scripts/gen_dataset_from_dem_v13.m` 取代。
- `scripts/validate_dataset.m` 和 `scripts/visualize_patches.m` 只理解上述旧目录结构。
- `scripts/analyze_wrap_count_distribution.m` 硬编码旧 `training_dataset_v13` 路径，
  不理解当前 split 和 grouped-MAT manifest。
- `helpers/` 中的 ALOS 批处理、裁 patch 和加噪脚本使用硬编码本机路径，并服务于
  更早的单基线工作流。

这些文件未删除，以便复现实验和追踪算法演变。它们不参与当前测试、Python 包、
Docker 镜像或正式数据集生成。若需恢复其中功能，应先改写为读取当前 JSON 配置、
固定 split manifest 和 grouped-MAT schema，而不是直接重新启用。
