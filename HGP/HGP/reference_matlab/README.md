# HGP MATLAB 参考脚本

本目录保存自 `cpp/HGP/MatlabScripts/HGP/` 归档的 **HGP 原版 `.m` 文件**（自 git 历史 `2210e719^` 恢复），仅供对照阅读与复现，**运行时不必全部加入 MATLAB path**。

## 文件列表

| 脚本 | C++ 迁移状态 | 说明 |
|:---|:---|:---|
| `HGP_iteration.m` | ❌ 仍在用 | SOCP 核心；运行副本亦在 `HGP/HGP_iteration.m` |
| `computeGrads.m` | ✅ 已迁 | → `HGP::computeGradientsInCpp` |
| `getFramesFromVectorField.m` | ✅ 已迁 | → `computeFramesFromVectorFieldInCpp` |
| `setFrames1.m` / `setFrames2.m` | ✅ 已迁 | → `updateFramesFromCurrentFzInCpp` |
| `frameFixABF.m` | ✅ 已迁 | → `FrameFixingClass::extractOneRingAngles` |
| `HGP_settings.m` | ⚠️ GUI 已删 | 参数现写死在 `HGP::run` |
| `HGP_report.m` | ⚠️ GUI 已删 | 结果写 `HGP.Result.*` workspace |

对照表见 [../MATLAB与C++对照表.md](../MATLAB与C++对照表.md)。

## 运行 HGP 混合管线时

仍需 MATLAB path 包含 **`HGP/`**（含 `HGP_iteration.m`），并安装 **CVX + MOSEK**。  
`reference_matlab/` 中的已迁移脚本**不需要**再加入 path，除非你要与 C++ 逐行对比或临时回退到纯 MATLAB 辅助逻辑。
