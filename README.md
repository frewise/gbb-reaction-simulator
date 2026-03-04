# GBB 反应虚拟模拟系统

**Groebke-Blackburn-Bienaymé Reaction Simulator**

## 📋 概述

GBB 反应是一种三组分反应 (Multicomponent Reaction, MCR)：

```
醛 (Aldehyde) + 2-氨基吡啶 (Amine) + 异腈 (Isocyanide)
              → 咪唑并 [1,2-a] 吡啶衍生物
```

该反应由 Groebke、Blackburn 和 Bienaymé 于 1998 年独立发现，广泛应用于药物化学和组合化学。

## 📁 文件结构

```
gbb-reaction/
├── README.md              # 本文件
├── gbb_lite.py            # 轻量级模拟器 (无需 RDKit)
├── gbb_simulation.py      # 完整模拟器 (需要 RDKit)
├── gbb_test_data.json     # 测试数据集 (15 个反应)
└── resources.md           # 相关资源列表
```

## 🚀 快速开始

### 方式 1: 轻量级模式 (无需安装)

```bash
python3 gbb_lite.py
```

生成 15 个测试反应，包含预测产率。

### 方式 2: 完整模式 (需要 RDKit)

```bash
# 安装 RDKit
conda install -c conda-forge rdkit

# 运行完整模拟
python3 gbb_simulation.py
```

## 📊 测试结果

运行 `gbb_lite.py` 生成的测试数据:

| 统计项 | 数值 |
|--------|------|
| 反应物组合总数 | 288 种 |
| 测试样本数 | 15 个 |
| 平均预测产率 | 72.8% |
| 最高预测产率 | 85.0% |
| 最低预测产率 | 62.0% |

### 示例反应

```
反应 ID: GBB-4-n-2-a-phe
反应物：4-nitrobenzaldehyde + 2-amino-5-bromopyridine + phenyl isocyanide
产物：5-bromo-3-4-nitrobenzyl-2-phenyl-imidazo[1,2-a]pyridine
预测产率：85.0%
```

## 🔬 反应机理

GBB 反应通过以下步骤进行:

1. **亚胺形成**: 醛 + 2-氨基吡啶 → 亚胺中间体
2. **异腈加成**: 亚胺 + 异腈 → 腈鎓离子
3. **环化**: 分子内亲核进攻 → 咪唑并吡啶环
4. **质子转移**: 生成最终产物

## 📦 反应物数据库

### 醛类 (8 种)
- benzaldehyde (苯甲醛)
- 4-methylbenzaldehyde (对甲基苯甲醛)
- 4-methoxybenzaldehyde (对甲氧基苯甲醛)
- 4-chlorobenzaldehyde (对氯苯甲醛)
- 4-nitrobenzaldehyde (对硝基苯甲醛)
- furfural (糠醛)
- acetaldehyde (乙醛)
- 3-pyridinecarboxaldehyde (烟碱醛)

### 胺类 (6 种)
- 2-aminopyridine (2-氨基吡啶)
- 2-amino-3-methylpyridine
- 2-amino-4-methylpyridine
- 2-amino-5-chloropyridine
- 2-amino-5-bromopyridine
- 2-amino-6-methylpyridine

### 异腈类 (6 种)
- tert-butyl isocyanide (叔丁基异腈)
- ethyl isocyanide (乙基异腈)
- methyl isocyanide (甲基异腈)
- phenyl isocyanide (苯基异腈)
- cyclohexyl isocyanide (环己基异腈)
- 4-methylphenyl isocyanide (对甲苯基异腈)

## 🎯 产率预测规则

基于经验规则的简化预测:

| 因素 | 影响 |
|------|------|
| 吸电子基团 (硝基、卤素) | +5~10% |
| 给电子基团 (甲氧基) | -5% |
| 大位阻基团 (叔丁基) | -5% |
| 芳基异腈 | -3% |
| 杂环醛 (糠醛) | +3% |

## 🔗 相关资源

### 代码库
- **RDKit**: https://github.com/rdkit/rdkit
- **aizynthfinder**: https://github.com/MolecularAI/aizynthfinder (逆合成规划)
- **awesome-cheminformatics**: https://github.com/hsiaoyi0504/awesome-cheminformatics

### 参考文献
- 原始论文：https://doi.org/10.3762/bjoc.20.143
- Beilstein J. Org. Chem. 2024, 20, 1498–1506

### 数据库
- **PubChem**: https://pubchem.ncbi.nlm.nih.gov/
- **ChEMBL**: https://www.ebi.ac.uk/chembl/
- **Reaxys**: https://www.reaxys.com/ (需要订阅)

## 📝 数据格式

测试数据 JSON 格式:

```json
{
  "reaction_id": "GBB-xxx-xxx-xxx",
  "aldehyde": "反应物名称",
  "amine": "反应物名称",
  "isocyanide": "反应物名称",
  "product_name": "产物名称",
  "product_smiles": "产物 SMILES 或反应物 SMILES 组合",
  "status": "predicted",
  "yield_predicted": 75.0
}
```

## ⚠️ 局限性

1. **简化模型**: 当前版本使用经验规则预测，非机器学习模型
2. **SMILES 生成**: 产物 SMILES 为近似表示，精确结构需要 RDKit
3. **立体化学**: 未考虑立体选择性
4. **溶剂/温度效应**: 未包含反应条件影响

## 🔮 下一步改进

1. 集成 RDKit 进行精确反应模拟
2. 添加机器学习产率预测模型
3. 扩展反应物数据库
4. 添加反应条件优化建议
5. 集成文献数据验证

## 📄 许可证

MIT License

---

*最后更新：2026-03-03*
