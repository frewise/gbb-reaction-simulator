# GBB 反应相关资源

## 🔍 已搜索的代码库

### 1. 化学信息学核心库

| 名称 | 链接 | 用途 |
|------|------|------|
| RDKit | https://github.com/rdkit/rdkit | 核心化学信息学库，支持反应处理 |
| Open Babel | https://github.com/openbabel/openbabel | 化学数据格式转换 |
| CDK | https://sourceforge.net/projects/cdk/ | Java 化学开发工具包 |
| Indigo | https://github.com/epam/Indigo | 通用分子工具包 |

### 2. 反应预测工具

| 名称 | 链接 | 用途 |
|------|------|------|
| aizynthfinder | https://github.com/MolecularAI/aizynthfinder | 逆合成规划工具 |
| pySiRC | https://github.com/jeffrichardchemistry/pySiRC | 反应速率常数预测 |
| RDchiral | https://github.com/connorcoley/rdchiral | RDKit 立体化学处理封装 |
| CGRtools | https://github.com/cimm-kzn/CGRtools | 反应凝缩图处理 |

### 3. 机器学习工具

| 名称 | 链接 | 用途 |
|------|------|------|
| DeepChem | https://github.com/deepchem/deepchem | 化学深度学习库 |
| Chemprop | https://github.com/chemprop/chemprop | 分子性质预测 |
| mol2vec | https://github.com/samoturk/mol2vec | 分子子结构向量表示 |
| Summit | https://github.com/sustainable-processes/summit | 化学反应优化 |

### 4. 数据资源

| 名称 | 链接 | 用途 |
|------|------|------|
| PubChem | https://pubchem.ncbi.nlm.nih.gov/ | 化合物数据库 |
| ChEMBL | https://www.ebi.ac.uk/chembl/ | 生物活性分子数据库 |
| Reaxys | https://www.reaxys.com/ | 化学反应数据库 (订阅) |

## 📚 GBB 反应文献

### 关键论文

1. **原始发现** (1998)
   - Groebke et al., Synlett 1998, 661-663
   - Blackburn et al., Tetrahedron Lett. 1998, 39, 5003-5006
   - Bienaymé et al., Eur. J. Org. Chem. 1998, 2329-2334

2. **综述文章**
   - Beilstein J. Org. Chem. 2024, 20, 1498–1506
   - DOI: https://doi.org/10.3762/bjoc.20.143

3. **应用进展**
   - 药物化学应用
   - 组合化学库构建
   - 天然产物合成

## 💻 安装指南

### RDKit 安装

**推荐方式 (Conda):**
```bash
conda create -n gbb python=3.10
conda activate gbb
conda install -c conda-forge rdkit
```

**pip 方式 (Linux/macOS):**
```bash
pip install rdkit-pypi
```

### 完整环境

```bash
# 创建环境
conda create -n gbb python=3.10 -y
conda activate gbb

# 安装核心包
conda install -c conda-forge rdkit -y
pip install pandas numpy matplotlib seaborn

# 可选：机器学习
pip install scikit-learn xgboost deepchem
```

## 🔬 验证方法

### 1. 结构验证
使用 RDKit 验证产物结构合理性:
```python
from rdkit import Chem
from rdkit.Chem import Descriptors

mol = Chem.MolFromSmiles(product_smiles)
mw = Descriptors.MolWt(mol)
logp = Descriptors.MolLogP(mol)
```

### 2. 文献对比
- 搜索 Reaxys/SciFinder 验证已知反应
- 对比实验产率与预测值

### 3. 计算验证
- DFT 计算反应能垒
- 过渡态搜索

## 📊 数据集扩展

### 建议的数据来源

1. **USPTO 专利数据库**
   - 提取 GBB 相关反应
   - 约 1000+ 已知反应

2. **Reaxys 提取**
   - 需要机构订阅
   - 高质量实验数据

3. **文献挖掘**
   - 使用 ChemDataExtractor
   - 自动提取反应条件

## 🎯 下一步行动

1. [ ] 配置 RDKit 环境
2. [ ] 运行完整反应模拟
3. [ ] 从文献收集实验数据
4. [ ] 训练产率预测模型
5. [ ] 开发 Web 界面

---

*资源整理时间：2026-03-03*
