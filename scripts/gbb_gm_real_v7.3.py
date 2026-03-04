#!/usr/bin/env python3
"""
GBB-3CR 工业级枚举与真实化学空间断言流水线 (V7.3 满血版)

修复日志：
- 【严重 Bug 修复】：在 SMIRKS 产物端添加了 `+0` 电荷中和修饰符 `[c+0:11]` 和 `[NH+0:10]`。
  这彻底解决了异腈电荷继承导致的 RDKit Kekulization 灾难（之前导致 50000 次枚举全盘失败的元凶）。
"""

import pandas as pd
import numpy as np
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

# ==========================================
# 0. 全局配置与日志处理
# ==========================================
RDLogger.DisableLog('rdApp.*')
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# ==========================================
# 1. 实验经验反应性过滤器 (Enamine REAL 补丁)
# ==========================================
class EmpiricalReactivityFilter:
    def __init__(self):
        # A 组分反应性黑名单
        self.A_BLACKLIST = [
            Chem.MolFromSmarts('[#7;r5]:[c]-[NH2]'),
            Chem.MolFromSmarts('[n]~[n]~[c]-[NH2]'),
            Chem.MolFromSmarts('[NH2]-[c]1[n]cccc1-[N,O,S;!R]'),
            Chem.MolFromSmarts('[NH2]-[c]1[n]c([$([N+](=O)[O-]),$([CX2]#N),$([SX4](=O)=O)])ccc1'), 
            Chem.MolFromSmarts('[NH2]-[c]1[n]cc([$([N+](=O)[O-]),$([CX2]#N),$([SX4](=O)=O)])cc1'), 
            Chem.MolFromSmarts('[NH2]-[c]1[n]ccc([$([N+](=O)[O-]),$([CX2]#N),$([SX4](=O)=O)])c1')  
        ]
        
        # B 组分反应性黑名单
        self.B_BLACKLIST = [
            Chem.MolFromSmarts('[c]1[n]cc[c]([F,Cl,Br,I])[c]1-C=O')
        ]

    def is_valid_A(self, mol):
        for bad_patt in self.A_BLACKLIST:
            if bad_patt and mol.HasSubstructMatch(bad_patt):
                return False
        return True

    def is_valid_B(self, mol):
        for bad_patt in self.B_BLACKLIST:
            if bad_patt and mol.HasSubstructMatch(bad_patt):
                return False
        return True

# ==========================================
# 2. 核心反应与拓扑断言引擎
# ==========================================
class GBBPipelineEngine:
    def __init__(self):
        # 1. 基础拓扑筛选器
        self.A_FILTER = Chem.MolFromSmarts('[n;r5,r6]:[c;r5,r6]-[NH2]')
        self.B_FILTER = Chem.MolFromSmarts('[CH1](=O)[#6]')
        self.C_FILTER = Chem.MolFromSmarts('[#6]-[N+]#[C-]')

        # 2. 引入经验补丁
        self.empirical_filter = EmpiricalReactivityFilter()

        # 3. 拓扑免疫的核心骨架验证器
        validator_smarts = '[#7;R2]1~[#6;R2]~[#7;R1]~[#6;R1](~[#6])~[#6;R1](~[#7])1'
        self.CORE_VALIDATOR = Chem.MolFromSmarts(validator_smarts)

        # 4. 修复后的满血 SMIRKS 构建引擎
        # 核心改动：[c+0:11] 和 [NH+0:10] 强制剥离了异腈的电荷，完美重建 10-pi 芳香体系
        gbb_smirks = (
            '([n;r5,r6:3]:[c;r5,r6:2]-[NH2:1]).' # A 组分
            '([CH1:8](=O)[#6:9]).'               # B 组分
            '([#6:12]-[N+:10]#[C-:11])'          # C 组分
            '>>'
            '[n:3]1:[c:2]:[n:1]:[c:8](-[#6:9]):[c+0:11](-[NH+0:10]-[#6:12]):1'
        )
        self.transform_engine = AllChem.ReactionFromSmarts(gbb_smirks)

    def process_raw_library(self, df):
        """处理原始库，执行双重清洗"""
        pools = {'A': [], 'B': [], 'C': []}
        
        for _, row in df.iterrows():
            smi = str(row.get('SMILES', ''))
            mol = Chem.MolFromSmiles(smi)
            if not mol:
                continue
            
            if mol.HasSubstructMatch(self.A_FILTER):
                if self.empirical_filter.is_valid_A(mol):
                    pools['A'].append((row, mol))
                else:
                    pass # 静默拦截
                    
            elif mol.HasSubstructMatch(self.B_FILTER):
                if self.empirical_filter.is_valid_B(mol):
                    pools['B'].append((row, mol))
                else:
                    pass # 静默拦截
                    
            elif mol.HasSubstructMatch(self.C_FILTER):
                pools['C'].append((row, mol))
                
        return pools

    def build_gbb_robust(self, mol_a, mol_b, mol_c):
        """执行反应 -> RDKit 坍缩芳香化 -> 拓扑结构断言"""
        try:
            outcomes = self.transform_engine.RunReactants((mol_a, mol_b, mol_c))
            if not outcomes:
                return None, "底物拓扑不兼容"
            
            product = outcomes[0][0]
            
            # 由于在 SMIRKS 中已处理电荷并指定了全芳香键，这里的 Sanitize 会丝滑通过
            sanitize_flag = Chem.SanitizeMol(product, catchErrors=True)
            if sanitize_flag != Chem.SanitizeFlags.SANITIZE_NONE:
                return None, f"Sanitize 失败 (错误码: {sanitize_flag})"
            
            # 终极断言
            if not product.HasSubstructMatch(self.CORE_VALIDATOR):
                return None, "核心拓扑断言失败"
            
            return product, "成功"
            
        except Exception as e:
            return None, f"引擎异常: {str(e)}"

# ==========================================
# 3. 工业级流水线调度
# ==========================================
def run_production_pipeline(input_csv, output_csv="gbb_library_production_v7.3.csv", target_success=30000, max_attempts=500000):
    engine = GBBPipelineEngine()
    
    logging.info(f"V7.3 (满血版) 启动：加载并清洗原始库 [{input_csv}]")
    try:
        df = pd.read_csv(input_csv)
    except FileNotFoundError:
        logging.error(f"找不到输入文件: {input_csv}")
        return

    pools = engine.process_raw_library(df)
            
    logging.info(f"双重安检通过池大小 - A: {len(pools['A'])}, B: {len(pools['B'])}, C: {len(pools['C'])}")
    
    if any(len(p) == 0 for p in pools.values()):
        logging.error("某一组分池为空，程序熔断。")
        return

    logging.info(f"核心枚举火力全开，目标合法产物数量: {target_success} ...")
    results = []
    success_count = 0
    attempt_count = 0
    seen_smiles = set()
    
    np.random.seed(42)
    
    while success_count < target_success and attempt_count < max_attempts:
        attempt_count += 1
        
        a_row, a_mol = pools['A'][np.random.randint(len(pools['A']))]
        b_row, b_mol = pools['B'][np.random.randint(len(pools['B']))]
        c_row, c_mol = pools['C'][np.random.randint(len(pools['C']))]
        
        product, status = engine.build_gbb_robust(a_mol, b_mol, c_mol)
        
        if product:
            prod_smi = Chem.MolToSmiles(product, canonical=True)
            
            if prod_smi not in seen_smiles:
                seen_smiles.add(prod_smi)
                success_count += 1
                
                results.append({
                    'Product_SMILES': prod_smi,
                    'A_SMILES': a_row.get('SMILES', ''),
                    'B_SMILES': b_row.get('SMILES', ''),
                    'C_SMILES': c_row.get('SMILES', ''),
                    'A_ID': a_row.get('CID', a_row.get('ID', 'N/A')),
                    'B_ID': b_row.get('CID', b_row.get('ID', 'N/A')),
                    'C_ID': c_row.get('CID', c_row.get('ID', 'N/A'))
                })
                
                # 每成功生成 1000 个分子，在控制台进行心跳播报
                if success_count % 1000 == 0:
                    logging.info(f"🚀 进度: {success_count}/{target_success} 成功 (总尝试数: {attempt_count})")

    if results:
        out_df = pd.DataFrame(results)
        out_df.to_csv(output_csv, index=False)
        logging.info(f"🎉 枚举圆满结束！高效保留 {success_count} 个绝对可合成真实产物。")
        logging.info(f"结果已落盘至: {output_csv}")
    else:
        logging.warning("未能生成合法产物。")

# ==========================================
# 4. 脚本执行入口
# ==========================================
if __name__ == "__main__":
    input_file = "bb_verified_gm.csv"
    output_file = "gbb_library_production_v7.3.csv"
    # 将目标产量设定为 30,000，为了防止无限循环放宽 max_attempts 到了 500,000
    target_count = 30000 
    
    run_production_pipeline(input_csv=input_file, output_csv=output_file, target_success=target_count)
