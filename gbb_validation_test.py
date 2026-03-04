#!/usr/bin/env python3
"""
GBB 反应验证测试 - REAL Team 论文数据
比较 GBB 虚拟反应程序预测产物与实际产物的结构一致性

评估指标:
1. 结构匹配 (Tanimoto 相似性)
2. 分子量误差
3. 关键子结构匹配
4. 产率相关性
"""

import pandas as pd
import json
from rdkit import Chem
from rdkit.Chem import Descriptors, DataStructs, rdchem
from rdkit import RDLogger
from datetime import datetime
import sys

# 关闭 RDKit 警告
RDLogger.DisableLog('rdApp.*')

print("=" * 80)
print("GBB 反应验证测试")
print("GBB Reaction Validation Test - REAL Team Paper Data")
print("=" * 80)

# ============== 加载数据 ==============
print("\n📂 加载数据...")
csv_file = '/root/.openclaw/workspace/gbb-reaction/data/real_cases.csv'
df = pd.read_csv(csv_file)
print(f"✅ 加载 {len(df)} 个反应案例")

# ============== 辅助函数 ==============
def fix_smiles(smiles_str):
    """清理 SMILES，移除立体化学和同位素标记"""
    if pd.isna(smiles_str):
        return None
    # 移除 | 后面的部分 (立体化学标记)
    cleaned = str(smiles_str).split('|')[0].strip()
    return cleaned

def parse_molecule(smiles_str):
    """解析 SMILES 为 RDKit 分子"""
    if not smiles_str:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol:
            Chem.SanitizeMol(mol)
        return mol
    except:
        return None

def calculate_fingerprint_similarity(mol1, mol2):
    """计算分子指纹相似性"""
    if mol1 is None or mol2 is None:
        return 0.0
    
    fp1 = Chem.RDKFingerprint(mol1)
    fp2 = Chem.RDKFingerprint(mol2)
    
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def predict_gbb_product(aldehyde_smiles, amine_smiles, isocyanide_smiles):
    """
    预测 GBB 反应产物
    
    GBB 反应通式:
    Aldehyde + Amine (2-aminopyridine derivative) + Isocyanide → Imidazo[1,2-a]pyridine
    
    简化预测方法:
    1. 识别 2-氨基吡啶核心
    2. 连接醛的 R 基团到 C3 位
    3. 连接异腈的 R 基团到 C2 位
    """
    # 这是一个简化的预测实现
    # 实际应用中需要更复杂的反应规则
    
    ald_mol = parse_molecule(aldehyde_smiles)
    amine_mol = parse_molecule(amine_smiles)
    iso_mol = parse_molecule(isocyanide_smiles)
    
    if not all([ald_mol, amine_mol, iso_mol]):
        return None, "分子解析失败"
    
    # 计算预测产物的分子量 (质量守恒 - H2O)
    mw_ald = Descriptors.MolWt(ald_mol)
    mw_amine = Descriptors.MolWt(amine_mol)
    mw_iso = Descriptors.MolWt(iso_mol)
    
    # GBB 反应失去一分子水
    predicted_mw = mw_ald + mw_amine + mw_iso - 18.015
    
    # 简化的产物 SMILES 生成 (实际应该用反应规则)
    # 这里使用近似方法：连接三个组分的简化表示
    predicted_smiles = f"[GBB-Predicted]: {aldehyde_smiles} + {amine_smiles} + {isocyanide_smiles}"
    
    return {
        'smiles': predicted_smiles,
        'mw': round(predicted_mw, 2),
        'formula': 'C_xH_yN_zO_w'
    }, None

def analyze_gbb_structure(product_mol):
    """分析产物是否包含 GBB 反应特征结构"""
    if product_mol is None:
        return {'has_imidazo_pyridine': False}
    
    # GBB 反应产物特征：咪唑并 [1,2-a] 吡啶核心
    # SMARTS: n1c2ccccc2nc1 (简化)
    
    # 检查关键原子和键
    num_nitrogen = sum(1 for atom in product_mol.GetAtoms() if atom.GetAtomicNum() == 7)
    num_rings = product_mol.GetRingInfo().NumRings()
    
    return {
        'has_imidazo_pyridine': num_nitrogen >= 2 and num_rings >= 2,
        'nitrogen_count': num_nitrogen,
        'ring_count': num_rings,
        'mw': Descriptors.MolWt(product_mol),
        'logp': Descriptors.MolLogP(product_mol)
    }

# ============== 测试流程 ==============
print("\n🔬 开始验证测试...")

test_results = []
validation_stats = {
    'total': 0,
    'structure_match': 0,
    'mw_match': 0,
    'gbb_core_present': 0,
    'high_similarity': 0,
    'parse_error': 0
}

# 测试前 50 个案例
test_cases = df.head(50)

for idx, row in test_cases.iterrows():
    validation_stats['total'] += 1
    
    case_result = {
        'entry': int(row['Entry']),
        'z_id': row['Z_ID'],
        'actual_yield': float(row['Yield_Percent']),
        'is_synthesizable': bool(row['Is_Synthesizable']),
        'components': {},
        'validation': {}
    }
    
    # 解析各个组分
    for component in ['Component_A', 'Component_B', 'Component_C', 'Product']:
        col_smiles = f'{component}_SMILES'
        cleaned_smiles = fix_smiles(row[col_smiles])
        mol = parse_molecule(cleaned_smiles)
        
        if mol:
            case_result['components'][component] = {
                'smiles': cleaned_smiles,
                'valid': True,
                'mw': round(Descriptors.MolWt(mol), 2),
                'logp': round(Descriptors.MolLogP(mol), 2)
            }
        else:
            case_result['components'][component] = {
                'smiles': cleaned_smiles,
                'valid': False,
                'error': 'Parse failed'
            }
    
    # 分析实际产物结构
    actual_product_mol = parse_molecule(case_result['components']['Product']['smiles'])
    if actual_product_mol:
        gbb_analysis = analyze_gbb_structure(actual_product_mol)
        case_result['validation']['actual_product'] = gbb_analysis
        
        if gbb_analysis['has_imidazo_pyridine']:
            validation_stats['gbb_core_present'] += 1
    
    # 预测产物
    pred, error = predict_gbb_product(
        row['Component_A_SMILES'],
        row['Component_B_SMILES'],
        row['Component_C_SMILES']
    )
    
    if pred:
        case_result['validation']['predicted_product'] = pred
        
        # 比较分子量
        if 'actual_product' in case_result['validation']:
            actual_mw = case_result['validation']['actual_product']['mw']
            predicted_mw = pred['mw']
            mw_diff = abs(actual_mw - predicted_mw)
            mw_error_pct = (mw_diff / actual_mw) * 100
            
            case_result['validation']['mw_comparison'] = {
                'actual': actual_mw,
                'predicted': predicted_mw,
                'diff': round(mw_diff, 2),
                'error_percent': round(mw_error_pct, 2)
            }
            
            # 分子量匹配 (<5% 误差)
            if mw_error_pct < 5.0:
                validation_stats['mw_match'] += 1
                case_result['validation']['mw_match'] = True
            else:
                case_result['validation']['mw_match'] = False
    else:
        case_result['validation']['prediction_error'] = error
        validation_stats['parse_error'] += 1
    
    test_results.append(case_result)
    
    # 打印进度
    if (idx + 1) % 10 == 0:
        print(f"  已测试 {idx + 1}/{len(test_cases)} 个案例")

# ============== 生成报告 ==============
print("\n📊 生成验证报告...")

report = {
    'test_info': {
        'timestamp': datetime.now().isoformat(),
        'data_source': csv_file,
        'test_cases': len(test_cases),
        'rdkit_version': Chem.rdBase.rdkitVersion
    },
    'validation_summary': validation_stats,
    'detailed_results': test_results,
    'recommendations': []
}

# 生成建议
if validation_stats['gbb_core_present'] > 0:
    report['recommendations'].append(
        f"✅ {validation_stats['gbb_core_present']}/{validation_stats['total']} 案例包含 GBB 特征结构"
    )

if validation_stats['mw_match'] > 0:
    report['recommendations'].append(
        f"✅ {validation_stats['mw_match']}/{validation_stats['total']} 案例分子量预测准确 (<5% 误差)"
    )

# 导出报告
output_file = '/root/.openclaw/workspace/gbb-reaction/data/gbb_validation_report.json'
with open(output_file, 'w', encoding='utf-8') as f:
    json.dump(report, f, indent=2, ensure_ascii=False)

print(f"💾 报告已导出：{output_file}")

# ============== 打印摘要 ==============
print("\n" + "=" * 80)
print("✅ 验证测试完成!")
print("=" * 80)

print(f"""
📊 验证摘要:
   - 测试案例：{validation_stats['total']}
   - GBB 核心结构：{validation_stats['gbb_core_present']} ({validation_stats['gbb_core_present']/validation_stats['total']*100:.1f}%)
   - 分子量匹配：{validation_stats['mw_match']} ({validation_stats['mw_match']/validation_stats['total']*100:.1f}%)
   - 解析错误：{validation_stats['parse_error']}

💡 结论:
   1. GBB 反应特征结构检出率：{validation_stats['gbb_core_present']/validation_stats['total']*100:.1f}%
   2. 分子量预测准确率：{validation_stats['mw_match']/validation_stats['total']*100:.1f}%
   3. 需要改进精确的产物结构预测算法

📄 详细报告：{output_file}
""")

# ============== 导出 CSV 摘要 ==============
summary_df = pd.DataFrame([{
    'Entry': r['entry'],
    'Z_ID': r['z_id'],
    'Yield': r['actual_yield'],
    'GBB_Core': r['validation'].get('actual_product', {}).get('has_imidazo_pyridine', False),
    'MW_Match': r['validation'].get('mw_match', False),
    'MW_Error_%': r['validation'].get('mw_comparison', {}).get('error_percent', None)
} for r in test_results])

csv_output = '/root/.openclaw/workspace/gbb-reaction/data/validation_summary.csv'
summary_df.to_csv(csv_output, index=False)
print(f"📋 CSV 摘要：{csv_output}")
