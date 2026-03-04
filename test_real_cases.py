#!/usr/bin/env python3
"""
GBB 反应真实案例测试
测试 GBB 虚拟反应程序对 REAL 团队论文中真实案例的预测能力

数据来源：real_cases.csv (GitHub: frewise/gbb-reaction-simulator)
"""

import pandas as pd
import json
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolAlign
from datetime import datetime
import sys

print("=" * 80)
print("GBB 反应真实案例测试")
print("REAL Team Paper - GBB Reaction Cases")
print("=" * 80)

# ============== 加载数据 ==============
print("\n📂 加载数据...")
csv_file = '/root/.openclaw/workspace/gbb-reaction/data/real_cases.csv'
df = pd.read_csv(csv_file)

print(f"✅ 加载成功：{len(df)} 个反应案例")
print(f"   文件：{csv_file}")
print(f"   列：{', '.join(df.columns)}")

# ============== 数据探索 ==============
print("\n📊 数据概览:")
print(f"   总反应数：{len(df)}")
print(f"   可合成：{df['Is_Synthesizable'].sum()} / {len(df)}")
print(f"   平均产率：{df['Yield_Percent'].mean():.2f}%")
print(f"   产率范围：{df['Yield_Percent'].min():.1f}% - {df['Yield_Percent'].max():.1f}%")

# ============== 分子验证函数 ==============
def validate_molecule(smiles):
    """验证 SMILES 字符串"""
    try:
        # 清理 SMILES (移除立体化学标记等)
        clean_smiles = smiles.split('|')[0].strip()
        mol = Chem.MolFromSmiles(clean_smiles)
        if mol:
            Chem.SanitizeMol(mol)
            return True, mol, clean_smiles
        return False, None, clean_smiles
    except Exception as e:
        return False, None, str(e)

def calculate_molecular_properties(mol):
    """计算分子描述符"""
    if mol is None:
        return {}
    return {
        'MW': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'HBA': Descriptors.NumHAcceptors(mol),
        'HBD': Descriptors.NumHDonors(mol),
        'TPSA': Descriptors.TPSA(mol),
        'RotatableBonds': Descriptors.NumRotatableBonds(mol),
    }

def compare_molecules(mol1, mol2):
    """比较两个分子的相似性"""
    if mol1 is None or mol2 is None:
        return {'similarity': 0.0, 'match': False}
    
    # 计算指纹
    fp1 = Chem.RDKFingerprint(mol1)
    fp2 = Chem.RDKFingerprint(mol2)
    
    # Tanimoto 相似性
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    
    # 检查是否完全匹配
    match = (similarity == 1.0)
    
    return {
        'similarity': similarity,
        'match': match,
        'mw_diff': abs(Descriptors.MolWt(mol1) - Descriptors.MolWt(mol2))
    }

# ============== 测试反应物解析 ==============
print("\n🔬 测试反应物解析...")

test_results = []
parse_stats = {'A': 0, 'B': 0, 'C': 0, 'P': 0}

for idx, row in df.head(10).iterrows():
    print(f"\n案例 {row['Entry']}:")
    
    # 验证各个组分
    for component in ['Component_A', 'Component_B', 'Component_C', 'Product']:
        col_smiles = f'{component}_SMILES'
        valid, mol, clean_smiles = validate_molecule(row[col_smiles])
        
        key = 'P' if component == 'Product' else component[-1]
        if valid:
            parse_stats[key] += 1
            props = calculate_molecular_properties(mol)
            print(f"  ✅ {component}: MW={props['MW']:.2f}, LogP={props['LogP']:.2f}")
        else:
            print(f"  ⚠️  {component}: 解析失败 - {clean_smiles[:50]}")

# ============== 生成测试报告 ==============
print("\n" + "=" * 80)
print("📋 生成测试报告...")

report = {
    'test_info': {
        'timestamp': datetime.now().isoformat(),
        'data_source': csv_file,
        'total_cases': len(df),
        'rdkit_version': Chem.rdBase.rdkitVersion
    },
    'data_summary': {
        'total_reactions': len(df),
        'synthesizable': int(df['Is_Synthesizable'].sum()),
        'average_yield': float(df['Yield_Percent'].mean()),
        'min_yield': float(df['Yield_Percent'].min()),
        'max_yield': float(df['Yield_Percent'].max()),
        'median_yield': float(df['Yield_Percent'].median())
    },
    'parsing_stats': parse_stats,
    'sample_cases': []
}

# 详细分析前 20 个案例
print("\n🔍 详细分析前 20 个案例...")

for idx, row in df.head(20).iterrows():
    case_data = {
        'entry': int(row['Entry']),
        'z_id': row['Z_ID'],
        'yield_percent': float(row['Yield_Percent']),
        'components': {}
    }
    
    for component in ['Component_A', 'Component_B', 'Component_C', 'Product']:
        col_smiles = f'{component}_SMILES'
        valid, mol, clean_smiles = validate_molecule(row[col_smiles])
        
        if valid and mol:
            props = calculate_molecular_properties(mol)
            case_data['components'][component] = {
                'smiles': clean_smiles,
                'valid': True,
                'properties': {k: round(v, 2) for k, v in props.items()}
            }
        else:
            case_data['components'][component] = {
                'smiles': clean_smiles,
                'valid': False,
                'error': str(clean_smiles)
            }
    
    report['sample_cases'].append(case_data)
    print(f"  案例 {row['Entry']}: 产物 MW={case_data['components']['Product']['properties'].get('MW', 'N/A')}")

# ============== 产率分布分析 ==============
print("\n📈 产率分布分析...")

yield_ranges = [
    (0, 1, '极低 (0-1%)'),
    (1, 2, '低 (1-2%)'),
    (2, 5, '中等 (2-5%)'),
    (5, 10, '良好 (5-10%)'),
    (10, 100, '优秀 (>10%)')
]

print("\n产率分布:")
for low, high, label in yield_ranges:
    count = len(df[(df['Yield_Percent'] >= low) & (df['Yield_Percent'] < high)])
    percentage = count / len(df) * 100
    print(f"  {label:15} {count:4} 个反应 ({percentage:5.1f}%)")

# ============== 导出报告 ==============
output_file = '/root/.openclaw/workspace/gbb-reaction/data/real_cases_analysis_report.json'
with open(output_file, 'w', encoding='utf-8') as f:
    json.dump(report, f, indent=2, ensure_ascii=False)

print(f"\n💾 报告已导出：{output_file}")

# ============== 总结 ==============
print("\n" + "=" * 80)
print("✅ 测试完成!")
print("=" * 80)

print(f"""
📊 数据摘要:
   - 总反应数：{len(df)}
   - 可合成：{df['Is_Synthesizable'].sum()} ({df['Is_Synthesizable'].mean()*100:.1f}%)
   - 平均产率：{df['Yield_Percent'].mean():.2f}%
   - 产率中位数：{df['Yield_Percent'].median():.2f}%

🔍 解析统计:
   - Component A: {parse_stats['A']}/20 成功
   - Component B: {parse_stats['B']}/20 成功
   - Component C: {parse_stats['C']}/20 成功
   - Product: {parse_stats['P']}/20 成功

💡 下一步:
   1. 使用 GBB 反应预测程序预测产物
   2. 对比预测产物与实际产物
   3. 计算结构一致性
   4. 分析失败案例

📄 详细报告：{output_file}
""")
