#!/usr/bin/env python3
"""
GBB v7.3 引擎验证测试
使用 REAL Team real_cases.csv 数据测试工业级 GBB 反应引擎
"""

import pandas as pd
import json
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdchem
from rdkit import RDLogger
from datetime import datetime
import sys

RDLogger.DisableLog('rdApp.*')

print("=" * 80)
print("GBB v7.3 引擎验证测试")
print("GBB v7.3 Engine Validation Test")
print("=" * 80)

# ============== v7.3 引擎 ==============
class GBBv73Engine:
    """GBB-3CR 工业级枚举引擎 V7.3"""
    
    def __init__(self):
        # SMIRKS 反应引擎 (修复电荷问题)
        gbb_smirks = (
            '([n;r5,r6:3]:[c;r5,r6:2]-[NH2:1]).'
            '([CH1:8](=O)[#6:9]).'
            '([#6:12]-[N+:10]#[C-:11])'
            '>>'
            '[n:3]1:[c:2]:[n:1]:[c:8](-[#6:9]):[c+0:11](-[NH+0:10]-[#6:12]):1'
        )
        self.transform_engine = AllChem.ReactionFromSmarts(gbb_smirks)
        
        # 核心骨架验证器
        validator_smarts = '[#7;R2]1~[#6;R2]~[#7;R1]~[#6;R1](~[#6])~[#6;R1](~[#7])1'
        self.CORE_VALIDATOR = Chem.MolFromSmarts(validator_smarts)
    
    def predict_product(self, smiles_a, smiles_b, smiles_c):
        """
        预测 GBB 反应产物
        
        参数:
            smiles_a: A 组分 (2-aminopyridine derivative)
            smiles_b: B 组分 (aldehyde)
            smiles_c: C 组分 (isocyanide)
        
        返回:
            (product_mol, status_message)
        """
        try:
            mol_a = Chem.MolFromSmiles(smiles_a)
            mol_b = Chem.MolFromSmiles(smiles_b)
            mol_c = Chem.MolFromSmiles(smiles_c)
            
            if not all([mol_a, mol_b, mol_c]):
                return None, "反应物解析失败"
            
            # 执行反应
            outcomes = self.transform_engine.RunReactants((mol_a, mol_b, mol_c))
            
            if not outcomes:
                return None, "反应引擎无产物"
            
            product = outcomes[0][0]
            
            # 清理和验证
            sanitize_flag = Chem.SanitizeMol(product, catchErrors=True)
            if sanitize_flag != Chem.SanitizeFlags.SANITIZE_NONE:
                return None, f"Sanitize 失败：{sanitize_flag}"
            
            # 核心骨架断言
            if not product.HasSubstructMatch(self.CORE_VALIDATOR):
                return None, "核心拓扑断言失败"
            
            return product, "成功"
            
        except Exception as e:
            return None, f"引擎异常：{str(e)}"

# ============== 加载数据 ==============
print("\n📂 加载数据...")
csv_file = '/root/.openclaw/workspace/gbb-reaction/data/real_cases.csv'
df = pd.read_csv(csv_file)
print(f"✅ 加载 {len(df)} 个反应案例")

# ============== 辅助函数 ==============
def clean_smiles(smiles_str):
    """清理 SMILES，移除立体化学标记"""
    if pd.isna(smiles_str):
        return None
    return str(smiles_str).split('|')[0].strip()

def calculate_similarity(mol1, mol2):
    """计算 Tanimoto 相似性"""
    if mol1 is None or mol2 is None:
        return 0.0
    
    fp1 = Chem.RDKFingerprint(mol1)
    fp2 = Chem.RDKFingerprint(mol2)
    
    from rdkit import DataStructs
    return DataStructs.TanimotoSimilarity(fp1, fp2)

# ============== 测试流程 ==============
print("\n🔬 开始 v7.3 引擎测试...")

engine = GBBv73Engine()
test_results = []

stats = {
    'total': 0,
    'success': 0,
    'exact_match': 0,
    'high_similarity': 0,
    'core_match': 0,
    'mw_match': 0,
    'failed': 0
}

# 测试全部 2052 个案例
test_cases = df
print(f"🔬 测试模式：全部 {len(df)} 个案例")

for idx, row in test_cases.iterrows():
    stats['total'] += 1
    
    case_result = {
        'entry': int(row['Entry']),
        'z_id': row['Z_ID'],
        'yield_percent': float(row['Yield_Percent']),
        'components': {},
        'prediction': {}
    }
    
    # 清理 SMILES
    smiles_a = clean_smiles(row['Component_A_SMILES'])
    smiles_b = clean_smiles(row['Component_B_SMILES'])
    smiles_c = clean_smiles(row['Component_C_SMILES'])
    actual_product = clean_smiles(row['Product_SMILES'])
    
    case_result['components'] = {
        'A': smiles_a,
        'B': smiles_b,
        'C': smiles_c,
        'Actual_Product': actual_product
    }
    
    # 使用 v7.3 引擎预测
    product_mol, status = engine.predict_product(smiles_a, smiles_b, smiles_c)
    
    if product_mol:
        stats['success'] += 1
        
        pred_smiles = Chem.MolToSmiles(product_mol, canonical=True)
        pred_mw = Descriptors.MolWt(product_mol)
        
        case_result['prediction'] = {
            'status': status,
            'smiles': pred_smiles,
            'mw': round(pred_mw, 2),
            'logp': round(Descriptors.MolLogP(product_mol), 2)
        }
        
        # 对比实际产物
        actual_mol = Chem.MolFromSmiles(actual_product)
        if actual_mol:
            actual_mw = Descriptors.MolWt(actual_mol)
            
            # 分子量对比
            mw_diff = abs(pred_mw - actual_mw)
            mw_error_pct = (mw_diff / actual_mw) * 100
            
            case_result['comparison'] = {
                'actual_mw': round(actual_mw, 2),
                'predicted_mw': round(pred_mw, 2),
                'mw_diff': round(mw_diff, 2),
                'mw_error_percent': round(mw_error_pct, 2)
            }
            
            # 分子量匹配 (<5% 误差)
            if mw_error_pct < 5.0:
                stats['mw_match'] += 1
                case_result['comparison']['mw_match'] = True
            else:
                case_result['comparison']['mw_match'] = False
            
            # 结构相似性
            similarity = calculate_similarity(product_mol, actual_mol)
            case_result['comparison']['tanimoto_similarity'] = round(similarity, 3)
            
            if similarity == 1.0:
                stats['exact_match'] += 1
                case_result['comparison']['match_type'] = 'Exact Match'
            elif similarity > 0.85:
                stats['high_similarity'] += 1
                case_result['comparison']['match_type'] = 'High Similarity'
            else:
                case_result['comparison']['match_type'] = 'Low Similarity'
            
            # 核心骨架验证
            if product_mol.HasSubstructMatch(engine.CORE_VALIDATOR):
                stats['core_match'] += 1
                case_result['prediction']['core_validated'] = True
            else:
                case_result['prediction']['core_validated'] = False
        else:
            case_result['comparison'] = {'error': 'Actual product parse failed'}
            stats['failed'] += 1
    else:
        case_result['prediction'] = {
            'status': status,
            'success': False
        }
        stats['failed'] += 1
    
    test_results.append(case_result)
    
    # 进度显示
    if (idx + 1) % 20 == 0:
        print(f"  已测试 {idx + 1}/{len(test_cases)} | 成功：{stats['success']} | 精确匹配：{stats['exact_match']}")

# ============== 生成报告 ==============
print("\n📊 生成测试报告...")

report = {
    'test_info': {
        'timestamp': datetime.now().isoformat(),
        'engine_version': 'v7.3',
        'data_source': csv_file,
        'test_cases': len(test_cases),
        'rdkit_version': Chem.rdBase.rdkitVersion
    },
    'summary_statistics': stats,
    'detailed_results': test_results
}

# 导出 JSON 报告
json_output = '/root/.openclaw/workspace/gbb-reaction/data/v73_engine_validation.json'
with open(json_output, 'w', encoding='utf-8') as f:
    json.dump(report, f, indent=2, ensure_ascii=False)

# 导出 CSV 摘要
summary_df = pd.DataFrame([{
    'Entry': r['entry'],
    'Z_ID': r['z_id'],
    'Yield': r['yield_percent'],
    'Prediction_Status': r['prediction'].get('status', 'Failed'),
    'Success': r['prediction'].get('status') == '成功',
    'MW_Match': r.get('comparison', {}).get('mw_match', False),
    'MW_Error_%': r.get('comparison', {}).get('mw_error_percent', None),
    'Tanimoto': r.get('comparison', {}).get('tanimoto_similarity', None),
    'Match_Type': r.get('comparison', {}).get('match_type', 'Failed')
} for r in test_results])

csv_output = '/root/.openclaw/workspace/gbb-reaction/data/v73_validation_summary.csv'
summary_df.to_csv(csv_output, index=False)

print(f"💾 JSON 报告：{json_output}")
print(f"📋 CSV 摘要：{csv_output}")

# ============== 打印摘要 ==============
print("\n" + "=" * 80)
print("✅ v7.3 引擎测试完成!")
print("=" * 80)

success_rate = stats['success'] / stats['total'] * 100 if stats['total'] > 0 else 0
exact_rate = stats['exact_match'] / stats['total'] * 100 if stats['total'] > 0 else 0
mw_match_rate = stats['mw_match'] / stats['total'] * 100 if stats['total'] > 0 else 0

print(f"""
📊 测试摘要:
   - 测试案例：{stats['total']}
   - 成功预测：{stats['success']} ({success_rate:.1f}%)
   - 精确匹配：{stats['exact_match']} ({exact_rate:.1f}%)
   - 高相似性：{stats['high_similarity']}
   - 核心验证：{stats['core_match']}
   - 分子量匹配：{stats['mw_match']} ({mw_match_rate:.1f}%)
   - 失败：{stats['failed']}

💡 结论:
   1. v7.3 引擎成功率：{success_rate:.1f}%
   2. 结构精确匹配率：{exact_rate:.1f}%
   3. 分子量预测准确率：{mw_match_rate:.1f}%

📄 详细报告:
   - JSON: {json_output}
   - CSV: {csv_output}
""")
