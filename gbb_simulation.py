#!/usr/bin/env python3
"""
GBB Reaction Virtual Simulation
Groebke-Blackburn-Bienaymé Reaction: Aldehyde + Amine + Isocyanide → Imidazo[1,2-a]pyridine

This script simulates the GBB multicomponent reaction using RDKit.
"""

import sys
from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions, Descriptors
from rdkit.Chem.Draw import MolToImage

# GBB 反应通式
# 醛 + 胺 (2-氨基吡啶) + 异腈 → 咪唑并 [1,2-a] 吡啶

def create_gbb_reaction():
    """
    创建 GBB 反应的反应式
    使用 SMARTS 模式定义三组分反应
    """
    # GBB 反应 SMARTS 模式
    # 醛：[C:1](=[O:2])[H] 
    # 2-氨基吡啶：[n:3]1ccccc1[N:4]
    # 异腈：[C:5]#[N:6]
    
    # 简化的反应模式 (需要 2-氨基吡啶衍生物)
    reaction_smarts = """
    [C:1](=[O:2])[H].[n:3]1c([N:4])cccc1.[C:5]#[N:6]>>
    [n:3]1c2c(cccc1)[nH][c:4](=[N+:5][C-:6])[c:1]2[O-:2]
    """
    
    # 更实用的版本 - 生成咪唑并吡啶核心
    gbb_smarts = """
    [C:1](=[O:2])[H].[nH:3]1c(N)cccc1.[C:4]#[N:5]>>
    [nH:3]1c2c(cccc1)nc([C:4])[c:1]2
    """
    
    try:
        rxn = rdChemReactions.ReactionFromSmarts(gbb_smarts)
        return rxn
    except Exception as e:
        print(f"创建反应式失败：{e}")
        return None

def get_test_reactants():
    """
    获取测试反应物数据集
    包含常见的 GBB 反应底物
    """
    reactants_db = {
        'aldehydes': [
            ('benzaldehyde', 'O=Cc1ccccc1'),
            ('4-methylbenzaldehyde', 'O=Cc1ccc(C)cc1'),
            ('4-methoxybenzaldehyde', 'O=Cc1ccc(OC)cc1'),
            ('4-chlorobenzaldehyde', 'O=Cc1ccc(Cl)cc1'),
            ('furfural', 'O=Cc1occc1'),
            ('acetaldehyde', 'CC=O'),
        ],
        'amines': [
            ('2-aminopyridine', 'Nc1ccccn1'),
            ('2-amino-3-methylpyridine', 'Cc1ncccc1N'),
            ('2-amino-5-chloropyridine', 'Nc1ccc(Cl)cn1'),
            ('2-amino-6-methylpyridine', 'Cc1cccc(N)n1'),
        ],
        'isocyanides': [
            ('tert-butyl isocyanide', 'CC(C)(C)[N+]#[C-]'),
            ('ethyl isocyanide', 'CC[N+]#[C-]'),
            ('phenyl isocyanide', 'c1ccccc1[N+]#[C-]'),
            ('cyclohexyl isocyanide', 'C1CCCCC1[N+]#[C-]'),
        ]
    }
    return reactants_db

def validate_molecule(smiles):
    """验证 SMILES 字符串"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            Chem.SanitizeMol(mol)
            return True, mol
        return False, None
    except:
        return False, None

def run_gbb_reaction(aldehyde_smiles, amine_smiles, isocyanide_smiles):
    """
    运行 GBB 反应模拟
    
    Args:
        aldehyde_smiles: 醛的 SMILES
        amine_smiles: 胺的 SMILES  
        isocyanide_smiles: 异腈的 SMILES
    
    Returns:
        产物列表
    """
    # 验证反应物
    valid_aldehyde, aldehyde_mol = validate_molecule(aldehyde_smiles)
    valid_amine, amine_mol = validate_molecule(amine_smiles)
    valid_isocyanide, isocyanide_mol = validate_molecule(isocyanide_smiles)
    
    if not all([valid_aldehyde, valid_amine, valid_isocyanide]):
        print("反应物验证失败")
        return None
    
    # 创建反应
    rxn = create_gbb_reaction()
    if not rxn:
        print("无法创建反应式")
        return None
    
    # 运行反应
    try:
        products = rxn.RunReactants((aldehyde_mol, amine_mol, isocyanide_mol))
        if products:
            return [Chem.MolToSmiles(p[0]) for p in products]
        return []
    except Exception as e:
        print(f"反应运行失败：{e}")
        return None

def generate_product_smiles(aldehyde_name, amine_name, isocyanide_name):
    """
    基于反应物名称生成预测产物 SMILES
    使用简化的连接规则
    """
    # 这是一个简化的预测，实际 GBB 反应更复杂
    # 产物核心：咪唑并 [1,2-a] 吡啶
    
    # 从反应物构建产物（简化版本）
    product_templates = {
        '2-aminopyridine': 'n1c2ccccc2nc1',  # 咪唑并吡啶核心
    }
    
    return product_templates.get(amine_name, None)

def main():
    print("=" * 60)
    print("GBB 反应虚拟模拟系统")
    print("Groebke-Blackburn-Bienaymé Reaction Simulator")
    print("=" * 60)
    
    # 获取测试数据
    reactants = get_test_reactants()
    
    print("\n📦 可用反应物数据库:")
    print(f"  醛类 (Aldehydes): {len(reactants['aldehydes'])} 种")
    print(f"  胺类 (Amines): {len(reactants['amines'])} 种")
    print(f"  异腈类 (Isocyanides): {len(reactants['isocyanides'])} 种")
    
    print("\n📋 测试反应组合:")
    print("-" * 60)
    
    # 测试几个组合
    test_cases = [
        ('benzaldehyde', '2-aminopyridine', 'tert-butyl isocyanide'),
        ('4-methylbenzaldehyde', '2-aminopyridine', 'ethyl isocyanide'),
        ('benzaldehyde', '2-amino-3-methylpyridine', 'phenyl isocyanide'),
    ]
    
    results = []
    for i, (ald, amine, iso) in enumerate(test_cases, 1):
        ald_smiles = next((s for n, s in reactants['aldehydes'] if n == ald), None)
        amine_smiles = next((s for n, s in reactants['amines'] if n == amine), None)
        iso_smiles = next((s for n, s in reactants['isocyanides'] if n == iso), None)
        
        print(f"\n测试 {i}: {ald} + {amine} + {iso}")
        print(f"  SMILES: {ald_smiles} + {amine_smiles} + {iso_smiles}")
        
        # 运行反应
        products = run_gbb_reaction(ald_smiles, amine_smiles, iso_smiles)
        
        if products:
            print(f"  ✅ 产物：{products[0][:50]}...")
            results.append({'reactants': (ald, amine, iso), 'product': products[0]})
        else:
            print(f"  ⚠️  无产物生成 (反应模式可能需要调整)")
            results.append({'reactants': (ald, amine, iso), 'product': None})
    
    print("\n" + "=" * 60)
    print("📊 模拟结果摘要:")
    print(f"  总测试数：{len(test_cases)}")
    print(f"  成功生成：{sum(1 for r in results if r['product'])}")
    print(f"  失败：{sum(1 for r in results if not r['product'])}")
    
    print("\n💡 说明:")
    print("  此模拟使用简化的反应模式。")
    print("  实际 GBB 反应需要更复杂的 SMARTS 模式和立体化学处理。")
    print("  推荐使用 RDKit 完整功能进行精确模拟。")
    
    return results

if __name__ == '__main__':
    try:
        main()
    except ImportError as e:
        print(f"❌ 错误：缺少依赖 - {e}")
        print("\n请安装 RDKit:")
        print("  pip install rdkit-pypi  # Linux/macOS")
        print("  conda install -c conda-forge rdkit  # 推荐")
