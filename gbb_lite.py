#!/usr/bin/env python3
"""
GBB 反应虚拟模拟系统 (无需 RDKit 版本)
Groebke-Blackburn-Bienaymé Reaction Simulator

本脚本演示 GBB 三组分反应的原理和数据结构，
可用于教学和数据准备。完整功能需要 RDKit。
"""

import json
from dataclasses import dataclass, asdict
from typing import List, Dict, Optional

@dataclass
class Molecule:
    """分子数据结构"""
    name: str
    smiles: str
    category: str  # aldehyde, amine, isocyanide, product
    mw: Optional[float] = None
    description: str = ""

@dataclass
class ReactionResult:
    """反应结果"""
    reaction_id: str
    aldehyde: str
    amine: str
    isocyanide: str
    product_name: str
    product_smiles: str
    status: str  # success, failed, predicted
    yield_predicted: Optional[float] = None

class GBBReactionDatabase:
    """GBB 反应数据库"""
    
    def __init__(self):
        self.aldehydes = [
            Molecule("benzaldehyde", "O=Cc1ccccc1", "aldehyde", 106.12, "苯甲醛"),
            Molecule("4-methylbenzaldehyde", "Cc1ccc(C=O)cc1", "aldehyde", 120.15, "对甲基苯甲醛"),
            Molecule("4-methoxybenzaldehyde", "COc1ccc(C=O)cc1", "aldehyde", 136.15, "对甲氧基苯甲醛"),
            Molecule("4-chlorobenzaldehyde", "Clc1ccc(C=O)cc1", "aldehyde", 140.57, "对氯苯甲醛"),
            Molecule("4-nitrobenzaldehyde", "[O-][N+](=O)c1ccc(C=O)cc1", "aldehyde", 151.12, "对硝基苯甲醛"),
            Molecule("furfural", "O=Cc1occc1", "aldehyde", 96.08, "糠醛"),
            Molecule("acetaldehyde", "CC=O", "aldehyde", 44.05, "乙醛"),
            Molecule("3-pyridinecarboxaldehyde", "O=Cc1cccnc1", "aldehyde", 107.11, "烟碱醛"),
        ]
        
        self.amines = [
            Molecule("2-aminopyridine", "Nc1ccccn1", "amine", 94.11, "2-氨基吡啶"),
            Molecule("2-amino-3-methylpyridine", "Cc1ncccc1N", "amine", 108.14, "2-氨基 -3-甲基吡啶"),
            Molecule("2-amino-5-chloropyridine", "Nc1ccc(Cl)cn1", "amine", 128.56, "2-氨基 -5-氯吡啶"),
            Molecule("2-amino-6-methylpyridine", "Cc1cccc(N)n1", "amine", 108.14, "2-氨基 -6-甲基吡啶"),
            Molecule("2-amino-4-methylpyridine", "Cc1cc(N)ccn1", "amine", 108.14, "2-氨基 -4-甲基吡啶"),
            Molecule("2-amino-5-bromopyridine", "Nc1ccc(Br)cn1", "amine", 173.01, "2-氨基 -5-溴吡啶"),
        ]
        
        self.isocyanides = [
            Molecule("tert-butyl isocyanide", "CC(C)(C)[N+]#[C-]", "isocyanide", 83.13, "叔丁基异腈"),
            Molecule("ethyl isocyanide", "CC[N+]#[C-]", "isocyanide", 57.09, "乙基异腈"),
            Molecule("phenyl isocyanide", "[N+]#[C-]c1ccccc1", "isocyanide", 103.12, "苯基异腈"),
            Molecule("cyclohexyl isocyanide", "[N+]#[C-]C1CCCCC1", "isocyanide", 109.17, "环己基异腈"),
            Molecule("methyl isocyanide", "C[N+]#[C-]", "isocyanide", 43.06, "甲基异腈"),
            Molecule("4-methylphenyl isocyanide", "Cc1ccc([N+]#[C-])cc1", "isocyanide", 117.15, "对甲苯基异腈"),
        ]
    
    def get_all_combinations(self) -> List[Dict]:
        """获取所有可能的反应组合"""
        combinations = []
        for ald in self.aldehydes:
            for amine in self.amines:
                for iso in self.isocyanides:
                    combinations.append({
                        'aldehyde': ald.name,
                        'amine': amine.name,
                        'isocyanide': iso.name,
                    })
        return combinations
    
    def predict_product(self, aldehyde_name: str, amine_name: str, isocyanide_name: str) -> Optional[ReactionResult]:
        """
        预测 GBB 反应产物
        
        GBB 反应通式:
        醛 + 2-氨基吡啶 + 异腈 → 咪唑并 [1,2-a] 吡啶衍生物
        
        产物结构预测规则:
        - 咪唑并 [1,2-a] 吡啶为核心骨架
        - 醛的 R 基团连接到 C3 位
        - 异腈的 R 基团连接到 C2 位
        - 胺的取代基保留在吡啶环上
        """
        ald = next((m for m in self.aldehydes if m.name == aldehyde_name), None)
        amine = next((m for m in self.amines if m.name == amine_name), None)
        iso = next((m for m in self.isocyanides if m.name == isocyanide_name), None)
        
        if not all([ald, amine, iso]):
            return None
        
        # 生成产物名称
        product_name = f"3-{ald.name.replace('aldehyde', 'yl').replace('e', '')}-2-{iso.name.replace(' isocyanide', '')}-imidazo[1,2-a]pyridine"
        if amine.name != "2-aminopyridine":
            substituent = amine.name.replace('2-amino-', '').replace('pyridine', '')
            product_name = f"{substituent}-{product_name}"
        
        # 生成简化 SMILES (实际结构需要 RDKit 精确计算)
        # 核心骨架：咪唑并 [1,2-a] 吡啶
        # 这是一个近似表示
        core = "c1ccc2nc(C)c3"  # 简化核心
        
        # 基于反应物构建近似 SMILES
        product_smiles = self._build_approximate_smiles(ald, amine, iso)
        
        # 预测产率 (基于电子效应和位阻的简化规则)
        predicted_yield = self._predict_yield(ald, amine, iso)
        
        return ReactionResult(
            reaction_id=f"GBB-{ald.name[:3]}-{amine.name[:3]}-{iso.name[:3]}",
            aldehyde=ald.name,
            amine=amine.name,
            isocyanide=iso.name,
            product_name=product_name,
            product_smiles=product_smiles,
            status="predicted",
            yield_predicted=predicted_yield
        )
    
    def _build_approximate_smiles(self, ald: Molecule, amine: Molecule, iso: Molecule) -> str:
        """构建近似产物 SMILES"""
        # 咪唑并 [1,2-a] 吡啶核心
        # 实际应用中应使用 RDKit 进行精确的反应物到产物转换
        
        # 简化规则：保留核心，添加取代基提示
        core = "c1ccc2[nH]c(R1)c(R2)2n1"
        
        # 这里返回一个标记 SMILES，表示需要进一步处理
        return f"[GBB-product]: {ald.smiles} + {amine.smiles} + {iso.smiles}"
    
    def _predict_yield(self, ald: Molecule, amine: Molecule, iso: Molecule) -> float:
        """
        基于经验规则预测产率
        
        影响因素:
        - 醛的电子效应：吸电子基团通常提高产率
        - 位阻效应：大位阻基团降低产率
        - 异腈类型：烷基异腈通常优于芳基异腈
        """
        base_yield = 75.0  # 基础产率
        
        # 醛的影响
        if 'nitro' in ald.name.lower():
            base_yield += 10  # 吸电子基团
        elif 'chloro' in ald.name.lower() or 'bromo' in ald.name.lower():
            base_yield += 5
        elif 'methoxy' in ald.name.lower():
            base_yield -= 5  # 给电子基团
        elif ald.name == 'furfural':
            base_yield += 3  # 杂环醛
        
        # 胺的影响
        if 'methyl' in amine.name.lower():
            base_yield -= 3  # 位阻
        elif 'chloro' in amine.name.lower() or 'bromo' in amine.name.lower():
            base_yield += 3
        
        # 异腈的影响
        if 'tert-butyl' in iso.name.lower():
            base_yield -= 5  # 大位阻
        elif 'phenyl' in iso.name.lower():
            base_yield -= 3  # 芳基异腈
        
        # 限制在合理范围
        return max(30.0, min(95.0, base_yield))

def generate_test_dataset(db: GBBReactionDatabase, num_samples: int = 10) -> List[ReactionResult]:
    """生成测试数据集"""
    import random
    
    results = []
    combinations = db.get_all_combinations()
    
    # 随机选择组合
    selected = random.sample(combinations, min(num_samples, len(combinations)))
    
    for combo in selected:
        result = db.predict_product(
            combo['aldehyde'],
            combo['amine'],
            combo['isocyanide']
        )
        if result:
            results.append(result)
    
    return results

def export_to_json(results: List[ReactionResult], filename: str):
    """导出结果为 JSON"""
    data = [asdict(r) for r in results]
    with open(filename, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"✅ 数据已导出到：{filename}")

def main():
    print("=" * 70)
    print("GBB 反应虚拟模拟系统")
    print("Groebke-Blackburn-Bienaymé Reaction Simulator")
    print("=" * 70)
    
    # 初始化数据库
    db = GBBReactionDatabase()
    
    print(f"\n📦 反应物数据库:")
    print(f"  醛类 (Aldehydes): {len(db.aldehydes)} 种")
    print(f"  胺类 (Amines): {len(db.amines)} 种")
    print(f"  异腈类 (Isocyanides): {len(db.isocyanides)} 种")
    print(f"  可能组合数：{len(db.get_all_combinations())} 种")
    
    # 生成测试数据集
    print("\n🔬 生成测试反应数据集...")
    results = generate_test_dataset(db, num_samples=15)
    
    print("\n" + "=" * 70)
    print("📊 反应预测结果:")
    print("=" * 70)
    
    for i, r in enumerate(results, 1):
        print(f"\n[{i}] {r.reaction_id}")
        print(f"    反应物：{r.aldehyde} + {r.amine} + {r.isocyanide}")
        print(f"    产物：{r.product_name}")
        print(f"    预测产率：{r.yield_predicted:.1f}%")
        print(f"    状态：{r.status}")
    
    # 导出到 JSON
    export_to_json(results, '/root/.openclaw/workspace/gbb-reaction/gbb_test_data.json')
    
    # 生成摘要统计
    print("\n" + "=" * 70)
    print("📈 统计摘要:")
    avg_yield = sum(r.yield_predicted for r in results) / len(results)
    print(f"  平均预测产率：{avg_yield:.1f}%")
    print(f"  最高预测产率：{max(r.yield_predicted for r in results):.1f}%")
    print(f"  最低预测产率：{min(r.yield_predicted for r in results):.1f}%")
    
    print("\n💡 下一步建议:")
    print("  1. 安装 RDKit 进行精确反应模拟:")
    print("     conda install -c conda-forge rdkit")
    print("  2. 运行 gbb_simulation.py 进行完整模拟")
    print("  3. 参考论文：https://doi.org/10.3762/bjoc.20.143")
    
    return results

if __name__ == '__main__':
    main()
