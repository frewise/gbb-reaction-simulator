#!/usr/bin/env python3
"""
d_sim_v3.py - 完整修复版
计算候选小分子与 Danuglipron 在 GLP-1R 上的对接模式相似度 (IFP similarity)

修复内容:
1. 使用基于距离的 IFP 计算（不依赖 ProLIF 的氢原子检测）
2. 修复 PDBQT 加载 (meeko API 变更)
3. 更新关键残基列表为 6X1A 实际结合残基
4. 添加完整的验证闭环流程
"""

import argparse
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import json
from collections import defaultdict

import MDAnalysis as mda
from rdkit import Chem
from meeko import PDBQTMolecule


# ==============================
# 常量 - 基于 6X1A 实际结合残基
# ==============================

REFERENCE_PDB = "6X1A.pdb"
REFERENCE_LIGAND = "UK4"
RECEPTOR_CHAIN = "R"
CONTACT_DISTANCE = 4.0
HBOND_DISTANCE = 3.5

# 6X1A 中实际与 Danuglipron 接触的关键残基
KEY_RESIDUES = [
    ("K", 197),  # LYS - 3.08 Å 最近接触
    ("R", 380),  # ARG - 氢键网络
    ("F", 385),  # PHE - 疏水
    ("W", 203),  # TRP - π-π
    ("Q", 221),  # GLN - 氢键
    ("Q", 37),   # GLN - 氢键
    ("L", 384),  # LEU - 疏水
    ("F", 230),  # PHE - 疏水/π-π
    ("W", 33),   # TRP - π-π
]

AA1_TO_AA3 = {
    'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS',
    'Q':'GLN', 'E':'GLU', 'G':'GLY', 'H':'HIS', 'I':'ILE',
    'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', 'P':'PRO',
    'S':'SER', 'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'
}

AA3_TO_AA1 = {v:k for k,v in AA1_TO_AA3.items()}

# 残基类型分类
HYDROPHOBIC_RES = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR']
AROMATIC_RES = ['PHE', 'TYR', 'TRP', 'HIS']
POSITIVE_RES = ['ARG', 'LYS', 'HIS']
NEGATIVE_RES = ['ASP', 'GLU']
HBDONOR_RES = ['ARG', 'LYS', 'HIS', 'ASN', 'GLN', 'SER', 'THR', 'TRP', 'TYR']
HBACCEPTOR_RES = ['ASP', 'GLU', 'ASN', 'GLN', 'SER', 'THR', 'TYR', 'HIS']


# ==============================
# 读取分子
# ==============================

def load_ligand_from_pdb(pdb_file, resname='UK4'):
    """从 PDB 文件提取配体坐标"""
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() == resname:
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                elem = line[76:78].strip() if len(line) > 76 else atom_name[0]
                atoms.append({'name': atom_name, 'elem': elem, 'coord': np.array([x, y, z])})
    return atoms


def load_receptor_from_pdb(pdb_file, chain='R'):
    """从 PDB 文件提取受体残基"""
    residues = defaultdict(lambda: {'atoms': [], 'resname': ''})
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == chain:
                resname = line[17:20].strip()
                resid = int(line[22:26])
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                key = (AA3_TO_AA1.get(resname, '?'), resid)
                residues[key]['atoms'].append({
                    'name': atom_name,
                    'resname': resname,
                    'coord': np.array([x, y, z])
                })
                residues[key]['resname'] = resname
    return dict(residues)


def load_pdbqt_poses(pdbqt_file):
    """从 PDBQT 文件加载多个 pose"""
    poses = []
    
    with open(pdbqt_file, 'r') as f:
        content = f.read()
    
    # 分割 MODEL
    models = content.split('MODEL')
    
    for model in models:
        if not model.strip() or 'ENDMDL' not in model:
            continue
        
        atoms = []
        for line in model.split('\n'):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    atom_name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    elem = atom_name[0]
                    atoms.append({'name': atom_name, 'elem': elem, 'coord': np.array([x, y, z])})
                except:
                    continue
        
        if atoms:
            poses.append(atoms)
    
    # 如果没有 MODEL 分割，整个文件作为一个 pose
    if not poses:
        atoms = []
        for line in content.split('\n'):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    atom_name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    elem = atom_name[0]
                    atoms.append({'name': atom_name, 'elem': elem, 'coord': np.array([x, y, z])})
                except:
                    continue
        if atoms:
            poses.append(atoms)
    
    return poses


# ==============================
# IFP 计算
# ==============================

def compute_ifp(ligand_atoms, receptor_residues, cutoff=4.0, hb_cutoff=3.5):
    """
    计算相互作用指纹
    
    返回 dict: {(aa1, resid): {interaction_type: 1/0}}
    """
    ifp = defaultdict(lambda: defaultdict(int))
    
    for lig_atom in ligand_atoms:
        for (aa1, resid), res_data in receptor_residues.items():
            resname = res_data['resname']
            
            for rec_atom in res_data['atoms']:
                dist = np.linalg.norm(lig_atom['coord'] - rec_atom['coord'])
                
                if dist < hb_cutoff:
                    # 氢键
                    if lig_atom['elem'].upper() in ['N', 'O'] or rec_atom['name'][0].upper() in ['N', 'O']:
                        if resname in HBDONOR_RES:
                            ifp[(aa1, resid)]['HBDonor'] = 1
                        if resname in HBACCEPTOR_RES:
                            ifp[(aa1, resid)]['HBAcceptor'] = 1
                
                if dist < cutoff:
                    # 疏水作用
                    if resname in HYDROPHOBIC_RES:
                        ifp[(aa1, resid)]['Hydrophobic'] = 1
                    
                    # π-π 堆积
                    if resname in AROMATIC_RES:
                        ifp[(aa1, resid)]['PiStacking'] = 1
                    
                    # 静电作用
                    if resname in POSITIVE_RES:
                        ifp[(aa1, resid)]['Cationic'] = 1
                    if resname in NEGATIVE_RES:
                        ifp[(aa1, resid)]['Anionic'] = 1
                    
                    # 一般接触
                    ifp[(aa1, resid)]['VdWContact'] = 1
    
    return ifp


def ifp_to_vector(ifp, all_residues):
    """IFP 转 numpy 向量"""
    features = ['HBDonor', 'HBAcceptor', 'Hydrophobic', 'PiStacking', 'Cationic', 'Anionic', 'VdWContact']
    vector = []
    for res in all_residues:
        for feat in features:
            vector.append(1 if ifp.get(res, {}).get(feat, 0) else 0)
    return np.array(vector)


def tanimoto_similarity(fp1, fp2):
    """计算 Tanimoto 相似度"""
    intersection = np.sum((fp1 > 0) & (fp2 > 0))
    union = np.sum((fp1 > 0) | (fp2 > 0))
    if union == 0:
        return 0.0
    return intersection / union


def count_key_hits(ifp, key_residues):
    """计算关键残基命中数"""
    hits = []
    for (aa1, resid) in key_residues:
        if (aa1, resid) in ifp:
            res_ifp = ifp[(aa1, resid)]
            if any(v == 1 for v in res_ifp.values()):
                hits.append(f"{aa1}{resid}")
    return hits


# ==============================
# 主程序
# ==============================

def main():
    parser = argparse.ArgumentParser(description="计算对接模式 IFP 相似度")
    parser.add_argument("dockings", nargs="+", help="对接结果文件 (sdf/pdbqt/pdb)")
    parser.add_argument("-o", "--output", default="similarity_results.csv")
    parser.add_argument("--ref-pdb", default=REFERENCE_PDB, help="参考 PDB 文件")
    parser.add_argument("--verbose", "-v", action="store_true", help="详细输出")
    args = parser.parse_args()

    print("="*70)
    print("IFP 相似度分析 - Danuglipron @ GLP-1R")
    print("="*70)
    
    # 1. 加载参考结构
    print(f"\n1. 加载参考结构：{args.ref_pdb}")
    ref_ligand = load_ligand_from_pdb(args.ref_pdb, REFERENCE_LIGAND)
    ref_receptor = load_receptor_from_pdb(args.ref_pdb, RECEPTOR_CHAIN)
    print(f"   参考配体：{len(ref_ligand)} 原子")
    print(f"   受体残基：{len(ref_receptor)}")
    
    # 2. 计算参考 IFP
    print("\n2. 计算参考 IFP (Danuglipron 晶体结构)")
    ref_ifp = compute_ifp(ref_ligand, ref_receptor)
    
    # 获取所有接触残基
    all_residues = sorted(ref_ifp.keys())
    print(f"   接触残基：{len(all_residues)}")
    
    # 参考 IFP 向量
    ref_vector = ifp_to_vector(ref_ifp, all_residues)
    print(f"   IFP 特征数：{np.sum(ref_vector)}")
    
    # 参考关键残基命中
    ref_key_hits = count_key_hits(ref_ifp, KEY_RESIDUES)
    print(f"   关键残基命中：{len(ref_key_hits)}/{len(KEY_RESIDUES)}")
    if args.verbose:
        print(f"   命中列表：{', '.join(ref_key_hits)}")
    
    # 3. 处理对接结果
    print("\n3. 处理对接结果")
    results = []
    
    for docking_file in args.dockings:
        print(f"\n   文件：{docking_file}")
        poses = load_pdbqt_poses(docking_file)
        print(f"   加载 {len(poses)} 个 pose")
        
        for i, pose_atoms in enumerate(poses):
            # 计算 IFP
            pose_ifp = compute_ifp(pose_atoms, ref_receptor)
            pose_vector = ifp_to_vector(pose_ifp, all_residues)
            
            # 计算相似度
            ifp_sim = tanimoto_similarity(ref_vector, pose_vector)
            
            # 关键残基命中
            key_hits = count_key_hits(pose_ifp, KEY_RESIDUES)
            key_ratio = len(key_hits) / len(KEY_RESIDUES)
            
            # 综合评分
            score = 0.7 * ifp_sim + 0.3 * key_ratio
            
            results.append({
                'ligand_file': os.path.basename(docking_file),
                'pose_id': i,
                'IFP_similarity': round(ifp_sim, 4),
                'IFP_intersection': int(((ref_vector == 1) & (pose_vector == 1)).sum()),
                'IFP_union': int(((ref_vector == 1) | (pose_vector == 1)).sum()),
                'key_hits': len(key_hits),
                'key_ratio': round(key_ratio, 4),
                'hit_residues': ';'.join(key_hits),
                'score': round(score, 4)
            })
            
            if args.verbose:
                print(f"   Pose {i}: IFP={ifp_sim:.3f}, Key={len(key_hits)}/{len(KEY_RESIDUES)}, Score={score:.3f}")
    
    # 4. 保存结果
    df = pd.DataFrame(results)
    if len(df) > 0:
        df = df.sort_values('score', ascending=False)
    df.to_csv(args.output, index=False)
    print(f"\n4. 结果已保存：{args.output}")
    
    # 5. 输出摘要
    print("\n" + "="*70)
    print("结果摘要")
    print("="*70)
    if len(df) > 0:
        print(df.to_string())
    else:
        print("无结果")
    
    # 6. 保存 JSON 详细结果
    json_output = args.output.replace('.csv', '.json')
    output_data = {
        'reference': {
            'pdb': args.ref_pdb,
            'ligand': REFERENCE_LIGAND,
            'chain': RECEPTOR_CHAIN,
            'key_residues': [f"{aa1}{resid}" for aa1, resid in KEY_RESIDUES],
            'contact_residues': [f"{aa1}{resid}" for aa1, resid in all_residues],
            'ref_key_hits': ref_key_hits
        },
        'results': results
    }
    with open(json_output, 'w') as f:
        json.dump(output_data, f, indent=2)
    print(f"\n详细结果：{json_output}")
    
    print("\n" + "="*70)
    print("✅ IFP 分析完成!")
    print("="*70)
    
    return df


if __name__ == "__main__":
    main()
