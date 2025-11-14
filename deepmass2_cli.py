# -*- coding: utf-8 -*-

# (本脚本用于替代 GUI 版本(DeepMASS2.py)，以命令行方式运行 DeepMASS2)
# Usage:
#   python deepmass2_cli.py <input_file> <output_dir>
#
# <input_file> : spectrum file (.mgf / .msp / .mat)
# <output_dir> : folder to save result CSVs
#

import os
import pickle
import pandas as pd
import numpy as np
from gensim.models import Word2Vec
import hnswlib

from core.importing.load_from_files import load_from_files
from core.main import identify_unknown, match_spectrum


# ==============================
# 加载基础资源（索引 + reference + word2vec）
# ==============================
def load_reference_index(index_path, dim):
    p = hnswlib.Index(space="l2", dim=dim)
    p.load_index(index_path)
    return p


def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)


# ==============================
# 主识别函数
# ==============================
def run_deepmass_identification(input_file, output_dir):
    # ---- 路径设置 ----
    index_pos = "data/references_index_positive_spec2vec.bin"
    index_neg = "data/references_index_negative_spec2vec.bin"
    ref_pos = "data/references_spectrums_positive.pickle"
    ref_neg = "data/references_spectrums_negative.pickle"
    database_file = "data/DeepMassStructureDB-v1.1.csv"

    # ---- load 数据库 ----
    database = pd.read_csv(database_file)

    # ---- 加载 hnsw 索引 ----
    p_pos = load_reference_index(index_pos, dim=300)
    p_neg = load_reference_index(index_neg, dim=300)

    # ---- 加载参考光谱 ----
    reference_positive = load_pickle(ref_pos)
    reference_negative = load_pickle(ref_neg)

    # ---- 加载 word2vec 模型 ----
    model_pos = Word2Vec.load("model/Ms2Vec_allGNPSpositive.hdf5")
    model_neg = Word2Vec.load("model/Ms2Vec_allGNPSnegative.hdf5")

    # ---- 加载未知光谱 ----
    spectrums = load_from_files([input_file])

    # ---- 输出目录 ----
    os.makedirs(output_dir, exist_ok=True)

    # ---- 依次识别 ----
    for i, s in enumerate(spectrums):
        print(f"[{i+1}/{len(spectrums)}] Identifying: {s.metadata.get('compound_name',f'spectrum_{i+1}')}")

        # 正负模式判断
        if s.metadata.get("ionmode") == "negative":
            sn = identify_unknown(
                s, p_neg, model_neg, reference_negative, database
            )
        else:
            sn = identify_unknown(
                s, p_pos, model_pos, reference_positive, database
            )

        # 生成输出文件名
        name = s.metadata.get("compound_name", f"spectrum_{i+1}")
        name = "".join(c for c in name if c.isalnum() or c == "_")

        # 识别表
        if sn.metadata.get("annotation") is not None:
            df = sn.metadata["annotation"]
        else:
            df = pd.DataFrame(columns=["Title", "MolecularFormula", "CanonicalSMILES", "InChIKey"])

        save_path = os.path.join(output_dir, f"{name}.csv")
        df.to_csv(save_path, index=True)

        print(f"  -> saved: {save_path}")

    print("\n全部识别完成！")


# ==============================
# main
# ==============================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="DeepMASS2 CLI")
    parser.add_argument("input_file", type=str, help="输入光谱文件")
    parser.add_argument("output_dir", type=str, help="输出 CSV 保存目录")

    args = parser.parse_args()

    run_deepmass_identification(args.input_file, args.output_dir)
