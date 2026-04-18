import pandas as pd
import re
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio import Align
import itertools
import sys
import os

CODON_TABLES = {
    'human': {'A': 'GCC', 'C': 'TGC', 'D': 'GAC', 'E': 'GAG', 'F': 'TTC', 'G': 'GGC', 'H': 'CAC', 'I': 'ATC', 'K': 'AAG', 'L': 'CTG', 'M': 'ATG', 'N': 'AAC', 'P': 'CCC', 'Q': 'CAG', 'R': 'CGC', 'S': 'AGC', 'T': 'ACC', 'V': 'GTG', 'W': 'TGG', 'Y': 'TAC', '*': 'TGA'},
    'ecoli': {'A': 'GCG', 'C': 'TGC', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGC', 'H': 'CAT', 'I': 'ATT', 'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAC', 'P': 'CCG', 'Q': 'CAG', 'R': 'CGC', 'S': 'AGC', 'T': 'ACC', 'V': 'GTG', 'W': 'TGG', 'Y': 'TAT', '*': 'TAA'},
    'yeast': {'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'K': 'AAG', 'L': 'TTG', 'M': 'ATG', 'N': 'AAC', 'P': 'CCA', 'Q': 'CAA', 'R': 'AGA', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT', '*': 'TAA'}
}

def parse_fasta(fasta_file):
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    return "".join([line.strip() for line in lines if not line.startswith(">")]).upper()

def get_all_candidates(full_seq, start_idx, direction):
    candidates = []
    for seq_len in range(15, 61):
        if direction == 1:
            end_idx = start_idx + seq_len
            if end_idx > len(full_seq): break
            binding = full_seq[start_idx: end_idx]
            has_gc_clamp = binding[-1].upper() in ['C', 'G']
        else:
            slice_start = start_idx - seq_len
            if slice_start < 0: break
            binding = full_seq[slice_start: start_idx]
            has_gc_clamp = binding[0].upper() in ['C', 'G']

        tm = mt.Tm_NN(binding, nn_table=mt.DNA_NN3, saltcorr=5, dnac1=250, dnac2=250, Na=50)
        candidates.append({'seq': binding, 'tm': tm, 'gc': has_gc_clamp, 'len': seq_len})
    return candidates

def get_best_primer_pair(fp_candidates, rp_candidates, min_overhang_len=9):
    fp_gc = [c for c in fp_candidates if c['gc']]
    rp_gc = [c for c in rp_candidates if c['gc']]
    fp_warn, rp_warn = "", ""

    if not fp_gc:
        if fp_candidates: fp_gc, fp_warn = [fp_candidates[len(fp_candidates) // 2]], "CRITICAL: No GC clamp for FP"
        else: return None
    if not rp_gc:
        if rp_candidates: rp_gc, rp_warn = [rp_candidates[len(rp_candidates) // 2]], "CRITICAL: No GC clamp for RP"
        else: return None

    pairs = []
    for f in fp_gc:
        for r in rp_gc:
            dtm = abs(f['tm'] - r['tm'])
            f_bounds = 24 <= f['len'] <= 36 and 58.0 <= f['tm'] <= 68.0
            r_bounds = 24 <= r['len'] <= 36 and 58.0 <= r['tm'] <= 68.0
            len_valid = (f['len'] + min_overhang_len <= 59) and (r['len'] + min_overhang_len <= 59)
            pairs.append({'fp': f, 'rp': r, 'dtm': dtm, 'dtm_valid': dtm <= 3.0, 'bounds_valid': f_bounds and r_bounds, 'len_valid': len_valid})

    def score(p):
        pen = sum([(24 - c['len']) if c['len'] < 24 else (c['len'] - 36) if c['len'] > 36 else 0 for c in [p['fp'], p['rp']]])
        pen += sum([(58.0 - c['tm']) if c['tm'] < 58.0 else (c['tm'] - 68.0) if c['tm'] > 68.0 else 0 for c in [p['fp'], p['rp']]])
        return pen

    # C0 最高优先级：优先过滤出预估总长度 <= 59bp 的候选组合
    valid_pool = [p for p in pairs if p['len_valid']]
    if not valid_pool:
        valid_pool = pairs  # 降级保底：如果由于 GC 极低等原因找不到 <= 59bp 的引物，退而求其次
        fp_warn = fp_warn + " | Total len > 59bp unavoidable" if fp_warn else "Total len > 59bp unavoidable"

    tier1 = [p for p in valid_pool if p['dtm_valid'] and p['bounds_valid']]
    if tier1: return min(tier1, key=lambda x: x['fp']['len'] + x['rp']['len']), fp_warn, rp_warn

    tier2 = [p for p in valid_pool if p['dtm_valid'] and not p['bounds_valid']]
    if tier2: return min(tier2, key=score), fp_warn, rp_warn

    tier3 = [p for p in valid_pool if not p['dtm_valid'] and p['bounds_valid']]
    if tier3: return min(tier3, key=lambda x: x['dtm']), fp_warn, rp_warn

    return min(valid_pool, key=lambda x: (x['dtm'], score(x))), fp_warn, rp_warn

def extract_mutation_info(variant, species):
    match = re.match(r'^([a-zA-Z])(\d+)([a-zA-Z\*])$', str(variant).strip())
    if not match: return None
    wt_aa, pos, mut_aa = match.group(1).upper(), int(match.group(2)), match.group(3).upper()
    return pos, CODON_TABLES[species.lower()][mut_aa]

def get_primer_footprint(fp_seq, rp_seq, full_seq):
    """通过局部序列比对，计算外部输入的引物对在野生型模板上的绝对覆盖范围(Footprint)"""
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    # 设置打分机制，容忍由突变带来的少量 Mismatch，严厉惩罚 Gap
    aligner.match_score = 1
    aligner.mismatch_score = -2
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    try:
        fp_aln = aligner.align(full_seq, fp_seq)[0]
        fp_start = fp_aln.aligned[0][0][0]
        fp_end = fp_aln.aligned[0][-1][-1]

        rp_rc = str(Seq(rp_seq).reverse_complement())
        rp_aln = aligner.align(full_seq, rp_rc)[0]
        rp_start = rp_aln.aligned[0][0][0]
        rp_end = rp_aln.aligned[0][-1][-1]
        return min(fp_start, rp_start), max(fp_end, rp_end)
    except Exception:
        return None, None

def design_core_primer(mut_seq_str, block_start, block_end, variant_name):
    """提取出的核心引物设计逻辑，支持单突变（3bp block）和相近双突变（长block）"""
    fp_cands = get_all_candidates(mut_seq_str, block_end, direction=1)
    rp_cands = get_all_candidates(mut_seq_str, block_start, direction=-1)

    # 提前计算引物悬垂区(Overhang)最小长度，用于评估 C0 最高优先级 (总长度 <= 59bp)
    # 单条引物分担的悬垂约等于 (突变跨度 + 最小Overlap(18)) / 2
    min_overhang_len = (block_end - block_start + 18) // 2

    pair_result = get_best_primer_pair(fp_cands, rp_cands, min_overhang_len)
    if not pair_result: return None

    best_pair, fp_warn, rp_warn = pair_result
    fp_bind_len, fp_tm = best_pair['fp']['len'], best_pair['fp']['tm']
    rp_bind_len, rp_tm = best_pair['rp']['len'], best_pair['rp']['tm']
    
    target_max_overlap_tm = min(fp_tm, rp_tm) - 5.0
    best_overlap, fallback_overlap, min_fallback_tm = None, None, float('inf')

    # 严格限制 Overlap 区域在 18-22bp，绝对不随突变跨度拉长而增加
    search_lengths = [18, 19, 20, 21, 22]
    all_valid_cands = []

    for L in search_lengths:
        # Overlap可以位于 [block_start - L, block_end] 之间的任何位置
        for o_start in range(block_start - L, block_end + 1):
            o_end = o_start + L
            if o_start < 0 or o_end > len(mut_seq_str): continue

            o_seq = mut_seq_str[o_start: o_end]
            tm = mt.Tm_NN(o_seq, nn_table=mt.DNA_NN3, saltcorr=5, dnac1=250, dnac2=250, Na=50)
            
            # 计算如果选择该 overlap 导致的两条引物分别的预计总长度
            fp_len_est = (block_end - o_start) + fp_bind_len
            rp_len_est = (o_end - block_start) + rp_bind_len
            max_len = max(fp_len_est, rp_len_est)

            cand = {'seq': o_seq, 'tm': tm, 'start': o_start, 'end': o_end, 'len': L, 'max_len': max_len}
            
            if tm < min_fallback_tm:
                min_fallback_tm = tm
                fallback_overlap = cand
                
            if tm <= target_max_overlap_tm:
                all_valid_cands.append(cand)

    if all_valid_cands:
        # 优先选择能让单条引物最大长度最小化的 Overlap 位置（即尽量居中），其次选择 Tm 最低的
        best_overlap = min(all_valid_cands, key=lambda x: (x['max_len'], x['tm']))

    warnings = [w for w in [fp_warn, rp_warn] if w]
    if not best_pair['dtm_valid']: warnings.append(f"Pair dTm > 3C")
    if not best_overlap:
        best_overlap = fallback_overlap
        if best_overlap: warnings.append(f"Overlap Tm NOT 5C lower than Binding Tms")
        else: return None # 彻底失败

    fp_seq = mut_seq_str[best_overlap['start']: block_end + fp_bind_len]
    rp_seq = str(Seq(mut_seq_str[block_start - rp_bind_len: best_overlap['end']]).reverse_complement())

    if len(fp_seq) > 59: warnings.append(f"LONG PRIMER: FP ({len(fp_seq)}bp) exceeds 59bp limit")
    if len(rp_seq) > 59: warnings.append(f"LONG PRIMER: RP ({len(rp_seq)}bp) exceeds 59bp limit")
    if len(fp_seq) > 90 or len(rp_seq) > 90: warnings.append("CRITICAL: Primer > 90bp, highly recommend PAGE purification or multi-fragment assembly")

    # 格式化引物名称，双突变默认用'+'连接，这里替换为'-'以符合命名规范
    primer_base_name = variant_name.replace('+', '-')

    return {
        'Variant': variant_name,
        'Forward_Name': f"{primer_base_name}-F",
        'Forward_Primer (5->3)': fp_seq, 'FP_Length': len(fp_seq), 'FP_Binding_Tm': round(fp_tm, 1),
        'Reverse_Name': f"{primer_base_name}-R",
        'Reverse_Primer (5->3)': rp_seq, 'RP_Length': len(rp_seq), 'RP_Binding_Tm': round(rp_tm, 1),
        'Pair_Avg_Tm': round((fp_tm + rp_tm) / 2, 1),
        'Pair_dTm': round(best_pair['dtm'], 1),
        'Overlap_Len': best_overlap['len'], 'Overlap_Tm': round(best_overlap['tm'], 1),
        'Warnings': " | ".join(warnings) if warnings else "None",
        'Action': 'Redesigned'
    }

def design_single_mutations(fasta_path, csv_path, species, overhang_len, output_path):
    full_seq = parse_fasta(fasta_path)
    mutations_df = pd.read_csv(csv_path, header=None, names=['Variant'])
    results = []

    for variant in mutations_df['Variant']:
        info = extract_mutation_info(variant, species)
        if not info: continue
        pos, optimal_mut_codon = info
        
        codon_start = overhang_len + (pos - 1) * 3
        codon_end = codon_start + 3
        mut_seq_str = full_seq[:codon_start] + optimal_mut_codon + full_seq[codon_end:]

        res = design_core_primer(mut_seq_str, codon_start, codon_end, variant)
        if res: results.append(res)

    pd.DataFrame(results).to_csv(output_path, index=False)
    print(f"✅ 成功生成 {len(results)} 对单突变引物，已保存至: {output_path}")

def design_double_mutations(fasta_path, csv_path, species, overhang_len, output_path):
    full_seq = parse_fasta(fasta_path)
    # 假设第一轮的结果包含Variant列，我们读取它来生成组合
    mutations_df = pd.read_csv(csv_path)
    if 'Variant' not in mutations_df.columns:
        print("错误: 输入的CSV文件必须包含 'Variant' 列。")
        return
        
    footprints = {}
    if 'Footprint_Start' in mutations_df.columns and 'Footprint_End' in mutations_df.columns:
        # 如果有缓存的 Footprint，直接使用
        footprints = {row['Variant']: (row['Footprint_Start'], row['Footprint_End']) for _, row in mutations_df.iterrows()}
    else:
        # 退级策略：如果没有 Footprint 列，但含有外部引物序列，则执行动态比对
        if 'Forward_Primer (5->3)' not in mutations_df.columns or 'Reverse_Primer (5->3)' not in mutations_df.columns:
            print("错误: 输入文件缺少 'Footprint' 数据，且没有完整的引物序列列('Forward_Primer (5->3)' 和 'Reverse_Primer (5->3)') 无法进行比对定位。")
            return
        
        print("⚠️ 发现未包含覆盖范围(Footprint)的输入文件。正在通过序列比对动态定位引物位置...")
        for idx, row in mutations_df.iterrows():
            variant = row['Variant']
            fp_seq = str(row['Forward_Primer (5->3)']).strip()
            rp_seq = str(row['Reverse_Primer (5->3)']).strip()
            
            if not fp_seq or not rp_seq or fp_seq == 'nan' or rp_seq == 'nan':
                print(f"⚠️ 警告: 变体 {variant} 缺少引物序列，跳过比对。")
                continue
                
            start, end = get_primer_footprint(fp_seq, rp_seq, full_seq)
            if start is not None and end is not None:
                footprints[variant] = (start, end)
            else:
                print(f"⚠️ 警告: 变体 {variant} 的引物在模板中匹配失败。")

    # 仅从成功定位到 Footprint 的列表中生成组合
    variants = [v for v in mutations_df['Variant'].tolist() if v in footprints]
    pairs = list(itertools.combinations(variants, 2))
    print(f"🔍 发现 {len(variants)} 个单突变，准备设计 {len(pairs)} 个双突变组合...")

    results = []

    for varA, varB in pairs:
        infoA = extract_mutation_info(varA, species)
        infoB = extract_mutation_info(varB, species)
        if not infoA or not infoB: continue

        posA, codonA = infoA
        posB, codonB = infoB

        # 计算在序列中的绝对坐标
        startA = overhang_len + (posA - 1) * 3
        startB = overhang_len + (posB - 1) * 3
        
        distance = abs(startA - startB)
        combined_name = f"{varA}+{varB}"

        fpA_start, fpA_end = footprints[varA]
        fpB_start, fpB_end = footprints[varB]

        # 动态碰撞判定：A突变是否在B引物范围内，或B突变是否在A引物范围内
        collision_A_in_B = (startA < fpB_end) and (startA + 3 > fpB_start)
        collision_B_in_A = (startB < fpA_end) and (startB + 3 > fpA_start)

        if not (collision_A_in_B or collision_B_in_A):
            # 互相不在对方引物范围内，无碰撞
            results.append({
                'Variant': combined_name,
                'Forward_Name': 'N/A', 'Reverse_Name': 'N/A',
                'Forward_Primer (5->3)': 'N/A', 'Reverse_Primer (5->3)': 'N/A',
                'Pair_Avg_Tm': 'N/A',
                'Warnings': f"No footprint collision (Dist: {distance}bp).",
                'Action': 'Use Single Primers Independently'
            })
        else:
            # 碰撞！需要将两个突变融合在同一个序列中并重新设计引物
            block_start = min(startA, startB)
            block_end = max(startA, startB) + 3
            
            # 生成包含双突变的序列
            mut_seq_list = list(full_seq)
            mut_seq_list[startA:startA+3] = list(codonA)
            mut_seq_list[startB:startB+3] = list(codonB)
            mut_seq_str = "".join(mut_seq_list)

            res = design_core_primer(mut_seq_str, block_start, block_end, combined_name)
            if res:
                res['Warnings'] += f" | Designed as single block (Dist: {distance}bp)"
                results.append(res)
            else:
                results.append({
                    'Variant': combined_name, 'Warnings': "FAILED to design bridging primer", 'Action': 'Failed'
                })

    # 调整列顺序并保存
    out_df = pd.DataFrame(results)
    out_df.to_csv(output_path, index=False)
    print(f"✅ 成功评估/设计 {len(results)} 对双突变组合，已保存至: {output_path}")

def get_user_file(prompt_message, file_description, base_dir):
    """Prompts the user for a filename, checks for its existence in base_dir, and returns the full path and filename."""
    print(prompt_message)
    while True:
        filename = input(f"请输入 {file_description} 文件名: ").strip()
        if not filename:
            print("错误: 文件名不能为空，请重新输入。")
            continue
        file_path = os.path.join(base_dir, filename)
        if os.path.exists(file_path):
            return file_path, filename
        else:
            print(f"❌ 错误: 在程序目录下找不到文件 '{filename}'。")
            print("请检查: 1. 文件是否已放置在本程序相同目录下。 2. 文件名是否输入正确（包含.fa或.csv等后缀）。\n")

def main():
    print("="*50)
    print("🧬 Gibson Mutagenesis Primer Designer v0.2.1")
    print("="*50)

    # Determine the base directory (where the script is located)
    if getattr(sys, 'frozen', False):
        base_dir = os.path.dirname(sys.executable)
    else:
        base_dir = os.path.dirname(os.path.abspath(__file__))
    
    print(f"ℹ️  请将所有输入文件放置于本程序所在目录:\n   {base_dir}\n")

    fasta_prompt = "第一步: 请提供您的 FASTA 模板序列文件。"
    fasta_path, _ = get_user_file(fasta_prompt, "FASTA (.fa/.fasta)", base_dir)
    
    species = "ecoli"
    overhang_len = 40

    print("\n" + "="*50)
    print("请选择功能:")
    print("1. 单突变引物设计 (Single Mutation)")
    print("2. 双突变组合引物设计 (Double Mutation Combinations)")
    
    choice = input("\n请输入数字 (1 或 2): ").strip()

    if choice == '1':
        csv_prompt = "\n第二步: 请提供您的突变位点清单文件。\n(文件应不含表头，每行一个突变，格式如: A123G)"
        csv_path, csv_filename = get_user_file(csv_prompt, "突变位点 CSV (.csv)", base_dir)
        base_name = os.path.splitext(csv_filename)[0]
        output_path = os.path.join(base_dir, f"{base_name}-single_mut_primers.csv")
        design_single_mutations(fasta_path, csv_path, species, overhang_len, output_path)
    elif choice == '2':
        csv_prompt = "\n第二步: 请提供【第一步生成的】单突变引物结果文件。"
        csv_path, csv_filename = get_user_file(csv_prompt, "单突变引物结果 CSV (.csv)", base_dir)
        base_name = os.path.splitext(csv_filename)[0].replace('-single_mut_primers', '')
        output_path = os.path.join(base_dir, f"{base_name}-double_mut_primers.csv")
        design_double_mutations(fasta_path, csv_path, species, overhang_len, output_path)
    else:
        print("无效的选择，程序已退出。")

if __name__ == "__main__":
    main()