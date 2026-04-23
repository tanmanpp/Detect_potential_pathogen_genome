#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
from pathlib import Path
import sys
import time
import shutil
import re

def run_cmd(cmd, log_fp: Path, cwd: Path | None = None, allow_fail: bool = False) -> bool:
    """執行 shell 指令，stdout/stderr 同步寫入 log，失敗就中止。"""
    cmd_str = " ".join(str(x) for x in cmd)
    with log_fp.open("a", encoding="utf-8") as f:
        f.write(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] CMD: {cmd_str}\n")
        f.flush()

        p = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=f,
            stderr=subprocess.STDOUT,
            text=True
        )
    if p.returncode != 0:
        if allow_fail:
            note(f"[WARN] 指令失敗但已允許繼續 (exit={p.returncode}): {cmd_str}", log_fp)
            return False
        raise RuntimeError(f"指令失敗 (exit={p.returncode}): {cmd_str}")

    return True

def ensure_exists(p: Path, what: str = "檔案/資料夾"):
    if not p.exists():
        raise FileNotFoundError(f"{what}不存在：{p}")
    
def note(msg: str, log_fp: Path | None = None):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)  # 直接吐到螢幕（很重要：flush=True）
    if log_fp is not None:
        with log_fp.open("a", encoding="utf-8") as f:
            f.write(line + "\n")


def extract_species_from_filename(fname: str) -> str:
    """
    檔名格式：mapping__<Category>__mapping_<Genus>_<species>_sorted.coverage
    取 Genus_species，即跳過開頭的 'mapping' 後的前兩節。
    """
    stem = fname.replace(".coverage", "")
    parts = stem.split("__")
    # parts = ["mapping", "Bacteria", "mapping_Bacillus_cereus_sorted"]

    if len(parts) >= 3:
        sub = parts[2].split("_")
        # sub = ["mapping", "Bacillus", "cereus", "sorted"]
        # 跳過第一個 "mapping"，取接下來兩節
        tokens = [t for t in sub if t.lower() != "mapping"]
        return "_".join(tokens[:2]) if len(tokens) >= 2 else "_".join(tokens)

    # fallback
    return stem

def parse_coverage_file(fp: Path) -> dict | None:
    """
    解析 samtools coverage 視覺化文字格式，抽出：
    Number of reads, Percent covered, Mean coverage
    """
    try:
        with fp.open("r", encoding="utf-8", errors="replace") as f:
            text = f.read()

        def find_val(pattern):
            m = re.search(pattern, text)
            return m.group(1).strip() if m else None

        # Number of reads（取括號前的數字，不含 filtered）
        num_reads_str = find_val(r"Number of reads:\s*([\d,]+)")
        # Percent covered
        pct_str = find_val(r"Percent covered:\s*([\d.]+)%")
        # Mean coverage（去掉結尾的 x）
        mean_str = find_val(r"Mean coverage:\s*([\d.]+)")

        num_reads = int(num_reads_str.replace(",", "")) if num_reads_str else None
        pct       = float(pct_str) if pct_str else None
        mean_cov  = float(mean_str) if mean_str else None

        if all(v is None for v in [num_reads, pct, mean_cov]):
            return None

        return {
            "Number of reads": num_reads,
            "Percent covered": pct,
            "Mean coverage":   mean_cov,
        }

    except Exception:
        return None
    
def parse_coverage_detailed(fp: Path) -> dict | None:
    """
    解析 samtools coverage 的 table 輸出（.coverage.detailed）
    並做 multi-contig weighted summary
    """
    import csv

    try:
        total_length = 0
        total_covbases = 0
        total_reads = 0
        weighted_depth_num = 0
        weighted_depth_den = 0
        contigs = []

        with fp.open("r", encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f, delimiter="\t")

            for row in reader:
                start = int(row["startpos"])
                end = int(row["endpos"])
                length = end - start + 1

                numreads = int(row["numreads"])
                covbases = int(row["covbases"])
                coverage = float(row["coverage"])
                meandepth = float(row["meandepth"])

                contigs.append({
                    "rname": row["#rname"],
                    "length": length,
                    "numreads": numreads,
                    "covbases": covbases,
                    "coverage": coverage,
                    "meandepth": meandepth,
                    "meanbaseq": float(row["meanbaseq"]) if row.get("meanbaseq") not in (None, "") else None,
                    "meanmapq": float(row["meanmapq"]) if row.get("meanmapq") not in (None, "") else None,
                })

                total_length += length
                total_covbases += covbases
                total_reads += numreads
                weighted_depth_num += meandepth * length
                weighted_depth_den += length

        if total_length == 0 or not contigs:
            return None

        overall_percent = total_covbases / total_length * 100.0
        overall_depth = weighted_depth_num / weighted_depth_den if weighted_depth_den > 0 else None
        best_contig = max(contigs, key=lambda x: x["coverage"])

        return {
            "Number of contigs": len(contigs),
            "Total length": total_length,
            "Covered bases": total_covbases,
            "Number of reads": total_reads,
            "Percent covered": overall_percent,
            "Mean coverage": overall_depth,
            "Best contig": best_contig["rname"],
            "Best contig coverage": best_contig["coverage"],
            "contigs": contigs,
        }

    except Exception:
        return None

def main():
    project_root = Path(__file__).resolve().parents[1]
    default_kraken_html_script = project_root / "tools" / "kraken2_html_excel_report.py"
    default_mapping_script = project_root / "tools" / "mapping_multi_genomes.py"
    default_blast_script = project_root / "tools" / "run_blastn_with_virulent_database.sh"
    default_report_script = project_root / "tools" / "generate_pipeline_summary_html.py"
    default_vfdb = project_root / "VFDB" / "VFDB_setA_nt_rename.fas"

    ap = argparse.ArgumentParser(
        description="Detect potential pathogen pipeline (len filter -> rm human -> flye -> kraken2/rcf -> mapping -> blastn)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    ap.add_argument("-i", "--reads", required=True, help="輸入 FASTQ（建議為 demux 後單一樣本 reads）")
    ap.add_argument("-sid", "--sample_id", required=True, help="樣本 ID（會用於輸出檔名/資料夾）")
    ap.add_argument("-o", "--outdir", default="output", help="輸出根目錄")
    ap.add_argument("-t", "--threads", type=int, default=8, help="threads")

    # Step 1: seqkit
    ap.add_argument("--min_len", type=int, default=200, help="SeqKit length filter 最小長度")

    # Step 2: human removal
    ap.add_argument("--human_index", required=True, help="minimap2 human index prefix（例如 hg38.mmi 或 minimap2 index prefix）")

    # Step 4: kraken2 + rcf
    ap.add_argument("--kraken_db", default="/home/yilun/YiLun/kraken_database/standard", help="Kraken2 DB 路徑")
    ap.add_argument("--rcf_env", default="DPPP", help="recentrifuge 的 conda env 名稱")
    ap.add_argument("--rcf_taxdb", default="/home/yilun/YiLun/kraken_database/standard", help="rcf -n 指向的 taxonomy/db 路徑")
    ap.add_argument("--kraken_sc", type=float, default=0.2, help="Kraken2 confidence score cutoff（k2 --confidence）")
    ap.add_argument("--kraken_html_script", default=str(default_kraken_html_script), help="生成 Kraken2 HTML 和 Excel 報告的 python 腳本路徑")
    # Step 5: mapping_multi_genomes.py
    ap.add_argument("--mapping_script", default=str(default_mapping_script), help="mapping_multi_genomes.py 路徑")
    ap.add_argument("--vg_bacteria", default=str(project_root / "Bacteria"), help="Bacteria genomes dir")
    ap.add_argument("--vg_virus", default=str(project_root / "Virus"), help="Virus genomes dir")
    ap.add_argument("--vg_parasite", default=str(project_root / "Parasite"), help="Parasite genomes dir")

    # Step 6: blastn
    ap.add_argument("--blast_script", default=str(default_blast_script),
                    help="run_blastn_with_virulent_database.sh 路徑")
    ap.add_argument("--vfdb", default=str(default_vfdb),
                    help="VFDB fasta（腳本的 -db 參數）")

    # Optional: step control
    ap.add_argument("--skip_len_filter", action="store_true", help="跳過 SeqKit 長度過濾（不建議）")
    ap.add_argument("--skip_rm_human", action="store_true", help="跳過去除人類序列（不建議）")
    ap.add_argument("--skip_flye", action="store_true", help="跳過 Flye（若你已經有組裝結果）")
    ap.add_argument("--skip_kraken", action="store_true", help="跳過 Kraken2/rcf")
    ap.add_argument("--skip_mapping", action="store_true", help="跳過 mapping_multi_genomes.py")
    ap.add_argument("--skip_blast", action="store_true", help="跳過 blastn")

    args = ap.parse_args()

    reads_in = Path(args.reads)
    ensure_exists(reads_in, "輸入 reads FASTQ")

    out_root = Path(args.outdir).resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    sid = args.sample_id

    # 輸出資料夾結構
    step1_dir = out_root / sid / "01_len_filter"
    step2_dir = out_root / sid / "02_rm_human"
    denovo_dir = out_root / sid / "03_denovo_flye"
    kraken_dir = out_root / sid / "04_kraken_rcf"
    mapping_dir = out_root / sid / "05_mapping_multi_genomes"
    blast_dir = out_root / sid / "06_blastn_vfdb"
    log_dir = out_root / sid / "logs"
    final_dir = out_root / sid / "final_results"
    final_dir.mkdir(parents=True, exist_ok=True)


    for d in [step1_dir, step2_dir, denovo_dir, kraken_dir, mapping_dir, blast_dir, log_dir]:
        d.mkdir(parents=True, exist_ok=True)

    log_fp = log_dir / f"{sid}.pipeline.log"

    # -------------------------
    # STEP 1) SeqKit length filter
    # -------------------------
    len200_fastq = step1_dir / f"{sid}.len{args.min_len}.fastq"
    
    if not args.skip_len_filter:
        note(f"STEP 1 START: SeqKit length filter (min_len={args.min_len})", log_fp)

        
        run_cmd(
            ["seqkit", "seq", "-m", str(args.min_len), str(reads_in), "-o", str(len200_fastq)],
            log_fp
        )
        ensure_exists(len200_fastq, "SeqKit 輸出 FASTQ")
        note(f"STEP 1 DONE: {len200_fastq}", log_fp)

    # -------------------------
    # STEP 2) Remove human reads (minimap2 + samtools extract unmapped)
    #   2.1 map to human -> BAM
    #   2.2 keep unmapped (-f 4) -> FASTQ
    # -------------------------
    human_bam = step2_dir / f"{sid}.map_human.bam"
    rmhum_fastq = step2_dir / f"{sid}.rmhum.fastq"
    rmhum_fasta = step2_dir / f"{sid}.rmhum.fasta"

    if not args.skip_rm_human:

        note("STEP 2 START: Remove human reads", log_fp)


        # minimap2 -> samtools sort -> bam
        # 注意：這裡用 map-ont（Nanopore）常見設定；你若要更嚴格可自行加 -p / -N 等參數
        run_cmd(
            ["bash", "-lc",
            f"minimap2 -t {args.threads} -ax map-ont {args.human_index} {len200_fastq} | "
            f"samtools sort -@ {args.threads} -o {human_bam} -"],
            log_fp
        )
        ensure_exists(human_bam, "Human mapping BAM")

        run_cmd(["samtools", "index", str(human_bam)], log_fp)

        # 抽出 unmapped reads
        cmd = f"""
        set -euo pipefail
        samtools fastq -f 4 {human_bam}> {rmhum_fastq} 
        """

        run_cmd(
            ["bash", "-lc", cmd],
            log_fp
        )

        run_cmd(
            ["seqkit", "fq2fa", str(rmhum_fastq), "-o", str(rmhum_fasta)],
        log_fp
        )

        ensure_exists(rmhum_fastq, "去人類後 reads FASTQ")
        note(f"STEP 2 DONE: {rmhum_fastq}", log_fp)

    # -------------------------
    # STEP 3) Flye de novo (meta)
    # -------------------------
    
    contigs_fa = denovo_dir / "assembly.fasta"  # Flye 常見輸出
    contigs_available = False
    if not args.skip_flye:
        note("STEP 3 START: Flye de novo assembly", log_fp)
        flye_ok = run_cmd(
            ["flye",
             "--nano-hq", str(rmhum_fastq),
             "--meta",
             "--threads", str(args.threads),
             "-o", str(denovo_dir)],
            log_fp,
            allow_fail=True
        )
        contigs_available = flye_ok and contigs_fa.exists() and contigs_fa.is_file() and contigs_fa.stat().st_size > 0
        if contigs_available:
            note(f"STEP 3 DONE: {contigs_fa}", log_fp)
        else:
            note("[WARN] Flye did not produce a usable assembly.fasta. The pipeline will continue with reads-based analyses and skip denovo-contig analyses.", log_fp)
    else:
        contigs_available = contigs_fa.exists() and contigs_fa.is_file() and contigs_fa.stat().st_size > 0
        if contigs_available:
            note(f"[INFO] --skip_flye enabled; using existing contigs: {contigs_fa}", log_fp)
        else:
            note("[INFO] --skip_flye enabled and no existing assembly.fasta was found. Denovo-contig analyses will be skipped.", log_fp)

    # -------------------------
    # STEP 4) Kraken2 + Recentrifuge + html & excel report
    # 會跑兩份：reads.rmhum.fastq & denovo contigs
    # -------------------------
    
    def kraken_and_rcf(input_fp: Path, tag: str):
        note("STEP 4 START: Kraken2 + Recentrifuge", log_fp)
        # Kraken output（分類逐行）+ report（彙總）
        k_out = kraken_dir / f"{sid}.{tag}.kraken.out"
       # k_report = kraken_dir / f"{sid}.{tag}.kraken.report"
        kraken_cmd = ["k2", "classify"] if shutil.which("k2") else ["kraken2"]

        run_cmd(
            [*kraken_cmd,
             "--db", args.kraken_db,
             "--threads", str(args.threads),
             "--output", str(k_out),
        #    "--report", str(k_report),
             "--confidence", str(args.kraken_sc),
             str(input_fp)],
            log_fp
        )
        ensure_exists(k_out, "Kraken2 output")
        #ensure_exists(k_report, "Kraken2 report")

        # Recentrifuge (rcf)
        rcf_out = kraken_dir / f"{sid}.{tag}"
        run_cmd(
            ["conda", "run", "-n", args.rcf_env,
             "rcf",
             "-n", args.rcf_taxdb,
             "-k", str(k_out),
             "-o", str(rcf_out)],
            log_fp
        )
        rcf_out_check = kraken_dir / f"{sid}.{tag}.rcf.html"
        ensure_exists(rcf_out_check, "Recentrifuge 輸出")

        # 生成 HTML 和 Excel python script
        ensure_exists(Path(args.kraken_html_script), "kraken2_html_excel_report.py")
        report_output_prefix = f"{sid}.{tag}.kraken_report"
        run_cmd(
            ["python", str(Path(args.kraken_html_script)),
             "-k", str(k_out),
             "--nodes", str(Path(args.kraken_db) / 'nodes.dmp'),
             "--names", str(Path(args.kraken_db) / 'names.dmp'),
             "--out_dir", str(kraken_dir),
             "-o", str(report_output_prefix)],
            log_fp
        )
        ensure_exists(kraken_dir / f"{report_output_prefix}.html", "Kraken2 HTML report")
        ensure_exists(kraken_dir / f"{report_output_prefix}.xlsx", "Kraken2 Excel report")

        note(f"STEP 4 DONE: Kraken2 + Recentrifuge & report for {tag}", log_fp)

    if not args.skip_kraken:
        kraken_and_rcf(rmhum_fastq, "rmhum_reads")
        if contigs_available:
            kraken_and_rcf(contigs_fa, "denovo_contigs")
        else:
            note("[INFO] STEP 4 SKIP: denovo_contigs Kraken2/Recentrifuge skipped because no usable Flye contigs were produced.", log_fp)
    
    

    # -------------------------
    # STEP 5) mapping_multi_genomes.py (跑三次：Bacteria/Virus/Parasite)
    # -------------------------
    if not args.skip_mapping:
        ensure_exists(Path(args.mapping_script), "mapping_multi_genomes.py")
        vg_sets = [
            ("Bacteria", Path(args.vg_bacteria)),
            ("Virus", Path(args.vg_virus)),
            ("Parasite", Path(args.vg_parasite)),
        ]
        for label, vgdir in vg_sets:
            note(f"STEP 5 START: mapping {label}", log_fp)
            ensure_exists(vgdir, f"{label} genome dir")
            out_d = mapping_dir / label
            out_d.mkdir(parents=True, exist_ok=True)

            run_cmd(
                ["python", str(Path(args.mapping_script)),
                 "-f1", str(rmhum_fastq),
                 "-VGdir", str(vgdir),
                 "-t", str(args.threads),
                 "-sid", str(sid),
                 "--keep_bam",
                 "-oD", str(out_d)],
                log_fp
            )
            note(f"STEP 5 DONE: mapping {label}", log_fp)

    # -------------------------
    # STEP 6) Blastn with virulence DB (reads & contigs 各跑一份)
    # -------------------------
    if not args.skip_blast:
        note("STEP 6 START: blastn with VFDB", log_fp)
        blast_script = Path(args.blast_script)
        ensure_exists(blast_script, "blast script")
        ensure_exists(Path(args.vfdb), "VFDB fasta")

        def run_blast(input_fp: Path, tag: str):
            out_fp_perfix = f"{sid}.{tag}.blastn"
            run_cmd(
                ["bash", str(blast_script),
                 "-db", str(Path(args.vfdb)),
                 "-i", str(input_fp),
                 "--prefix", str(out_fp_perfix),
                 "--outdir", str(blast_dir)],
                log_fp
            )
            out_fp = blast_dir / f"{out_fp_perfix}.blast.tsv"
            ensure_exists(out_fp, "blastn result")

        run_blast(rmhum_fasta, "rmhum_reads")
        if contigs_available:
            run_blast(contigs_fa, "denovo_contigs")
        else:
            note("[INFO] STEP 6 SKIP: denovo_contigs blastn skipped because no usable Flye contigs were produced.", log_fp)
        note("STEP 6 DONE: blastn with VFDB", log_fp)


    def _copy_to_final(src: Path, dst_dir: Path, new_name: str | None = None) -> Path:
        """複製檔案到 final，保留 metadata。回傳目的地路徑。"""
        dst_dir.mkdir(parents=True, exist_ok=True)
        dst = dst_dir / (new_name if new_name else src.name)
        shutil.copy2(src, dst)
        return dst

    def _is_nonempty_file(p: Path) -> bool:
        return p.exists() and p.is_file() and p.stat().st_size > 0

    def _tsv_has_second_line(tsv: Path) -> bool:
        """判斷 TSV 是否至少有 2 行（header + 至少 1 筆命中）。"""
        if not _is_nonempty_file(tsv):
            return False
        try:
            with tsv.open("r", encoding="utf-8", errors="replace") as f:
                _ = f.readline()          # line 1
                line2 = f.readline()      # line 2
                return bool(line2.strip())
        except Exception:
            return False

    # -------------------------
    # FINAL) Collect results to final_results
    # 規則：
    # - 03_denovo_flye/assembly.fasta -> final/{sid}_denovo.fasta
    # - 04_kraken_rcf/*.html -> final/
    # - 05_mapping_multi_genomes/ 內所有 size>0 的檔案 -> final/（扁平化命名）
    # - 06_blastn_vfdb/*.tsv 若至少有第二行 -> final/
    # 並輸出 denovo/kraken/mapping/blast 是否有結果
    # -------------------------
    note("FINAL START: Collect results", log_fp)

    denovo_has_result = False
    kraken_has_result = False
    mapping_has_result = False
    blast_has_result = False

    # 1) denovo
    if _is_nonempty_file(contigs_fa):
        _copy_to_final(contigs_fa, final_dir, new_name=f"{sid}_denovo.fasta")
        denovo_has_result = True
    else:
        note(f"[WARN] denovo 沒有可用 contigs：{contigs_fa}", log_fp)

    # 2) kraken/rcf html + report.py html/excel
    htmls = list(kraken_dir.glob("*.html"))
    htmls = [x for x in htmls if _is_nonempty_file(x)]
    excels = list(kraken_dir.glob("*.xlsx"))
    excels = [x for x in excels if _is_nonempty_file(x)]
    if htmls or excels:
        for h in htmls:
            _copy_to_final(h, final_dir)
        for e in excels:
            _copy_to_final(e, final_dir)
        kraken_has_result = True
    else:
        note(f"[WARN] kraken/rcf + report.py 沒有可收集的 html/excel：{kraken_dir}", log_fp)

    

    # 3) mapping：收集所有 size>0 的檔案（避免撞名：用相對路徑扁平化）以及刪除 0 size 檔案
    if mapping_dir.exists():
        copied_any = False
        for f in mapping_dir.rglob("*"):
            if not f.is_file():
                continue

            # 刪除所有 0 byte 檔案（不分副檔名）
            if f.stat().st_size == 0 and f.name.endswith(".coverage"):
                run_cmd(["rm", "-f", str(f)], log_fp, allow_fail=True)
                continue

            # 只收集 .coverage 結尾
            if not f.name.endswith(".coverage"):
                continue

            rel = f.relative_to(mapping_dir).as_posix()
            flat_name = rel.replace("/", "__").replace(" ", "_")
            _copy_to_final(f, final_dir, new_name=f"mapping__{flat_name}")
            copied_any = True

        mapping_has_result = copied_any
        if not copied_any:
            note(f"[WARN] mapping 目錄存在但沒有 size>0 的 .coverage 檔案可收集：{mapping_dir}", log_fp)

        # ---- 生成 mapping_summary.csv（以 final_dir 中實際保留的 .coverage 為準，
        #      回 mapping_dir 找對應的 .sorted.coverage.detailed 來整理） ----
        import csv

        summary_rows_mapping = []

        # 先找 final_results 中真正被保留的 mapping coverage 檔
        final_mapping_cov_files = sorted(final_dir.glob("mapping__*.coverage"))

        for final_cov in final_mapping_cov_files:
            # final_cov.name 例如：
            # mapping__Bacteria__subdir__mapping_Escherichia_coli_sorted.coverage

            # 還原它在 mapping_dir 中的相對路徑
            # 先去掉前綴 mapping__
            rel_flat = final_cov.name[len("mapping__"):]
            # 還原成原始相對路徑
            rel_path_str = rel_flat.replace("__", "/")

            # 原本 coverage 在 mapping_dir 內的位置
            original_cov = mapping_dir / rel_path_str

            # 對應的 detailed 檔：
            # xxx.coverage -> xxx.coverage.detailed
            original_detailed = Path(str(original_cov) + ".detailed")

            if not original_detailed.exists() or not original_detailed.is_file():
                note(f"[WARN] 找不到對應的 detailed coverage 檔：{original_detailed}", log_fp)
                continue

            if original_detailed.stat().st_size == 0:
                note(f"[WARN] detailed coverage 檔為空：{original_detailed}", log_fp)
                continue

            parsed = parse_coverage_detailed(original_detailed)
            if parsed is None:
                note(f"[WARN] 無法解析 detailed coverage 檔：{original_detailed.name}", log_fp)
                continue

            # species 名稱仍沿用原本 final coverage 檔名來抽
            species = extract_species_from_filename(final_cov.name)

            # 分類（Bacteria / Virus / Parasite）
            category = ""
            m = re.match(r"mapping__([^_]+)__", final_cov.name)
            if m:
                category = m.group(1)

            summary_rows_mapping.append({
                "Category": category,
                "Species": species,
                "Number of contigs": parsed["Number of contigs"],
                "Total length": parsed["Total length"],
                "Covered bases": parsed["Covered bases"],
                "Number of reads": parsed["Number of reads"],
                "Percent covered": round(parsed["Percent covered"], 4) if parsed["Percent covered"] is not None else None,
                "Mean coverage": round(parsed["Mean coverage"], 4) if parsed["Mean coverage"] is not None else None,
                "Best contig": parsed["Best contig"],
                "Best contig coverage": round(parsed["Best contig coverage"], 4) if parsed["Best contig coverage"] is not None else None,
                "Final coverage file": final_cov.name,
                "Detailed source file": str(original_detailed.relative_to(mapping_dir)),
            })

        if summary_rows_mapping:
            csv_out = final_dir / "mapping_summary.csv"
            fieldnames = [
                "Category",
                "Species",
                "Number of contigs",
                "Total length",
                "Covered bases",
                "Number of reads",
                "Percent covered",
                "Mean coverage",
                "Best contig",
                "Best contig coverage",
                "Final coverage file",
                "Detailed source file",
            ]
            with csv_out.open("w", newline="", encoding="utf-8") as csvf:
                writer = csv.DictWriter(csvf, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(summary_rows_mapping)

            note(f"[INFO] mapping_summary.csv 已生成：{csv_out}（{len(summary_rows_mapping)} 筆）", log_fp)
        else:
            note("[WARN] 沒有可匯整的 detailed coverage 資料，mapping_summary.csv 未生成", log_fp)

    else:
        note(f"[WARN] mapping 目錄不存在：{mapping_dir}", log_fp)

    # 4) blast tsv：至少有第二行才算有結果
    tsvs = list(blast_dir.glob("*.tsv"))
    copied_any = False
    for tsv in tsvs:
        if _tsv_has_second_line(tsv):
            _copy_to_final(tsv, final_dir)
            copied_any = True
        else:
            # 有檔但沒命中（只有header或空）
            if tsv.exists():
                note(f"[INFO] blast tsv 無第二行（可能無命中），不收集：{tsv.name}", log_fp)

    blast_has_result = copied_any
    if not copied_any:
        note(f"[WARN] blast 沒有可收集的有效 tsv（需至少兩行）：{blast_dir}", log_fp)


    # -------------------------
    # FINAL HTML REPORT
    # -------------------------
    report_script = default_report_script
    ensure_exists(report_script, "generate_pipeline_summary_html.py")

    summary_html = final_dir / f"{sid}.summary_report.html"

    run_cmd(
        ["python", str(report_script),
        "--final_dir", str(final_dir),
        "--sid", str(sid),
        "--mapping_reads_dir", str(mapping_dir),
        "--fastq", str(rmhum_fastq),
        "--out", str(summary_html)],
        log_fp
    )
    ensure_exists(summary_html, "summary html report")
    note(f"[INFO] summary html 已生成：{summary_html}", log_fp)


    # summary（你要的 4 個步驟是否有結果）
    summary = (
        f"[SUMMARY] results collected to: {final_dir}\n"
        f"  - denovo : {'YES' if denovo_has_result else 'NO'}\n"
        f"  - kraken : {'YES' if kraken_has_result else 'NO'}\n"
        f"  - mapping: {'YES' if mapping_has_result else 'NO'}\n"
        f"  - blastn : {'YES' if blast_has_result else 'NO'}"
    )
    note(summary, log_fp)

    print(f"[OK] Pipeline finished. Output: {out_root / sid}")
    print(f"[LOG] {log_fp}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
