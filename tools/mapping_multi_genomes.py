import argparse
import os
import subprocess
from pathlib import Path
from colorama import Fore, Back, Style

def find_genomes(genome_dir: str):
    p = Path(genome_dir)
    if not p.exists() or not p.is_dir():
        raise FileNotFoundError(f"Genome dir not found: {genome_dir}")

    exts = {".fa", ".fna", ".fasta", ".fas"}
    genomes = [x for x in p.iterdir() if x.is_file() and x.suffix.lower() in exts]

    if not genomes:
        raise FileNotFoundError(f"No genome fasta files found in: {genome_dir}")
    return sorted(genomes)

def safe_name_from_fasta(path: Path) -> str:
    return path.stem

def run(cmd: str):
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f"Command failed (exit={ret}): {cmd}")

def run_and_capture(cmd: str) -> str:
    result = subprocess.run(
        cmd,
        shell=True,
        check=True,
        text=True,
        capture_output=True
    )
    return result.stdout.strip()

def remove_if_exists(path: str):
    if os.path.exists(path):
        os.remove(path)

def main():
    parser = argparse.ArgumentParser(
        description="Map ONT reads to ALL virus genomes in a folder; output only to coverage."
    )

    parser.add_argument("-f1", required=True, help="fastq R1 file (ONT single-end fastq)")
    parser.add_argument("-VGdir", required=True, help="Folder containing multiple virus genome FASTA files")
    parser.add_argument("-mt", default="minimap2", type=str, help="Mapping method (default: minimap2)")
    parser.add_argument("-t", default=20, type=int, help="Number of threads (default: 20)")
    parser.add_argument("-sid", required=True, type=str, help="Sample id (kept for compatibility / logging)")
    parser.add_argument("-oD", required=True, help="Output directory for results")

    parser.add_argument(
        "--keep_bam",
        action="store_true",
        help="Keep intermediate BAM and BAM index (.bai). Default: delete them after coverage is generated."
    )

    args = parser.parse_args()

    fastq_1 = os.path.abspath(args.f1)
    genome_dir = os.path.abspath(args.VGdir)
    output_dir = os.path.abspath(args.oD)
    threads = args.t
    method = args.mt
    keep_bam = args.keep_bam

    if not os.path.exists(fastq_1):
        raise FileNotFoundError(f"FASTQ not found: {fastq_1}")

    os.makedirs(output_dir, exist_ok=True)

    genomes = find_genomes(genome_dir)

    print(Back.GREEN + f"Mapping method is {method}." + Style.RESET_ALL)
    print(Back.GREEN + f"Found {len(genomes)} genome(s) in: {genome_dir}" + Style.RESET_ALL)

    if method != "minimap2":
        print(Fore.RED + f"Unsupported mapping method: {method}" + Style.RESET_ALL)
        raise SystemExit(1)

    for genome_path in genomes:
        genome_fa = str(genome_path.resolve())
        genome_name = safe_name_from_fasta(genome_path)

        sam_path = os.path.join(output_dir, f"mapping_{genome_name}.sam")
        bam_path = os.path.join(output_dir, f"mapping_{genome_name}_sorted.bam")
        bai_path = bam_path + ".bai"
        cov_path = os.path.join(output_dir, f"mapping_{genome_name}_sorted.coverage")
        cov_path_detailed = os.path.join(output_dir, f"mapping_{genome_name}_sorted.coverage.detailed")

        print(Back.GREEN + f"\n[Genome] {genome_name}  ({genome_path.name})" + Style.RESET_ALL)

        print(Back.GREEN + "Running minimap2..." + Style.RESET_ALL)
        run(f"minimap2 -ax map-ont -t {threads} {genome_fa} {fastq_1} > {sam_path}")

        print(Back.GREEN + "Doing samtools sort..." + Style.RESET_ALL)
        run(f"samtools sort -@ {threads} -o {bam_path} {sam_path}")
        run(f"rm -f {sam_path}")

        print(Back.GREEN + "Doing samtools index/coverage..." + Style.RESET_ALL)
        run(f"samtools index {bam_path}")
        run(f"samtools coverage {bam_path} > {cov_path_detailed}")
        run(f"samtools coverage -m {bam_path} > {cov_path}")

        # 檢查 mapped reads 數量
        mapped_reads = int(run_and_capture(f"samtools view -c -F 4 {bam_path}"))
        print(Back.GREEN + f"Mapped reads = {mapped_reads}" + Style.RESET_ALL)

        # 若完全沒有 mapping，刪除該 genome 的相關檔案
        if mapped_reads == 0:
            print(Fore.YELLOW + f"No mapped reads for {genome_name}, removing related files..." + Style.RESET_ALL)
            remove_if_exists(sam_path)
            remove_if_exists(cov_path)
            remove_if_exists(cov_path_detailed)
            remove_if_exists(bam_path)
            remove_if_exists(bai_path)
            continue

        # 若有 mapping，但使用者沒指定 keep_bam，就只刪除 bam/bai
        if not keep_bam:
            print(Back.GREEN + "Cleaning intermediate BAM/BAI..." + Style.RESET_ALL)
            remove_if_exists(bam_path)
            remove_if_exists(bai_path)

    print(Back.GREEN + "\nDone" + Style.RESET_ALL)

if __name__ == "__main__":
    main()