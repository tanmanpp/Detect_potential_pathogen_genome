#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import html
import os
import re
import json
from datetime import datetime
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote, quote_plus
try:
    import pysam
except ImportError:
    pysam = None


def read_text(fp: Path) -> str:
    return fp.read_text(encoding="utf-8", errors="replace")


def read_fastq_ids(fastq_path: Path) -> set:
    """從 fastq / fastq.gz 讀取所有 read ID（@ 行第一個 token）。"""
    import gzip
    ids = set()
    opener = gzip.open if fastq_path.suffix in {".gz", ".gzip"} else open
    with opener(str(fastq_path), "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0 and line.startswith("@"):
                ids.add(line[1:].split()[0].rstrip())
    return ids


def guess_file_label(fp: Path) -> str:
    name = fp.name.lower()
    if "rmhum" in name or "rm_human" in name:
        return "rm_human"
    if "denovo" in name or "contig" in name:
        return "contigs"
    return "other"


def extract_species_from_mapping_filename(fname: str) -> str:
    stem = fname.replace(".coverage", "")
    parts = stem.split("__")
    if len(parts) >= 3:
        sub = parts[2].split("_")
        tokens = [t for t in sub if t.lower() != "mapping"]
        return "_".join(tokens[:2]) if len(tokens) >= 2 else "_".join(tokens)
    return stem


def parse_coverage_metrics(text: str) -> dict:
    def find_val(pattern: str):
        m = re.search(pattern, text)
        return m.group(1).strip() if m else None

    num_reads_str = find_val(r"Number of reads:\s*([\d,]+)")
    pct_str = find_val(r"Percent covered:\s*([\d.]+)%")
    mean_str = find_val(r"Mean coverage:\s*([\d.]+)")

    return {
        "Number of reads": int(num_reads_str.replace(",", "")) if num_reads_str else None,
        "Percent covered": float(pct_str) if pct_str else None,
        "Mean coverage": float(mean_str) if mean_str else None,
    }


def parse_tsv(fp: Path, max_rows: int = 200):
    rows = []
    with fp.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for i, row in enumerate(reader):
            rows.append(row)
            if i + 1 >= max_rows:
                break

    if not rows:
        return [], [], 0

    header = rows[0]
    body = rows[1:]
    total_rows = 0
    with fp.open("r", encoding="utf-8", errors="replace", newline="") as f:
        for total_rows, _ in enumerate(f, start=1):
            pass
    data_rows = max(total_rows - 1, 0)
    return header, body, data_rows


def file_link_text(fp: Path, base_dir: Path) -> str:
    rel = fp.relative_to(base_dir).as_posix()
    return quote(rel)


def species_google_url(species: str) -> str:
    # Replace underscores with spaces for better search results
    query = species.replace("_", " ")
    return f"https://www.google.com/search?q={quote_plus(query)}"


def html_table_with_species_links(headers, rows, max_preview=None):
    """Render table where the first column (Species) is a Google search link."""
    if not headers and not rows:
        return "<p class='muted'>無資料</p>"

    preview_rows = rows if max_preview is None else rows[:max_preview]
    out = ["<div class='table-wrap'><table>"]
    out.append("<thead><tr>")
    for h in headers:
        out.append(f"<th>{html.escape(str(h))}</th>")
    out.append("</tr></thead><tbody>")

    for row in preview_rows:
        out.append("<tr>")
        for i, cell in enumerate(row):
            if i == 0:
                # Species column — add Google search link
                display = html.escape(str(cell))
                url = species_google_url(str(cell))
                out.append(f"<td><a class='species-link' href='{url}' target='_blank' rel='noopener'>{display}</a></td>")
            else:
                out.append(f"<td>{html.escape(str(cell))}</td>")
        out.append("</tr>")

    out.append("</tbody></table></div>")
    return "".join(out)


def html_table(headers, rows, max_preview=None):
    if not headers and not rows:
        return "<p class='muted'>無資料</p>"

    preview_rows = rows if max_preview is None else rows[:max_preview]
    out = ["<div class='table-wrap'><table>"]
    out.append("<thead><tr>")
    for h in headers:
        out.append(f"<th>{html.escape(str(h))}</th>")
    out.append("</tr></thead><tbody>")

    for row in preview_rows:
        out.append("<tr>")
        for cell in row:
            out.append(f"<td>{html.escape(str(cell))}</td>")
        out.append("</tr>")

    out.append("</tbody></table></div>")
    return "".join(out)


def html_table_with_auto_species_links(headers, rows, max_preview=None):
    """Render table and auto-link any column whose header contains 'species' (case-insensitive)."""
    if not headers and not rows:
        return "<p class='muted'>無資料</p>"

    # 找出所有 header 含 "species" 的欄位 index
    species_cols = {i for i, h in enumerate(headers) if "species" in str(h).lower()}

    preview_rows = rows if max_preview is None else rows[:max_preview]
    out = ["<div class='table-wrap'><table>"]
    out.append("<thead><tr>")
    for h in headers:
        out.append(f"<th>{html.escape(str(h))}</th>")
    out.append("</tr></thead><tbody>")

    for row in preview_rows:
        out.append("<tr>")
        for i, cell in enumerate(row):
            if i in species_cols and str(cell).strip():
                display = html.escape(str(cell))
                url = species_google_url(str(cell))
                out.append(f"<td><a class='species-link' href='{url}' target='_blank' rel='noopener'>{display}</a></td>")
            else:
                out.append(f"<td>{html.escape(str(cell))}</td>")
        out.append("</tr>")

    out.append("</tbody></table></div>")
    return "".join(out)



def collect_mapping(final_dir: Path):
    cats = {"Bacteria": [], "Parasite": [], "Virus": []}

    csv_fp = final_dir / "mapping_summary.csv"
    if not csv_fp.exists():
        return cats

    with csv_fp.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            category = row.get("Category", "").strip()
            if category not in cats:
                continue

            species = row.get("Species", "").strip()

            metrics = {
                "Number of reads": int(float(row["Number of reads"])) if row.get("Number of reads") else None,
                "Percent covered": float(row["Percent covered"]) if row.get("Percent covered") else None,
                "Mean coverage": float(row["Mean coverage"]) if row.get("Mean coverage") else None,
            }

            final_cov_name = row.get("Final coverage file", "").strip()
            final_cov_fp = final_dir / final_cov_name if final_cov_name else None

            # 保留 .coverage 的完整文字內容，給 HTML 展開顯示
            raw_text = ""
            if final_cov_fp and final_cov_fp.exists() and final_cov_fp.is_file():
                raw_text = read_text(final_cov_fp)
            else:
                raw_text = (
                    f"[WARN] 找不到對應的 .coverage 檔\n"
                    f"Final coverage file: {final_cov_name}"
                )

            cats[category].append({
                "file": final_cov_fp if final_cov_fp else Path(f"{species}.coverage"),
                "species": species,
                "metrics": metrics,
                "raw_text": raw_text,
                "image": None,
            })

    # 若未來有把 mapping 圖另外複製進 final_dir，仍可自動掛上
    image_suffixes = {".png", ".jpg", ".jpeg", ".svg", ".webp"}
    mapping_images = {}
    for fp in final_dir.iterdir():
        if fp.is_file() and fp.suffix.lower() in image_suffixes and fp.name.startswith("mapping__"):
            mapping_images[fp.stem] = fp

    for cat, items in cats.items():
        for item in items:
            stem = item["file"].stem
            item["image"] = mapping_images.get(stem)

    return cats


def collect_mapping_bam_pie_data(mapping_reads_dir: Path, fastq_ids: set = None):
    """
    mapping_reads_dir/
      Bacteria/ Parasite/ Virus/  →  mapping_*_sorted.bam

    fastq_ids: set of read IDs from input fastq (optional).
               If provided, only reads present in fastq are counted,
               and unmapped_count is computed.

    Returns:
    {
        "sets":           {label: set_of_read_ids},
        "genome_meta":    {label: {"category":..., "genome_name":...}},
        "fastq_total":    int,
        "unmapped_count": int,
        "total_mapped_unique_reads": int,
    }
    """
    if pysam is None:
        raise ImportError("需要安裝 pysam 才能讀取 BAM：pip install pysam")

    categories = ["Bacteria", "Parasite", "Virus"]
    genome_sets = {}   # label -> set of primary-mapped read IDs
    genome_meta = {}   # label -> {category, genome_name}

    for category in categories:
        cat_dir = mapping_reads_dir / category
        if not cat_dir.exists() or not cat_dir.is_dir():
            continue
        for bam_fp in sorted(cat_dir.glob("mapping_*_sorted.bam")):
            gname = re.sub(r"^mapping_", "", bam_fp.name)
            gname = re.sub(r"_sorted\.bam$", "", gname)
            label = f"{category} | {gname}"
            read_ids = set()
            with pysam.AlignmentFile(str(bam_fp), "rb") as bam:
                for aln in bam.fetch(until_eof=True):
                    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                        continue
                    if aln.query_name:
                        read_ids.add(aln.query_name)
            # Intersect with fastq universe if provided
            if fastq_ids is not None:
                read_ids &= fastq_ids
            genome_sets[label] = read_ids
            genome_meta[label] = {"category": category, "genome_name": gname}

    all_mapped: set = set()
    for s in genome_sets.values():
        all_mapped |= s
    total_mapped = len(all_mapped)

    fastq_total    = len(fastq_ids) if fastq_ids else 0
    unmapped_count = len(fastq_ids - all_mapped) if fastq_ids else 0

    return {
        "sets":           genome_sets,
        "genome_meta":    genome_meta,
        "fastq_total":    fastq_total,
        "unmapped_count": unmapped_count,
        "total_mapped_unique_reads": total_mapped,
    }

def collect_blast(final_dir: Path):
    sections = {
        "rm_human": {"blastn_best": [], "species_counts": [], "other_tsv": []},
        "contigs": {"blastn_best": [], "species_counts": [], "other_tsv": []},
    }

    for fp in sorted(final_dir.glob("*.tsv")):
        label = guess_file_label(fp)
        if label not in sections:
            continue

        name = fp.name.lower()
        if "blastn.best.tsv" in name:
            sections[label]["blastn_best"].append(fp)
        elif "species_counts.tsv" in name:
            sections[label]["species_counts"].append(fp)
        else:
            sections[label]["other_tsv"].append(fp)

    return sections


# Category badge colours
CAT_COLORS = {
    "Bacteria": ("#e8f5e9", "#2e7d32", "#43a047"),
    "Parasite": ("#fff3e0", "#e65100", "#fb8c00"),
    "Virus":    ("#fce4ec", "#880e4f", "#e91e63"),
}


def render_mapping_section(mapping_data, final_dir: Path) -> str:
    out = ["<section id='mapping'><h1>🧬 Mapping 結果</h1>"]
    out.append("<p class='muted'>依 Bacteria / Parasite / Virus 分類，並統整 final_results 中的 mapping 檔案。</p>")

    for cat in ["Virus", "Bacteria", "Parasite"]:
        items = sorted(
            mapping_data.get(cat, []),
            key=lambda x: x["metrics"]["Number of reads"] or 0,
            reverse=True,
        )
        bg, fg, accent = CAT_COLORS.get(cat, ("#f4f7fb", "#172033", "#2457d6"))
        out.append(
            f"<div id='mapping-{cat}' class='cat-block' style='--cat-bg:{bg};--cat-fg:{fg};--cat-accent:{accent}'>"
            f"<div class='cat-header'><span class='cat-badge'>{cat}</span>"
            f"<span class='cat-count'>{len(items)} 筆</span></div>"
        )

        if not items:
            out.append("<p class='muted' style='padding:12px 0'>此分類沒有 mapping 結果。</p></div>")
            continue

        # Summary table for the category
        rows = []
        for item in items:
            m = item["metrics"]
            rows.append([
                item["species"],
                "" if m["Number of reads"] is None else m["Number of reads"],
                "" if m["Percent covered"] is None else f'{m["Percent covered"]:.2f}',
                "" if m["Mean coverage"] is None else f'{m["Mean coverage"]:.2f}',
                item["file"].name,
            ])

        out.append("<div class='summary-table-wrap'>")
        out.append(html_table_with_species_links(
            ["Species", "Number of reads", "Percent covered (%)", "Mean coverage", "Source file"],
            rows
        ))
        out.append("</div>")

        # Individual collapsible cards — one per species
        for item in items:
            m = item["metrics"]
            species_display = html.escape(item["species"])
            google_url = species_google_url(item["species"])

            reads_val = m["Number of reads"] if m["Number of reads"] is not None else "NA"
            covered_val = f'{m["Percent covered"]:.2f}%' if m["Percent covered"] is not None else "NA"
            mean_val = f'{m["Mean coverage"]:.2f}' if m["Mean coverage"] is not None else "NA"

            out.append("<details class='species-card'>")
            out.append(
                f"<summary class='species-summary'>"
                f"<span class='species-name'>"
                f"<a class='species-link' href='{google_url}' target='_blank' rel='noopener' "
                f"   onclick='event.stopPropagation()'>{species_display}</a>"
                f"</span>"
                f"<span class='summary-pills'>"
                f"<span class='pill pill-reads'>reads: {reads_val}</span>"
                f"<span class='pill pill-covered'>covered: {covered_val}</span>"
                f"<span class='pill pill-mean'>mean cov: {mean_val}</span>"
                f"</span>"
                f"<span class='expand-icon'>▸</span>"
                f"</summary>"
            )

            out.append("<div class='species-body'>")

            if item.get("image"):
                rel = file_link_text(item["image"], final_dir)
                out.append(f"<img class='plot' src='{rel}' alt='mapping plot for {species_display}'>")

            out.append(
                "<div class='coverage-pre-wrap'>"
                "<div class='coverage-pre-label'>📄 .coverage 內容</div>"
                f"<pre class='coverage-pre'>{html.escape(item['raw_text'])}</pre>"
                "</div>"
            )

            out.append("</div></details>")  # species-body / species-card

        out.append("</div>")  # cat-block

    out.append("</section>")
    return "".join(out)


def render_blast_section(blast_data, final_dir: Path) -> str:
    out = ["<section id='blast'><h1>🔬 Blast virulence 結果</h1>"]
    out.append("<p class='muted'>依 rm_human 與 contigs 兩類整理，並顯示 blastn.best.tsv 與 species_counts.tsv。</p>")

    for label in ["rm_human", "contigs"]:
        bucket = blast_data.get(label, {})
        out.append(f"<div id='blast-{label}' class='blast-block'><div class='blast-label'>{html.escape(label)}</div>")

        for kind, title in [
            ("species_counts", "species_counts.tsv"),
            ("blastn_best", "blastn.best.tsv"),
        ]:
            files = bucket.get(kind, [])
            out.append(f"<h3 class='blast-sub'>{title}</h3>")
            if not files:
                out.append("<p class='muted'>沒有找到此類檔案。</p>")
                continue

            for fp in files:
                header, rows, total = parse_tsv(fp, max_rows=60)
                rel = file_link_text(fp, final_dir)
                if kind == "blastn_best":
                    # 預設收起，點標題列才展開
                    out.append("<details class='blast-card'>")
                    out.append(
                        f"<summary class='blast-card-summary'>"
                        f"<span class='blast-card-title'>📋 {html.escape(fp.name)}</span>"
                        f"<span class='muted'>資料列數：{total}</span>"
                        f"<a class='blast-card-link' href='{rel}' target='_blank' "
                        f"   onclick='event.stopPropagation()'>開啟原始檔 ↗</a>"
                        f"<span class='expand-icon'>▸</span>"
                        f"</summary>"
                    )
                    out.append("<div class='blast-card-body'>")
                    out.append(html_table(header, rows, max_preview=30))
                    if total > 30:
                        out.append(f"<p class='muted' style='margin-top:8px'>僅預覽前 30 筆，完整內容請開啟原始檔。</p>")
                    out.append("</div></details>")
                else:
                    out.append("<div class='card'>")
                    out.append(
                        f"<div class='card-meta'>"
                        f"<strong>{html.escape(fp.name)}</strong>"
                        f"<span class='muted'>　資料列數：{total}</span>"
                        f"　<a href='{rel}'>開啟原始檔 ↗</a>"
                        f"</div>"
                    )
                    out.append(html_table_with_auto_species_links(header, rows, max_preview=30))
                    if total > 30:
                        out.append(f"<p class='muted' style='margin-top:8px'>僅預覽前 30 筆，完整內容請開啟原始檔。</p>")
                    out.append("</div>")

        other = bucket.get("other_tsv", [])
        if other:
            out.append("<details class='raw-details'><summary>其他 TSV 檔案</summary><ul class='other-list'>")
            for fp in other:
                rel = file_link_text(fp, final_dir)
                out.append(f"<li><a href='{rel}'>{html.escape(fp.name)}</a></li>")
            out.append("</ul></details>")

        out.append("</div>")  # blast-block

    out.append("</section>")
    return "".join(out)


def build_sample_info(sid: str, final_dir: Path, fastq_path: Path = None, upset_data: dict = None) -> list[tuple[str, str]]:
    info = [
        ("Sample ID", sid),
        ("Report generated", datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        ("Final results", str(final_dir)),
    ]
    if fastq_path is not None:
        info.append(("Input FASTQ", fastq_path.name))
        info.append(("Input FASTQ path", str(fastq_path)))
    if upset_data and upset_data.get("fastq_total") is not None:
        info.append(("Input reads", f"{upset_data['fastq_total']:,}"))
    return info


def render_index(mapping_data, blast_data, sample_info) -> str:
    mapping_counts = {k: len(v) for k, v in mapping_data.items()}
    blast_counts = {
        k: {
            "blastn_best": len(v["blastn_best"]),
            "species_counts": len(v["species_counts"]),
        }
        for k, v in blast_data.items()
    }
    sample_info_html = "".join(
        f"""
        <div class="sample-info-item">
          <div class="sample-info-label">{html.escape(label)}</div>
          <div class="sample-info-value">{html.escape(value)}</div>
        </div>
        """
        for label, value in sample_info
    )

    return f"""
    <section class="hero">
      <div class="hero-title">
        <div>
          <h1 style="margin:0 0 6px">Pathogen Pipeline Summary</h1>
          <p style="margin:0;opacity:.8">此報告整理 <strong>Mapping</strong> 與 <strong>Blast virulence</strong> 結果。</p>
        </div>
      </div>
      <div class="sample-info">
        {sample_info_html}
      </div>
      <div class="grid">
        <a class="stat virus"    href="#mapping-Virus"><div class="label">🦠 Virus mapping</div><div class="value">{mapping_counts.get("Virus", 0)}</div></a>
        <a class="stat bacteria" href="#mapping-Bacteria"><div class="label">🧫 Bacteria mapping</div><div class="value">{mapping_counts.get("Bacteria", 0)}</div></a>
        <a class="stat parasite" href="#mapping-Parasite"><div class="label">🪱 Parasite mapping</div><div class="value">{mapping_counts.get("Parasite", 0)}</div></a>
        <a class="stat"          href="#blast-rm_human"><div class="label">rm_human blastn.best</div><div class="value">{blast_counts.get("rm_human", {}).get("blastn_best", 0)}</div></a>
        <a class="stat"          href="#blast-contigs"><div class="label">contigs blastn.best</div><div class="value">{blast_counts.get("contigs", {}).get("blastn_best", 0)}</div></a>
        <a class="stat"          href="#blast"><div class="label">species_counts files</div><div class="value">{blast_counts.get("rm_human", {}).get("species_counts", 0) + blast_counts.get("contigs", {}).get("species_counts", 0)}</div></a>
      </div>
      <nav class="topnav">
        <a href="#mapping">🧬 Mapping</a>
        <a href="#blast">🔬 Blast virulence</a>
      </nav>
    </section>
    """

def render_upset_section(upset_data) -> str:
    if not upset_data or not upset_data.get("sets"):
        return """<section id='mapping-upset'>
          <h1>📊 Mapping read overlap (UpSet Plot)</h1>
          <p class='muted'>沒有找到可用的 BAM 資料。</p>
        </section>"""

    sets       = upset_data["sets"]
    gm         = upset_data.get("genome_meta", {})
    fastq_total   = upset_data.get("fastq_total", 0)
    unmapped_count = upset_data.get("unmapped_count", 0)
    total_mapped  = upset_data.get("total_mapped_unique_reads", 0)

    from itertools import combinations

    labels = sorted(sets.keys(), key=lambda lbl: len(sets[lbl]), reverse=True)
    n = len(labels)

    # ── Compute all non-empty intersections ──────────────────────────
    # intersection_map: frozenset(indices) -> count
    # We want EXACT intersections (only those indices, not subsets).
    # i.e., reads that map to exactly this combination.
    read_to_genomes: dict = {}
    for idx, lbl in enumerate(labels):
        for rid in sets[lbl]:
            read_to_genomes.setdefault(rid, set()).add(idx)

    combo_counts: dict = {}  # frozenset of indices -> count
    for rid, idxs in read_to_genomes.items():
        key = frozenset(idxs)
        combo_counts[key] = combo_counts.get(key, 0) + 1

    # Per-genome total (may include overlap reads)
    genome_totals = [len(sets[lbl]) for lbl in labels]

    # Per-genome exclusive reads (mapped to this genome only)
    genome_exclusive = [
        combo_counts.get(frozenset([i]), 0)
        for i in range(n)
    ]

    # Category colours
    CAT_COLOR = {"Bacteria": "#43a047", "Parasite": "#fb8c00", "Virus": "#e91e63"}
    genome_colors = [
        CAT_COLOR.get(gm.get(lbl, {}).get("category", ""), "#5c7cfa")
        for lbl in labels
    ]

    # Build intersections list sorted by size desc
    intersections = []
    for combo, cnt in sorted(combo_counts.items(), key=lambda x: -x[1]):
        intersections.append({
            "indices": sorted(combo),
            "count": cnt,
        })

    upset_json = json.dumps({
        "labels":        labels,
        "colors":        genome_colors,
        "totals":        genome_totals,
        "exclusive":     genome_exclusive,
        "intersections": intersections,
        "fastq_total":   fastq_total,
        "unmapped":      unmapped_count,
        "total_mapped":  total_mapped,
    })

    if fastq_total > 0:
        stats_html = (
            f"Input fastq: <strong>{fastq_total:,}</strong> reads　"
            f"| Mapped (unique): <strong>{total_mapped:,}</strong>　"
            f"| Unmapped: <strong>{unmapped_count:,}</strong>"
        )
    else:
        stats_html = f"Total primary-mapped (unique): <strong>{total_mapped:,}</strong>"

    out = ["<section id='mapping-upset'><h1>📊 Mapping read overlap (UpSet Plot)</h1>"]
    out.append(
        "<p class='muted'>"
        "上方直條：每種 intersection 的 read 數（只計入 <em>恰好</em> 比對到該組合的 reads）。"
        "左側橫條：每個 genome 的 primary mapped reads 總數。"
        "中間點矩陣：實心點連線代表參與該 intersection 的 genome。"
        "<br>滑鼠移至直條可查看詳細資訊。"
        "</p>"
    )
    out.append(
        f"""<div class="card" style="overflow-x:auto;">
  <div id="upset-chart"></div>
  <script id="upset-data" type="application/json">{upset_json}</script>
  <p class="muted" style="margin-top:8px;text-align:center">{stats_html}</p>
</div>"""
    )
    out.append("</section>")
    return "".join(out)

_UPSET_JS = """

<script>
(function () {
  const tag = document.getElementById("upset-data");
  const wrap = document.getElementById("upset-chart");
  if (!tag || !wrap || typeof d3 === "undefined") return;

  const D = JSON.parse(tag.textContent);
  const labels     = D.labels;
  const colors     = D.colors;
  const totals     = D.totals;
  const inters     = D.intersections;
  const fastqTotal = D.fastq_total;
  const totalMapped = D.total_mapped;
  const n          = labels.length;

  // ── Layout: dynamically sized, never squished ──────────────────────
  // Left panel width is computed to fit the longest label + bar + text.
  // Minimum COL_W ensures bars are readable even with many intersections.
  const BAR_H   = 32;   // row height — taller for readability
  const DOT_R   = 10;
  const PAD_TOP = 190;  // top panel height (intersection bars + axis)
  const PAD_BOT = 28;
  const PAD_R   = 32;
  const COL_W   = 52;   // wider columns to avoid overlap

  // Left panel: label area (fixed) + bar area (fixed max width)
  const LABEL_W = 260;  // px reserved for genome label text
  const BAR_MAX = 160;  // max length of set-size bar
  const NUM_W   = 90;   // px for "N (xx.x%)" text to the left of bar
  const PAD_L   = NUM_W + BAR_MAX + LABEL_W + 12;  // total left width

  const totalW  = PAD_L + inters.length * COL_W + PAD_R;
  const totalH  = PAD_TOP + n * BAR_H + PAD_BOT + 24;

  // Make SVG scrollable horizontally without squishing
  d3.select(wrap).style("overflow-x","auto").style("overflow-y","visible");
  const svg = d3.select(wrap).append("svg")
    .attr("width", totalW)
    .attr("height", totalH)
    .style("font-family","'Segoe UI',Arial,sans-serif")
    .style("font-size","12px")
    .style("display","block");

  const maxInter  = d3.max(inters, d => d.count) || 1;
  const maxTotal  = d3.max(totals) || 1;
  const yBarScale = d3.scaleLinear().domain([0, maxInter]).range([PAD_TOP - 14, 10]);
  const xTotScale = d3.scaleLinear().domain([0, maxTotal]).range([BAR_MAX, 0]);

  // Denominator for percentage (prefer fastq total if available)
  const pctBase = fastqTotal > 0 ? fastqTotal : totalMapped;

  // ── Tooltip ──────────────────────────────────────────────────────
  const tip = d3.select("body").append("div")
    .style("position","absolute").style("display","none")
    .style("background","rgba(15,22,50,0.93)").style("color","#e5ecff")
    .style("padding","10px 14px").style("border-radius","10px")
    .style("font-size","13px").style("pointer-events","none")
    .style("max-width","380px").style("line-height","1.75")
    .style("box-shadow","0 4px 20px rgba(0,0,0,.35)").style("z-index","9999");

  const exclusive = D.exclusive;

  // ── Left panel: set-size bars ─────────────────────────────────────
  // Layout (left → right): [NUM_W: "N (xx%)"] [BAR_MAX: bar] [LABEL_W: label]
  // Scale is reversed (large→left, 0→right), so bar grows leftward from right edge of BAR_MAX zone.
  const lG = svg.append("g").attr("transform","translate(0,"+PAD_TOP+")");

  labels.forEach(function(lbl, i) {
    const y     = i * BAR_H + BAR_H / 2;
    const pct   = pctBase > 0 ? (totals[i] / pctBase * 100).toFixed(1) : "–";
    const excl  = exclusive[i];
    const share = totals[i] - excl;
    const exclPct  = pctBase > 0 ? (excl  / pctBase * 100).toFixed(1) : "–";

    // Bar right edge = NUM_W + BAR_MAX; grows leftward.
    // Full bar spans from xTotScale(totals[i]) to BAR_MAX (within zone offset by NUM_W).
    const fullLeft  = NUM_W + xTotScale(totals[i]);
    const fullRight = NUM_W + BAR_MAX;               // = xTotScale(0)
    const fullW     = fullRight - fullLeft;

    // Exclusive segment: leftmost portion (deeper colour)
    const exclLeft  = NUM_W + xTotScale(excl);
    const exclW     = fullRight - exclLeft;           // exclusive bar width

    // Shared segment fills the rest (from fullLeft to exclLeft)
    const sharedW   = exclLeft - fullLeft;

    // Shared portion (lighter = 40% opacity of base colour)
    if (sharedW > 0) {
      lG.append("rect")
        .attr("x", fullLeft).attr("y", y - BAR_H * 0.32)
        .attr("width", sharedW).attr("height", BAR_H * 0.64)
        .attr("fill", colors[i]).attr("rx", 3).attr("opacity", 0.35);
    }

    // Exclusive portion (full opacity, slightly darker feel)
    if (exclW > 0) {
      lG.append("rect")
        .attr("x", exclLeft).attr("y", y - BAR_H * 0.32)
        .attr("width", exclW).attr("height", BAR_H * 0.64)
        .attr("fill", colors[i]).attr("rx", 3).attr("opacity", 0.9);
    }

    // Invisible hover target over entire bar
    const hoverRect = lG.append("rect")
      .attr("x", fullLeft).attr("y", y - BAR_H * 0.5)
      .attr("width", Math.max(fullW, 4)).attr("height", BAR_H)
      .attr("fill", "transparent").style("cursor","pointer");

    const tipHtml =
      "<b>" + lbl + "</b><br>" +
      "Total mapped: <b>" + totals[i].toLocaleString() + "</b> (" + pct + "%)<br>" +
      "Exclusive: <b>" + excl.toLocaleString() + "</b> (" + exclPct + "%)<br>" +
      "Shared (overlap): <b>" + share.toLocaleString() + "</b>" +
      (pctBase > 0 ? " (" + (share/pctBase*100).toFixed(1) + "%)" : "") +
      "<br><span style='font-size:11px;color:#a0b8e0'>■ dark = exclusive &nbsp; □ light = shared</span>";

    hoverRect
      .on("mousemove", function(event) {
        tip.style("display","block")
           .style("left",(event.pageX+16)+"px").style("top",(event.pageY-12)+"px")
           .html(tipHtml);
      })
      .on("mouseleave", function(){ tip.style("display","none"); });

    // Count + percent (left of bar, right-aligned to NUM_W)
    lG.append("text")
      .attr("x", NUM_W - 6).attr("y", y + 4)
      .attr("text-anchor","end").attr("font-size","11px").attr("fill","#3a4f7a")
      .text(totals[i].toLocaleString() + " (" + pct + "%)");

    // Genome label (right of bar zone)
    const short = lbl.length > 36 ? lbl.slice(0,34) + "…" : lbl;
    const t = lG.append("text")
      .attr("x", NUM_W + BAR_MAX + 10).attr("y", y + 4)
      .attr("text-anchor","start").attr("font-size","12px")
      .attr("fill", colors[i]).attr("font-weight","600").text(short);
    t.append("title").text(lbl);
  });

  // Axis above bar zone (offset by NUM_W)
  const lAxis = d3.axisTop(xTotScale).ticks(4).tickFormat(d3.format(".2s"));
  lG.append("g")
    .attr("transform","translate("+NUM_W+",0)")
    .call(lAxis)
    .call(function(g){ g.select(".domain").remove(); })
    .call(function(g){ g.selectAll("text").style("fill","#5c6b8d").style("font-size","10px"); });

  // "Set size" label centred over bar zone
  lG.append("text")
    .attr("x", NUM_W + BAR_MAX / 2).attr("y", -36)
    .attr("text-anchor","middle").attr("font-size","11px").attr("fill","#5c6b8d")
    .text("Set size  (% of " + (fastqTotal > 0 ? "all reads" : "mapped reads") + ")");



  // ── Dot matrix ────────────────────────────────────────────────────
  const mG = svg.append("g").attr("transform","translate("+PAD_L+","+PAD_TOP+")");
  labels.forEach(function(_, i) {
    mG.append("line")
      .attr("x1", 0).attr("x2", inters.length * COL_W)
      .attr("y1", i*BAR_H + BAR_H/2).attr("y2", i*BAR_H + BAR_H/2)
      .attr("stroke","#e2e9f5").attr("stroke-width", 1);
  });
  inters.forEach(function(inter, ci) {
    const cx = ci * COL_W + COL_W/2;
    if (inter.indices.length > 1) {
      const ys = inter.indices.map(function(i){ return i*BAR_H + BAR_H/2; });
      mG.append("line")
        .attr("x1",cx).attr("x2",cx)
        .attr("y1",d3.min(ys)).attr("y2",d3.max(ys))
        .attr("stroke","#3a4f7a").attr("stroke-width",3)
        .attr("stroke-linecap","round");
    }
    labels.forEach(function(lbl, i) {
      const active = inter.indices.indexOf(i) >= 0;
      const circle = mG.append("circle")
        .attr("cx", cx).attr("cy", i*BAR_H + BAR_H/2)
        .attr("r", DOT_R)
        .attr("fill",   active ? colors[i] : "#dde4f0")
        .attr("stroke", active ? colors[i] : "#c5cfe0")
        .attr("stroke-width", 1.5)
        .style("cursor", active ? "pointer" : "default");

      if (active) {
        // Tooltip shows: this genome, all genomes in intersection, read count + pct
        const allNames = inter.indices.map(function(j){ return labels[j]; });
        const interLabel = allNames.length === 1
          ? "<b>" + lbl + "</b>"
          : allNames.map(function(n){ return n === lbl ? "<b><u>"+n+"</u></b>" : n; }).join(" ∩ ");
        const pct = pctBase > 0 ? " ("+(inter.count/pctBase*100).toFixed(2)+"%)" : "";
        const base = fastqTotal > 0 ? " of all reads" : " of mapped";
        const tipHtml = interLabel
          + "<br>Intersection reads: <b>" + inter.count.toLocaleString() + "</b>" + pct + base
          + (allNames.length > 1 ? "<br><span style='font-size:11px;color:#a0b8e0'>Underlined = this genome</span>" : "");

        circle
          .on("mousemove", function(event) {
            tip.style("display","block")
               .style("left",(event.pageX+16)+"px").style("top",(event.pageY-12)+"px")
               .html(tipHtml);
          })
          .on("mouseleave", function(){ tip.style("display","none"); })
          .on("mouseover", function(){
            d3.select(this).attr("r", DOT_R + 3).attr("stroke-width", 2.5);
          })
          .on("mouseout", function(){
            d3.select(this).attr("r", DOT_R).attr("stroke-width", 1.5);
          });
      }
    });
  });

  // ── Top bars: intersection sizes ──────────────────────────────────
  const tG = svg.append("g").attr("transform","translate("+PAD_L+",0)");
  inters.forEach(function(inter, ci) {
    const cx  = ci * COL_W + COL_W/2;
    const by  = yBarScale(inter.count);
    const bh  = PAD_TOP - 14 - by;
    const col = inter.indices.length === 1 ? colors[inter.indices[0]] : "#5c7cfa";

    const bar = tG.append("rect")
      .attr("x", cx - COL_W * 0.38).attr("y", by)
      .attr("width", COL_W * 0.76).attr("height", bh)
      .attr("fill", col).attr("rx", 3).attr("opacity", 0.85)
      .style("cursor","pointer");

    if (bh > 16) {
      tG.append("text")
        .attr("x", cx).attr("y", by - 3)
        .attr("text-anchor","middle").attr("font-size","10px").attr("fill","#3a4f7a")
        .text(inter.count.toLocaleString());
    }

    bar.on("mousemove", function(event) {
        const names = inter.indices.map(function(i){ return "<b>"+labels[i]+"</b>"; }).join(" ∩ ");
        const pct   = pctBase > 0 ? " ("+(inter.count/pctBase*100).toFixed(2)+"%)" : "";
        const base  = fastqTotal > 0 ? " of all reads" : " of mapped";
        tip.style("display","block")
           .style("left",(event.pageX+16)+"px").style("top",(event.pageY-12)+"px")
           .html(names+"<br>Reads: <b>"+inter.count.toLocaleString()+"</b>"+pct+base);
      })
      .on("mouseleave", function(){ tip.style("display","none"); })
      .on("click", function() {
        const sel = d3.select(this).classed("hi");
        tG.selectAll("rect").classed("hi",false).attr("opacity",0.85);
        if (!sel) d3.select(this).classed("hi",true).attr("opacity",1);
      });
  });

  const tAxis = d3.axisLeft(yBarScale).ticks(5).tickFormat(d3.format(".2s"));
  tG.append("g").call(tAxis)
    .call(function(g){ g.select(".domain").remove(); })
    .call(function(g){ g.selectAll("text").style("fill","#5c6b8d").style("font-size","10px"); });
  tG.append("text")
    .attr("transform","rotate(-90)").attr("x",-PAD_TOP/2).attr("y",-42)
    .attr("text-anchor","middle").attr("font-size","11px").attr("fill","#5c6b8d")
    .text("Intersection size");
})();
</script>
"""

def build_html(final_dir: Path, sid: str, mapping_reads_dir: Path = None, fastq_path: Path = None) -> str:
    mapping_data = collect_mapping(final_dir)
    blast_data   = collect_blast(final_dir)

    upset_data = None
    if mapping_reads_dir is not None:
        fastq_ids  = read_fastq_ids(fastq_path) if fastq_path else None
        upset_data = collect_mapping_bam_pie_data(mapping_reads_dir, fastq_ids=fastq_ids)

    sample_info = build_sample_info(sid, final_dir, fastq_path=fastq_path, upset_data=upset_data)

    body = []
    body.append(render_index(mapping_data, blast_data, sample_info))
    if upset_data is not None:
        body.append(render_upset_section(upset_data))
    body.append(render_mapping_section(mapping_data, final_dir))
    body.append(render_blast_section(blast_data, final_dir))

    upset_js_block = _UPSET_JS
    return f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(sid)} pipeline summary</title>
<script src="https://cdn.jsdelivr.net/npm/d3@7/dist/d3.min.js"></script>
<style>
/* ── Design tokens ── */
:root {{
  --bg: #0b1020;
  --panel: #141b34;
  --panel2: #1b2445;
  --text: #edf2ff;
  --muted: #b8c0e0;
  --line: #30406e;
  --accent: #8fb3ff;

  --radius-sm: 10px;
  --radius-md: 16px;
  --radius-lg: 20px;

  --shadow-sm: 0 2px 8px rgba(23,32,51,0.06);
  --shadow-md: 0 4px 20px rgba(23,32,51,0.10);
}}

/* ── Smooth scroll ── */
html {{ scroll-behavior: smooth; }}

/* ── Reset ── */
*, *::before, *::after {{ box-sizing: border-box; }}
body {{
  margin: 0;
  font-family: "Segoe UI", Arial, "Noto Sans TC", sans-serif;
  background: #eef2f9;
  color: #172033;
  line-height: 1.6;
  font-size: 14px;
}}
h1, h2, h3 {{ margin: 0 0 .5em; font-weight: 700; }}
a {{ color: #2457d6; text-decoration: none; }}
a:hover {{ text-decoration: underline; }}

/* ── Layout ── */
.wrapper {{
  max-width: 1320px;
  margin: 0 auto;
  padding: 28px 20px;
}}

/* ── Hero ── */
.hero {{
  background: linear-gradient(135deg, #0b1020 0%, #1b2a5c 60%, #0f3460 100%);
  color: var(--text);
  padding: 32px 28px 24px;
  border-radius: var(--radius-lg);
  margin-bottom: 28px;
  box-shadow: var(--shadow-md);
}}
.hero-title {{ margin-bottom: 20px; }}
.hero h1 {{ font-size: 22px; letter-spacing: .3px; }}

.sample-info {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
  gap: 10px;
  margin: 0 0 20px;
}}
.sample-info-item {{
  background: rgba(255,255,255,0.08);
  border: 1px solid rgba(255,255,255,0.12);
  padding: 12px 14px;
  border-radius: var(--radius-md);
}}
.sample-info-label {{
  font-size: 11px;
  color: #c5d0f0;
  margin-bottom: 4px;
  text-transform: uppercase;
  letter-spacing: .3px;
}}
.sample-info-value {{
  font-size: 14px;
  font-weight: 600;
  word-break: break-word;
}}

.grid {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(155px, 1fr));
  gap: 12px;
  margin: 0 0 20px;
}}
.stat {{
  background: rgba(255,255,255,0.08);
  border: 1px solid rgba(255,255,255,0.12);
  padding: 14px 16px;
  border-radius: var(--radius-md);
  backdrop-filter: blur(6px);
  transition: background .2s;
}}
.stat:hover {{ background: rgba(255,255,255,0.13); }}
a.stat {{
  display: block;
  text-decoration: none;
  color: inherit;
  cursor: pointer;
  outline: none;
}}
a.stat:hover {{
  background: rgba(255,255,255,0.18);
  transform: translateY(-2px);
  box-shadow: 0 6px 20px rgba(0,0,0,0.2);
}}
a.stat {{ transition: background .2s, transform .15s, box-shadow .15s; }}
.stat.bacteria {{ border-left: 3px solid #66bb6a; }}
.stat.parasite  {{ border-left: 3px solid #ffa726; }}
.stat.virus     {{ border-left: 3px solid #ef5350; }}
.stat .label {{ font-size: 11px; color: #c5d0f0; margin-bottom: 4px; }}
.stat .value {{ font-size: 30px; font-weight: 800; letter-spacing: -1px; }}

.topnav {{ display: flex; gap: 6px; flex-wrap: wrap; }}
.topnav a {{
  color: white;
  background: rgba(255,255,255,0.12);
  border: 1px solid rgba(255,255,255,0.2);
  padding: 7px 18px;
  border-radius: 999px;
  font-size: 13px;
  font-weight: 600;
  transition: background .2s;
}}
.topnav a:hover {{ background: rgba(255,255,255,0.22); text-decoration: none; }}

/* ── Section titles ── */
section {{ margin-bottom: 36px; }}
section > h1 {{
  font-size: 18px;
  border-left: 4px solid #4a7aff;
  padding-left: 12px;
  margin-bottom: 16px;
  color: #172033;
}}

/* ── Category block (Bacteria / Parasite / Virus) ── */
.cat-block {{
  background: white;
  border: 1px solid #dce5f5;
  border-radius: var(--radius-lg);
  padding: 20px 22px 16px;
  margin-bottom: 20px;
  box-shadow: var(--shadow-sm);
  scroll-margin-top: 20px;
}}
.cat-header {{
  display: flex;
  align-items: center;
  gap: 12px;
  margin-bottom: 16px;
}}
.cat-badge {{
  background: var(--cat-bg, #eef4ff);
  color: var(--cat-fg, #172033);
  border: 1.5px solid var(--cat-accent, #4a7aff);
  font-weight: 700;
  font-size: 15px;
  padding: 4px 16px;
  border-radius: 999px;
}}
.cat-count {{
  font-size: 13px;
  color: #5c6b8d;
}}

/* ── Summary table ── */
.summary-table-wrap {{
  margin-bottom: 18px;
}}

/* ── Species collapsible card ── */
.species-card {{
  border: 1px solid #e2e9f5;
  border-radius: var(--radius-md);
  margin-bottom: 10px;
  overflow: hidden;
  background: #f8faff;
  transition: box-shadow .2s;
}}
.species-card[open] {{
  box-shadow: var(--shadow-md);
  background: white;
}}
.species-summary {{
  display: flex;
  align-items: center;
  flex-wrap: wrap;
  gap: 10px;
  padding: 13px 16px;
  cursor: pointer;
  list-style: none;
  user-select: none;
  background: inherit;
  transition: background .15s;
}}
.species-summary::-webkit-details-marker {{ display: none; }}
.species-summary:hover {{ background: #eef4ff; }}
.species-card[open] > .species-summary {{
  background: #eef4ff;
  border-bottom: 1px solid #d8e4fa;
}}

.species-name {{
  font-weight: 700;
  font-size: 14px;
  flex: 0 0 auto;
}}
.species-link {{
  color: #1a44b8;
  text-decoration: none;
  border-bottom: 1px dashed #8fb3ff;
}}
.species-link:hover {{ border-bottom-style: solid; }}

.summary-pills {{
  display: flex;
  flex-wrap: wrap;
  gap: 6px;
  flex: 1;
}}
.pill {{
  display: inline-flex;
  align-items: center;
  border-radius: 999px;
  padding: 3px 11px;
  font-size: 12px;
  font-weight: 600;
  white-space: nowrap;
}}
.pill-reads   {{ background: #e8f0fe; color: #1a44b8; }}
.pill-covered {{ background: #e6f4ea; color: #1e6823; }}
.pill-mean    {{ background: #fff3e0; color: #8c4a00; }}

.expand-icon {{
  margin-left: auto;
  font-size: 13px;
  color: #8fa3c8;
  transition: transform .2s;
  flex: 0 0 auto;
}}
.species-card[open] .expand-icon {{ transform: rotate(90deg); }}

.species-body {{
  padding: 16px 18px;
}}

/* ── Plot image ── */
.plot {{
  width: 100%;
  max-width: 980px;
  display: block;
  border: 1px solid #d8e0ef;
  border-radius: var(--radius-md);
  margin-bottom: 12px;
}}

/* ── Coverage pre block ── */
.coverage-pre-wrap {{
  margin-top: 12px;
  border: 1px solid #2a3a5e;
  border-radius: var(--radius-sm);
  overflow: hidden;
}}
.coverage-pre-label {{
  background: #1a2540;
  color: #8fb3ff;
  font-size: 11px;
  font-weight: 700;
  letter-spacing: .6px;
  text-transform: uppercase;
  padding: 6px 14px;
}}
.coverage-pre {{
  margin: 0;
  border-radius: 0;
  border-top: 1px solid #2a3a5e;
}}

pre {{
  white-space: pre-wrap;
  word-break: break-word;
  background: #0f172a;
  color: #e5ecff;
  padding: 16px;
  border-radius: 0;
  overflow-x: auto;
  font-size: 12.5px;
  line-height: 1.55;
  margin: 0;
}}

/* ── Table ── */
.table-wrap {{
  overflow-x: auto;
  background: white;
  border: 1px solid #d8e0ef;
  border-radius: var(--radius-md);
}}
table {{
  border-collapse: collapse;
  width: 100%;
  min-width: 680px;
}}
th, td {{
  border-bottom: 1px solid #e8eef8;
  padding: 9px 12px;
  text-align: left;
  vertical-align: top;
  font-size: 13px;
}}
th {{
  background: #eef4ff;
  font-size: 12px;
  color: #3a4f7a;
  letter-spacing: .3px;
  text-transform: uppercase;
  position: sticky;
  top: 0;
}}
tr:last-child td {{ border-bottom: none; }}
tr:hover td {{ background: #f7faff; }}
a.species-link {{ font-weight: 600; }}

/* ── Blast section ── */
.blast-block {{
  background: white;
  border: 1px solid #dce5f5;
  border-radius: var(--radius-lg);
  padding: 20px 22px 16px;
  margin-bottom: 20px;
  box-shadow: var(--shadow-sm);
  scroll-margin-top: 20px;
}}
.blast-label {{
  display: inline-block;
  background: #e8edfa;
  color: #223388;
  font-weight: 700;
  font-size: 14px;
  padding: 4px 16px;
  border-radius: 999px;
  margin-bottom: 14px;
}}
.blast-sub {{
  font-size: 14px;
  color: #3a4f7a;
  margin: 16px 0 10px;
  padding-left: 10px;
  border-left: 3px solid #b0c4f0;
}}

/* ── Blastn best collapsible card ── */
.blast-card {{
  border: 1px solid #d0daf0;
  border-radius: var(--radius-md);
  margin-bottom: 12px;
  overflow: hidden;
  background: #f8faff;
  transition: box-shadow .2s;
}}
.blast-card[open] {{
  box-shadow: var(--shadow-md);
  background: white;
}}
.blast-card-summary {{
  display: flex;
  align-items: center;
  gap: 12px;
  flex-wrap: wrap;
  padding: 12px 16px;
  cursor: pointer;
  list-style: none;
  user-select: none;
  transition: background .15s;
}}
.blast-card-summary::-webkit-details-marker {{ display: none; }}
.blast-card-summary:hover {{ background: #e8eef8; }}
.blast-card[open] > .blast-card-summary {{
  background: #eef4ff;
  border-bottom: 1px solid #d0daf0;
}}
.blast-card-title {{
  font-weight: 700;
  font-size: 13px;
  color: #1a2f6a;
  flex: 1;
}}
.blast-card-link {{
  font-size: 12px;
  color: #2457d6;
  white-space: nowrap;
}}
.blast-card-body {{
  padding: 14px 16px;
}}

/* ── Generic card ── */
.card {{
  background: #f8faff;
  border: 1px solid #dce5f5;
  border-radius: var(--radius-md);
  padding: 16px 18px;
  margin-bottom: 14px;
}}
.card-meta {{
  margin-bottom: 10px;
  font-size: 13px;
}}

/* ── Other list ── */
.other-list {{ padding-left: 20px; margin: 10px 0 4px; }}
.other-list li {{ margin-bottom: 4px; font-size: 13px; }}

/* ── Muted ── */
.muted {{ color: #5c6b8d; font-size: 13px; }}
</style>
</head>
<body>
<div class="wrapper">
{''.join(body)}
</div>

{upset_js_block}
</body>
</html>"""


def main():
    ap = argparse.ArgumentParser(description="Generate summary HTML for pathogen pipeline final_results.")
    ap.add_argument("--final_dir",        required=True,  help="final_results directory")
    ap.add_argument("--sid",              required=True,  help="sample id")
    ap.add_argument("--out",              default=None,   help="output html path")
    ap.add_argument("--mapping_reads_dir",default=None,
                    help="directory with Bacteria/Parasite/Virus sub-folders containing sorted BAM files")
    ap.add_argument("--fastq",            default=None,
                    help="input fastq or fastq.gz; provides read universe for UpSet plot")
    args = ap.parse_args()

    final_dir = Path(args.final_dir).resolve()
    if not final_dir.exists():
        raise FileNotFoundError(f"final_dir does not exist: {final_dir}")

    out = Path(args.out) if args.out else (final_dir / f"{args.sid}.summary_report.html")

    mapping_reads_dir = Path(args.mapping_reads_dir).resolve() if args.mapping_reads_dir else None
    if mapping_reads_dir and not mapping_reads_dir.exists():
        raise FileNotFoundError(f"mapping_reads_dir does not exist: {mapping_reads_dir}")

    fastq_path = Path(args.fastq).resolve() if args.fastq else None
    if fastq_path and not fastq_path.exists():
        raise FileNotFoundError(f"fastq file does not exist: {fastq_path}")

    html_text = build_html(final_dir, args.sid,
                           mapping_reads_dir=mapping_reads_dir,
                           fastq_path=fastq_path)
    out.write_text(html_text, encoding="utf-8")
    print(f"[OK] summary html written to: {out}")


if __name__ == "__main__":
    main()
