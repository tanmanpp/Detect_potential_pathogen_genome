#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import html
import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


MAJOR_RANKS = [
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

DISPLAY_RANK_PRIORITY = {r: i for i, r in enumerate(MAJOR_RANKS)}
DISPLAY_RANK_PRIORITY["unclassified"] = 999
DISPLAY_RANK_PRIORITY["no rank"] = 998
DISPLAY_RANK_PRIORITY["root"] = -1


def open_text_auto(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def sanitize_seq_id(seq_id: str) -> str:
    return seq_id.strip().split()[0]


def parse_kraken_output(kraken_file):
    """
    Kraken2 standard output format (common):
    status \t seq_id \t taxid \t length \t lca_info

    Returns a list of dicts.
    """
    records = []

    with open_text_auto(kraken_file) as f:
        for lineno, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                print(f"[warning] Skip line {lineno}: fewer than 4 columns", file=sys.stderr)
                continue

            status = parts[0].strip()
            seq_id = sanitize_seq_id(parts[1])
            taxid = str(parts[2]).strip()

            raw_len = parts[3].strip()
            if "|" in raw_len:
                try:
                    length = sum(int(x) for x in raw_len.split("|") if x.isdigit())
                except Exception:
                    length = raw_len
            else:
                try:
                    length = int(raw_len)
                except Exception:
                    length = raw_len

            records.append({
                "status": status,
                "seq_id": seq_id,
                "taxid": taxid,
                "length": length,
            })

    return records


def load_taxonomy(nodes_dmp, names_dmp):
    parent = {}
    rank_map = {}
    name_map = {}

    with open_text_auto(nodes_dmp) as f:
        for line in f:
            parts = [x.strip() for x in line.split("|")]
            if len(parts) < 3:
                continue
            taxid = parts[0]
            parent_taxid = parts[1]
            tax_rank = parts[2]
            parent[taxid] = parent_taxid
            rank_map[taxid] = tax_rank

    with open_text_auto(names_dmp) as f:
        for line in f:
            parts = [x.strip() for x in line.split("|")]
            if len(parts) < 4:
                continue
            taxid = parts[0]
            name_txt = parts[1]
            name_class = parts[3]
            if name_class == "scientific name":
                name_map[taxid] = name_txt

    return parent, rank_map, name_map


def lineage_of_taxid(taxid, parent):
    out = []
    seen = set()
    cur = str(taxid)

    while cur and cur not in seen:
        seen.add(cur)
        out.append(cur)
        p = parent.get(cur)
        if p is None or p == cur:
            break
        cur = p

    out.reverse()
    return out


def new_tree_node(name="root", taxid="0", rank="root"):
    return {
        "name": name,
        "taxid": taxid,
        "rank": rank,
        "children": {},
        "seqs": [],
        "subtree_count": 0,
        "direct_count": 0,
        "pct_of_all": 0.0,
        "pct_of_classified": 0.0,
        "path_names": "",
        "path_taxids": "",
    }


def build_tree(records, parent, rank_map, name_map):
    root = new_tree_node()

    for rec in records:
        seq_id = rec["seq_id"]
        taxid = rec["taxid"]
        status = rec["status"]
        length = rec["length"]

        if status != "C" or taxid in ("0", "", "unclassified"):
            if "unclassified" not in root["children"]:
                root["children"]["unclassified"] = new_tree_node(
                    name="Unclassified",
                    taxid="0",
                    rank="unclassified"
                )
            root["children"]["unclassified"]["seqs"].append({
                "seq_id": seq_id,
                "length": length,
                "taxid": taxid,
                "status": status,
            })
            continue

        filtered = extract_major_lineage_with_domain(
            taxid=taxid,
            parent=parent,
            rank_map=rank_map,
            name_map=name_map
        )

        node = root
        for tid, n, r in filtered:
            if tid not in node["children"]:
                node["children"][tid] = new_tree_node(name=n, taxid=tid, rank=r)
            node = node["children"][tid]

        node["seqs"].append({
            "seq_id": seq_id,
            "length": length,
            "taxid": taxid,
            "status": status,
        })

    return root


def compute_metrics(node, total_all, total_classified, parent_path_names=None, parent_path_taxids=None):
    if parent_path_names is None:
        parent_path_names = []
    if parent_path_taxids is None:
        parent_path_taxids = []

    current_path_names = parent_path_names.copy()
    current_path_taxids = parent_path_taxids.copy()

    if node["rank"] != "root":
        current_path_names.append(node["name"])
        current_path_taxids.append(node["taxid"])

    node["path_names"] = " > ".join(current_path_names)
    node["path_taxids"] = " > ".join(current_path_taxids)
    node["direct_count"] = len(node["seqs"])

    total = node["direct_count"]
    for child in node["children"].values():
        total += compute_metrics(
            child,
            total_all=total_all,
            total_classified=total_classified,
            parent_path_names=current_path_names,
            parent_path_taxids=current_path_taxids,
        )

    node["subtree_count"] = total
    node["pct_of_all"] = (100.0 * total / total_all) if total_all else 0.0

    if node["rank"] == "unclassified":
        node["pct_of_classified"] = 0.0
    else:
        node["pct_of_classified"] = (100.0 * total / total_classified) if total_classified else 0.0

    return total


def sort_tree(node):
    for child in node["children"].values():
        sort_tree(child)

    node["children"] = dict(sorted(
        node["children"].items(),
        key=lambda kv: (
            DISPLAY_RANK_PRIORITY.get(kv[1]["rank"], 500),
            -kv[1]["subtree_count"],
            kv[1]["name"].lower()
        )
    ))

    node["seqs"].sort(
        key=lambda x: (
            -(x["length"] if isinstance(x["length"], int) else -1),
            x["seq_id"]
        )
    )


def flatten_tree(node, rows=None, depth=0, order_counter=None):
    if rows is None:
        rows = []
    if order_counter is None:
        order_counter = {"n": 0}

    if node["rank"] != "root":
        order_counter["n"] += 1
        rows.append({
            "tree_order": order_counter["n"],
            "depth": depth,
            "taxid": node["taxid"],
            "name": node["name"],
            "rank": node["rank"],
            "path_names": node["path_names"],
            "path_taxids": node["path_taxids"],
            "direct_count": node["direct_count"],
            "subtree_count": node["subtree_count"],
            "pct_of_all": round(node["pct_of_all"], 6),
            "pct_of_classified": round(node["pct_of_classified"], 6),
        })

    for child in node["children"].values():
        flatten_tree(child, rows, depth + 1, order_counter)

    return rows


def collect_seq_assignments(records, parent, rank_map, name_map):
    rows = []

    for rec in records:
        seq_id = rec["seq_id"]
        taxid = rec["taxid"]
        status = rec["status"]
        length = rec["length"]

        if status != "C" or taxid in ("0", "", "unclassified"):
            rows.append({
                "seq_id": seq_id,
                "status": status,
                "length": length,
                "assigned_taxid": taxid,
                "assigned_name": "Unclassified",
                "assigned_rank": "unclassified",
                "path_names": "Unclassified",
                "path_taxids": "0",
            })
            continue

        filtered = extract_major_lineage_with_domain(
            taxid=taxid,
            parent=parent,
            rank_map=rank_map,
            name_map=name_map
        )

        final_name = name_map.get(taxid, taxid)
        final_rank = rank_map.get(taxid, "no rank")

        filtered_names = [x[1] for x in filtered]
        filtered_taxids = [x[0] for x in filtered]

        rows.append({
            "seq_id": seq_id,
            "status": status,
            "length": length,
            "assigned_taxid": taxid,
            "assigned_name": final_name,
            "assigned_rank": final_rank,
            "path_names": " > ".join(filtered_names),
            "path_taxids": " > ".join(filtered_taxids),
        })

    return rows


def tree_to_json_for_html(node):
    return {
        "name": node["name"],
        "taxid": node["taxid"],
        "rank": node["rank"],
        "subtree_count": node["subtree_count"],
        "direct_count": node["direct_count"],
        "pct_of_all": round(node["pct_of_all"], 4),
        "pct_of_classified": round(node["pct_of_classified"], 4),
        "children": [tree_to_json_for_html(ch) for ch in node["children"].values()],
        "seqs": [
            {
                "seq_id": x["seq_id"],
                "length": x["length"],
                "taxid": x["taxid"],
                "status": x["status"],
            }
            for x in node["seqs"]
        ]
    }

def extract_major_lineage_with_domain(taxid, parent, rank_map, name_map):
    """
    回傳固定從 superkingdom(domain) 開始的 major lineage。
    即使中間某些 rank 缺失，也至少會先有 domain。
    回傳格式:
      [(taxid, name, rank), ...]
    """
    lineage = lineage_of_taxid(taxid, parent)

    domain_node = None
    major_nodes = []

    for tid in lineage:
        r = rank_map.get(tid, "no rank")
        n = name_map.get(tid, tid)

        if r == "superkingdom":
            domain_node = (tid, n, "superkingdom")

        if r in MAJOR_RANKS:
            major_nodes.append((tid, n, r))

    # 若有找到 domain，強制放最前面
    if domain_node is not None:
        major_nodes_no_domain = [x for x in major_nodes if x[2] != "superkingdom"]
        return [domain_node] + major_nodes_no_domain

    # 若完全找不到 superkingdom，退而求其次：
    # 仍然回傳現有 major lineage；若連 major 都沒有，就放自己
    if major_nodes:
        return major_nodes

    return [(taxid, name_map.get(taxid, taxid), rank_map.get(taxid, "no rank"))]

def generate_html(tree, total_all, total_classified, total_unclassified, html_out):
    tree_json = tree_to_json_for_html(tree)
    tree_json_str = json.dumps(tree_json, ensure_ascii=False)

    pct_classified = (100.0 * total_classified / total_all) if total_all else 0.0
    pct_unclassified = (100.0 * total_unclassified / total_all) if total_all else 0.0

    html_text = f"""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Kraken2 Taxonomy Browser</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<style>
:root {{
  --bg: #0d1117;
  --bg2: #161b22;
  --bg3: #1c2128;
  --border: #30363d;
  --border2: #21262d;
  --text: #e6edf3;
  --text2: #8b949e;
  --text3: #6e7681;
  --accent: #58a6ff;
  --accent2: #3fb950;
  --accent3: #d2a8ff;
  --warn: #f85149;
  --warn2: #ff7b72;
  --gold: #e3b341;
  --teal: #39d353;
  --rank-superkingdom: #ff7b72;
  --rank-kingdom: #ffa657;
  --rank-phylum: #e3b341;
  --rank-class: #3fb950;
  --rank-order: #39c5cf;
  --rank-family: #58a6ff;
  --rank-genus: #bc8cff;
  --rank-species: #ff79c6;
  --rank-other: #8b949e;
  --radius: 8px;
  --radius-lg: 12px;
}}
* {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{
  font-family: 'Inter', "Microsoft JhengHei", sans-serif;
  background: var(--bg);
  color: var(--text);
  min-height: 100vh;
  font-size: 14px;
  line-height: 1.6;
}}
/* ── header ── */
.site-header {{
  background: var(--bg2);
  border-bottom: 1px solid var(--border);
  padding: 16px 28px;
  display: flex;
  align-items: center;
  gap: 14px;
  position: sticky;
  top: 0;
  z-index: 100;
  backdrop-filter: blur(8px);
}}
.site-header .logo {{
  font-family: 'JetBrains Mono', monospace;
  font-size: 18px;
  font-weight: 600;
  color: var(--accent);
  letter-spacing: -0.5px;
}}
.site-header .logo span {{ color: var(--text2); font-weight: 400; }}
.site-header .badge {{
  font-family: 'JetBrains Mono', monospace;
  font-size: 11px;
  background: var(--bg3);
  border: 1px solid var(--border);
  color: var(--text2);
  padding: 2px 8px;
  border-radius: 20px;
}}
.main-wrap {{
  max-width: 1400px;
  margin: 0 auto;
  padding: 24px 28px;
}}
/* ── cards ── */
.card {{
  background: var(--bg2);
  border: 1px solid var(--border);
  border-radius: var(--radius-lg);
  padding: 20px 24px;
  margin-bottom: 16px;
}}
.card-title {{
  font-size: 11px;
  font-weight: 600;
  text-transform: uppercase;
  letter-spacing: 1px;
  color: var(--text3);
  margin-bottom: 14px;
}}
/* ── stats ── */
.stats-grid {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
  gap: 12px;
}}
.stat-card {{
  background: var(--bg3);
  border: 1px solid var(--border2);
  border-radius: var(--radius);
  padding: 16px 20px;
  position: relative;
  overflow: hidden;
}}
.stat-card::before {{
  content: '';
  position: absolute;
  top: 0; left: 0; right: 0;
  height: 2px;
}}
.stat-card.total::before {{ background: var(--accent); }}
.stat-card.classified::before {{ background: var(--accent2); }}
.stat-card.unclassified::before {{ background: var(--warn); }}
.stat-label {{
  font-size: 11px;
  color: var(--text3);
  text-transform: uppercase;
  letter-spacing: 0.8px;
  font-weight: 600;
  margin-bottom: 6px;
}}
.stat-val {{
  font-family: 'JetBrains Mono', monospace;
  font-size: 26px;
  font-weight: 600;
  color: var(--text);
  line-height: 1.1;
}}
.stat-val.total {{ color: var(--accent); }}
.stat-val.classified {{ color: var(--accent2); }}
.stat-val.unclassified {{ color: var(--warn); }}
.stat-sub {{
  font-size: 12px;
  color: var(--text3);
  margin-top: 3px;
  font-family: 'JetBrains Mono', monospace;
}}
/* ── controls ── */
.controls-row {{
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
  align-items: center;
}}
.file-label {{
  display: inline-flex;
  align-items: center;
  gap: 8px;
  background: var(--bg3);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: 7px 14px;
  cursor: pointer;
  font-size: 13px;
  color: var(--text2);
  transition: border-color 0.15s, color 0.15s;
}}
.file-label:hover {{ border-color: var(--accent); color: var(--accent); }}
.file-label svg {{ width: 15px; height: 15px; flex-shrink: 0; }}
#seqFileInput {{ display: none; }}
.search-wrap {{
  display: flex;
  align-items: center;
  background: var(--bg3);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: 0 10px;
  gap: 8px;
  flex: 1;
  min-width: 200px;
  max-width: 380px;
  transition: border-color 0.15s;
}}
.search-wrap:focus-within {{ border-color: var(--accent); }}
.search-wrap svg {{ width: 15px; height: 15px; color: var(--text3); flex-shrink: 0; }}
#searchBox {{
  background: none;
  border: none;
  outline: none;
  color: var(--text);
  font-size: 13px;
  padding: 8px 0;
  width: 100%;
  font-family: 'Inter', sans-serif;
}}
#searchBox::placeholder {{ color: var(--text3); }}
.btn {{
  display: inline-flex;
  align-items: center;
  gap: 6px;
  padding: 8px 14px;
  border-radius: var(--radius);
  border: 1px solid var(--border);
  background: var(--bg3);
  color: var(--text2);
  cursor: pointer;
  font-size: 13px;
  font-family: 'Inter', sans-serif;
  transition: background 0.15s, border-color 0.15s, color 0.15s;
}}
.btn:hover {{ background: var(--bg); border-color: var(--accent); color: var(--accent); }}
.btn svg {{ width: 14px; height: 14px; }}
.file-status {{
  margin-top: 12px;
  font-size: 13px;
  color: var(--text3);
  display: flex;
  align-items: center;
  gap: 8px;
}}
.ok {{ color: var(--accent2); font-weight: 500; }}
.warn {{ color: var(--warn); font-weight: 500; }}
/* ── tree ── */
#taxonomyTree {{
  padding: 4px 0;
}}
details.tax-node {{
  margin-left: 20px;
  border-left: 1px solid var(--border2);
  margin-top: 2px;
  margin-bottom: 2px;
  transition: border-color 0.2s;
}}
details.tax-node[open] > summary {{
  background: rgba(88, 166, 255, 0.05);
  border-radius: var(--radius) var(--radius) 0 0;
}}
details.tax-node:hover > summary {{
  background: rgba(255,255,255,0.03);
}}
summary {{
  cursor: pointer;
  padding: 6px 10px 6px 14px;
  list-style: none;
  display: flex;
  align-items: center;
  flex-wrap: wrap;
  gap: 4px;
  border-radius: var(--radius);
  user-select: none;
  position: relative;
}}
summary::-webkit-details-marker {{ display: none; }}
summary::before {{
  content: '▶';
  font-size: 9px;
  color: var(--text3);
  margin-right: 4px;
  transition: transform 0.15s;
  flex-shrink: 0;
}}
details[open] > summary::before {{ transform: rotate(90deg); }}
.tax-name {{
  font-weight: 600;
  color: var(--text);
  font-size: 14px;
  text-decoration: none;
  transition: color 0.15s;
}}
a.tax-name:hover {{
  color: var(--accent);
  text-decoration: underline;
  text-underline-offset: 3px;
}}
.rank-badge {{
  font-family: 'JetBrains Mono', monospace;
  font-size: 10px;
  padding: 1px 6px;
  border-radius: 4px;
  font-weight: 600;
  text-transform: uppercase;
  letter-spacing: 0.5px;
}}
.rank-superkingdom {{ background: rgba(255,123,114,0.15); color: var(--rank-superkingdom); }}
.rank-kingdom {{ background: rgba(255,166,87,0.15); color: var(--rank-kingdom); }}
.rank-phylum {{ background: rgba(227,179,65,0.15); color: var(--rank-phylum); }}
.rank-class {{ background: rgba(63,185,80,0.15); color: var(--rank-class); }}
.rank-order {{ background: rgba(57,197,207,0.15); color: var(--rank-order); }}
.rank-family {{ background: rgba(88,166,255,0.15); color: var(--rank-family); }}
.rank-genus {{ background: rgba(188,140,255,0.15); color: var(--rank-genus); }}
.rank-species {{ background: rgba(255,121,198,0.15); color: var(--rank-species); }}
.rank-unclassified {{ background: rgba(139,148,158,0.15); color: var(--rank-other); }}
.rank-other {{ background: rgba(139,148,158,0.1); color: var(--rank-other); }}
.tax-stats {{
  display: flex;
  align-items: center;
  gap: 10px;
  margin-left: 4px;
  flex-wrap: wrap;
}}
.stat-chip {{
  font-family: 'JetBrains Mono', monospace;
  font-size: 11px;
  padding: 1px 7px;
  border-radius: 4px;
  white-space: nowrap;
}}
.chip-count {{
  background: rgba(248,81,73,0.12);
  color: #ff7b72;
  border: 1px solid rgba(248,81,73,0.2);
}}
.chip-pct-all {{
  background: rgba(88,166,255,0.1);
  color: #79c0ff;
  border: 1px solid rgba(88,166,255,0.2);
}}
.chip-pct-cls {{
  background: rgba(63,185,80,0.1);
  color: #56d364;
  border: 1px solid rgba(63,185,80,0.2);
}}
.chip-direct {{
  background: rgba(188,140,255,0.1);
  color: #d2a8ff;
  border: 1px solid rgba(188,140,255,0.15);
}}
.tax-taxid {{
  font-family: 'JetBrains Mono', monospace;
  font-size: 10px;
  color: var(--text3);
  margin-left: 2px;
}}
/* ── seq list ── */
.seq-section {{
  margin: 8px 12px 10px 28px;
  padding: 10px 14px;
  background: var(--bg3);
  border: 1px solid var(--border2);
  border-radius: var(--radius);
}}
.seq-section-header {{
  font-size: 11px;
  font-weight: 600;
  color: var(--text3);
  text-transform: uppercase;
  letter-spacing: 0.8px;
  margin-bottom: 8px;
  padding-bottom: 6px;
  border-bottom: 1px solid var(--border2);
  display: flex;
  align-items: center;
  gap: 6px;
}}
.seq-section-header .seq-count-badge {{
  font-family: 'JetBrains Mono', monospace;
  background: rgba(88,166,255,0.1);
  color: var(--accent);
  border: 1px solid rgba(88,166,255,0.2);
  padding: 1px 6px;
  border-radius: 4px;
  font-size: 10px;
}}
.seq-list {{
  list-style: none;
  display: flex;
  flex-direction: column;
  gap: 3px;
}}
.seq-list li {{
  display: flex;
  align-items: center;
  gap: 0;
  padding: 4px 8px;
  border-radius: 6px;
  background: transparent;
  transition: background 0.1s;
}}
.seq-list li:hover {{ background: rgba(255,255,255,0.04); }}
.seq-idx {{
  font-family: 'JetBrains Mono', monospace;
  font-size: 10px;
  color: var(--text3);
  width: 32px;
  flex-shrink: 0;
  text-align: right;
  margin-right: 10px;
}}
.seq-id-link {{
  font-family: 'JetBrains Mono', monospace;
  font-size: 12px;
  color: var(--accent);
  text-decoration: none;
  font-weight: 500;
  flex: 1;
  min-width: 0;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}}
.seq-id-link:hover {{ color: #a5c8ff; text-decoration: underline; }}
.seq-len {{
  font-family: 'JetBrains Mono', monospace;
  font-size: 11px;
  color: var(--text3);
  background: var(--bg2);
  border: 1px solid var(--border2);
  padding: 1px 6px;
  border-radius: 4px;
  margin-left: 8px;
  white-space: nowrap;
  flex-shrink: 0;
}}
.seq-taxid {{
  font-family: 'JetBrains Mono', monospace;
  font-size: 10px;
  color: var(--text3);
  margin-left: 6px;
  white-space: nowrap;
  flex-shrink: 0;
}}
/* ── more-seq ── */
details.more-seq {{
  margin-top: 6px;
}}
details.more-seq > summary {{
  cursor: pointer;
  display: inline-flex;
  align-items: center;
  gap: 6px;
  font-size: 12px;
  color: var(--accent);
  padding: 4px 8px;
  border-radius: 6px;
  border: 1px dashed rgba(88,166,255,0.3);
  list-style: none;
  user-select: none;
  transition: background 0.15s, border-color 0.15s;
  background: rgba(88,166,255,0.05);
}}
details.more-seq > summary::-webkit-details-marker {{ display: none; }}
details.more-seq > summary::before {{
  content: '▶';
  font-size: 9px;
  transition: transform 0.15s;
}}
details.more-seq[open] > summary::before {{ transform: rotate(90deg); }}
details.more-seq > summary:hover {{
  background: rgba(88,166,255,0.1);
  border-color: rgba(88,166,255,0.5);
}}
details.more-seq > .seq-list {{
  margin-top: 4px;
  border-top: 1px solid var(--border2);
  padding-top: 4px;
}}
/* ── note box ── */
.info-box {{
  display: flex;
  gap: 12px;
  background: rgba(88,166,255,0.07);
  border: 1px solid rgba(88,166,255,0.2);
  border-radius: var(--radius);
  padding: 12px 16px;
  font-size: 13px;
  color: var(--text2);
  line-height: 1.7;
}}
.info-box svg {{ flex-shrink: 0; margin-top: 2px; color: var(--accent); }}
/* ── rank legend ── */
.rank-legend {{
  display: flex;
  flex-wrap: wrap;
  gap: 6px;
  margin-top: 10px;
}}
.hidden {{ display: none !important; }}
/* ── drop zone ── */
.drop-zone {{
  border: 2px dashed var(--border);
  border-radius: var(--radius);
  padding: 24px 20px 16px;
  text-align: center;
  cursor: pointer;
  transition: border-color 0.2s, background 0.2s;
  position: relative;
  user-select: none;
}}
.drop-zone:hover {{
  border-color: var(--accent);
  background: rgba(88,166,255,0.04);
}}
.drop-zone.drag-over {{
  border-color: var(--accent);
  background: rgba(88,166,255,0.08);
  box-shadow: 0 0 0 3px rgba(88,166,255,0.15);
}}
.drop-zone.drag-over .drop-icon svg {{
  stroke: var(--accent);
  transform: translateY(-4px) scale(1.1);
}}
.drop-icon {{
  margin-bottom: 8px;
}}
.drop-icon svg {{
  width: 36px;
  height: 36px;
  stroke: var(--text3);
  transition: stroke 0.2s, transform 0.2s;
}}
.drop-text {{
  font-size: 14px;
  font-weight: 600;
  color: var(--text2);
  margin-bottom: 4px;
}}
.drop-subtext {{
  font-size: 12px;
  color: var(--text3);
  margin-bottom: 10px;
  font-family: 'JetBrains Mono', monospace;
}}
.drop-status {{
  font-size: 13px;
  color: var(--text3);
  padding-top: 8px;
  border-top: 1px solid var(--border2);
}}
/* scrollbar */
::-webkit-scrollbar {{ width: 8px; height: 8px; }}
::-webkit-scrollbar-track {{ background: var(--bg2); }}
::-webkit-scrollbar-thumb {{ background: var(--border); border-radius: 4px; }}
::-webkit-scrollbar-thumb:hover {{ background: var(--text3); }}
</style>
</head>
<body>

<header class="site-header">
  <div class="logo">kraken2<span>://</span>browser</div>
  <span class="badge">v2.0</span>
</header>

<div class="main-wrap">

<div class="card">
  <div class="card-title">Classification Summary</div>
  <div class="stats-grid">
    <div class="stat-card total">
      <div class="stat-label">Total Reads</div>
      <div class="stat-val total">{total_all:,}</div>
      <div class="stat-sub">100.00%</div>
    </div>
    <div class="stat-card classified">
      <div class="stat-label">Classified</div>
      <div class="stat-val classified">{total_classified:,}</div>
      <div class="stat-sub">{pct_classified:.2f}% of total</div>
    </div>
    <div class="stat-card unclassified">
      <div class="stat-label">Unclassified</div>
      <div class="stat-val unclassified">{total_unclassified:,}</div>
      <div class="stat-sub">{pct_unclassified:.2f}% of total</div>
    </div>
  </div>
</div>

<div class="card">
  <div class="card-title">How to Use</div>
  <div class="info-box">
    <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" viewBox="0 0 16 16">
      <path d="M8 15A7 7 0 1 1 8 1a7 7 0 0 1 0 14zm0 1A8 8 0 1 0 8 0a8 8 0 0 0 0 16z"/>
      <path d="m8.93 6.588-2.29.287-.082.38.45.083c.294.07.352.176.288.469l-.738 3.468c-.194.897.105 1.319.808 1.319.545 0 1.178-.252 1.465-.598l.088-.416c-.2.176-.492.246-.686.246-.275 0-.375-.193-.304-.533L8.93 6.588zM9 4.5a1 1 0 1 1-2 0 1 1 0 0 1 2 0z"/>
    </svg>
    <div>
      載入 <strong>FASTA / FASTQ</strong> 檔案後，點擊序列 ID 即可直接送出 NCBI BLAST 查詢。<br>
      每個分類節點顯示：<strong>subtree reads</strong>（含子節點）、<strong>direct reads</strong>（僅此節點）、<strong>% of all</strong>（佔全部序列）、<strong>% of classified</strong>（佔已分類序列）。
    </div>
  </div>
  <div class="rank-legend" style="margin-top:14px;">
    <span class="rank-badge rank-superkingdom">superkingdom</span>
    <span class="rank-badge rank-kingdom">kingdom</span>
    <span class="rank-badge rank-phylum">phylum</span>
    <span class="rank-badge rank-class">class</span>
    <span class="rank-badge rank-order">order</span>
    <span class="rank-badge rank-family">family</span>
    <span class="rank-badge rank-genus">genus</span>
    <span class="rank-badge rank-species">species</span>
    <span class="rank-badge rank-other">other</span>
  </div>
</div>

<div class="card">
  <!-- Drop Zone -->
  <div id="dropZone" class="drop-zone" onclick="document.getElementById('seqFileInput').click()">
    <input type="file" id="seqFileInput" accept=".fasta,.fa,.fq,.fastq,.gz">
    <div class="drop-icon">
      <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="1.5">
        <path stroke-linecap="round" stroke-linejoin="round" d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5m-13.5-9L12 3m0 0l4.5 4.5M12 3v13.5"/>
      </svg>
    </div>
    <div class="drop-text">拖曳 FASTA / FASTQ 到此處</div>
    <div class="drop-subtext">或點擊選擇檔案 &nbsp;·&nbsp; 支援 .fasta .fa .fastq .fq .gz</div>
    <div id="fileStatus" class="drop-status">尚未載入序列檔。</div>
  </div>

  <!-- Controls row -->
  <div class="controls-row" style="margin-top: 14px;">
    <div class="search-wrap">
      <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2">
        <circle cx="11" cy="11" r="8"/><path stroke-linecap="round" d="m21 21-4.35-4.35"/>
      </svg>
      <input type="text" id="searchBox" placeholder="搜尋 taxon name 或 seq ID…">
    </div>
    <button class="btn" onclick="expandAll()">
      <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2">
        <path stroke-linecap="round" stroke-linejoin="round" d="M4 8V4m0 0h4M4 4l5 5m11-1V4m0 0h-4m4 0l-5 5M4 16v4m0 0h4m-4 0l5-5m11 5l-5-5m5 5v-4m0 4h-4"/>
      </svg>
      全部展開
    </button>
    <button class="btn" onclick="collapseAll()">
      <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2">
        <path stroke-linecap="round" stroke-linejoin="round" d="M9 9V4.5M9 9H4.5M9 9 3.75 3.75M9 15v4.5M9 15H4.5M9 15l-5.25 5.25M15 9h4.5M15 9V4.5M15 9l5.25-5.25M15 15h4.5M15 15v4.5m0-4.5 5.25 5.25"/>
      </svg>
      全部收合
    </button>
    <button class="btn" onclick="clearSearch()">
      <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2">
        <path stroke-linecap="round" stroke-linejoin="round" d="M6 18L18 6M6 6l12 12"/>
      </svg>
      清除搜尋
    </button>
  </div>
</div>

<div class="card" style="padding: 16px 20px;">
  <div id="taxonomyTree"></div>
</div>

</div>

<script>
const TREE_DATA = {tree_json_str};
let SEQ_MAP = {{}};

function escapeHtml(text) {{
  return String(text)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#039;");
}}

function formatPct(x) {{
  return Number(x).toFixed(2) + "%";
}}

function formatNum(n) {{
  return Number(n).toLocaleString();
}}

function openBlastBySeqId(seqId) {{
  if (!(seqId in SEQ_MAP)) {{
    alert("找不到序列 ID: " + seqId + "\\n請先確認是否已載入正確的 FASTA/FASTQ 檔案。");
    return;
  }}
  const seq = SEQ_MAP[seqId];
  const form = document.createElement("form");
  form.method = "POST";
  form.action = "https://blast.ncbi.nlm.nih.gov/Blast.cgi";
  form.target = "_blank";
  const params = {{
    PAGE_TYPE: "BlastSearch",
    PROGRAM: "blastn",
    DATABASE: "core_nt",
    QUERY: seq
  }};
  for (const key in params) {{
    const input = document.createElement("input");
    input.type = "hidden";
    input.name = key;
    input.value = params[key];
    form.appendChild(input);
  }}
  document.body.appendChild(form);
  form.submit();
  document.body.removeChild(form);
}}

function getRankClass(rank) {{
  const r = (rank || "").toLowerCase();
  const valid = ["superkingdom","kingdom","phylum","class","order","family","genus","species","unclassified"];
  return valid.includes(r) ? "rank-" + r : "rank-other";
}}

// ── BUG FIX: renderSeqList now uses a SINGLE flat list with data-in-more-seq attribute
// so querySelectorAll(".seq-list li") correctly finds ALL li items including hidden ones.
function renderSeqList(seqs) {{
  if (!seqs || seqs.length === 0) return "";

  const MAX_SHOW = 20;
  const total = seqs.length;

  let html = `<div class="seq-section">`;
  html += `<div class="seq-section-header">
    Reads assigned here
    <span class="seq-count-badge">${{total}}</span>
  </div>`;

  // All visible items (first MAX_SHOW)
  const first = seqs.slice(0, MAX_SHOW);
  html += '<ul class="seq-list">';
  first.forEach((rec, idx) => {{
    const seqId = escapeHtml(rec.seq_id);
    const length = escapeHtml(String(rec.length));
    const taxid = escapeHtml(String(rec.taxid));
    html += `<li data-searchtext="${{seqId.toLowerCase()}}">
      <span class="seq-idx">${{idx + 1}}</span>
      <a class="seq-id-link" href="javascript:void(0)" onclick="openBlastBySeqId('${{seqId}}')">${{seqId}}</a>
      <span class="seq-len">${{length}} bp</span>
      <span class="seq-taxid">taxid:${{taxid}}</span>
    </li>`;
  }});
  html += '</ul>';

  // More items in a <details> — IMPORTANT: the inner <ul> also has class "seq-list"
  // so that querySelectorAll(".seq-list li") covers ALL li items for search
  if (total > MAX_SHOW) {{
    const rest = seqs.slice(MAX_SHOW);
    html += `<details class="more-seq">
      <summary>顯示剩餘 ${{rest.length}} 條 reads ▾</summary>
      <ul class="seq-list">`;
    rest.forEach((rec, idx) => {{
      const seqId = escapeHtml(rec.seq_id);
      const length = escapeHtml(String(rec.length));
      const taxid = escapeHtml(String(rec.taxid));
      html += `<li data-searchtext="${{seqId.toLowerCase()}}">
        <span class="seq-idx">${{MAX_SHOW + idx + 1}}</span>
        <a class="seq-id-link" href="javascript:void(0)" onclick="openBlastBySeqId('${{seqId}}')">${{seqId}}</a>
        <span class="seq-len">${{length}} bp</span>
        <span class="seq-taxid">taxid:${{taxid}}</span>
      </li>`;
    }});
    html += `</ul></details>`;
  }}

  html += `</div>`;
  return html;
}}

function renderNode(node) {{
  if (node.rank === "root") {{
    return node.children.map(ch => renderNode(ch)).join("\\n");
  }}

  const name = escapeHtml(node.name);
  const rank = escapeHtml(node.rank);
  const taxid = escapeHtml(node.taxid);
  const rankClass = getRankClass(node.rank);
  const pctAll = formatPct(node.pct_of_all);
  const pctClassified = formatPct(node.pct_of_classified);
  const searchText = (node.name + " " + node.rank + " " + node.taxid).toLowerCase();

  let childHtml = "";
  if (node.children && node.children.length > 0) {{
    childHtml += node.children.map(ch => renderNode(ch)).join("\\n");
  }}
  childHtml += renderSeqList(node.seqs);

  const directLabel = node.direct_count > 0
    ? `<span class="stat-chip chip-direct">direct: ${{formatNum(node.direct_count)}}</span>`
    : "";

  const googleUrl = "https://www.google.com/search?q=" + encodeURIComponent(node.name);

  return `
    <details class="tax-node" data-searchtext="${{escapeHtml(searchText)}}">
      <summary>
        <a class="tax-name" href="${{googleUrl}}" target="_blank" rel="noopener" onclick="event.stopPropagation()">${{name}}</a>
        <span class="rank-badge ${{rankClass}}">${{rank}}</span>
        <span class="tax-stats">
          <span class="stat-chip chip-count">reads: ${{formatNum(node.subtree_count)}}</span>
          ${{directLabel}}
          <span class="stat-chip chip-pct-all">all: ${{pctAll}}</span>
          <span class="stat-chip chip-pct-cls">cls: ${{pctClassified}}</span>
        </span>
        <span class="tax-taxid">txid${{taxid}}</span>
      </summary>
      ${{childHtml}}
    </details>
  `;
}}

function renderTree() {{
  document.getElementById("taxonomyTree").innerHTML = renderNode(TREE_DATA);
}}

function expandAll() {{
  document.querySelectorAll("details.tax-node, details.more-seq").forEach(d => d.open = true);
}}

function collapseAll() {{
  document.querySelectorAll("details.tax-node, details.more-seq").forEach(d => d.open = false);
}}

function clearSearch() {{
  document.getElementById("searchBox").value = "";
  // 移除所有 hidden，並把搜尋時強制打開的 details 全部關回去
  document.querySelectorAll("details.tax-node").forEach(el => {{
    el.classList.remove("hidden");
    el.open = false;
  }});
  document.querySelectorAll("details.more-seq").forEach(el => {{
    el.classList.remove("hidden");
    el.open = false;
  }});
  document.querySelectorAll(".seq-list li, .seq-section").forEach(el => {{
    el.classList.remove("hidden");
  }});
}}

// ── BUG FIX: Unified search that correctly handles all li items (including those
// inside more-seq details). The key fix is:
// 1. We query ALL ".seq-list li" — this now includes li inside more-seq because
//    more-seq inner ul also has class "seq-list"
// 2. When a li inside a more-seq is hit, we walk up and open the more-seq details
// 3. We never rely on .more-seq-wrap (removed)
function doSearch(keyword) {{
  const kw = keyword.trim().toLowerCase();

  const allTaxNodes = document.querySelectorAll("details.tax-node");
  const allSeqItems = document.querySelectorAll(".seq-list li");  // covers ALL li now
  const allMoreSeqDetails = document.querySelectorAll("details.more-seq");
  const allSeqSections = document.querySelectorAll(".seq-section");

  if (!kw) {{
    allTaxNodes.forEach(el => el.classList.remove("hidden"));
    allSeqItems.forEach(el => el.classList.remove("hidden"));
    allMoreSeqDetails.forEach(el => el.classList.remove("hidden"));
    allSeqSections.forEach(el => el.classList.remove("hidden"));
    return;
  }}

  // Hide everything first
  allTaxNodes.forEach(el => el.classList.add("hidden"));
  allSeqItems.forEach(el => el.classList.add("hidden"));
  allMoreSeqDetails.forEach(el => el.classList.add("hidden"));
  allSeqSections.forEach(el => el.classList.add("hidden"));

  function showAncestors(el) {{
    let p = el.parentElement;
    while (p) {{
      if (p.classList) p.classList.remove("hidden");
      if (p.tagName && p.tagName.toLowerCase() === "details") p.open = true;
      p = p.parentElement;
    }}
  }}

  // 1. Tax node hits: show the node + all its seq items + all more-seq inside
  allTaxNodes.forEach(node => {{
    const txt = (node.dataset.searchtext || "").toLowerCase();
    if (txt.includes(kw)) {{
      node.classList.remove("hidden");
      node.open = true;
      // Show all seq items and sections inside this node
      node.querySelectorAll(".seq-list li").forEach(li => li.classList.remove("hidden"));
      node.querySelectorAll(".seq-section").forEach(s => s.classList.remove("hidden"));
      node.querySelectorAll("details.more-seq").forEach(d => {{
        d.classList.remove("hidden");
        d.open = true;
      }});
      showAncestors(node);
    }}
  }});

  // 2. Seq ID hits: show the matching li and walk up to reveal parents
  allSeqItems.forEach(li => {{
    const txt = (li.dataset.searchtext || "").toLowerCase();
    if (txt.includes(kw)) {{
      li.classList.remove("hidden");
      // The seq-section containing this li must also be visible
      const section = li.closest(".seq-section");
      if (section) section.classList.remove("hidden");
      // If this li lives inside a more-seq details, open it
      const moreSeq = li.closest("details.more-seq");
      if (moreSeq) {{
        moreSeq.classList.remove("hidden");
        moreSeq.open = true;
      }}
      showAncestors(li);
    }}
  }});
}}

function parseFasta(text) {{
  const map = {{}};
  let currentId = null;
  let chunks = [];
  const lines = text.split(/\\r?\\n/);

  for (const line of lines) {{
    if (!line) continue;
    if (line.startsWith(">")) {{
      if (currentId !== null) {{
        map[currentId] = chunks.join("");
      }}
      currentId = line.slice(1).trim().split(/\\s+/)[0];
      chunks = [];
    }} else {{
      chunks.push(line.trim());
    }}
  }}
  if (currentId !== null) {{
    map[currentId] = chunks.join("");
  }}
  return map;
}}

function parseFastq(text) {{
  const map = {{}};
  const lines = text.split(/\\r?\\n/);
  for (let i = 0; i + 3 < lines.length; i += 4) {{
    const h = lines[i];
    const s = lines[i + 1];
    if (h && h.startsWith("@")) {{
      const seqId = h.slice(1).trim().split(/\\s+/)[0];
      map[seqId] = (s || "").trim();
    }}
  }}
  return map;
}}

function detectSeqFormat(text) {{
  const lines = text.split(/\\r?\\n/);
  for (const line of lines) {{
    if (!line.trim()) continue;
    if (line.startsWith(">")) return "fasta";
    if (line.startsWith("@")) return "fastq";
  }}
  return "unknown";
}}

function processFile(file) {{
  if (!file) return;

  const reader = new FileReader();
  reader.onload = function(e) {{
    const text = e.target.result;
    const fmt = detectSeqFormat(text);

    if (fmt === "fasta") {{
      SEQ_MAP = parseFasta(text);
    }} else if (fmt === "fastq") {{
      SEQ_MAP = parseFastq(text);
    }} else {{
      SEQ_MAP = {{}};
      document.getElementById("fileStatus").innerHTML =
        '<span class="warn">無法辨識檔案格式，請載入 FASTA 或 FASTQ。</span>';
      return;
    }}

    const n = Object.keys(SEQ_MAP).length;
    document.getElementById("fileStatus").innerHTML =
      '<span class="ok">✓ ' + escapeHtml(file.name) + '</span> &nbsp;·&nbsp; 共 ' + n.toLocaleString() + ' 條序列';
  }};
  reader.readAsText(file);
}}

// Click-to-browse
document.getElementById("seqFileInput").addEventListener("change", function(event) {{
  processFile(event.target.files[0]);
}});

// Drag & Drop
const dropZone = document.getElementById("dropZone");

dropZone.addEventListener("dragenter", function(e) {{
  e.preventDefault();
  e.stopPropagation();
  dropZone.classList.add("drag-over");
}});

dropZone.addEventListener("dragover", function(e) {{
  e.preventDefault();
  e.stopPropagation();
  dropZone.classList.add("drag-over");
}});

dropZone.addEventListener("dragleave", function(e) {{
  e.preventDefault();
  e.stopPropagation();
  // only remove when truly leaving the zone (not entering a child element)
  if (!dropZone.contains(e.relatedTarget)) {{
    dropZone.classList.remove("drag-over");
  }}
}});

dropZone.addEventListener("drop", function(e) {{
  e.preventDefault();
  e.stopPropagation();
  dropZone.classList.remove("drag-over");
  const file = e.dataTransfer.files[0];
  if (file) {{
    processFile(file);
  }}
}});

document.getElementById("searchBox").addEventListener("input", function() {{
  doSearch(this.value);
}});

renderTree();
</script>
</body>
</html>
"""

    with open(html_out, "w", encoding="utf-8") as f:
        f.write(html_text)


def write_excel(summary_rows, excel_out):
    df_summary = pd.DataFrame(summary_rows)

    # 依照樹形順序，不再另外排序
    # 只保留 taxonomy_summary
    with pd.ExcelWriter(excel_out, engine="openpyxl") as writer:
        df_summary.to_excel(writer, sheet_name="taxonomy_summary", index=False)

        ws1 = writer.sheets["taxonomy_summary"]
        ws1.freeze_panes = "A2"

        width_map_1 = {
            "A": 12,  # tree_order
            "B": 8,   # depth
            "C": 14,  # taxid
            "D": 32,  # name
            "E": 16,  # rank
            "F": 60,  # path_names
            "G": 50,  # path_taxids
            "H": 14,  # direct_count
            "I": 14,  # subtree_count
            "J": 14,  # pct_of_all
            "K": 18,  # pct_of_classified
        }
        for col, width in width_map_1.items():
            ws1.column_dimensions[col].width = width


def main():
    parser = argparse.ArgumentParser(
        description="Generate single HTML taxonomy browser + Excel report from Kraken2 output."
    )
    parser.add_argument("-k", "--kraken", required=True, help="Kraken2 output file")
    parser.add_argument("--nodes", required=True, help="NCBI taxonomy nodes.dmp")
    parser.add_argument("--names", required=True, help="NCBI taxonomy names.dmp")
    parser.add_argument("-o", "--out_prefix", required=True, help="Output prefix")
    parser.add_argument("-d", "--out_dir", default=".", help="Output directory")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    html_out = out_dir / f"{args.out_prefix}.html"
    excel_out = out_dir / f"{args.out_prefix}.xlsx"

    print("[1/6] Reading Kraken2 output...", file=sys.stderr)
    records = parse_kraken_output(args.kraken)

    total_all = len(records)
    total_classified = sum(1 for x in records if x["status"] == "C" and x["taxid"] not in ("0", "", "unclassified"))
    total_unclassified = total_all - total_classified

    print("[2/6] Loading taxonomy...", file=sys.stderr)
    parent, rank_map, name_map = load_taxonomy(args.nodes, args.names)

    print("[3/6] Building taxonomy tree...", file=sys.stderr)
    tree = build_tree(records, parent, rank_map, name_map)

    print("[4/6] Computing counts and percentages...", file=sys.stderr)
    compute_metrics(tree, total_all=total_all, total_classified=total_classified)
    sort_tree(tree)

    print("[5/6] Preparing Excel tables...", file=sys.stderr)
    summary_rows = flatten_tree(tree)

    print("[6/6] Writing outputs...", file=sys.stderr)
    generate_html(tree, total_all, total_classified, total_unclassified, html_out)
    write_excel(summary_rows, excel_out)

    print(f"HTML written to: {html_out}", file=sys.stderr)
    print(f"Excel written to: {excel_out}", file=sys.stderr)


if __name__ == "__main__":
    main()