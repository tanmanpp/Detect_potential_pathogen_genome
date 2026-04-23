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
<style>
body {{
    font-family: Arial, "Microsoft JhengHei", sans-serif;
    margin: 20px;
    background: #f7f7f7;
    color: #222;
}}
h1 {{
    margin-bottom: 8px;
}}
.panel {{
    background: white;
    border: 1px solid #ddd;
    border-radius: 12px;
    padding: 14px 18px;
    margin-bottom: 16px;
}}
.note {{
    color: #555;
    line-height: 1.6;
}}
.summary-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
    gap: 10px;
}}
.summary-card {{
    background: #fafafa;
    border: 1px solid #ececec;
    border-radius: 10px;
    padding: 10px 12px;
}}
.controls {{
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    align-items: center;
}}
input[type="file"], input[type="text"] {{
    padding: 6px;
}}
button {{
    padding: 7px 12px;
    border-radius: 8px;
    border: 1px solid #bbb;
    background: white;
    cursor: pointer;
}}
button:hover {{
    background: #f0f0f0;
}}
details.tax-node {{
    margin-left: 16px;
    margin-top: 4px;
    margin-bottom: 4px;
}}
summary {{
    cursor: pointer;
    padding: 4px 0;
    line-height: 1.5;
}}
.tax-name {{
    font-weight: 600;
}}
.tax-rank {{
    color: #444;
    margin-left: 6px;
}}
.tax-count {{
    color: #b00020;
    margin-left: 8px;
    font-weight: 600;
}}
.tax-pct {{
    color: #0b57d0;
    margin-left: 8px;
    font-weight: 600;
}}
.tax-taxid {{
    color: #777;
    margin-left: 8px;
    font-size: 0.9em;
}}
.seq-list {{
    margin-top: 6px;
    margin-bottom: 8px;
}}
.seq-list li {{
    margin: 3px 0;
}}
.seq-len {{
    color: #666;
}}
.seq-meta {{
    color: #777;
    margin-left: 8px;
    font-size: 0.9em;
}}
a {{
    text-decoration: none;
}}
a:hover {{
    text-decoration: underline;
}}
.warn {{
    color: #b00020;
    font-weight: 600;
}}
.ok {{
    color: #136c2e;
    font-weight: 600;
}}
.hidden {{
    display: none;
}}
</style>
</head>
<body>
<h1>Kraken2 Taxonomy Browser</h1>

<div class="panel">
  <div class="summary-grid">
    <div class="summary-card"><b>總序列數</b><br>{total_all}</div>
    <div class="summary-card"><b>已分類序列數</b><br>{total_classified} ({pct_classified:.2f}%)</div>
    <div class="summary-card"><b>未分類序列數</b><br>{total_unclassified} ({pct_unclassified:.2f}%)</div>
  </div>
</div>

<div class="panel">
  <div class="note">
    若要點擊序列 ID 後直接跳到 NCBI BLAST，請先在下方載入你的 FASTA 或 FASTQ 檔案。<br>
    每個分類節點旁會顯示：
    <br>1. 此節點底下的序列總數
    <br>2. 占全部序列百分比 (% of all)
    <br>3. 占已分類序列百分比 (% of classified)
  </div>
</div>

<div class="panel">
  <div class="controls">
    <label><b>載入 FASTA/FASTQ：</b></label>
    <input type="file" id="seqFileInput">
    <input type="text" id="searchBox" placeholder="搜尋 taxon name 或 seq ID">
    <button onclick="expandAll()">全部展開</button>
    <button onclick="collapseAll()">全部收合</button>
    <button onclick="clearSearch()">清除搜尋</button>
  </div>
  <div id="fileStatus" class="note" style="margin-top:10px;">尚未載入序列檔。</div>
</div>

<div class="panel">
  <div id="taxonomyTree"></div>
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

function openBlastBySeqId(seqId) {{
  if (!(seqId in SEQ_MAP)) {{
    alert("找不到序列 ID: " + seqId + "\\n請先確認是否已載入正確的 FASTA/FASTQ 檔案。");
    return;
  }}

  const seq = SEQ_MAP[seqId];

  // 建立 POST form
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

function renderSeqList(seqs) {{
  if (!seqs || seqs.length === 0) return "";

  const MAX_SHOW = 20;

  let html = '<ul class="seq-list">';

  // ===== 前 20 條 =====
  const first = seqs.slice(0, MAX_SHOW);
  for (const rec of first) {{
    const seqId = escapeHtml(rec.seq_id);
    const length = escapeHtml(String(rec.length));
    const taxid = escapeHtml(String(rec.taxid));

    html += `
      <li data-searchtext="${{seqId.toLowerCase()}}">
        <a href="javascript:void(0)" onclick="openBlastBySeqId('${{seqId}}')">${{seqId}}</a>
        <span class="seq-len">(${{length}} bp)</span>
        <span class="seq-meta">taxid:${{taxid}}</span>
      </li>
    `;
  }}

  // ===== 剩下的 =====
  if (seqs.length > MAX_SHOW) {{
    html += `
      <details class="more-seq">
        <summary>顯示更多 (${{seqs.length - MAX_SHOW}} 條)</summary>
        <ul>
    `;

    const rest = seqs.slice(MAX_SHOW);
    for (const rec of rest) {{
      const seqId = escapeHtml(rec.seq_id);
      const length = escapeHtml(String(rec.length));
      const taxid = escapeHtml(String(rec.taxid));

      html += `
        <li data-searchtext="${{seqId.toLowerCase()}}">
          <a href="javascript:void(0)" onclick="openBlastBySeqId('${{seqId}}')">${{seqId}}</a>
          <span class="seq-len">(${{length}} bp)</span>
          <span class="seq-meta">taxid:${{taxid}}</span>
        </li>
      `;
    }}

    html += `
        </ul>
      </details>
    `;
  }}

  html += '</ul>';
  return html;
}}

function renderNode(node) {{
  if (node.rank === "root") {{
    return node.children.map(ch => renderNode(ch)).join("\\n");
  }}

  const name = escapeHtml(node.name);
  const rank = escapeHtml(node.rank);
  const taxid = escapeHtml(node.taxid);
  const count = escapeHtml(String(node.subtree_count));
  const pctAll = formatPct(node.pct_of_all);
  const pctClassified = formatPct(node.pct_of_classified);
  const searchText = (node.name + " " + node.rank + " " + node.taxid).toLowerCase();

  let childHtml = "";
  if (node.children && node.children.length > 0) {{
    childHtml += node.children.map(ch => renderNode(ch)).join("\\n");
  }}
  childHtml += renderSeqList(node.seqs);

  return `
    <details class="tax-node" data-searchtext="${{escapeHtml(searchText)}}">
      <summary>
        <span class="tax-name">${{name}}</span>
        <span class="tax-rank">[${{rank}}]</span>
        <span class="tax-count">count: ${{count}};</span>
        <span class="tax-pct">% of all: ${{pctAll}},</span>
        <span class="tax-pct">% of classified: ${{pctClassified}}</span>
        <span class="tax-taxid">taxid:${{taxid}}</span>
      </summary>
      ${{childHtml}}
    </details>
  `;
}}

function renderTree() {{
  document.getElementById("taxonomyTree").innerHTML = renderNode(TREE_DATA);
}}

function expandAll() {{
  document.querySelectorAll("details.tax-node").forEach(d => d.open = true);
}}

function collapseAll() {{
  document.querySelectorAll("details.tax-node").forEach(d => d.open = false);
}}

function clearSearch() {{
  document.getElementById("searchBox").value = "";

  document.querySelectorAll("details.tax-node, .seq-list li, .more-seq-wrap, details.more-seq").forEach(el => {{
    el.classList.remove("hidden");
  }});
}}

function doSearch(keyword) {{
  const kw = keyword.trim().toLowerCase();

  const allTaxNodes = document.querySelectorAll("details.tax-node");
  const allSeqItems = document.querySelectorAll(".seq-list li");
  const allMoreSeqWraps = document.querySelectorAll(".more-seq-wrap");
  const allMoreSeqDetails = document.querySelectorAll("details.more-seq");

  // 沒輸入時全部恢復
  if (!kw) {{
    allTaxNodes.forEach(node => {{
      node.classList.remove("hidden");
    }});
    allSeqItems.forEach(li => {{
      li.classList.remove("hidden");
    }});
    allMoreSeqWraps.forEach(wrap => {{
      wrap.classList.remove("hidden");
    }});
    allMoreSeqDetails.forEach(d => {{
      d.classList.remove("hidden");
    }});
    return;
  }}

  // 先全部隱藏
  allTaxNodes.forEach(node => {{
    node.classList.add("hidden");
  }});
  allSeqItems.forEach(li => {{
    li.classList.add("hidden");
  }});
  allMoreSeqWraps.forEach(wrap => {{
    wrap.classList.add("hidden");
  }});
  allMoreSeqDetails.forEach(d => {{
    d.classList.add("hidden");
  }});

  // 先處理 taxon node 命中
  allTaxNodes.forEach(node => {{
    const txt = (node.dataset.searchtext || "").toLowerCase();
    if (txt.includes(kw)) {{
      node.classList.remove("hidden");
      node.open = true;

      
      // 命中的 taxon node，順便顯示底下所有 seq list 項目與 more-seq
      node.querySelectorAll(".seq-list li").forEach(li => li.classList.remove("hidden"));
      node.querySelectorAll("details.more-seq").forEach(d => {{
        d.classList.remove("hidden");
        d.open = true;
      }});


      let p = node.parentElement;
      while (p) {{
        if (p.classList) {{
          p.classList.remove("hidden");
        }}
        if (p.tagName && p.tagName.toLowerCase() === "details") {{
          p.open = true;
        }}
        p = p.parentElement;
      }}
    }}
  }});

  // 再處理 seq id 命中
  allSeqItems.forEach(li => {{
    const txt = (li.dataset.searchtext || "").toLowerCase();
    if (txt.includes(kw)) {{
      li.classList.remove("hidden");

      let p = li.parentElement;
      while (p) {{
        if (p.classList) {{
          p.classList.remove("hidden");
        }}

        if (p.tagName && p.tagName.toLowerCase() === "details") {{
          p.open = true;
        }}

        p = p.parentElement;
      }}
    }}
  }});

  // 若某個 more-seq 裡面有命中的 li，顯示該 more-seq-wrap 與 more-seq
  allMoreSeqDetails.forEach(detail => {{
    const visibleHit = detail.querySelector("li:not(.hidden)");
    if (visibleHit) {{
      detail.classList.remove("hidden");
      detail.open = true;

      const wrap = detail.closest(".more-seq-wrap");
      if (wrap) {{
        wrap.classList.remove("hidden");
      }}

      let p = detail.parentElement;
      while (p) {{
        if (p.classList) {{
          p.classList.remove("hidden");
        }}
        if (p.tagName && p.tagName.toLowerCase() === "details") {{
          p.open = true;
        }}
        p = p.parentElement;
      }}
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

document.getElementById("seqFileInput").addEventListener("change", function(event) {{
  const file = event.target.files[0];
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
      '<span class="ok">已載入 ' + escapeHtml(file.name) + '</span>，共 ' + n + ' 條序列。';
  }};
  reader.readAsText(file);
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