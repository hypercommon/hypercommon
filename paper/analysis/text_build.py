"""Build sentence networks from long English Wikipedia featured articles, with the
articles' sections and subsections (nested) as ground-truth groups.

The recipe is the Hi-C top-k recipe (hic_build.py) moved to text: nodes are units
along a 1D order (sentences in reading order), each node keeps its k most similar
partners among nearby units, and the groups are contiguous, nested segments
(sections), like contact domains along a chromosome.

Source data (fetched as plain JSON from the public MediaWiki API, never executed):
  * the member list of Category:Featured articles (namespace 0), and
  * for each article tried, its plain-text extract with section headings
    (action=query, prop=extracts, explaintext=1, exsectionformat=wiki) and the
    revision id it was rendered from (prop=revisions).
Both are cached as JSON under RAW_DIR (default D:/tmp/text_raw), so a rebuild
does not fetch again. Text is CC BY-SA 4.0 (see data_text/README.md).

Usage (run from paper/ with an isolated interpreter; the cache is untrusted data):
  python -I analysis/text_build.py build --raw D:/tmp/text_raw --out analysis/data_text
  python -I analysis/text_build.py table --out analysis/data_text      # README tables

`build` writes, for each selected article and each k, a graph folder
<name>_k<k>/ and its position-only control <name>_k<k>_null/, plus
selection_log.tsv and summary.json.  Everything below is fixed before any
detection result was seen.

Selection rule:
  1. All namespace-0 members of Category:Featured articles (fetched once, cached),
     sorted by title, shuffled with random.Random(SEED).
  2. Walk the shuffled list; fetch each article; build its units and groups
     (below); accept it if N_MIN <= n <= N_MAX and it has >= MIN_GROUPS groups.
     Stop at N_ARTICLES accepted articles.

Units (nodes):
  * The extract is split into sections at its '== heading ==' lines (any level);
    text before the first heading is the lead.  Sections whose heading is one of
    END_MATTER (references, notes, see also, ...) are dropped with all their
    subsections.
  * Each non-empty line (paragraph or list item) is split into sentences with
    the fixed regex splitter split_sentences() (no NLP model).
  * A sentence is kept only if it has >= MIN_TOKENS content tokens (scikit-learn
    default token pattern, lower-cased, English stop words removed).  Kept
    sentences are numbered 0..n-1 in reading order.

Similarity: TF-IDF fitted on the article's own kept sentences (scikit-learn
  TfidfVectorizer: lowercase, stop_words='english', default token pattern,
  smooth idf, l2 norm); similarity = cosine = dot product of the rows.

Edges: for each node i, the candidates are nodes j with 0 < |i-j| <= W and
  sim(i, j) > 0; keep the k candidates with the largest sim (ties: smaller |i-j|,
  then smaller j).  Undirected, unweighted edge set = union over i.

Groups: the lead and every kept section at every heading level, each as the set
  of kept sentences inside it including its subsections.  Groups with fewer than
  3 nodes, exact duplicates of an earlier group (a section whose only content is
  one subsection) and any group equal to the whole article are dropped.

Control (<name>_null): same nodes and groups; each node keeps its k nearest
  nodes by position (|i-j|, ties to the smaller id), union over i - a pure band.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import re
import sys
import time

import numpy as np

API = "https://en.wikipedia.org/w/api.php"
UA = "seq-graph-builder/0.1 (research script; python-requests)"
CATEGORY = "Category:Featured articles"
SEED = 2026
N_ARTICLES = 15
N_MIN, N_MAX = 300, 3000
MIN_GROUPS = 5
MIN_TOKENS = 3
W = 40
KS = (8, 16)
MIN_GROUP = 3
MAX_TRIES = 400
SLEEP = 2.0

END_MATTER = {
    "see also", "notes", "references", "sources", "citations", "bibliography",
    "further reading", "external links", "footnotes", "works cited",
    "explanatory notes", "general references", "cited sources", "cited works",
    "notes and references", "references and notes", "notes and citations",
    "primary sources", "secondary sources", "select bibliography",
    "selected bibliography", "general and cited sources", "general bibliography",
    "citations and notes", "content notes", "informational notes",
    "reference notes", "endnotes", "web sources", "book sources",
    "journal sources", "news sources", "online sources", "other sources",
}

ABBREV = {
    "mr", "mrs", "ms", "dr", "st", "jr", "sr", "mt", "ft", "lt", "col", "gen",
    "capt", "sgt", "rev", "prof", "no", "nos", "vol", "vols", "vs", "etc", "e.g",
    "i.e", "c", "ca", "approx", "inc", "ltd", "co", "corp", "bros", "jan", "feb",
    "mar", "apr", "jun", "jul", "aug", "sep", "sept", "oct", "nov", "dec",
    "u.s", "u.k", "u.n", "a.d", "b.c", "p", "pp", "fig", "op", "ed", "eds",
    "hon", "gov", "sen", "rep", "adm", "cmdr", "maj", "brig", "pvt", "esq",
    "dept", "univ", "assn", "est", "al", "cf", "viz", "lit",
}

HEAD_RE = re.compile(r"^(={2,6})\s*(.*?)\s*\1\s*$")
# candidate break: [.!?] + optional closing quotes/brackets + whitespace,
# followed by an uppercase letter, a digit or an opening quote/bracket
BREAK_RE = re.compile(r"([.!?][\"'\u201d\u2019)\]]*)\s+(?=[\"'\u201c\u2018(\[]?[A-Z0-9])")
TOKEN_RE = re.compile(r"(?u)\b\w\w+\b")


# ----------------------------------------------------------------------------
# fetching (polite: one request at a time, SLEEP s between requests)
# ----------------------------------------------------------------------------
_session = None
_last = [0.0]


def api(params):
    import requests
    global _session
    if _session is None:
        _session = requests.Session()
        _session.headers["User-Agent"] = UA
    p = dict(params, format="json", formatversion="2", maxlag="5")
    for attempt in range(6):
        wait = _last[0] + SLEEP - time.time()
        if wait > 0:
            time.sleep(wait)
        try:
            r = _session.get(API, params=p, timeout=60)
            _last[0] = time.time()
            if r.status_code == 200:
                d = r.json()
                if "error" in d and d["error"].get("code") == "maxlag":
                    raise IOError("maxlag")
                if "error" in d:
                    raise ValueError(d["error"])
                return d
            ra = r.headers.get("Retry-After", "")
            raise IOError(f"HTTP {r.status_code}", int(ra) if ra.isdigit() else 0)
        except (IOError, OSError) as e:
            _last[0] = time.time()
            wait = max(5 * 2 ** attempt, e.args[1] if len(e.args) > 1 else 0)
            print(f"  request failed ({e.args[0]}); retry in {wait}s", file=sys.stderr)
            time.sleep(wait)
    raise IOError("giving up")


def featured_list(raw):
    path = os.path.join(raw, "featured_articles.json")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    members, cont = [], {}
    while True:
        d = api(dict(action="query", list="categorymembers", cmtitle=CATEGORY,
                     cmnamespace="0", cmlimit="500", cmprop="ids|title", **cont))
        members += [dict(pageid=m["pageid"], title=m["title"])
                    for m in d["query"]["categorymembers"]]
        if "continue" not in d:
            break
        cont = {"cmcontinue": d["continue"]["cmcontinue"], "continue": d["continue"]["continue"]}
    out = dict(fetched_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               category=CATEGORY, members=members)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=0)
    print(f"category: {len(members)} members", file=sys.stderr)
    return out


def article(raw, pageid):
    d_dir = os.path.join(raw, "articles")
    os.makedirs(d_dir, exist_ok=True)
    path = os.path.join(d_dir, f"{pageid}.json")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    d = api(dict(action="query", prop="extracts|revisions", pageids=str(pageid),
                 explaintext="1", exsectionformat="wiki", rvprop="ids|timestamp"))
    pg = d["query"]["pages"][0]
    rev = pg["revisions"][0]
    out = dict(pageid=pg["pageid"], title=pg["title"], revid=rev["revid"],
               rev_timestamp=rev["timestamp"],
               fetched_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               extract=pg.get("extract", ""))
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, ensure_ascii=False)
    return out


# ----------------------------------------------------------------------------
# text -> sections -> sentences
# ----------------------------------------------------------------------------
def split_sentences(line):
    """Fixed regex sentence splitter; no break after a known abbreviation or a
    single-letter initial."""
    out, start = [], 0
    for m in BREAK_RE.finditer(line):
        end = m.end(1)
        before = line[start:m.start(1) + 1]
        words = before[:-1].split()
        last = words[-1].lower().strip("(\"'\u201c\u2018") if words else ""
        if line[m.start(1)] == "." and (last in ABBREV or
                                        (len(last) == 1 and last.isalpha())):
            continue
        s = line[start:end].strip()
        if s:
            out.append(s)
        start = m.end()
    s = line[start:].strip()
    if s:
        out.append(s)
    return out


def parse_sections(extract):
    """Return list of (path tuple, level, [lines]).  path () is the lead."""
    secs = [((), 1, [])]
    stack = []  # (level, title)
    dropped_at = None  # level of the end-matter heading being skipped
    for raw_line in extract.split("\n"):
        line = raw_line.strip()
        m = HEAD_RE.match(line)
        if m:
            level = len(m.group(1))
            title = m.group(2).strip()
            if dropped_at is not None and level > dropped_at:
                continue
            dropped_at = None
            if title.lower() in END_MATTER:
                dropped_at = level
                continue
            while stack and stack[-1][0] >= level:
                stack.pop()
            stack.append((level, title))
            secs.append((tuple(t for _, t in stack), level, []))
            continue
        if dropped_at is not None or not line:
            continue
        secs[-1][2].append(line)
    return secs


def units_and_groups(extract, stop):
    secs = parse_sections(extract)
    units = []  # dict(text, path, sent_index)
    n_sent = 0
    for path, _, lines in secs:
        for line in lines:
            for s in split_sentences(line):
                toks = [t for t in TOKEN_RE.findall(s.lower()) if t not in stop]
                if len(toks) >= MIN_TOKENS:
                    units.append(dict(text=s, path=path, sent_index=n_sent))
                n_sent += 1
    n = len(units)
    # groups: lead + every section path prefix, in document order
    order, members = [], {}
    for i, u in enumerate(units):
        keys = [("lead",)] if u["path"] == () else \
               [("sec",) + u["path"][:d] for d in range(1, len(u["path"]) + 1)]
        for key in keys:
            if key not in members:
                members[key] = []
                order.append(key)
            members[key].append(i)
    groups, labels, seen = [], [], set()
    n_small = n_dup = n_whole = 0
    for key in order:
        g = tuple(members[key])
        if len(g) < MIN_GROUP:
            n_small += 1
            continue
        if len(g) == n:
            n_whole += 1
            continue
        if g in seen:
            n_dup += 1
            continue
        seen.add(g)
        groups.append(g)
        labels.append("(lead)" if key == ("lead",) else " > ".join(key[1:]))
    info = dict(n_sentences_all=n_sent, n_sections=len(order),
                n_groups_dropped_small=n_small, n_groups_dropped_dup=n_dup,
                n_groups_dropped_whole=n_whole)
    return units, groups, labels, info


# ----------------------------------------------------------------------------
# graphs
# ----------------------------------------------------------------------------
def tfidf_sim(texts):
    from sklearn.feature_extraction.text import TfidfVectorizer
    X = TfidfVectorizer(lowercase=True, stop_words="english").fit_transform(texts)
    return (X @ X.T).toarray()


def topk_edges(n, k, score):
    """score(i, cand) -> array of scores; keep score > 0; ties by |i-j|, then j."""
    E = set()
    for i in range(n):
        cand = np.array([j for j in range(max(0, i - W), min(n, i + W + 1)) if j != i])
        s = score(i, cand)
        keep = s > 0
        cand, s = cand[keep], s[keep]
        if len(cand) == 0:
            continue
        order = np.lexsort((cand, np.abs(cand - i), -s))
        for j in cand[order[:k]]:
            j = int(j)
            E.add((min(i, j), max(i, j)))
    return sorted(E)


def graph_stats(n, edges, groups):
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    memb = [set() for _ in range(n)]
    for gi, g in enumerate(groups):
        for u in g:
            memb[u].add(gi)
    # deepest (smallest) group of each node
    leaf = [min(ms, key=lambda gi: (len(groups[gi]), gi)) if ms else None for ms in memb]
    d = np.array([b - a for a, b in edges])
    return dict(m=len(edges), mean_degree=2 * len(edges) / n,
                clustering=nx.average_clustering(G),
                n_components=nx.number_connected_components(G),
                n_isolated=sum(1 for v in G if G.degree(v) == 0),
                median_dist=float(np.median(d)),
                frac_dist_le5=float(np.mean(d <= 5)),
                frac_edges_within_group=float(np.mean([bool(memb[a] & memb[b])
                                                       for a, b in edges])),
                frac_edges_same_leaf=float(np.mean([leaf[a] is not None and leaf[a] == leaf[b]
                                                    for a, b in edges])))


def write_graph(od, edges, groups, nodes_rows):
    os.makedirs(od, exist_ok=True)
    with open(os.path.join(od, "edges.tsv"), "w", newline="\n") as fh:
        fh.writelines(f"{a}\t{b}\n" for a, b in edges)
    with open(os.path.join(od, "groups.txt"), "w", newline="\n") as fh:
        fh.writelines(" ".join(map(str, g)) + "\n" for g in groups)
    with open(os.path.join(od, "nodes.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("id\tposition\tsentence_index\tsection_path\n")
        fh.writelines(r for r in nodes_rows)


def slug(title):
    s = re.sub(r"[^a-z0-9]+", "_", title.lower()).strip("_")
    return s[:40].rstrip("_")


def cmd_build(args):
    from sklearn.feature_extraction.text import ENGLISH_STOP_WORDS
    stop = set(ENGLISH_STOP_WORDS)
    os.makedirs(args.out, exist_ok=True)
    fl = featured_list(args.raw)
    titles = sorted(fl["members"], key=lambda m: m["title"])
    rng = random.Random(SEED)
    rng.shuffle(titles)
    log, accepted = [], []
    for tries, m in enumerate(titles[:MAX_TRIES], 1):
        a = article(args.raw, m["pageid"])
        units, groups, labels, info = units_and_groups(a["extract"], stop)
        n = len(units)
        if not (N_MIN <= n <= N_MAX):
            why = f"n={n} outside {N_MIN}-{N_MAX}"
        elif len(groups) < MIN_GROUPS:
            why = f"groups={len(groups)} < {MIN_GROUPS}"
        else:
            why = "accepted"
        log.append(dict(order=tries, pageid=a["pageid"], title=a["title"],
                        revid=a["revid"], n=n, n_groups=len(groups), result=why))
        print(f"{tries:4d} {a['title'][:50]:50s} n={n:5d} groups={len(groups):3d} {why}",
              file=sys.stderr)
        if why == "accepted":
            accepted.append((a, units, groups, labels, info))
            if len(accepted) == N_ARTICLES:
                break
    with open(os.path.join(args.out, "selection_log.tsv"), "w", encoding="utf-8",
              newline="\n") as fh:
        fh.write("order\tpageid\ttitle\trevid\tn\tn_groups\tresult\n")
        for r in log:
            fh.write("\t".join(str(r[c]) for c in
                               ("order", "pageid", "title", "revid", "n", "n_groups",
                                "result")) + "\n")
    summary = []
    for a, units, groups, labels, info in accepted:
        n = len(units)
        sim = tfidf_sim([u["text"] for u in units])
        sizes = [len(g) for g in groups]
        covered = set().union(*map(set, groups))
        cnt = np.zeros(n, int)
        for g in groups:
            cnt[list(g)] += 1
        depth = max(len(u["path"]) for u in units)
        rows = [f"{i}\t{i}\t{u['sent_index']}\t"
                f"{'(lead)' if u['path'] == () else ' > '.join(u['path'])}\n"
                for i, u in enumerate(units)]
        for k in KS:
            name = f"wiki_{slug(a['title'])}_k{k}"
            real = topk_edges(n, k, lambda i, c: sim[i, c])
            null = topk_edges(n, k, lambda i, c: 1.0 / (1.0 + np.abs(c - i)))
            write_graph(os.path.join(args.out, name), real, groups, rows)
            write_graph(os.path.join(args.out, name + "_null"), null, groups, rows)
            sr, sn = graph_stats(n, real, groups), graph_stats(n, null, groups)
            row = dict(graph=name, title=a["title"], pageid=a["pageid"],
                       revid=a["revid"], rev_timestamp=a["rev_timestamp"],
                       url=f"https://en.wikipedia.org/w/index.php?oldid={a['revid']}",
                       k=k, W=W, n=n, **info, max_heading_depth=depth,
                       **sr,
                       n_groups=len(groups), group_size_min=min(sizes),
                       group_size_median=float(np.median(sizes)),
                       group_size_max=max(sizes),
                       frac_nodes_in_group=len(covered) / n,
                       n_nodes_in_multiple_groups=int((cnt > 1).sum()),
                       **{"null_" + key: v for key, v in sn.items()})
            summary.append(row)
            print(json.dumps({c: row[c] for c in ("graph", "n", "m", "n_groups",
                                                  "frac_edges_within_group",
                                                  "null_frac_edges_within_group")}))
    with open(os.path.join(args.out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)


def cmd_table(args):
    with open(os.path.join(args.out, "summary.json"), encoding="utf-8") as fh:
        S = json.load(fh)
    print("| graph | n | m | mean degree | clustering (control) | components (isolated) | "
          "median edge length (control) | groups | group size min-median-max | "
          "within | within, control | same leaf | same leaf, control |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in S:
        print(f"| {r['graph']} | {r['n']} | {r['m']} | {r['mean_degree']:.2f} | "
              f"{r['clustering']:.3f} ({r['null_clustering']:.3f}) | "
              f"{r['n_components']} ({r['n_isolated']}) | "
              f"{r['median_dist']:g} ({r['null_median_dist']:g}) | {r['n_groups']} | "
              f"{r['group_size_min']}-{r['group_size_median']:g}-{r['group_size_max']} | "
              f"{r['frac_edges_within_group']:.3f} | {r['null_frac_edges_within_group']:.3f} | "
              f"{r['frac_edges_same_leaf']:.3f} | {r['null_frac_edges_same_leaf']:.3f} |")
    for k in KS:
        R = [r for r in S if r["k"] == k]
        med = lambda c: float(np.median([r[c] for r in R]))
        print(f"k={k}: median within {med('frac_edges_within_group'):.3f} vs control "
              f"{med('null_frac_edges_within_group'):.3f}; same leaf "
              f"{med('frac_edges_same_leaf'):.3f} vs {med('null_frac_edges_same_leaf'):.3f}; "
              f"median clustering {med('clustering'):.3f} vs {med('null_clustering'):.3f}")


def main():
    for st in (sys.stdout, sys.stderr):
        st.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("build")
    a.add_argument("--raw", default="D:/tmp/text_raw")
    a.add_argument("--out", required=True)
    a = sub.add_parser("table")
    a.add_argument("--out", required=True)
    args = ap.parse_args()
    {"build": cmd_build, "table": cmd_table}[args.cmd](args)


if __name__ == "__main__":
    main()
