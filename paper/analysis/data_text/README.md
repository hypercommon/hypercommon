# Wikipedia sentence graphs with nested sections as ground truth

30 graphs and 30 position-only controls (60 folders): 15 English Wikipedia
featured articles, each at k = 8 and k = 16. Nodes are sentences in reading
order; each sentence is joined to its k most similar sentences (TF-IDF cosine)
among the 40 sentences before and after it; the ground-truth groups are the
article's sections and subsections, nested like Hi-C contact domains. Built by
`paper/analysis/text_build.py` with the same top-k-along-a-1D-order recipe as
`hic_build.py`.

## Sources

| What | URL |
|---|---|
| Article list | `Category:Featured articles`, namespace 0, via `https://en.wikipedia.org/w/api.php?action=query&list=categorymembers&cmtitle=Category:Featured_articles&cmnamespace=0` (7,031 members, fetched 2026-11-13 00:05 UTC) |
| Article text | `https://en.wikipedia.org/w/api.php?action=query&prop=extracts\|revisions&explaintext=1&exsectionformat=wiki&rvprop=ids\|timestamp&pageids=<id>` (TextExtracts plain text with `== heading ==` lines; the revision id comes from the same request) |

All requests went to the public MediaWiki API, one at a time with at least 2 s
between them, under a generic descriptive User-Agent with no personal details;
no account or login was used. The API rate-limited us a few times (HTTP 429);
the script waits (Retry-After or exponential back-off) and retries. 92
successful requests in total: 15 for the category list (500 titles per page)
and 77 articles (one each, listed in `selection_log.tsv`).
The raw JSON is cached under `D:/tmp/text_raw/` (`featured_articles.json`,
`articles/<pageid>.json`) and is not in the repository.

### Licence and attribution

Wikipedia text is © the Wikipedia contributors and licensed under the Creative
Commons Attribution-ShareAlike 4.0 International licence (CC BY-SA 4.0),
https://creativecommons.org/licenses/by-sa/4.0/ (see
https://en.wikipedia.org/wiki/Wikipedia:Copyrights). The graphs are derived
from the revisions listed below; the authors of each article are listed in its
page history (`https://en.wikipedia.org/w/index.php?title=<title>&action=history`).
The repository files hold no sentence text, only node ids, positions and
section headings (in `nodes.tsv`); as adaptations they are shared under the same
licence, CC BY-SA 4.0. Changes made: text split into sentences, end matter
removed, sentences turned into a similarity graph.

| Article | page id | revision id (permanent link: `https://en.wikipedia.org/w/index.php?oldid=<revid>`) | revision date | sentences | n (kept) | groups | heading depth | folder stem |
|---|---|---|---|---|---|---|---|---|
| Don Bradman | 87021 | 1385116950 | 2026-11-11 | 588 | 578 | 29 | 3 | wiki_don_bradman |
| Juan Manuel de Rosas | 490071 | 1369849474 | 2026-08-17 | 338 | 336 | 20 | 2 | wiki_juan_manuel_de_rosas |
| Parallel computing | 145162 | 1378851310 | 2026-10-06 | 326 | 306 | 34 | 4 | wiki_parallel_computing |
| Imperial Trans-Antarctic Expedition | 2183356 | 1378298421 | 2026-10-03 | 330 | 328 | 18 | 3 | wiki_imperial_trans_antarctic_expedition |
| Ulysses S. Grant | 31752 | 1384335961 | 2026-11-07 | 762 | 759 | 44 | 3 | wiki_ulysses_s_grant |
| Lesley J. McNair | 1248450 | 1383755385 | 2026-11-04 | 380 | 375 | 37 | 4 | wiki_lesley_j_mcnair |
| PlayStation (console) | 23964 | 1385233566 | 2026-11-12 | 344 | 342 | 21 | 3 | wiki_playstation_console |
| British logistics in the Normandy campaign | 56060616 | 1339840786 | 2026-02-22 | 428 | 427 | 20 | 2 | wiki_british_logistics_in_the_normandy_campai |
| Raymond Leane | 6307138 | 1379854800 | 2026-10-12 | 384 | 382 | 15 | 4 | wiki_raymond_leane |
| First-move advantage in chess | 16365141 | 1358527793 | 2026-06-09 | 439 | 397 | 19 | 4 | wiki_first_move_advantage_in_chess |
| Sexuality after spinal cord injury | 46895003 | 1368365059 | 2026-08-08 | 306 | 305 | 21 | 3 | wiki_sexuality_after_spinal_cord_injury |
| Draining and development of the Everglades | 17375601 | 1374787786 | 2026-09-14 | 330 | 328 | 17 | 2 | wiki_draining_and_development_of_the_everglad |
| Ludwigsburg Palace | 849697 | 1374061164 | 2026-09-09 | 300 | 300 | 13 | 2 | wiki_ludwigsburg_palace |
| Peking opera | 618649 | 1381868948 | 2026-10-25 | 393 | 390 | 24 | 2 | wiki_peking_opera |
| Battle of Saipan | 1129826 | 1375586871 | 2026-09-18 | 333 | 330 | 25 | 3 | wiki_battle_of_saipan |

"sentences" = all sentences after end matter is removed; "n (kept)" = those
with at least 3 content tokens. Folder names are `<stem>_k8`, `<stem>_k16` and
the controls `<stem>_k8_null`, `<stem>_k16_null` (stems are cut at 40
characters).

## Selection rule (fixed before building any graph)

1. Take all namespace-0 members of `Category:Featured articles` (fetched once
   and cached), sort them by title, and shuffle them with Python
   `random.Random(2026)`.
2. Walk the shuffled list. For each article, build its nodes and groups
   (below) and accept it if 300 <= n <= 3000 and it has >= 5 groups. Stop at
   15 accepted articles.
3. Result: 77 articles were tried; 15 accepted, 62 rejected, all for n < 300
   (most featured articles have 80-300 kept sentences). `selection_log.tsv`
   lists every article tried, with its revision id, n, group count and the
   reason.

The rule looks only at length and section count, never at any method's
score. No community detection was run while building these graphs.

## Construction

1. **Sections.** The extract is cut at its heading lines
   (`== H2 ==` ... `====== H6 ======`). Text before the first heading is the
   lead. A section whose heading (case-insensitive) is end matter (`See also`,
   `Notes`, `References`, `Sources`, `Citations`, `Bibliography`,
   `Further reading`, `External links`, `Footnotes`, `Works cited` and similar;
   the full list is `END_MATTER` in the script) is dropped together with its
   subsections. TextExtracts already leaves out tables, infoboxes, image
   captions and footnote markers.
2. **Sentences.** Each non-empty line (paragraph or list item) is split by a
   fixed regex: break after `.`, `!` or `?` (plus any closing quotes or
   brackets) followed by white space and an upper-case letter, digit or opening
   quote/bracket; no break after a single-letter initial (`Ulysses S. Grant`)
   or a listed abbreviation (`Mr.`, `Dr.`, `St.`, `U.S.`, `c.`, month names,
   ...). No NLP model is used.
3. **Nodes.** A sentence is kept if it has >= 3 content tokens (lower-cased,
   scikit-learn default token pattern `\b\w\w+\b`, scikit-learn English stop
   words removed); 5,883 of 5,981 sentences (98.4%) were kept over the 15
   articles. Kept sentences are numbered 0..n-1 in reading order. Why: shorter
   fragments ("He was 13.", list stubs) have almost empty TF-IDF vectors and
   would join to arbitrary neighbours.
4. **Similarity.** TF-IDF fitted on the article's own kept sentences
   (`TfidfVectorizer(lowercase=True, stop_words="english")`, smooth idf,
   l2-normalised rows); similarity = cosine.
5. **Edges (exact rule).** For each node i, the candidates are the nodes j with
   0 < |i - j| <= W = 40 and similarity > 0. Keep the k candidates with the
   largest similarity (ties: smaller |i - j|, then smaller j). The undirected,
   unweighted edge set is the union of these choices over all i (degree >= k
   unless fewer than k candidates have similarity > 0). Each edge written once
   as `u<TAB>v`, u < v; no self-loops. k = 8 and k = 16.
6. **Groups.** The lead and every kept section at every heading level, each as
   the set of kept sentences inside it, including its subsections (so groups
   nest). Dropped: groups with < 3 nodes (7 over all articles, short leads
   or subsections), exact duplicates of an earlier group (a section whose only
   content is one subsection; 1) and a group equal to the whole
   article (none). Every top-level section is a group, so groups cover all
   nodes except the sentences of a dropped short lead (frac_nodes_in_group
   0.995-1.000).
7. **Control (`_null`).** Same nodes, same `groups.txt` and `nodes.tsv`;
   each node keeps its k nearest nodes by position (|i - j|, ties to the
   smaller id), union over i. A pure band with no content signal, as in
   `hic_null.py`.

Why these choices: sentences give n of a few hundred per article, inside the
300-3000 target, so no windowing was needed. The +-40 window plays the role of
the Hi-C distance decay: without it the TF-IDF neighbours of a sentence would
be scattered over the whole article (the lead in particular repeats every
section's key words), and the graph would no longer be ordered along a line.
40 is about two median-sized sections. k = 8 and 16 match the Hi-C k values
that keep m small. All of these were set before any graph was built.

## Files

Each folder has:

- `edges.tsv` - one undirected edge per line, `u<TAB>v`, ids 0..n-1, u < v, no header.
- `groups.txt` - one group per line, space-separated sorted node ids, size >= 3,
  in document order (lead first, then each section before its subsections).
- `nodes.tsv` - header `id position sentence_index section_path`; `position`
  = id (rank in reading order among kept sentences), `sentence_index` = rank
  among all sentences before the 3-token filter, `section_path` = headings from
  top level down joined by ` > `, or `(lead)`.

Also `summary.json` (one row per graph: n, m, mean degree, clustering,
n_components, n_isolated, median edge length in positions, fraction of edges
<= 5 positions, n_groups and group sizes, frac_nodes_in_group,
frac_edges_within_group, frac_edges_same_leaf, group-drop counts, and the same
graph statistics for the control prefixed `null_`) and `selection_log.tsv`.

## Reproduce

Python 3.12.10 in a separate venv with numpy 2.5.4, scipy 1.18.1, scikit-learn
1.9.1, networkx 3.7, requests 2.34.2. From `paper/`:

    python -I analysis/text_build.py build --raw D:/tmp/text_raw --out analysis/data_text
    python -I analysis/text_build.py table --out analysis/data_text   # table below

With the cache present nothing is fetched and the output is byte-identical
(checked). Without the cache the script fetches the category list and the
articles again; both change over time (new featured articles, new revisions),
so the selection could differ. The cached revision ids above pin the text.

## Statistics

"within" = fraction of edges whose two ends share at least one group (in
practice: lie in the same top-level section, or both in the lead). "same leaf"
= both ends have the same deepest group. Edge length = |i - j| in positions.
Clustering is networkx `average_clustering`. "components (isolated)" counts
connected components and, in brackets, isolated nodes; all controls are
connected.

| graph | n | m | mean degree | clustering (control) | components (isolated) | median edge length (control) | groups | group size min-median-max | within | within, control | same leaf | same leaf, control |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| wiki_don_bradman_k8 | 578 | 2774 | 9.60 | 0.266 (0.650) | 14 (13) | 16 (3) | 29 | 3-22-316 | 0.800 | 0.966 | 0.443 | 0.880 |
| wiki_don_bradman_k16 | 578 | 4960 | 17.16 | 0.388 (0.711) | 14 (13) | 17 (5) | 29 | 3-22-316 | 0.779 | 0.938 | 0.398 | 0.795 |
| wiki_juan_manuel_de_rosas_k8 | 336 | 1446 | 8.61 | 0.301 (0.655) | 9 (8) | 16 (3) | 20 | 13-23.5-83 | 0.648 | 0.948 | 0.433 | 0.897 |
| wiki_juan_manuel_de_rosas_k16 | 336 | 2469 | 14.70 | 0.414 (0.719) | 9 (8) | 17 (5) | 20 | 13-23.5-83 | 0.620 | 0.907 | 0.384 | 0.815 |
| wiki_parallel_computing_k8 | 306 | 1400 | 9.15 | 0.349 (0.657) | 4 (3) | 11 (3) | 34 | 3-9.5-97 | 0.750 | 0.922 | 0.446 | 0.750 |
| wiki_parallel_computing_k16 | 306 | 2232 | 14.59 | 0.417 (0.721) | 4 (3) | 14 (5) | 34 | 3-9.5-97 | 0.720 | 0.855 | 0.362 | 0.589 |
| wiki_imperial_trans_antarctic_expedition_k8 | 328 | 1553 | 9.47 | 0.282 (0.656) | 4 (3) | 15 (3) | 18 | 9-21.5-220 | 0.804 | 0.970 | 0.446 | 0.894 |
| wiki_imperial_trans_antarctic_expedition_k16 | 328 | 2630 | 16.04 | 0.390 (0.720) | 4 (3) | 16 (5) | 18 | 9-21.5-220 | 0.796 | 0.946 | 0.418 | 0.811 |
| wiki_ulysses_s_grant_k8 | 759 | 3711 | 9.78 | 0.304 (0.648) | 4 (3) | 13 (3) | 44 | 6-21-284 | 0.869 | 0.974 | 0.515 | 0.862 |
| wiki_ulysses_s_grant_k16 | 759 | 6732 | 17.74 | 0.412 (0.709) | 4 (3) | 16 (5) | 44 | 6-21-284 | 0.844 | 0.953 | 0.437 | 0.758 |
| wiki_lesley_j_mcnair_k8 | 375 | 1866 | 9.95 | 0.314 (0.654) | 5 (4) | 11 (3) | 37 | 3-13-122 | 0.700 | 0.922 | 0.434 | 0.785 |
| wiki_lesley_j_mcnair_k16 | 375 | 3388 | 18.07 | 0.394 (0.717) | 5 (4) | 15 (5) | 37 | 3-13-122 | 0.659 | 0.854 | 0.328 | 0.628 |
| wiki_playstation_console_k8 | 342 | 1739 | 10.17 | 0.260 (0.655) | 3 (2) | 14 (3) | 21 | 5-17-166 | 0.818 | 0.964 | 0.496 | 0.869 |
| wiki_playstation_console_k16 | 342 | 3185 | 18.63 | 0.362 (0.719) | 3 (2) | 16 (5) | 21 | 5-17-166 | 0.795 | 0.932 | 0.437 | 0.754 |
| wiki_british_logistics_in_the_normandy_campai_k8 | 427 | 2063 | 9.66 | 0.310 (0.653) | 5 (4) | 13 (3) | 20 | 4-23.5-201 | 0.795 | 0.965 | 0.531 | 0.889 |
| wiki_british_logistics_in_the_normandy_campai_k16 | 427 | 3260 | 15.27 | 0.389 (0.715) | 5 (4) | 15 (5) | 20 | 4-23.5-201 | 0.772 | 0.928 | 0.484 | 0.799 |
| wiki_raymond_leane_k8 | 382 | 1865 | 9.76 | 0.257 (0.654) | 3 (2) | 14 (3) | 15 | 8-28-259 | 0.875 | 0.974 | 0.565 | 0.909 |
| wiki_raymond_leane_k16 | 382 | 3339 | 17.48 | 0.356 (0.717) | 3 (2) | 17 (5) | 15 | 8-28-259 | 0.872 | 0.951 | 0.518 | 0.835 |
| wiki_first_move_advantage_in_chess_k8 | 397 | 1930 | 9.72 | 0.275 (0.653) | 5 (3) | 15 (3) | 19 | 4-17-213 | 0.817 | 0.950 | 0.585 | 0.889 |
| wiki_first_move_advantage_in_chess_k16 | 397 | 3474 | 17.50 | 0.372 (0.716) | 5 (3) | 17 (5) | 19 | 4-17-213 | 0.778 | 0.914 | 0.518 | 0.814 |
| wiki_sexuality_after_spinal_cord_injury_k8 | 305 | 1559 | 10.22 | 0.277 (0.657) | 1 (0) | 12 (3) | 21 | 4-19-114 | 0.700 | 0.943 | 0.442 | 0.854 |
| wiki_sexuality_after_spinal_cord_injury_k16 | 305 | 2759 | 18.09 | 0.374 (0.721) | 1 (0) | 14 (5) | 21 | 4-19-114 | 0.662 | 0.898 | 0.381 | 0.742 |
| wiki_draining_and_development_of_the_everglad_k8 | 328 | 1491 | 9.09 | 0.287 (0.656) | 4 (3) | 14 (3) | 17 | 9-23-86 | 0.702 | 0.947 | 0.508 | 0.879 |
| wiki_draining_and_development_of_the_everglad_k16 | 328 | 2328 | 14.20 | 0.396 (0.720) | 4 (3) | 15 (5) | 17 | 9-23-86 | 0.679 | 0.905 | 0.459 | 0.789 |
| wiki_ludwigsburg_palace_k8 | 300 | 1505 | 10.03 | 0.273 (0.657) | 1 (0) | 12 (3) | 13 | 7-27-129 | 0.882 | 0.975 | 0.524 | 0.901 |
| wiki_ludwigsburg_palace_k16 | 300 | 2553 | 17.02 | 0.394 (0.722) | 1 (0) | 15 (5) | 13 | 7-27-129 | 0.864 | 0.956 | 0.469 | 0.812 |
| wiki_peking_opera_k8 | 390 | 1816 | 9.31 | 0.342 (0.654) | 3 (2) | 12 (3) | 24 | 6-20-97 | 0.816 | 0.943 | 0.519 | 0.867 |
| wiki_peking_opera_k16 | 390 | 2791 | 14.31 | 0.434 (0.717) | 3 (2) | 13 (5) | 24 | 6-20-97 | 0.794 | 0.899 | 0.470 | 0.766 |
| wiki_battle_of_saipan_k8 | 330 | 1618 | 9.81 | 0.274 (0.656) | 5 (4) | 14 (3) | 25 | 4-19-157 | 0.745 | 0.932 | 0.409 | 0.842 |
| wiki_battle_of_saipan_k16 | 330 | 2775 | 16.82 | 0.382 (0.720) | 5 (4) | 15 (5) | 25 | 4-19-157 | 0.720 | 0.873 | 0.366 | 0.716 |

Medians over the 15 articles:

| | within | within, control | same leaf | same leaf, control | clustering | clustering, control |
|---|---|---|---|---|---|---|
| k = 8 | 0.800 | 0.950 | 0.496 | 0.879 | 0.282 | 0.655 |
| k = 16 | 0.778 | 0.914 | 0.437 | 0.789 | 0.390 | 0.719 |

## Caveats

- **Content edges cross section borders more than the band does.** Unlike
  Hi-C (where real edges fell inside domains 2-4 points more often than the
  band null), here the real graphs have fewer within-group edges than their
  controls (median 0.80 vs 0.95 at k = 8). TF-IDF neighbours are spread over
  the whole +-40 window (median edge length 11-17 positions against 3-5 for
  the band), so many reach into the next section; prose is also repetitive in
  names (a biography's subject appears in every section). The graphs are much
  less lattice-like than Hi-C or the band: clustering 0.26-0.43 against
  0.65-0.72. This is a property of the fixed rule, reported as found; nothing
  was tuned.
- **Isolated nodes.** 0-13 nodes per graph (same at k = 8 and 16) have no
  candidate with similarity > 0 inside the window and are chosen by no one
  (e.g. statistical lines in "Don Bradman"). They are kept, as Hi-C keeps
  every covered bin; most graphs therefore have a few components, nearly all
  of them single nodes.
- **Groups are nested and unequal.** Top-level sections hold up to 316
  sentences (the "Test career" section of "Don Bradman"), while
  subsections can be 3 sentences; 168-629 nodes per article are in more than
  one group. Section boundaries are editorial, not topical ground truth, and
  the lead summarises the whole article.
- **Text extraction.** TextExtracts drops tables, infoboxes and most lists'
  markup, but keeps some list items as lines; the regex splitter can miss or
  add a break at unusual abbreviations. 1.6% of sentences were dropped by the
  3-token filter.
- **Selection.** "Long" here means >= 300 kept sentences, which only about one
  featured article in five reaches (15 of 77 tried); the sample is therefore
  biased towards long biographies and histories. One language, one snapshot.

## Music sub-family

Not built. The folder budget (at most about 60 including controls) is used by
the 15 text articles, and the music corpus was optional.
