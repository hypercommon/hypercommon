"""
Sequence-only control for the protein graphs: on the same residues, residue i is joined to
the residues within h positions along the chain, with h chosen so that the mean degree
(about 2h) matches the contact graph's; the groups are the same CATH domains. Written as
<chain>_null/ next to each 8 A graph (data_protein/).

    python protein_null.py
"""
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data_protein")


def main():
    for name in sorted(os.listdir(D)):
        src = os.path.join(D, name)
        if "_" in name or not os.path.isfile(os.path.join(src, "edges.tsv")):
            continue
        E = [l.split() for l in open(os.path.join(src, "edges.tsv")) if l.strip()]
        E = [(int(a), int(b)) for a, b in E if a.isdigit()]
        n = 1 + max(max(e) for e in E)
        h = max(1, round(len(E) / n))            # mean degree 2m/n, so about 2h = 2m/n
        dst = os.path.join(D, name + "_null")
        os.makedirs(dst, exist_ok=True)
        with open(os.path.join(dst, "edges.tsv"), "w", newline="\n") as fh:
            for i in range(n):
                for j in range(i + 1, min(n, i + h + 1)):
                    fh.write(f"{i}\t{j}\n")
        for f in ("groups.txt", "nodes.tsv"):
            shutil.copy(os.path.join(src, f), os.path.join(dst, f))
        print(name, n, len(E), h)


if __name__ == "__main__":
    main()
