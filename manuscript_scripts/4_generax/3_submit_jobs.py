import subprocess
import os

all_genes = os.listdir("families/")
mem = "500"

complete = []
with open("out/done/done.txt") as f:
    for line in f:
        toks = line.strip().split()
        complete.append(toks[-1].replace(".o",""))


for g in all_genes:
    gene_name = g.split(".")[0]
    if gene_name in complete:
        continue
    lsf_prefix = ["bsub", "-q", "small", "-J", gene_name, "-G", \
        "team216", "-o", gene_name + ".o", "-e", gene_name + ".e", \
        '-R"select[mem>' + mem + '] rusage[mem=' + mem + '] span[hosts=1]"',\
         '-M' + mem]
    command = ["/nfs/users/nfs_g/gh11/GeneRax/build/bin/generax", "-f", "families/" + g, "-s", "tree.nwk", "-p", "results/" + gene_name]
    subprocess.call(lsf_prefix + command)
