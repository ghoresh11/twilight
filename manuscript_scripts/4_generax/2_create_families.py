import os

all_genes = os.listdir("links/")

for g in all_genes:
    gene_name = g.split(".")[0]
    with open("families/" + gene_name + ".txt", "w") as out:
        out.write("[FAMILIES]\n")
        out.write("- " + gene_name + "\n")
        out.write("starting_gene_tree = /lustre/scratch118/infgen/team216/gh11/gene_distances/trees/" + gene_name + ".treefile\n")
        out.write("alignment = /lustre/scratch118/infgen/team216/gh11/gene_distances/msas/" + gene_name + "_msa.fa\n")
        out.write("mapping = /lustre/scratch118/infgen/team216/gh11/gene_distances/test_generax/links/" + gene_name + ".link\n")
        out.write("subst_model = GTR+G\n")
