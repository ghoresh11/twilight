import networkx as nx
import os

## step 1: for each gene, say which clusters:
clusters_per_gene = {}
with open("/Users/gh11/poppunk_pangenome/4_pairwise_roary/231019_corrected/freqs.csv") as f:
    for line in f:
        toks = line.strip().split(",")
        if toks[0] == "Gene":
            clusters = toks[1:]
            continue
        curr_gene = toks[0]
        for c in ("(",")","'","/"):
            curr_gene = curr_gene.replace(c, "_")
        curr_freqs = list(map(float,toks[1:]))
        clusters_per_gene[curr_gene] = {}
        for i in range(0, len(curr_freqs)):
            if curr_freqs[i] == 0:
                continue
            clusters_per_gene[curr_gene][clusters[i]] = curr_freqs[i] ## save frequency of a gene in a cluster

## Step 2: create a network with all the clusters and their metadata
cluster_details = {}
with open("/Users/gh11/Submissions/my_thesis/Chapter4/figures/cluster_graphics.csv") as f:
    for line in f:
        toks = line.strip().split(",")
        cluster_details[toks[0]] = {"Phylogroup":toks[1]}
with open("/Users/gh11/Submissions/my_thesis/Chapter4/figures_thesis/ecoli_cluster_summary.csv", encoding='utf-8-sig') as f:
    for line in f:
        if line.startswith("Cluster"):
            continue
        toks = line.strip().split(",")
        cluster_details[toks[0]]["ST"] = toks[1]
        cluster_details[toks[0]]["MDR"] = toks[2]
        cluster_details[toks[0]]["hypervirulent"] = toks[3]
        cluster_details[toks[0]]["Pathotype"] = toks[4]

G = nx.Graph()
for c in cluster_details:
    if c in ["43","49","21","Cluster"]:
        continue
    G.add_node(c, ST = cluster_details[c]["ST"], Pathotype = cluster_details[c]["Pathotype"],\
    MDR = cluster_details[c]["MDR"], Phylogroup = cluster_details[c]["Phylogroup"],\
    hypervirulent = cluster_details[c]["hypervirulent"])

## find distance between every two clusters
distances = {}
with open("/Users/gh11/Submissions/my_thesis/Chapter4/figures/phylo_dist.csv") as f:
    for line in f:
        if line.startswith("iso"):
            continue
        toks = line.strip().split("\t")
        distances[toks[3]] = float(toks[2])


gene_type_graphs = {}
## in addition to making a graph for each gene type, make one
gene_type_graphs["all"] = G.copy()
gene_type_graphs["distant"] = G.copy()

## Step 3: add an edge for each type of gene class with the number of genes two cluster share
with open("/Users/gh11/poppunk_pangenome/9_gene_properties/classification_w_props_v2.csv") as f:
    for line in f:
        if line.startswith("gene"):
            continue
        toks = line.strip().split("\t")
        curr_gene = toks[0]
        gene_type = toks[6]
        gene_type2 = toks[8]
        num_members = toks[4]
        mean_freq = toks[9]


        if gene_type not in gene_type_graphs:  ## make multiple graphs, for all and for
            gene_type_graphs[gene_type] = G.copy()

        if gene_type2 not in gene_type_graphs:  ## make multiple graphs, for all and for
            gene_type_graphs[gene_type2] = G.copy()

        curr_clusters = list(clusters_per_gene[curr_gene].keys())


        for i in range(0,len(curr_clusters)-1):
            for j in range(i+1,len(curr_clusters)):
                key_dist = curr_clusters[j]+"-"+curr_clusters[i]
                if key_dist not in distances:
                    key_dist = curr_clusters[i]+"-"+curr_clusters[j]
                for key in (gene_type, gene_type2, "all", "distant"):
                    if not gene_type_graphs[key].has_edge(curr_clusters[i], curr_clusters[j]):
                        gene_type_graphs[key].add_edge(curr_clusters[i], curr_clusters[j], weight = 0)
                    if key == "distant":
                        if gene_type == "Varied" and distances[key_dist] > 0.3: ## add to distant network only in these cases
                            gene_type_graphs[key][curr_clusters[i]][curr_clusters[j]]["weight"] += 1
                    else:
                        gene_type_graphs[key][curr_clusters[i]][curr_clusters[j]]["weight"] += 1 ## count how many genes they share, frequency doesnt

## add edges with 0 to the networks
for g in gene_type_graphs:
    for c1 in cluster_details:
        for c2 in cluster_details:
            if c1 == c2 or c1 in ["21","43","49","Cluster"] or c2 in ["21","43","49","Cluster"]:
                continue
            if not gene_type_graphs[g].has_edge(c1, c2):
                gene_type_graphs[g].add_edge(c1, c2, weight = 0)

## Step 4: Write the graph to a file
out = open("num_shared_gene_v2.csv","w")
out.write("ClusterA\tClusterB\tGene_class\tCount\n")
for g in gene_type_graphs:
    nx.write_gml(gene_type_graphs[g], os.path.join("graphs_v2", g.replace(" ","_")) + ".gml")
    for e in gene_type_graphs[g].edges(data = True):
        out.write("\t".join([e[0],e[1],g,str(e[2]["weight"])]) + "\n")
out.close()
