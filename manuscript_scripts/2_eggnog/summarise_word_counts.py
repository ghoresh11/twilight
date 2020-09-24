import os
import re
from itertools import combinations

def read_words_to_ignore():
    ignore = []
    with open("common_words.txt") as f:
        for line in f:
            toks = line.strip().split(",")
            ignore += toks
    return ignore

def read_classification_file(infile):
    print("Reading classification file...")
    gene_to_genetype = {}
    with open(infile) as f:
        for line in f:
            if line.startswith("gene"):
                continue
            toks = line.strip().split("\t")
            gene_to_genetype[toks[0]] = toks[8]
    return gene_to_genetype

def read_eggnog_annotation(eggnog_file, ignore, gene_to_genetype):
    ''' read an eggnot annotation file and count all the words for each
    COG category'''
    print("Reading eggnog file...")
    word_counts = {}
    for value in set(gene_to_genetype.values()):
        word_counts[value] = {}

    with open(eggnog_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            if len(toks) < 12:
                continue
            if toks[0] not in gene_to_genetype:
                print(toks[0])
                continue
           
            curr_word_counts = word_counts[gene_to_genetype[toks[0]]]
            
           
            cogs_cat = toks[-2].split(",")

            words_temp = toks[-1].lower()
            words_temp = re.sub('[^a-zA-Z0-9\- \n\.]', '', words_temp)
            print(words_temp)
            words_temp = words_temp.split()
            words = []
            for w in words_temp: ## remove common words
                if w not in ignore:
                    words.append(w)
            
            for c in cogs_cat:
                c = c.replace(" ","")
                if c not in curr_word_counts:
                    curr_word_counts[c] = {}


                ## take all word combinations
                for start, end in combinations(range(len(words)), 2):
                    phrase = " ".join(words[start:end+1])
                    if phrase not in curr_word_counts[c]:
                        curr_word_counts[c][phrase] = 0
                    curr_word_counts[c][phrase] += 1

    for key in word_counts:
        clean(word_counts[key])
    return word_counts

def clean(word_counts):
    '''clean up the word word_counts
    alot of words are substrings of others, when the counts are the
    same only keep the longer one.
    Also, remove anything that only appears 10 times or less'''
    
    for c in word_counts:
        to_delete = set()
        curr_words = list(word_counts[c].keys())
        for i in range(0, len(curr_words)-1):
            ## already been deleted
            if curr_words[i] in to_delete:
                continue
            ## ignore uncommon words
            # if word_counts[c][curr_words[i]] < 3:
            #     to_delete.add(curr_words[i])
            #     continue
            for j in range(i+1, len(curr_words)):
                ## already been removed
                if curr_words[i] in to_delete:
                    continue
                if curr_words[j] in to_delete:
                    continue
                # if word_counts[c][curr_words[j]] < 3: ## ignore uncommon words
                #     to_delete.add(curr_words[j])
                #     continue
                ## they have the same count almost exactly (+-3)
                if abs(word_counts[c][curr_words[i]] - word_counts[c][curr_words[j]]) < 3:
                    ## they're substrings of each other
                    if curr_words[i] in curr_words[j]:
                        to_delete.add(curr_words[i])
                    elif curr_words[j] in curr_words[i]:
                        to_delete.add(curr_words[j])
        for item in to_delete:
            del word_counts[c][item]
    return


def generate_outputs(word_counts):
    print("Generating outputs....")
    for gene_type in word_counts:
        curr_word_counts = word_counts[gene_type]
        out = open("words/" + gene_type.replace(" ","_") + "_words.csv", "w")
        out.write("Class\tCOG\tWord\tLabel\tCount\n")
        for c in curr_word_counts:
            for w in curr_word_counts[c]:
                label = w
                if len(w) > 30:
                    label = w[:30]
                out.write("\t".join(map(str,[gene_type.replace("_"," "), c, w, label, curr_word_counts[c][w]])) + "\n")
        out.close()
    return

def run():
    ignore = read_words_to_ignore()
    gene_to_genetype = read_classification_file("/Users/gh11/poppunk_pangenome/5_classify_genes/classification_v2.csv")
    word_counts = read_eggnog_annotation("eggnog_results.emapper.annotations", ignore, gene_to_genetype)
    generate_outputs(word_counts)
    return

if __name__ == "__main__":
    run()
