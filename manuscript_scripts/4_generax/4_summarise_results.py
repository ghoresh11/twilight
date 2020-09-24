import os

''' given the output from running generax, count how many events each gene had'''

all_directories = os.listdir("results/")

# out = open("all_stats.csv","w")
# out.write("Name,Likelihood,Duplication,Loss,Transfer\n")
#
# i = 1
# for d in all_directories:
#     print(i)
#     i += 1
#     filepath = os.path.join("results",d,"results", d, "stats.txt")
#     if not os.path.exists(filepath):
#         continue
#     row = 1
#     with open(filepath) as f:
#         for line in f:
#             if row == 1:
#                 likelihood = line.strip().split()[-1]
#                 row += 1
#                 continue
#             toks = line.strip().split()[3:]
#             out.write(",".join([d, likelihood] + toks) + "\n")
#
# out.close()


out = open("all_transfers.csv", "w")
out.write("Name,Node,Type\n")

i=0
for d in all_directories:
    print(i)
    i+=1
    filepath = os.path.join("results",d,"reconciliations", d + "_transfers.txt")
    if not os.path.exists(filepath):
        continue
    with open(filepath) as f:
        for line in f:
            toks = line.strip().split()
            out.write(",".join([d,toks[0],"donor"]) + "\n")
            out.write(",".join([d,toks[1],"reciever"]) + "\n")


out.close()
