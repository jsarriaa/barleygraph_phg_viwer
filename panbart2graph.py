#!/usr/bin/env python3

input_file = "/agave/compbio/jsarria/barleygraph_phg_viwer/panbart_anno.tsv"
output_file = "/agave/compbio/jsarria/barleygraph_phg_viwer/panbart_graph.tsv"

with open(input_file) as fin, open(output_file, "w") as fout:
    for line in fin:
        line = line.rstrip()

        if not line or line.startswith("#"):
            continue

        fields = line.split("\t")

        if len(fields) < 5:
            continue

        gene = fields[0]
        graph_chr = fields[1]
        graph_start = fields[2]
        graph_end = fields[3]
        strand = fields[4]

        fout.write(
            f"{gene}\t{graph_chr}\t{graph_start}\t{graph_end}\t{strand}\n"
        )
