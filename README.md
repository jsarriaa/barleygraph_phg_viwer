# barleygraph_phg_viwer
Hosting an online browser tool to visualize the global barley pangenome with Practical Haplotype Graph

# This for the gene indexing
```
scp jsarria@161.111.227.25:~contrera/barleymap/datasets/pgsb_HC/pgsb_HC.morexv3 .
scp jsarria@161.111.227.25:~contrera/barleymap/datasets/pgsb_LC/pgsb_LC.morexv3 .
cat pgsb_LC.morexv3 pgsb_HC.morexv3 > pgsb_genes.tsv
```

# To get the json files from the hapIDranges.tsv
```
python3 ranges2tubemap.py Pan76-mmap_pro/output/hapIDranges.tsv Pan76-mmap_pro/output/ref_ranges.bed /agave/compbio/jsarria/Pan76_tube_map.json
```
