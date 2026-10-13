# barleygraph_phg_viwer

Static browser viewer for exploring pangenome graphs and their chromosomes.

## Available graphs and data layout

Choose a graph in the viewer before selecting a chromosome:

- **Pan76 (Barley)** — chromosomes `chr1H`–`chr7H` and `chrUn`
- **Wheat30 (Bread wheat)** — chromosomes `chr1A`–`chr7A` and `chr1B`–`chr7B`

Chromosome JSON files are organized by graph under `chromosomes/`:

```text
chromosomes/
├── Pan76/
│   └── Pan76_tube_map_chr*.json.gz
└── Wheat30/
    └── Wheat30_mmap_pro_chr*.json.gz
```

Pan76 gene searches use the PGSB and PanBaRT indexes. Wheat30 gene searches
use the annotations in the currently selected chromosome.

# This for the gene indexing
```
scp jsarria@161.111.227.25:~contrera/barleymap/datasets/pgsb_HC/pgsb_HC.morexv3 .
scp jsarria@161.111.227.25:~contrera/barleymap/datasets/pgsb_LC/pgsb_LC.morexv3 .
cat pgsb_LC.morexv3 pgsb_HC.morexv3 > pgsb_genes.tsv
```

# To get the json files from the hapIDranges.tsv
```
python3 ranges2tubemap.py /scratch/PHG_barleymap/Pan76/Pan76-mmap_pro/output/hapIDranges.tsv /scratch/PHG_barleymap/Pan76/Pan76-mmap_pro/output/ref_ranges.bed /agave/compbio/jsarria/Pan76_tube_map /scratch/PHG_barleymap/Pan76/Pan76-mmap_pro/vcf_dbs/hvcf_files/
```

> Enter the website:
> https://jsarriaa.github.io/barleygraph_phg_viwer/
