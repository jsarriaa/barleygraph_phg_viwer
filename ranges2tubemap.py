import sys
import json
import time
import glob
import os
from collections import defaultdict
import gzip

def format_time(seconds):
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    return f"{h}h {m}m {s}s" if h > 0 else f"{m}m {s}s"

def save_chromosome(chrom, data, prefix):
    out_file = f"{prefix}_{chrom}.json.gz"
    print(f"\nCompressing and saving {chrom} to {out_file}...")
    with gzip.open(out_file, 'wt', encoding='utf-8') as f:
        json.dump(data, f)
    print(f"Saved {chrom} successfully.")

def main():
    if len(sys.argv) < 5:
        print("Usage: python split_tsv_to_rails.py <input.tsv> <ref_ranges.bed> <out_prefix> <bed_folder>")
        print("Example: python split_tsv_to_rails.py hapIDranges.tsv ref_ranges.bed out_map /scratch/.../hvcf_files/")
        sys.exit(1)

    tsv_file = sys.argv[1]
    bed_file = sys.argv[2]
    out_prefix = sys.argv[3]
    bed_folder = sys.argv[4]

    # 1. LOAD ALL STRAND DATA FROM PANGENOME BED FILES
    print(f"Scanning genome BED files in {bed_folder} for strand data...")
    hash_to_strand = {}
    bed_files = glob.glob(os.path.join(bed_folder, "*.bed"))
    
    if not bed_files:
        print(f"WARNING: No .bed files found in {bed_folder}")
    
    for bf in bed_files:
        with open(bf, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split()
                if len(parts) >= 5:
                    # parts[3] = strand (+ or -)
                    # parts[4] = local checksum / hash
                    hash_to_strand[parts[4]] = parts[3]
                    
    print(f"Loaded strand information for {len(hash_to_strand)} unique haplotype blocks.")

    # 2. LOAD REFERENCE ANNOTATIONS
    print("Loading BED annotations (Interval mapping)...")
    annotations = defaultdict(list)
    with open(bed_file, 'r') as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= 4:
                chrom = p[0]
                start = int(p[1])
                end = int(p[2])
                
                # Immediately filter out 'intergenic' noise and only keep the actual genes
                genes = [g for g in p[3].split(',') if 'intergenic' not in g.lower()]
                if genes:
                    annotations[chrom].append((start, end, ",".join(genes)))
    
    # Sort for fast overlapping checks
    for chrom in annotations:
        annotations[chrom].sort(key=lambda x: x[0])

    def get_annotation(chrom, q_start, q_end):
        q_start, q_end = int(q_start), int(q_end)
        matches = set()
        for (s, e, anno) in annotations.get(chrom, []):
            if s < q_end and e > q_start: # Overlaps!
                matches.add(anno)
            if s > q_end: # Past the region, stop looking
                break
        return ",".join(matches)

    # 3. STREAM TSV AND BUILD JSON
    print("Streaming TSV file...")
    with open(tsv_file, 'r') as f:
        header = f.readline().strip().split('\t')
        samples = header[3:]
        
        current_chrom = None
        data = None
        prev_state = {s: None for s in samples}
        col_idx = 0
        start_time = time.time()

        for line_num, line in enumerate(f):
            if line.startswith('##'): continue
            cols = line.strip().split('\t')
            if len(cols) < 3: continue
            
            chrom, start, end = cols[0], cols[1], cols[2]

            if chrom != current_chrom:
                if current_chrom is not None:
                    save_chromosome(current_chrom, data, out_prefix)
                
                current_chrom = chrom
                col_idx = 0
                data = {"samples": samples, "columns": [], "segments": [], "transitions": []}
                prev_state = {s: None for s in samples}

            anno_text = get_annotation(chrom, start, end)
            
            # Group hashes for this column
            hash_groups = defaultdict(list)
            for s_idx, h in enumerate(cols[3:]):
                if h != '.': hash_groups[h].append(samples[s_idx])
            
            # Determine the consensus strand for this column (majority rule)
            col_strands = [hash_to_strand.get(h.strip('<>'), '+') for h in hash_groups.keys() if h != '.']
            consensus_strand = max(set(col_strands), key=col_strands.count) if col_strands else '+'
            
            data["columns"].append({
                "col": col_idx, 
                "chrom": chrom, 
                "start": start, 
                "end": end, 
                "anno": anno_text,
                "strand": consensus_strand  # Added to the column for the UI arrows
            })
            
            for h, samps in hash_groups.items():
                clean_h = h.strip('<>') # Remove brackets for dictionary lookup
                data["segments"].append({
                    "id": f"c{col_idx}_{h}", 
                    "col": col_idx, 
                    "hash": h, 
                    "samples": samps,
                    "strand": hash_to_strand.get(clean_h, "?")  # Now it will find the match!
                })

            if col_idx > 0:
                trans_groups = defaultdict(list)
                for s in samples:
                    prev_h = prev_state[s]
                    curr_h = cols[3+samples.index(s)]
                    if prev_h and prev_h != '.' and curr_h != '.':
                        trans_groups[(prev_h, curr_h)].append(s)
                
                for (ph, ch), samps in trans_groups.items():
                    data["transitions"].append({
                        "source": f"c{col_idx-1}_{ph}", "target": f"c{col_idx}_{ch}", "samples": samps
                    })

            for s_idx, s in enumerate(samples): prev_state[s] = cols[3+s_idx]
            col_idx += 1

            if line_num % 500 == 0:
                elapsed = time.time() - start_time
                sys.stdout.write(f"\rReading line {line_num}... | Chromosome: {current_chrom} | Elapsed: {format_time(elapsed)}")
                sys.stdout.flush()

        if current_chrom is not None:
            save_chromosome(current_chrom, data, out_prefix)

    print("\nAll chromosomes processed successfully!")

if __name__ == "__main__":
    main()