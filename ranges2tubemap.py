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
    sample_hash_to_strand = {}
    bed_files = glob.glob(os.path.join(bed_folder, "*.bed"))
    
    if not bed_files:
        print(f"\nCRITICAL WARNING: No .bed files found in {bed_folder}")
    else:
        print(f"Found {len(bed_files)} BED files. Parsing...")
        
    for bf in bed_files:
        with open(bf, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split()
                if len(parts) >= 10:
                    strand = parts[3]
                    local_hash = parts[4]
                    sample = parts[5]  # Column 6 is the genome name (e.g., HID357)
                    ref_hash = parts[9]
                    
                    # Map the strand specific to THIS genome
                    sample_hash_to_strand[(sample, local_hash)] = strand
                    sample_hash_to_strand[(sample, ref_hash)] = strand
                    
    print(f"Loaded strand information for {len(sample_hash_to_strand)} Genome-Hash pairs.")

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
            
            # Determine the consensus strand for the COLUMN (for the UI arrow)
            # We just take a peek at the first available strand in this column
            consensus_strand = "+"
            for h in hash_groups.keys():
                if h != '.':
                    clean_h = h.strip('<>')
                    # Look up the first sample's strand for this hash to guess column direction
                    first_samp = hash_groups[h][0]
                    consensus_strand = sample_hash_to_strand.get((first_samp, clean_h), "+")
                    break

            data["columns"].append({
                "col": col_idx, 
                "chrom": chrom, 
                "start": start, 
                "end": end, 
                "anno": anno_text,
                "strand": consensus_strand
            })
            
            for h, samps in hash_groups.items():
                clean_h = h.strip('<>')
                
                # Build a dictionary of exactly which strand each genome is on
                strand_dict = {}
                for s in samps:
                    strand_dict[s] = sample_hash_to_strand.get((s, clean_h), "?")
                    
                data["segments"].append({
                    "id": f"c{col_idx}_{h}", 
                    "col": col_idx, 
                    "hash": h, 
                    "samples": samps,
                    "strands": strand_dict  # This replaces the single strand string
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