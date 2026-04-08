import os
import subprocess
import sqlite3
import yaml
import argparse
import shutil
import pandas as pd
from pathlib import Path

# === Parse command-line arguments ===
parser = argparse.ArgumentParser(description="RNA Processing Pipeline")
parser.add_argument("--config", required=True, help="Path to the config.yaml file")
parser.add_argument("--sample", required=True, help="Sample name to process")
args = parser.parse_args()

# === Load Configuration ===
with open(args.config, "r") as config_file:
    config = yaml.safe_load(config_file)

# Helper to get tool parameters with fallback to sample and global
sample_config = config.get("samples", {}).get(args.sample, {})

def get_param(section_path, key, default=None):
    def deep_get(obj, path):
        for part in path.split('.'):
            obj = obj.get(part, {})
        return obj if isinstance(obj, dict) else {}

    # Try sample-specific override
    sample_val = deep_get(sample_config.get("params", {}), section_path).get(key)
    if sample_val is not None:
        return sample_val

    # Try global/tools config
    global_val = deep_get(config, section_path).get(key)
    if global_val is not None:
        return global_val

    # Use default or raise
    if default is not None:
        return default
    raise KeyError(f"[CONFIG ERROR] Missing: {section_path}.{key}")

# === Core Parameters ===
OUTPUT_DIR = config["samples"][args.sample]["output_dir"]
THREADS = config["global"].get("threads", 4)

TRIMMOMATIC_ADAPTERS = get_param("tools.trimmomatic", "adapter_file")
SORTMERNA_DB = get_param("tools.sortmerna", "db_dir")
DIAMOND_DB = get_param("tools.diamond", "db")
SQLITE_DB = get_param("taxonomy", "sqlite_db")
NODES_DMP = get_param("taxonomy", "nodes_dmp")
NAMES_DMP = get_param("taxonomy", "names_dmp")

# Validate taxonomy paths
for param_name, value in {
    "sqlite_db": SQLITE_DB,
    "nodes_dmp": NODES_DMP,
    "names_dmp": NAMES_DMP,
}.items():
    if not value:
        raise ValueError(f"[CONFIG ERROR] taxonomy.{param_name} is missing in the configuration.")

# ?? Debug Path Checks
import os
print("[DEBUG] SQLITE_DB path:", SQLITE_DB)
print("[DEBUG] Exists:", os.path.exists(SQLITE_DB))
print("[DEBUG] Readable:", os.access(SQLITE_DB, os.R_OK))
print("[DEBUG] Current Working Directory:", os.getcwd())

# === Sample Info ===
if args.sample not in config.get("samples", {}):
    raise ValueError(f"[CONFIG ERROR] Sample '{args.sample}' not found under 'samples:' section in config!")

SAMPLE_KEY = args.sample
R1_PATH = config["samples"][SAMPLE_KEY]["r1"]
R2_PATH = config["samples"][SAMPLE_KEY]["r2"]
READ_BASENAME = os.path.basename(R1_PATH).split("_")[0]

# === Taxonomy Functions ===
def load_taxdump():
    nodes, names = {}, {}
    with open(NODES_DMP, "r") as f:
        for line in f:
            fields = line.strip().split("|")
            taxid = fields[0].strip()
            parent = fields[1].strip()
            rank = fields[2].strip()
            nodes[taxid] = {"parent": parent, "rank": rank}
    with open(NAMES_DMP, "r") as f:
        for line in f:
            fields = line.strip().split("|")
            taxid = fields[0].strip()
            name = fields[1].strip()
            name_class = fields[3].strip()
            if name_class == "scientific name":
                names[taxid] = name
    return nodes, names

def get_full_taxonomy(taxid, nodes, names):
    lineage = []
    while taxid != "1" and taxid in nodes:
        rank = nodes[taxid]["rank"]
        name = names.get(taxid, f"Unknown_{taxid}")
        lineage.append(f"{rank}:{name}")
        taxid = nodes[taxid]["parent"]
    return ";".join(lineage[::-1])

# === Step 1: Preprocess Reads ===
def preprocess_reads():
    print(f" Preprocessing sample: {READ_BASENAME}")
    sample_output_dir = os.path.join(OUTPUT_DIR, READ_BASENAME)
    os.makedirs(sample_output_dir, exist_ok=True)

    for sub in ["trimmed_reads", "sortmerna_output", "merged_reads", "fasta_files", "diamond_output", "annotated_tsv", "full_annotated_tsv", "logs", "qc_reports"]:
        os.makedirs(os.path.join(sample_output_dir, sub), exist_ok=True)

    # Raw QC
    raw_qc_dir = os.path.join(sample_output_dir, "qc_reports", "raw_reads")
    os.makedirs(raw_qc_dir, exist_ok=True)
    subprocess.run(f"fastqc {R1_PATH} {R2_PATH} -o {raw_qc_dir}", shell=True, check=True)

    trimmed_r1 = os.path.join(sample_output_dir, "trimmed_reads", "trimmed_1.fastq.gz")
    trimmed_r2 = os.path.join(sample_output_dir, "trimmed_reads", "trimmed_2.fastq.gz")
    trimmomatic_log = os.path.join(sample_output_dir, "logs", "trimmomatic.log")

    trimmomatic_cmd = (
        f"trimmomatic PE -threads {THREADS} {R1_PATH} {R2_PATH} "
        f"{trimmed_r1} {sample_output_dir}/trimmed_reads/unpaired_1.fastq.gz "
        f"{trimmed_r2} {sample_output_dir}/trimmed_reads/unpaired_2.fastq.gz "
        f"ILLUMINACLIP:{TRIMMOMATIC_ADAPTERS}:2:30:10 "
        f"LEADING:{get_param('trimmomatic', 'leading', 3)} "
        f"TRAILING:{get_param('trimmomatic', 'trailing', 3)} "
        f"SLIDINGWINDOW:{get_param('trimmomatic', 'slidingwindow', '4:15')} "
        f"MINLEN:{get_param('trimmomatic', 'minlen', 36)} > {trimmomatic_log} 2>&1"
    )
    subprocess.run(trimmomatic_cmd, shell=True, check=True)

    # Trimmed QC
    trimmed_qc_dir = os.path.join(sample_output_dir, "qc_reports", "trimmed_reads")
    os.makedirs(trimmed_qc_dir, exist_ok=True)
    subprocess.run(f"fastqc {trimmed_r1} {trimmed_r2} -o {trimmed_qc_dir}", shell=True, check=True)

    # SortMeRNA
    non_rrna_fwd = os.path.join(sample_output_dir, "sortmerna_output", "non_rRNA_fwd.fq.gz")
    non_rrna_rev = os.path.join(sample_output_dir, "sortmerna_output", "non_rRNA_rev.fq.gz")
    sortmerna_log = os.path.join(sample_output_dir, "logs", "sortmerna.log")
    kvdb_path = os.path.join(sample_output_dir, "sortmerna_output", "kvdb")
    if os.path.exists(kvdb_path):
        shutil.rmtree(kvdb_path)

    sortmerna_cmd = (
        f"sortmerna --ref {SORTMERNA_DB} --reads {trimmed_r1} --reads {trimmed_r2} "
        f"--workdir {os.path.join(sample_output_dir, 'sortmerna_output')} --fastx --paired_in "
        f"--out2 --other {sample_output_dir}/sortmerna_output/non_rRNA "
        f"--aligned {sample_output_dir}/sortmerna_output/rRNA --threads {THREADS} > {sortmerna_log} 2>&1"
    )
    subprocess.run(sortmerna_cmd, shell=True, check=True)
    return sample_output_dir, non_rrna_fwd, non_rrna_rev

# === Step 2: Merge Reads ===
def merge_reads(sample_output_dir, non_rrna_fwd, non_rrna_rev):
    print(f" Merging paired reads...")
    merged_reads_prefix = os.path.join(sample_output_dir, "merged_reads", "merged_reads")
    merged_reads = f"{merged_reads_prefix}.assembled.fastq"
    pear_log = os.path.join(sample_output_dir, "logs", "pear.log")

    pear_cmd = f"pear -f {non_rrna_fwd} -r {non_rrna_rev} -o {merged_reads_prefix} > {pear_log} 2>&1"
    subprocess.run(pear_cmd, shell=True, check=True)
    
    return merged_reads

# === Step 2.5: Deduplicate Merged Reads ===
def deduplicate_reads(sample_output_dir, merged_reads):
    print(" Deduplicating merged reads using Clumpify...")

    dedup_dir = os.path.join(sample_output_dir, "deduplicated_reads")
    os.makedirs(dedup_dir, exist_ok=True)

    intermediate_out = os.path.join(dedup_dir, "merged_deduplicated_raw.fastq")  # temporary file
    final_out = os.path.join(dedup_dir, "merged_deduplicated.fastq")  # used downstream

    clumpify_log = os.path.join(sample_output_dir, "logs", "clumpify.log")

    # Step 1: Clumpify deduplication
    clumpify_cmd = f"clumpify.sh in={merged_reads} out={intermediate_out} dedupe >> {clumpify_log} 2>&1"
    subprocess.run(clumpify_cmd, shell=True, check=True)

    print(" Cleaning deduplicated reads with Fastp (trimming + overrep removal)...")

    fastp_json = os.path.join(sample_output_dir, "qc_reports", "dedup_reads_fastp.json")
    fastp_html = os.path.join(sample_output_dir, "qc_reports", "dedup_reads_fastp.html")
    fastp_log = os.path.join(sample_output_dir, "logs", "fastp_dedup.log")

    # Step 2: Fastp cleaning
    fastp_cmd = (
        f"fastp -i {intermediate_out} -o {final_out} "
        f"--trim_poly_g --trim_poly_x --detect_adapter_for_pe "
        f"--thread {THREADS} "
        f"--json {fastp_json} --html {fastp_html} > {fastp_log} 2>&1"
    )
    subprocess.run(fastp_cmd, shell=True, check=True)

    # FastQC for MultiQC inclusion
    print(" Running FastQC on cleaned deduplicated reads...")
    fastqc_outdir = os.path.join(sample_output_dir, "qc_reports", "fastqc_dedup")
    os.makedirs(fastqc_outdir, exist_ok=True)

    fastqc_cmd = f"fastqc {final_out} -o {fastqc_outdir} --threads {THREADS}"
    subprocess.run(fastqc_cmd, shell=True, check=True)

    return final_out


# === Step 3: Convert to FASTA ===
def convert_to_fasta(sample_output_dir, dedup_reads):
    print(f" Converting to FASTA...")
    merged_fasta = os.path.join(sample_output_dir, "fasta_files", f"{READ_BASENAME}.fasta")
    with open(merged_fasta, "w") as out_f:
        subprocess.run(["seqtk", "seq", "-A", dedup_reads], stdout=out_f, check=True)
    return merged_fasta


# === Step 4: Run DIAMOND ===
def run_diamond(sample_output_dir, fasta_file):
    print(f"Running DIAMOND...")

    diamond_dir = os.path.join(sample_output_dir, "diamond_output")
    os.makedirs(diamond_dir, exist_ok=True)

    # Output paths
    daa_file = os.path.join(diamond_dir, f"{READ_BASENAME}.daa")
    m8_file = os.path.join(diamond_dir, f"{READ_BASENAME}.m8")
    tsv_file = os.path.join(diamond_dir, f"{READ_BASENAME}.tsv")

    # Step 4.1: Run DIAMOND to generate .daa file
    diamond_daa_cmd = (
        f"diamond blastx --db {DIAMOND_DB} --query {fasta_file} --out {daa_file} "
        f"--evalue {config['tools']['diamond']['evalue']} "
        f"--max-target-seqs {config['tools']['diamond']['max_targets']} "
        f"--threads {THREADS} --outfmt 100"
    )
    subprocess.run(diamond_daa_cmd, shell=True, check=True)

    # Step 4.2: Convert .daa to .m8 (tabular output)
    diamond_view_cmd = (
        f"diamond view -a {daa_file} -o {m8_file} "
        f"-f 6 qseqid sseqid pident length evalue bitscore"
    )
    subprocess.run(diamond_view_cmd, shell=True, check=True)

    # Step 4.3: Create a .tsv copy of .m8
    shutil.copyfile(m8_file, tsv_file)

    return daa_file, m8_file, tsv_file

# === Step 5: Annotate Results + Taxonomy ===
def annotate_results(sample_output_dir, diamond_out):
    print(f" Annotating with taxonomy...")

    nodes, names = load_taxdump()
    conn = sqlite3.connect(SQLITE_DB)
    cursor = conn.cursor()

    df = pd.read_csv(diamond_out, sep="\t", header=None)
    df.columns = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]

    df["accession"] = df["sseqid"].str.split(".").str[0]
    taxid_map = {}
    for acc in df["accession"].unique():
        cursor.execute("SELECT taxid FROM accession_taxid WHERE accession=?", (acc,))
        result = cursor.fetchone()
        taxid_map[acc] = str(result[0]) if result else "Unknown_TaxID"
    df["taxid"] = df["accession"].map(taxid_map)

    annotated_file = os.path.join(sample_output_dir, "annotated_tsv", f"{READ_BASENAME}_annotated.tsv")
    df.to_csv(annotated_file, sep="\t", index=False, header=False)

    df["full_taxonomy"] = df["taxid"].apply(
        lambda t: get_full_taxonomy(t, nodes, names) if t != "Unknown_TaxID" else "Unknown_Taxonomy"
    )
    full_annotated = os.path.join(sample_output_dir, "full_annotated_tsv", f"{READ_BASENAME}_annotated_full_taxonomy.tsv")
    df.to_csv(full_annotated, sep="\t", index=False, header=False)
    conn.close()

    return full_annotated

# === Step 6: Filter RNA Viruses ===
def filter_rna_viruses(full_annotated_path):
    print(" Filtering RNA viruses and RNA bacteriophages...")

    df = pd.read_csv(full_annotated_path, sep="\t", header=None)

    rna_keywords = ["ssRNA", "dsRNA", "Riboviria", "Orthornavirae",
                    "Kitrinoviricota", "Pisuviricota", "Duplornaviricota",
                    "Leviviridae", "Cystoviridae", "Bacteriophage"]

    dna_exclude_keywords = ["dsDNA", "ssDNA", "DNA virus", "Caudovirales"]

    filtered_df = df[df.iloc[:, 8].str.contains('|'.join(rna_keywords), case=False, na=False)]
    filtered_rna_df = filtered_df[~filtered_df.iloc[:, 8].str.contains('|'.join(dna_exclude_keywords), case=False, na=False)]

    out_path = os.path.join(os.path.dirname(full_annotated_path), "filtered_rna_viruses.tsv")
    filtered_rna_df.to_csv(out_path, sep="\t", index=False, header=False)

    print(" RNA virus filtering complete.")
    return out_path

# === Run MultiQC ===
def run_multiqc(base_output_dir, sample_name):
    print(f" Running MultiQC to summarize QC metrics for {sample_name}...")

    multiqc_output_dir = os.path.join(base_output_dir, f"multiqc_{sample_name}")
    os.makedirs(multiqc_output_dir, exist_ok=True)

    sample_dir = os.path.join(base_output_dir, sample_name)  # Only this sample
    multiqc_report_name = f"multiqc_report_{sample_name}.html"
    multiqc_cmd = f"multiqc {sample_dir} -o {multiqc_output_dir} -n {multiqc_report_name}"
    subprocess.run(multiqc_cmd, shell=True, check=True)

    print(f" MultiQC report generated at: {multiqc_output_dir}/{multiqc_report_name}")

# === Main Runner ===
def main():
    sample_output_dir, non_rrna_fwd, non_rrna_rev = preprocess_reads()

    merged_reads = merge_reads(sample_output_dir, non_rrna_fwd, non_rrna_rev)
    dedup_reads = deduplicate_reads(sample_output_dir, merged_reads)

    fasta_file = convert_to_fasta(sample_output_dir, dedup_reads)

    # Unpack only what you need
    _, _, tsv_file = run_diamond(sample_output_dir, fasta_file)

    full_annotated = annotate_results(sample_output_dir, tsv_file)
    filter_rna_viruses(full_annotated)
    run_multiqc(OUTPUT_DIR, os.path.basename(sample_output_dir))

    print(f" Sample {os.path.basename(sample_output_dir)} processed successfully!")

if __name__ == "__main__":
    main()
