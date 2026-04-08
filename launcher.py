#!/usr/bin/env python3

import os
import yaml
import subprocess
import argparse
import re
from collections import defaultdict
from pathlib import Path
import tempfile

def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def run_pipeline(script_path, project_id, sample_id, r1_path, r2_path, base_config):
    output_dir = str(Path(base_config['global']['output_base']) / project_id / sample_id)

    # Merge base config with dynamic sample block
    full_config = dict(base_config)  # Shallow copy
    full_config.setdefault("samples", {})[sample_id] = {
        'r1': r1_path,
        'r2': r2_path,
        'output_dir': output_dir
    }

    with tempfile.NamedTemporaryFile(mode='w+', suffix='.yaml', delete=False) as temp_yaml:
        yaml.dump(full_config, temp_yaml)
        temp_yaml_path = temp_yaml.name

    print(f"[INFO] Running sample: {sample_id}")
    print(f"[DEBUG] Temp config file path: {temp_yaml_path}")

    with open(temp_yaml_path) as f:
        print("[DEBUG] Temp config contents:")
        print(f.read())

    try:
        subprocess.run(
            ['python3', base_config['global']['script_path'], '--config', temp_yaml_path, '--sample', sample_id],
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Processing failed for {sample_id}: {e}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True, help='Path to config YAML')
    args = parser.parse_args()

    base_config = load_config(args.config)
    print("Config loaded successfully")

    input_base = Path(base_config['global']['input_base'])
    print(f"[DEBUG] input_base: {input_base}")
    
    if not input_base.exists():
        print(f"[ERROR] Input base directory does not exist: {input_base}")
        return

    found_any = False

    for project_dir in input_base.iterdir():
        print(f"[DEBUG] Found in input_base: {project_dir}")
        if not project_dir.is_dir():
            continue
        project_id = project_dir.name

        # Group files by sample_id
        samples = defaultdict(dict)
        for file in project_dir.iterdir():
            match = re.match(r"(SRR\d+)_([12])\.fastq\.gz", file.name)
            if match:
                sample_id, pair = match.groups()
                samples[sample_id][pair] = file

        for sample_id, files in samples.items():
            r1_path = files.get("1")
            r2_path = files.get("2")

            if not r1_path or not r2_path:
                print(f"[WARNING] Missing R1 or R2 for sample {sample_id}")
                continue

            found_any = True
            run_pipeline(
                script_path=base_config['global']['script_path'],
                project_id=project_id,
                sample_id=sample_id,
                r1_path=str(r1_path),
                r2_path=str(r2_path),
                base_config=base_config
            )

    if not found_any:
        print("[WARNING] No valid samples found. Check your folder structure.")

if __name__ == '__main__':
    main()
