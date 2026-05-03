#!/usr/bin/env python3
import os
import subprocess
import platform
import shutil
import sys

# Configuration
DB = "dataset/db"
INSPECT = os.path.join(DB, "inspect.txt")
INSPECT_SKIP = os.path.join(DB, "inspect.txt.skip")
BARENA = "./bin/barena"
MAMBA = ""
UNAME = platform.system()
THREADS = 8 # Adjust if needed
ENV_CANDIDATES = ("kraken", "base")
MAMBA_CANDIDATES = ("mamba", "micromamba")

# Dataset files
ZYMO_GZ = "dataset/Zymo-nanopore.fastq.gz"
ZYMO_FQ = "dataset/Zymo-nanopore.fastq"
ZYMO_KRAKEN = "dataset/zymo.kraken"

SRR_R1_GZ = "dataset/SRR19995508_R1.fastq.gz"
SRR_R2_GZ = "dataset/SRR19995508_R2.fastq.gz"
SRR_R1_FQ = "dataset/SRR19995508_R1.fastq"
SRR_R2_FQ = "dataset/SRR19995508_R2.fastq"
PAIRED_KRAKEN = "dataset/paired.kraken"

def tool_cmd(program, args=""):
    cmd = f"{MAMBA} {program}" if MAMBA else program
    if args:
        cmd = f"{cmd} {args}"
    return cmd

def can_run(cmd):
    return subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    ).returncode == 0

def resolve_mamba_prefix():
    runners = [runner for runner in MAMBA_CANDIDATES if shutil.which(runner)]
    for env in ENV_CANDIDATES:
        for runner in runners:
            prefix = f"{runner} run -n {env}"
            if (
                can_run(f"{prefix} extract_kraken_reads.py --help")
                and can_run(f"{prefix} kractor --version")
            ):
                return prefix
    return ""

def verify_startup_tools():
    global MAMBA
    MAMBA = resolve_mamba_prefix()
    if MAMBA:
        print(f"[tools] Using external tools via: {MAMBA}")
    else:
        print("[tools] No working mamba/micromamba env found; using external tools from PATH")

    checks = [
        ("extract_kraken_reads.py", tool_cmd("extract_kraken_reads.py", "--help")),
        ("barena", f"{BARENA} --version"),
        ("kractor", tool_cmd("kractor", "--version")),
    ]
    for label, cmd in checks:
        print(f"[tools] Checking {label}: {cmd}")
        if not can_run(cmd):
            raise RuntimeError(f"Cannot run {label}: {cmd}")

def run_cmd(cmd):
    print(f">>> {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def prepare():
    print("--- Preparing Database and Files ---")
    
    # 1. Generate inspect.txt if missing
    if not os.path.exists(INSPECT):
        run_cmd(f"{tool_cmd('kraken2-inspect')} --db {DB} > {INSPECT}")
    
    # 2. Prepare kractor-compatible inspect file (no header)
    run_cmd(f"grep -v '^#' {INSPECT} > {INSPECT_SKIP}")

    # 3. Unzip files for Benchmark 1 and 3 (only if not already there)
    # Note: These are huge, we check for space or existence
    if not os.path.exists(ZYMO_FQ):
        print(f"Unzipping {ZYMO_GZ}...")
        run_cmd(f"gunzip -c {ZYMO_GZ} > {ZYMO_FQ}")
    
    if not os.path.exists(SRR_R1_FQ):
        print(f"Unzipping {SRR_R1_GZ}...")
        run_cmd(f"gunzip -c {SRR_R1_GZ} > {SRR_R1_FQ}")
    if not os.path.exists(SRR_R2_FQ):
        print(f"Unzipping {SRR_R2_GZ}...")
        run_cmd(f"gunzip -c {SRR_R2_GZ} > {SRR_R2_FQ}")

    # 4. Generate Kraken output files
    if not os.path.exists(ZYMO_KRAKEN):
        print("Generating Kraken output for Zymo...")
        run_cmd(f"{tool_cmd('kraken2')} --db {DB} {ZYMO_GZ} --output {ZYMO_KRAKEN} --threads {THREADS}")

    if not os.path.exists(PAIRED_KRAKEN):
        print("Generating Kraken output for Paired SRR...")
        run_cmd(f"{tool_cmd('kraken2')} --db {DB} --paired {SRR_R1_GZ} {SRR_R2_GZ} --output {PAIRED_KRAKEN} --threads {THREADS}")

def run_benchmark(bench_id, name, fastq_input, kraken_input, taxon, is_paired=False):
    out_base = f"tests/kractor-benchmark/{bench_id}_{UNAME}"
    os.makedirs(out_base, exist_ok=True)
    
    print(f"\n=== Running Benchmark {bench_id}: {name} ===")
    
    # Barena output directory
    barena_outdir = os.path.join(out_base, "barena_out")
    os.makedirs(barena_outdir, exist_ok=True)
    
    # Fastq inputs string
    if is_paired:
        r1, r2 = fastq_input
        input_args = f"-1 {r1} -2 {r2}"
        barena_output_args = f"-c {barena_outdir}/class_#.fq"
        krakentools_args = f"-s {r1} -s2 {r2} -o {out_base}/kt_R1.fq -o2 {out_base}/kt_R2.fq"
        kractor_args = f"-i {r1} -i {r2} -o {out_base}/kr_R1.fq -o {out_base}/kr_R2.fq"
        # For barena stream
        kraken_cmd = f"{tool_cmd('kraken2')} --db {DB} --threads {THREADS} --paired {r1} {r2}"
    else:
        r1 = fastq_input
        input_args = f"-1 {r1}"
        barena_output_args = f"-c {barena_outdir}/class.fq"
        krakentools_args = f"-s {r1} -o {out_base}/kt.fq"
        kractor_args = f"-i {r1} -o {out_base}/kr.fq"
        kraken_cmd = f"{tool_cmd('kraken2')} --db {DB} --threads {THREADS} {r1}"

    # Commands for hyperfine
    cmds = [
        # Barena modes
        f'"{BARENA} {input_args} -k {kraken_input} -d {INSPECT} -t {taxon} {barena_output_args}"',
        # KrakenTools
        f'"{tool_cmd("extract_kraken_reads.py")} -k {kraken_input} -t {taxon} {krakentools_args} --fastq-output"',
        # Kractor
        f'"{tool_cmd("kractor")} -k {kraken_input} -t {taxon} -r {INSPECT_SKIP} {kractor_args}"'
     ]
    
    # Tool names
    names = [
        "barena_classified",
        "krakentools",
        "kractor"
    ]
    
    hyperfine_cmd = [
        "hyperfine",
        "--warmup 1",
        "-r 5",
        f"--prepare 'rm -f {out_base}/kr*.fq {out_base}/kt*.fq {barena_outdir}/class*.fq'",
        f"--export-json {out_base}/results.json",
        f"--export-markdown {out_base}/results.md",
        f"--export-csv {out_base}/results.csv"
    ]
    
    for n, c in zip(names, cmds):
        hyperfine_cmd.append(f'-n "{n}" {c}')
    
    full_hyperfine = " ".join(hyperfine_cmd)
    run_cmd(full_hyperfine)

def main():
    try:
        verify_startup_tools()

        # Check barena binary
        if not os.path.exists(BARENA):
            print(f"Error: {BARENA} not found. Please build it first.")
            # sys.exit(1) # We'll try to continue or user can run it after building
        
        prepare()
        
        # Benchmark 1: ONT, .fastq, taxid 96241
        run_benchmark(1, "ONT_Unpaired_Fastq", ZYMO_FQ, ZYMO_KRAKEN, 96241)
        
        # Benchmark 2: ONT, .fastq.gz, taxid 96241
        run_benchmark(2, "ONT_Unpaired_FastqGz", ZYMO_GZ, ZYMO_KRAKEN, 96241)
        
        # Benchmark 3: Paired, .fastq, taxid 590
        run_benchmark(3, "Illumina_Paired_Fastq", (SRR_R1_FQ, SRR_R2_FQ), PAIRED_KRAKEN, 590, is_paired=True)
        
        # Benchmark 4: Paired, .fastq.gz, taxid 2
        run_benchmark(4, "Illumina_Paired_FastqGz", (SRR_R1_GZ, SRR_R2_GZ), PAIRED_KRAKEN, 2, is_paired=True)
        
    except Exception as e:
        print(f"Error during benchmark: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
