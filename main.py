#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path

def run_snakemake(args):
    """Constructs and executes the Snakemake command."""
    
    cmd = ["snakemake"]
    
    # Core options
    cmd.extend(["--use-conda", "--use-singularity"])
    cmd.extend(["--cores", str(args.cores)])
    
    # Config file
    if args.config:
        cmd.extend(["--configfile", args.config])
        
    # Execution mode
    if args.mode == "slurm":
        cmd.extend(["--executor", "slurm"])
        cmd.extend(["--latency-wait", "60"])
        cmd.extend(["--jobs", str(args.jobs)])
        # Add singularity args for GPU if needed
        cmd.extend(["--singularity-args", "'--nv'"])
    elif args.mode == "local":
        pass # Default
        
    # Dry run
    if args.dry_run:
        cmd.append("--dry-run")
        
    # Target rules/files
    if args.targets:
        cmd.extend(args.targets)

    print(f"Running command: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Snakemake workflow failed with exit code {e.returncode}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        print("Error: 'snakemake' command not found. Please ensure it is installed and in your PATH.")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="HaloArchaea Metagenomics Workflow CLI")
    
    parser.add_argument(
        "--mode", 
        choices=["local", "slurm"], 
        default="local", 
        help="Execution mode (default: local)"
    )
    parser.add_argument(
        "--cores", 
        type=int, 
        default=1, 
        help="Number of cores to use (local mode) or max jobs (slurm mode)"
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=50,
        help="Max concurrent jobs (SLURM mode only)"
    )
    parser.add_argument(
        "--config", 
        help="Path to custom config file (default: uses Snakemake defaults)"
    )
    parser.add_argument(
        "-n", "--dry-run", 
        action="store_true", 
        help="Perform a dry run"
    )
    parser.add_argument(
        "targets", 
        nargs="*", 
        help="Specific rules or files to build"
    )

    args = parser.parse_args()
    
    run_snakemake(args)

if __name__ == "__main__":
    main()