""" 
This scripts runs the step_06_getStats.R script 
"""
import os
import sys
import subprocess
import argparse

from aspring import __version__


def test_rscript_availability():
    try:
        cmd = ["Rscript", "--version"]
        subprocess.run(cmd,
                       check=True,
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        raise Exception(
            "Rscript not found. Please install R and make sure Rscript is in your PATH."
        )


def run_r_script(gene, path_data):
    test_rscript_availability()
    aspring_path = os.path.dirname(os.path.abspath(__file__))
    if '.tox' not in aspring_path:
        package_path = os.path.abspath(os.path.join(aspring_path, "..", ".."))
        r_script_file = os.path.join(package_path, "src/aspring/R_script",
                                    "step_06_getStats.R")
        # Initialize the R environment
        subprocess.run(
            ["Rscript", "-e", "if (!require('renv')) install.packages('renv')"],
            check=True)
        subprocess.run(
            ["Rscript", "-e", "renv::restore(project='" + package_path + "')"],
            check=True)
    else:
        r_script_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 
            "R_script", "step_06_getStats.R")

    # Run the R script
    try:
        subprocess.run(
            ["Rscript", r_script_file, "--gene", gene, "--path_data", path_data],
            check=True)
    except subprocess.CalledProcessError as err:
        # reraise the exception with a more informative message
        raise Exception("R script failed.") from err


def parse_args(args):
    parser = argparse.ArgumentParser(description="Get stats for a gene")
    parser.add_argument('--gene', required=True, help='Gene name')
    parser.add_argument('--path_data',
                        required=True,
                        help='Path to dir containing Thoraxe outputs')
    parser.add_argument('--version',
                        action='version',
                        version=f"aspring {__version__}")
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])
    gene = args.gene
    path_data = args.path_data
    run_r_script(gene, path_data)


if __name__ == '__main__':
    run()
