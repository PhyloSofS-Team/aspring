""" 
This scripts runs the step_06_getStats.R script 
"""         
import os
import sys
import subprocess

def parse_args(args):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gene', required=True, help='gene name')
    parser.add_argument('--out', required=True, help='output directory')
    parser.add_argument('--rscript', required=True, help='Rscript path')

    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    gene = args.gene
    out = args.out
    rscript = args.rscript

    # run R script
    cmd = [rscript, 'step_06_getStats.R', gene, out]
    subprocess.call(cmd)
    
if __name__ == '__main__':
    main(sys.argv[1:])



