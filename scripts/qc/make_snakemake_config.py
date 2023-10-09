import sys
import os
import json
import argparse
import re
import math
import csv
# this script makes a config file for snakemake in JSON
# snakefile for germline/somantic variant calling pipeline

def main():

    args = parse_args()

    config = {}
    
    ## parameters for project path
    config["project_dir"] = args.project_dir

    # reading the CSV file for metadata
    Sequence_lengths = {}
    with open(args.metadata, mode ='r') as file:
        
        csvFile = csv.reader(file,delimiter=',', quotechar='"')
            
        # read other lines after header
        for line in csvFile:
            Bioproject = line[0]
            Sequence_length = line[1]
            Sequence_lengths[Bioproject] = Sequence_length


    ## parameters for inputs
    config["in"] = {
        "Bioproject": args.Bioproject,
        "Normal_Run": args.Normal_Run,
        "Tumor_Run": args.Tumor_Run,
        "CaseName": args.CaseName,
        "Cancer_Type": Sequence_lengths[args.Bioproject]
    }

    ## paths for outputs
    config["out"] = {
        "outdir": args.outdir
    }

    ## parameters for cluster job resources/envs
    mem_int = int(re.sub("G", "", args.memory))
    mem_2parallel = math.floor(mem_int * 0.5)
    mem_4parallel = math.floor(mem_int * 0.25)
    config["resources"] = {
        "threads": args.threads,
        "mem": args.memory,
        "threads_2parallel": math.floor(args.threads * 0.5),
        "mem_2parallel": str(mem_2parallel)+"G",
        "threads_4parallel": math.floor(args.threads * 0.25),
        "mem_4parallel": str(mem_4parallel)+"G"
    }

    ## envs
    config["envs"] = {
        "wesenv": "wes_env",
        "annovarenv": "annovar_env",
        "mutect2_env": "mutect2_env",
        "strelka_env": "strelka_env",
        "java17": "java17"
    }

    ## script folders
    config["scripts"] = {
        "annovar_dir": args.annovar_dir,
        "share_dir": args.share_dir
    }
  
    
    with open(args.out, "w") as conf:
        json.dump(config, conf, indent=4)





def parse_args():
    parser = argparse.ArgumentParser(prog='make_snakemake_config.py', description="Create config file in json format for snakemake pipeline.")

    ## required ##
    essential_args = parser.add_argument_group("Required")
    essential_args.add_argument("--project_dir", type=str, help="Path to local project code repo. Required.", required=True)
    essential_args.add_argument("--out", type=str, help="File name of the output config file in json. Required.", required=True)
    essential_args.add_argument("--outdir", type=str, help="Output (head) directory for all results. Required.", required=True)
    essential_args.add_argument("--Bioproject", type=str, help="Bioproject source of the data.", required=True)
    essential_args.add_argument("--Normal_Run", type=str, help="Run accession for normal sample. Deliminated with '-' if merged run.", required=True)
    essential_args.add_argument("--Tumor_Run", type=str, help="Run accession for tumor sample. Deliminated with '-' if merged run.", required=True)
    essential_args.add_argument("--CaseName", type=str, help="Case name for the dog.", required=True)
    essential_args.add_argument("--metadata", type=str, help="Metadata table in CSV format. Required columns are 'Case_ID', 'Sample_ID', 'Bioproject', and 'Status'.", required=True)

    ## optional ##
    optional_args = parser.add_argument_group('Optional')
    # optional_args.add_argument("--coverageseq", type=str, help="[McClintock option] Sequences of the TEs for coverage module.", required=False)
    # optional_args.add_argument("--prefix", type=str, help="Prefix of output.", default="empirical", required=False)    
    optional_args.add_argument("--annovar_dir", type=str, help="Annovar script folder.", required=False)
    optional_args.add_argument("--share_dir", type=str, help="Share script folder", required=False)
    optional_args.add_argument("--readlength", type=str, help="Metadata table in CSV format for read length", required=False)

    ## resources options ##
    resources_args = parser.add_argument_group('Optional resources')
    resources_args.add_argument("--threads", type=int, help="The number of processors to use for individual cluster jobs. [default = 8]", default=8, required=False)
    resources_args.add_argument("--memory", type=str, help="The number of memory in 'G' to use for individual cluster jobs. [default = '60G']", default="60G", required=False)

    args = parser.parse_args()

    ## required ##
    # parse project repo path
    args.project_dir = os.path.abspath(args.project_dir)

    # parse output path for the json config
    args.out = os.path.abspath(args.out)

    # parse head output folder of results
    args.outdir = os.path.abspath(args.outdir)

    # parse path for metadata table
    args.metadata = os.path.abspath(args.metadata)

    ## optional ##
    # parse annovar_dir
    if args.annovar_dir is not None:
        args.annovar_dir = os.path.abspath(args.annovar_dir)
    else:
        args.annovar_dir = os.path.abspath("/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf")

    # parse share_dir
    if args.share_dir is not None:
        args.share_dir = os.path.abspath(args.share_dir)
    else:
        args.share_dir = os.path.abspath("/work/szlab/kh31516_Lab_Share_script")

    # parse readlength
    if args.readlength is not None:
        args.readlength = os.path.abspath(args.readlength)
    else:
        args.readlength = os.path.abspath(args.project_dir+"/metadata/data_new_readlength.csv")

    # # parse samples
    # for i in args.samples:
    #     if re.search("(^SRR)|(^ERR)|(^DRR)", i) is None:
    #         sys.stderr.write("Provide available SRA accession numbers")
    #         sys.exit(1)


    # parse threads
    args.threads = int(args.threads)
    
    # parse memory
    if re.search("G$", args.memory) is None:
        sys.stderr.write("Only memory in 'G' is allowed for --memory option, i.e, '20G'")
        sys.exit(1)


    return args


if __name__ == "__main__":
    main()