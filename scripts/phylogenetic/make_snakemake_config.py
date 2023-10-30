import sys
import os
import json
import argparse
import re
# this script makes a config file for snakemake in JSON

def main():

    args = parse_args()

    config = {}
    
    ## parameters for inputs
    vcffiles = []
    with open(args.vcffilelist, "r") as file:
        for line in file:
            if not line.startswith("#"):
                line = line.strip()
                line = line.split("\t")
                vcffiles.append(line[1])



    config["in"] = {
        "vcffiles" : vcffiles,
        "metatable": args.metadata,
        "vcffilelist": args.vcffilelist,
        "breedSpecific": args.breedSpecific,
        "somaticMutation": args.somaticMutation
    }

    ## parameters for project path
    config["project_dir"] = args.project_dir
    config["breed_dir"] = args.breed_dir

    ## paths for outputs
    config["out"] = {
        "outdir": args.outdir
    }

    ## parameters for individual cluster job resources/envs
    config["resources"] = {
        "threads": args.threads,
        "mem": args.memory
    }

    ## envs
    config["envs"] = {
        "phylogenetics": "phylogenetics"
    }

    ## scripts
    config["scripts"] = {
        "mask_vcf": args.project_dir+"/scripts/phylogenetic/mask_vcf.py"
    }
   
  
    
    with open(args.out, "w") as conf:
        json.dump(config, conf, indent=4)





def parse_args():
    parser = argparse.ArgumentParser(prog='make_snakemake_config.py', description="Create config file in json format for snakemake pipeline.")

    ## required ##
    essential_args = parser.add_argument_group("Required")
    essential_args.add_argument("--project_dir", type=str, help="Path to local project code repo. Required.", required=True)
    essential_args.add_argument("--breed_dir", type=str, help="Path to local project code repo for breed prediction pipeline. Required.", required=True)
    essential_args.add_argument("--out", type=str, help="File name of the output config file in json. Required.", required=True)
    essential_args.add_argument("--outdir", type=str, help="Output (head) directory for all results. Required.", required=True)
    essential_args.add_argument("--metadata", type=str, help="Metadata table in CSV format. Required columns are 'Case_ID', 'Sample_ID', 'Bioproject', and 'Status'.", required=True)
    essential_args.add_argument("--vcffilelist", type=str, help="VCF file list used in breed prediction pipeline", required=True)
    essential_args.add_argument("--breedSpecific", type=str, help="breed specific sites created in breed prediction pipeline", required=True)
    essential_args.add_argument("--somaticMutation", type=str, help="Known somatic mutations.", required=True)

    ## resources options ##
    resources_args = parser.add_argument_group('Optional resources')
    resources_args.add_argument("--threads", type=int, help="The number of processors to use for individual cluster jobs. [default = 8]", default=8, required=False)
    resources_args.add_argument("--memory", type=str, help="The number of memory in 'G' to use for individual cluster jobs. [default = '60G']", default="60G", required=False)

    args = parser.parse_args()

    ## required ##
    # parse project repo path
    args.project_dir = os.path.abspath(args.project_dir)

    # parse breed repo path
    args.breed_dir = os.path.abspath(args.breed_dir)

    # parse output path for the json config
    args.out = os.path.abspath(args.out)

    # parse head output folder of results
    args.outdir = os.path.abspath(args.outdir)

    # parse path for metadata table
    args.metadata = os.path.abspath(args.metadata)

    # parse path for vcf file list 
    args.vcffilelist = os.path.abspath(args.vcffilelist)

    # parse path for breedSpecific file 
    args.breedSpecific = os.path.abspath(args.breedSpecific)

    # parse path for somaticMutation file 
    args.somaticMutation = os.path.abspath(args.somaticMutation)

    ## resources options ##
    # parse threads
    args.threads = int(args.threads)
    
    # parse memory
    if re.search("G$", args.memory) is None:
        sys.stderr.write("Only memory in 'G' is allowed for --memory option, e.g., '20G'")
        sys.exit(1)


    return args


if __name__ == "__main__":
    main()