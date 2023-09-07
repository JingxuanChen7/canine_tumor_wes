import sys
import os
import json
import argparse
import re
import csv
# this script makes a config file for snakemake in JSON
# summary snakefile for submitting individual jobs

def main():

    args = parse_args()

    config = {}
    
    ## parameters for inputs
    inputs = {}

    # reading the CSV file for metadata
    with open(args.metadata, mode ='r') as file:
        
        csvFile = csv.reader(file,delimiter=',', quotechar='"')
        
        # read the header line
        header = next(csvFile)
        index = {}
        for info in ['Case_ID', 'Sample_ID', 'Status', 'Bioproject']:        
            if info in header:
                index[info] = header.index(info)
            else:
                sys.stderr.write(info+" column is not in provided metatable.")
                sys.exit(1)
            
        # read other lines after header
        for line in csvFile:
            Case_ID = line[index['Case_ID']]
            Sample_ID = line[index['Sample_ID']]
            Status = line[index['Status']]
            Bioproject = line[index['Bioproject']]
            
            # initialize dictionary if not present
            if inputs.get(Bioproject) == None:
                inputs[Bioproject] = {}
                
            if inputs[Bioproject].get(Case_ID) == None:
                inputs[Bioproject][Case_ID] = ['NA', 'NA']
            
            if Status == "Normal":
                inputs[Bioproject][Case_ID][0] = Sample_ID
            elif Status == "Tumor":
                inputs[Bioproject][Case_ID][1] = Sample_ID

    config["in"] = inputs


    ## parameters for project path
    config["project_dir"] = args.project_dir

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
        "wesenv": "wes_env"
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
    essential_args.add_argument("--metadata", type=str, help="Metadata table in CSV format. Required columns are 'Case_ID', 'Sample_ID', 'Bioproject', and 'Status'.", required=True)

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