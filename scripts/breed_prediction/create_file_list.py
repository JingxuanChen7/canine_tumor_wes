#!/usr/bin/env python3

import argparse
import gzip
import sys
import os
import csv

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metatable",
        action = "store",
        dest = "metatable",
        required = True,
        help = "metadata table in qouted csv format")
    parser.add_argument("--filelist",
        action = "store",
        dest = "filelist",
        required = True,
        help = "file list for annovar, one file per line")
    parser.add_argument("--output",
        action = "store",
        dest = "output",
        required = True,
        help = "output of formatted annovar file list")
    parser.add_argument("--vcflist",
        action = "store",
        dest = "vcflist",
        required = True,
        help = "file list for original vcfs, one file per line")
    parser.add_argument("--depthlist",
        action = "store",
        dest = "depthlist",
        required = True,
        help = "file list for DepthOfCoverage BED files, one file per line")
    parser.add_argument("--vcfout",
        action = "store",
        dest = "vcfout",
        required = True,
        help = "output of formatted vcf file and depthofcoverage list")


    args = parser.parse_args()

    if os.path.exists(args.metatable):
        args.metatable = os.path.abspath(args.metatable)
    else:
        print("metatable input missing\n")
        sys.exit()
    
    if os.path.exists(args.filelist): 
        args.filelist = os.path.abspath(args.filelist)
    else:
        print("filelist input missing\n")
        sys.exit()

    if os.path.exists(args.vcflist): 
        args.vcflist = os.path.abspath(args.vcflist)
    else:
        print("vcflist input missing\n")
        sys.exit()

    if os.path.exists(args.depthlist): 
        args.depthlist = os.path.abspath(args.depthlist)
    else:
        print("depthlist input missing\n")
        sys.exit()

    args.output = os.path.abspath(args.output)
    args.vcfout = os.path.abspath(args.vcfout)


    # reading the CSV file for caseID and sampleID match
    sample2case = {}
    with open(args.metatable, mode ='r') as file:
        
        csvFile = csv.reader(file, delimiter=',', quotechar='"')
        
        # read the header line
        header = next(csvFile)
        index = {}
        for info in ['Case_ID', 'Sample_ID']:        
            if info in header:
                index[info] = header.index(info)
            else:
                sys.stderr.write(info+" column is not in provided metatable.")
                sys.exit(1)
            
        # read other lines after header
        for line in csvFile:
            Case_ID = line[index['Case_ID']]
            Sample_ID = line[index['Sample_ID']]
            sample2case[Sample_ID] = Case_ID

    # output annovar list file
    output = open(args.output, "w")
    output.write("#Sample_id\tSampleName\tGermlineAnnovarPath\n")

    with open(args.filelist, mode ='r') as file:
        for line in file:
            line = line.strip("\n")
            #initialize caseID
            caseID = None
            # get SRR sample ID from file name
            sampleID = line.split("/")[-1].split("_")[0]
            # few datasets contains complicated pairs
            datasets = [
                "/project/szlab/Kun_Lin/Pan_Cancer/Melanoma/Germline",
                "/project/szlab/Kun_Lin/Pan_Cancer/ValidationData/HSA/Germline/",
                "PRJNA680382"
            ]
            if any(dataset in line for dataset in datasets):
                # get case ID from path
                caseID = line.split("/")[-2]
            else:
                if sampleID in sample2case.keys():
                    # get corresponding case ID from meta
                    caseID = sample2case[sampleID]
            if caseID is not None:
                # write to output
                output.write(sampleID+"\t"+caseID+"\t"+line+"\n")

    output.close()

    # output vcf list file
    output2 = open(args.vcfout, "w")
    output2.write("#Sample_id\tvcf_file_path\tdepth_of_coverage_file_path\n")
    
    # create a dict of depth of coverage files
    depthfiles = {}
    with open(args.depthlist, mode ='r') as file:
        for line in file:
            line = line.strip("\n")
            # get SRR sample ID from file name
            sampleID = line.split("/")[-1].split("_")[0]
            # only keep first hit if any replicates
            if depthfiles.get(sampleID) == None:
                depthfiles[sampleID] = line

    with open(args.vcflist, mode ='r') as file:
        for line in file:
            line = line.strip("\n")
            # get SRR sample ID from file name
            sampleID = line.split("/")[-1].split("_")[0]
            # make sure the sample is in provided metatable
            # previously, the metatable can be pre-filtered
            if sampleID in depthfiles.keys() and sampleID in sample2case.keys():
                output2.write(sampleID+"\t"+line+"\t"+depthfiles[sampleID]+"\n")


    output2.close()

        
    

if __name__ == "__main__":
    main()
