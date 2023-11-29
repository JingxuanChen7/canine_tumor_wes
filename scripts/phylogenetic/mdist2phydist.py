#!/usr/bin/env python3

import argparse
import gzip
import sys
import os

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mdist",
        action = "store",
        dest = "mdist",
        required = True,
        help = "mdist file from PLINK, can be gzipped")
    parser.add_argument("--id",
        action = "store",
        dest = "id",
        required = True,
        help = "id file from PLINK")
    parser.add_argument("--phydist",
        action = "store",
        dest = "phydist",
        required = True,
        help = "output distance matrix for phylip neighbor")
    parser.add_argument("--keeptaxa",
        action = "store_true",
        dest = "keeptaxa",
        default = False,
        required = False,
        help = "output distance matrix for phylip neighbor")
    args = parser.parse_args()

    if os.path.exists(args.mdist):
        args.mdist = os.path.abspath(args.mdist)
    else:
        print("\nmdist input missing\n")
        sys.exit()
    
    if os.path.exists(args.id): 
        args.id = os.path.abspath(args.id)
    else:
        print("\id input missing\n")
        sys.exit()

    args.phydist = os.path.abspath(args.phydist)


    # save sample id into a vector
    samples = []
    nSamples = 0
    with open(args.id, "r") as ids:
        for line in ids:
            line = line.strip("\n")
            record = line.split("\t")
            samples.append(record[0])
            nSamples += 1

    # output file
    output = open(args.phydist, "w")
    output.write("\t"+str(nSamples)+"\n")
    

    # parse distance matrix
    if args.mdist.lower().endswith(".gz"):
        opener = gzip.open
    else:
        opener = open

    count = 1
    # output.write(samples[0]+"     \n")
    with opener(args.mdist, "rb") as mdist:
        for line in mdist:
            line = line.decode()
            line = line.strip("\n")
            value = line.split("\t")
            #value2Float = ["{:.4f}".format(float(i)) for i in value]
            # value2Float2value = [str(i) for i in value2Float]
            if args.keeptaxa:
                value.insert(0, samples[count-1])
            else:
                value.insert(0, str(count))
            formatted_row = "".join([str(column)[:9].ljust(10) for column in value])
            output.write(formatted_row+"\n")
            count += 1
        
    output.close()

if __name__ == "__main__":
    main()
