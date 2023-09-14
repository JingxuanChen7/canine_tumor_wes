#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import random
import sys
from pathlib import Path

# Dictionary of IUPAC ambiguities for nucleotides
# '*' is a deletion in GATK, deletions are ignored in consensus, lowercase consensus is used when an
# 'N' or '*' is part of the genotype. Capitalization is used by some software but ignored by Geneious
# for example
AMBIG = {
    "A"    :"A", "C"    :"C", "G"    :"G", "N"    :"N", "T"    :"T",
    "*A"   :"a", "*C"   :"c", "*G"   :"g", "*N"   :"n", "*T"   :"t",
    "AC"   :"M", "AG"   :"R", "AN"   :"a", "AT"   :"W", "CG"   :"S",
    "CN"   :"c", "CT"   :"Y", "GN"   :"g", "GT"   :"K", "NT"   :"t",
    "*AC"  :"m", "*AG"  :"r", "*AN"  :"a", "*AT"  :"w", "*CG"  :"s",
    "*CN"  :"c", "*CT"  :"y", "*GN"  :"g", "*GT"  :"k", "*NT"  :"t",
    "ACG"  :"V", "ACN"  :"m", "ACT"  :"H", "AGN"  :"r", "AGT"  :"D",
    "ANT"  :"w", "CGN"  :"s", "CGT"  :"B", "CNT"  :"y", "GNT"  :"k",
    "*ACG" :"v", "*ACN" :"m", "*ACT" :"h", "*AGN" :"r", "*AGT" :"d",
    "*ANT" :"w", "*CGN" :"s", "*CGT" :"b", "*CNT" :"y", "*GNT" :"k",
    "ACGN" :"v", "ACGT" :"N", "ACNT" :"h", "AGNT" :"d", "CGNT" :"b",
    "*ACGN":"v", "*ACGT":"N", "*ACNT":"h", "*AGNT":"d", "*CGNT":"b",
    "*"    :"-", "*ACGNT":"N",
}

# define the number of information lines in the input matrix
num_info_col = 7

def extract_sample_names(vcf_file):
    """
    Extract sample names from VAF matrix
    """
    if vcf_file.lower().endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
    sample_names = []
    with opener(vcf_file, "rt") as vcf:
        for line in vcf:
            line = line.strip("\n")
            if line.startswith("Gene"):
                record = line.split("\t")
                sample_names = [record[i].replace("./", "") for i in range(num_info_col, len(record))]
                break
    return sample_names

def is_anomalous(record, num_samples):
    """
    Determine if the number of samples in current record corresponds to number of samples described
    in the line '#CHROM'
    """
    return bool(len(record) != num_samples + num_info_col)

def num_genotypes(record, num_samples):
    """
    Get number of genotypes in VCF record, total number of samples - missing genotypes
    """
    missing = 0
    for i in range(num_info_col, num_samples + num_info_col):
        if record[i].startswith("NA"):
            missing += 1
    return num_samples - missing

def is_snp(record):
    """
    Determine if current VCF record is a SNP (single nucleotide polymorphism) as opposed to MNP
    (multinucleotide polymorphism)
    """
    return bool(len(record[3]) == 1 and len(record[4]) == 1)

def get_matrix_column(record, num_samples, resolve_IUPAC, seed):
    """
    Transform a VCF record into a phylogenetic matrix column with nucleotides instead of numbers
    """
    # thresholds for homozygous sites
    homo_alt = 0.8
    homo_ref = 0.2
    # seed for random generator
    random.seed(seed)

    nt_dict = {str(0): record[3].replace("-","*").upper(), "NA": "N"}
    alt = record[4].replace("-", "*") # always 1 element
    alt = alt.split(",") # always 1 element
    for n in range(len(alt)):
        nt_dict[str(n+1)] = alt[n] # always 3 elements in dict
    column = ""
    for i in range(num_info_col, num_samples + num_info_col):
        vaf = record[i] 
        try:
            if vaf == "NA":
                geno_nuc = nt_dict["NA"]
            elif float(vaf) > homo_alt: # homozygous alt
                geno_nuc = nt_dict[str(1)]
            elif float(vaf) < homo_ref: # homozygous ref
                geno_nuc = nt_dict[str(0)]
            else:
                geno_nuc = "".join(sorted(set([nt_dict[j] for j in [str(1),str(0)]])))
        except KeyError:
            return "malformed"
        
        if resolve_IUPAC is False:
            column += AMBIG[geno_nuc]
        else:
            if vaf == "NA": 
                column += AMBIG[geno_nuc]
            else:
                # random select ref or alt based on VAF
                for index in random.choices(population = [str(1),str(0)], weights = [float(vaf), 1 - float(vaf)], k=1):
                    column += AMBIG[nt_dict[index]]
                # column += AMBIG[nt_dict[random.choices(population = [str(1),str(0)], weights = [float(vaf), 1 - float(vaf)], k=1)]]
    return column

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input",
        action = "store",
        dest = "filename",
        required = True,
        help = "Name of the input VCF file, can be gzipped")
    parser.add_argument("--output-folder",
        action = "store",
        dest = "folder",
        default = "./",
        help = "Output folder name, it will be created if it does not exist (same folder as input by "
               "default)")
    parser.add_argument("--output-prefix",
        action = "store",
        dest = "prefix",
        help = "Prefix for output filenames (same as the input VCF filename without the extension by "
               "default)")
    parser.add_argument("-m", "--min-samples-locus",
        action = "store",
        dest = "min_samples_locus",
        type = int,
        default = 4,
        help = "Minimum of samples required to be present at a locus (default=4)")
    parser.add_argument("-o", "--outgroup",
        action = "store",
        dest = "outgroup",
        default = "",
        help = "Name of the outgroup in the matrix. Sequence will be written as first taxon in the "
               "alignment.")
    parser.add_argument("-r", "--resolve-IUPAC",
        action = "store_true",
        dest = "resolve_IUPAC",
        help = "Randomly resolve heterozygous genotypes to avoid IUPAC ambiguities in the matrices "
               "(disabled by default)")
    parser.add_argument("-w", "--write-used-sites",
        action = "store_true",
        dest = "write_used",
        help = "Save the list of coordinates that passed the filters and were used in the alignments "
               "(disabled by default)")
    parser.add_argument("-s", "--seed",
        type = int,
        action = "store",
        dest = "seed",
        default = 123,
        help = "Seed for random generators. (Default: 123)")
    parser.add_argument("--select_sites",
        type = str,
        action = "store",
        dest = "select_sites",
        help = "Text file for selecting (breed-specific) sites.")
    
    args = parser.parse_args()

    outgroup = args.outgroup.split(",")[0].split(";")[0]

    # Get samples names and number of samples in VCF
    if Path(args.filename).exists():
        sample_names = extract_sample_names(args.filename)
    else:
        print("\nInput VCF file not found, please verify the provided path")
        sys.exit()
    num_samples = len(sample_names)
    if num_samples == 0:
        print("\nSample names not found in VCF, your file may be corrupt or missing the header.\n")
        sys.exit()
    print("\nConverting file '{}':\n".format(args.filename))
    print("Number of samples in VCF: {:d}".format(num_samples))

    # If the 'min_samples_locus' is larger than the actual number of samples in VCF readjust it
    args.min_samples_locus = min(num_samples, args.min_samples_locus)

    # Output filename will be the same as input file, indicating the minimum of samples specified
    if not args.prefix:
        parts = Path(args.filename).name.split(".")
        args.prefix = []
        for p in parts:
            if p.lower() == "vcf":
                break
            else:
                args.prefix.append(p)
        args.prefix = ".".join(args.prefix)
    args.prefix += ".min" + str(args.min_samples_locus)

    # Check if outfolder exists, create it if it doesn't
    if not Path(args.folder).exists():
        Path(args.folder).mkdir(parents=True)

    outfile = str(Path(args.folder, args.prefix))

    # We need to create an intermediate file to hold the sequence data vertically and then transpose
    # it to create the matrices
    # only support fasta output format
    temporal = open(outfile+".tmp", "w")

    # Process file for selected sites
    list_sites = []
    if args.select_sites:
        if Path(args.select_sites).exists():
            with open(args.select_sites, "r") as sites:
                for line in sites:
                    line = line.strip("\n")
                    if not line.startswith("Gene"):
                        record = line.split("\t")
                        chromosome, pos = record[1].split(":")
                        ref, alt = record[2].split(">")
                        list_sites.append([chromosome, pos, ref, alt])
        else:
            print("\nInput selecting-sites file not found, please verify the provided path")
            sys.exit()

    ##########################
    # PROCESS GENOTYPES IN VCF

    # if args.write_used:
    #     used_sites = open(outfile+".used_sites.tsv", "w")
    #     used_sites.write("#CHROM\tPOS\tNUM_SAMPLES\n")

    if args.filename.lower().endswith(".gz"):
        opener = gzip.open
    else:
        opener = open

    with opener(args.filename, "rt") as vcf:
        # Initialize line counter
        snp_num = 0
        snp_accepted = 0
        snp_shallow = 0
        mnp_num = 0
        snp_biallelic = 0

        while 1:
            # Load large chunks of file into memory
            vcf_chunk = vcf.readlines(50000)
            if not vcf_chunk:
                break

            for line in vcf_chunk:
                line = line.strip()

                if line and not line.startswith("Gene"): # skip the header line
                    # Split line into columns
                    record = line.split("\t")
                    # check if selected sites, skip if not selected
                    if args.select_sites:
                        if record[1:5] not in list_sites:
                            continue
                    # Keep track of number of genotypes processed
                    snp_num += 1
                    # Print progress every 500000 lines
                    if snp_num % 500000 == 0:
                        print("{:d} genotypes processed.".format(snp_num))
                    if is_anomalous(record, num_samples):
                        print("Skipping malformed line:\n{}".format(line))
                        continue
                    else:
                        # Check if the SNP has the minimum number of samples required
                        num_samples_locus = num_genotypes(record, num_samples)
                        if  num_samples_locus < args.min_samples_locus:
                            # Keep track of loci rejected due to exceeded missing data
                            snp_shallow += 1
                            continue
                        else:
                            # Check that neither REF nor ALT contain MNPs
                            if is_snp(record):
                                # Uncomment for debugging
                                # print(record)
                                # Transform VCF record into an alignment column
                                site_tmp = get_matrix_column(record, num_samples,
                                                                args.resolve_IUPAC, args.seed)
                                # Uncomment for debugging
                                # print(site_tmp)
                                # Write entire row of single nucleotide genotypes to temp file
                                if site_tmp == "malformed":
                                    print("Skipping malformed line:\n{}".format(line))
                                    continue
                                else:
                                    # Add to running sum of accepted SNPs
                                    snp_accepted += 1
                                    temporal.write(site_tmp+"\n")
                                    # if args.write_used:
                                    #     used_sites.write(record[0] + "\t"
                                    #                         + record[1] + "\t"
                                    #                         + str(num_samples_locus) + "\n")
                            else:
                                # Keep track of loci rejected due to multinucleotide genotypes
                                mnp_num += 1

        # Print useful information about filtering of SNPs
        print("Total of genotypes processed: {:d}".format(snp_num))
        print("Genotypes excluded because they exceeded the amount "
              "of missing data allowed: {:d}".format(snp_shallow))
        print("Genotypes that passed missing data filter but were "
              "excluded for being MNPs: {:d}".format(mnp_num))
        print("SNPs that passed the filters: {:d}".format(snp_accepted))

    # if args.write_used:
    #     print("Used sites saved to: '" + outfile + ".used_sites.tsv'")
    #     used_sites.close()
    print("")

    temporal.close()

    #######################
    # WRITE OUTPUT MATRICES

    output_fas = open(outfile+".fasta", "w")


    # Get length of longest sequence name
    len_longest_name = 0
    for name in sample_names:
        if len(name) > len_longest_name:
            len_longest_name = len(name)

    # Write outgroup as first sequence in alignment if the name is specified
    idx_outgroup = None
    if outgroup in sample_names:
        idx_outgroup = sample_names.index(outgroup)

        with open(outfile+".tmp") as tmp_seq:
            seqout = ""

            # This is where the transposing happens
            for line in tmp_seq:
                seqout += line[idx_outgroup]

            # Write FASTA line
            output_fas.write(">"+sample_names[idx_outgroup]+"\n"+seqout+"\n")

            # Print current progress
            print("Outgroup, '{}', added to the matrix(ces).".format(outgroup))

    # Write sequences of the ingroup
    for s in range(0, len(sample_names)):
        if s != idx_outgroup:
            with open(outfile+".tmp") as tmp_seq:
                seqout = ""

                # This is where the transposing happens
                for line in tmp_seq:
                    seqout += line[s]

                # Write FASTA line
                output_fas.write(">"+sample_names[s]+"\n"+seqout+"\n")

                # Print current progress
                print("Sample {:d} of {:d}, '{}', added to the nucleotide matrix(ces).".format(
                                                        s+1, len(sample_names), sample_names[s]))

    print()
    print("FASTA matrix saved to: " + outfile+".fasta")
    output_fas.close()
    Path(outfile+".tmp").unlink()

    print( "\nDone!\n")

if __name__ == "__main__":
    main()