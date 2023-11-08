import sys
from cyvcf2 import VCF, Writer
import csv
import resource # limit memory usage

def memory_limit(memgb):
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (memgb * 1000 * 1000 * 1000, hard))


def main():
    depthFile = sys.argv[1]
    vcfFile = sys.argv[2]
    fname = sys.argv[3]
    minCoverage = int(sys.argv[4])
    nthread = int(sys.argv[5])
    mem = str(sys.argv[6])

    # limit memory usage
    memgb = int(mem.replace('G', ''))
    memory_limit(memgb)

    # read vcf file
    vcf = VCF(vcfFile, threads = nthread)
    # define output vcf file
    w = Writer(fname, vcf, mode = 'wz')
    
    with open(depthFile, mode ='r') as file:
        csvFile = csv.reader(file, delimiter=',', quoting=csv.QUOTE_NONE)

        for depthLine,variant in zip(csvFile, vcf):
            depthVector = []
            pos = ":".join([variant.CHROM, str(variant.POS)])
            # print(str(len(depthLine)))
            # print(str(len(variant.genotypes)))
            # double check if the position is the same
            if depthLine[0] == pos:
                for col in range(1, len(depthLine)):
                    # print(str(col)+":"+str(depthLine[col]))
                    try:
                        int(depthLine[col])
                    except ValueError:
                        depthLine[col] = 0
                    if int(depthLine[col]) < minCoverage:
                        variant.genotypes[col - 1] = [-1, -1, False]
            else:
                print("Error: position in depth file and vcf file are not the same")
                sys.exit(1)
            variant.genotypes = variant.genotypes
            w.write_record(variant)

    w.close()
    vcf.close()

    # # create a matrix of read depth
    # # the other of sites (rows) and samples (columns)
    # # are identical as vcf input
    # depthMatrix = {}
    # with open(depthFile, mode ='r') as file:
    #     csvFile = csv.reader(file, delimiter=',', quoting=csv.QUOTE_NONE)
    #     ncol = 0
    #     # nrow = 0
    #     for line in csvFile:
    #         pos = line[0]
    #         depthMatrix[pos] = []
    #         for col in range(1, len(line)):
    #             depthMatrix[pos].append(line[col])
    #         if ncol == 0:
    #             ncol = len(line) - 1
    #         # nrow += 1
    #         # depthMatrix.append(line)


    # # vcf = VCF('breedSample_breedSpecific_merged_germline_variants_SNP.vcf.gz', threads = 10)
    # vcf = VCF(vcfFile, threads = nthread)

    # # for variant in vcf('chr12:220157-220157'):
    # #     print(variant.gt_types)
    # # minCoverage = 10

    # # fname = "test_out.vcf.gz"
    # w = Writer(fname, vcf, mode = 'wz')
    # rowCounter = 0
    # for variant in vcf:
    #     pos = ":".join([variant.CHROM, str(variant.POS)])
    #     if pos in depthMatrix:
    #         for colCounter in range(ncol):
    #             # if depthMatrix[pos][colCounter] == "":
    #             #     depthMatrix[pos][colCounter] = 0
    #             if int(depthMatrix[pos][colCounter]) < minCoverage:
    #                 variant.genotypes[colCounter] = [-1, -1, False]
    #             variant.genotypes = variant.genotypes
    #             w.write_record(variant)
    # w.close()

if __name__ == "__main__":
    main()