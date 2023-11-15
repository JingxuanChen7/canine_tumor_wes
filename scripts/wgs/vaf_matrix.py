import sys
from cyvcf2 import VCF
import resource # limit memory usage
import gzip

def memory_limit(memgb):
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (memgb * 1000 * 1000 * 1000, hard))

def main():
    vcfFile = sys.argv[1]
    fname = sys.argv[2]
    nthread = int(sys.argv[3])
    mem = str(sys.argv[4])

    # limit memory usage
    memgb = int(mem.replace('G', ''))
    memory_limit(memgb)

    # read vcf file
    vcf = VCF(vcfFile, threads = nthread)
    # define output vaf matrix file

    with gzip.open(fname, mode ='wt') as file:
        # write header
        header = ["Chromosome", "Position", "Ref", "Alt"]
        header.extend(vcf.samples)
        file.write("\t".join(header)+"\n")
        # deal with multiple ALT alleles
        for variant in vcf:
            for alt_idx in range(0, len(variant.ALT)):
                outline = []
                outline.extend([variant.CHROM, str(variant.POS), variant.REF])
                outline.append(variant.ALT[alt_idx])  
                for sample_idx in range(0, len(variant.genotypes)):
                    vaf = str(variant.format('VAF')[sample_idx][alt_idx])
                    if variant.format('DP')[sample_idx] < 10:
                        vaf = "NA"
                    if vaf == "0.0":
                        vaf = "0"
                    outline.append(vaf)
                file.write("\t".join(outline)+"\n")
        
    vcf.close()

if __name__ == "__main__":
    main()