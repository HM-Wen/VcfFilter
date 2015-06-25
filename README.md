# VcfFilter

Filtering of polymorphisms according to thresholds of quality score, read depth and frequency of base calls.
------------------------------------------------------------------------------------------------------------

# Reference: "Polymorphisms and minihaplotypes in the VvNAC26 gene associate with berry size
#             variation in grapevine" (Tello J et al., BMC Plant Biology 2015)

The script tar_FilteringVariants.pl takes and extended VCF file and filter SNP and Indels according to thresholds of quality score, 
read depth and frequency of base calls. These are the parameters and thresholds for the filtering:

	1.  SNP_HOMOZ_MIN_RATIO (minimum RATIO for polimorphism of type "homozygous SNP") = 0.942,
	2.  SNP_HOMOZ_MIN_DEPTH (minimum DEPTH for polimorphism of type "homozygous SNP") = 10,
	3.  SNP_HOMOZ_MIN_QUAL  (minimum QUAL for polimorphism of type "homozygous SNP") = 76,
	4.  SNP_HETER_MIN_RATIO (minimum RATIO for polimorphism of type "heterozygous SNP") = 0.38,
	5.  SNP_HETER_MIN_DEPTH (minimum DEPTH for polimorphism of type "heterozygous SNP") = 20,
	6.  SNP_HETER_MIN_QUAL (minimum QUAL for polimorphism of type "heterozygous SNP") =  57,
	7.  INDEL_HOMOZ_MIN_RATIO (minimum RATIO for polimorphism of type "homozygous Indel") = 0.83,
	8.  INDEL_HOMOZ_MIN_DEPTH (minimum DEPTH for polimorphism of type "homozygous Indel") = 14,
	9.  INDEL_HOMOZ_MIN_QUAL (minimum QUAL for polimorphism of type "homozygous Indel") = 43,
	10. INDEL_HETER_MIN_RATIO (minimum RATIO for polimorphism of type "homozygous Indel") = 0.28,
	11. INDEL_HETER_MIN_DEPTH (minimum DEPTH for polimorphism of type "homozygous Indel") = 38,
	12. INDEL_HETER_MIN_QUAL (minimum QUAL for polimorphism of type "homozygous Indel") = 33

with RATIO = ratio "calls_alternative_allele/calls_reference_allele", 
     DEPTH = coverage (num of reads) for the position of the polimorphism
     QUAL = quality of the polimorphism calling (determined by the variant caller, samtools)

- The thresholds 1 to 12 are the mean less one standard deviation value for the distribution of the corresponding parameter in the entire set of samples.
For example, taking into account all the depth values in the set of heterozygous SNPs, the mean value is 80 and the standard deviation is 60, so the threshold for the parameter SNP_HETER_MIN_DEPTH is 20.
