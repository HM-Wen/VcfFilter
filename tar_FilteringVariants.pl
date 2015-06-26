#/usr/bin/perl -w
#
# tar_filteringVariants.pl
# ========================
#
# Author: Rafael Torres (ICVV) 
#
# Description: It takes an extended VCF file and filter SNP and Indels according to thresholds of quality score, 
# read depth and frequency of base calls.
#
# Usage: 
#	      perl tar_filteringVariants.pl param1 param2
#		   param1: vcf extended file with the vcf data for the gene(s) of reference in all the varieties.
#                  param2: output folder (it will host all the generated files)
#
# Reference: "Polymorphisms and minihaplotypes in the VvNAC26 gene associate with berry size
#             variation in grapevine" (Tello J et al., BMC Plant Biology 2015)
#

use strict;

# THRESHOLDS:
use constant {
	SNP_HOMOZ_MIN_RATIO => 0.942,
	SNP_HOMOZ_MIN_DEPTH => 10,
	SNP_HOMOZ_MIN_QUAL => 76,
	SNP_HETER_MIN_RATIO => 0.38,
	SNP_HETER_MIN_DEPTH => 20,
	SNP_HETER_MIN_QUAL => 57,
	INDEL_HOMOZ_MIN_RATIO => 0.83,
	INDEL_HOMOZ_MIN_DEPTH => 14,
	INDEL_HOMOZ_MIN_QUAL => 43,
	INDEL_HETER_MIN_RATIO => 0.28,
	INDEL_HETER_MIN_DEPTH => 38,
	INDEL_HETER_MIN_QUAL => 33,
	RESCUE_THRESHOLD_RATIO => .7,
	RESCUE_THRESHOLD_SAMPLES => 10
};
#

###### time #########
my $time = time();
 
my ($sec, $min, $hours, $mday, $month, $year, $wday, $yday, $isdst) = localtime($time);
 
my @day = qw (Sunday Monday Tuesday Wednesday Thursday Friday Saturday);
 
my @month_name = qw (January February March April May June July August September October November December);
 
#print STATS "init: ", $day[$wday], " ",$month_name[$month], " ",$mday," ", $hours,":",$min,":",$sec, " ",1900+$year,"\n";
print "init time script $0: ", $day[$wday], " ",$month_name[$month], " ",$mday," ", $hours,":",$min,":",$sec, " ",1900+$year,"\n";
##################

# vcf extended file of all varieties for the gene(s) of reference
my $VCFPLUS = $ARGV[0];
print "VCFPLUS (arg 1): $VCFPLUS\n";
open (VCFPLUS, $VCFPLUS) || die "Can't open $VCFPLUS";

# output folder
my $tasselFormatGenerateDir=$ARGV[1];

# output file
my $st=$tasselFormatGenerateDir.'/stats_filtering.txt';
open (STATS,">$st") or die "Cannot open $st\n";

# output file
my $flt=$tasselFormatGenerateDir.'/filtered.txt';
open (FILTERED,">$flt") or die "Cannot open $flt\n";

# output file
my $svd=$tasselFormatGenerateDir.'/saved.txt';
open (SAVED,">$svd") or die "Cannot open $svd\n";

# output file
my $fs=$tasselFormatGenerateDir.'/f&s.txt';
open (FILTANDS,">$fs") or die "Cannot open $fs\n";


my %hs_gene_sample=();
my %hs_gene_pos=();
my %hs_gene_type_check=();

my %hs_sample_genotype=();
my %hs_discarded_genotype=();
my %hs_discarded_genotype_count=();
my %hs_discarded_genotype_temp=();
my %hs_evaluated_genotype=();
my %hs_anysample_genotype=();

my %hs_auxiliar=();


while(<VCFPLUS>)
 {
# Format:
#
#    GENE
#	 SNPID
#	 SAMPLE
#    TYPE
#	  CHROM
#     POS
#     ID
#     REF
#     ALT
#     QUAL
#     FILTER
#     INFO (DP: combined depth across samples, e.g. DP=154)
#     FORMAT
#     VALUES
#    STRAND (+ or -)
#    DEPTH (coverage)
#    REFCNT (num. counts equal to reference allele)
#    ALTCNT (num. counts equal to alternative allele)
#    RATIO % (ratio alternative/reference * 100)
#    HOMOHETERO (homozygous or heterozygous mutation)
#
# Example of record:	
#VIT_00s0179g00150	chrUn_7415097_G_T	01VALD41-04-C4	SNP	chrUn	-4311	.	G	T	19.1	.	DP=28;VDB=0.0200;AF1=0.5;AC1=1;DP4=20,3,5,0;MQ=56;FQ=22;PV4=1,0.00087,3.6e-05,1	GT:PL:GQ	0/1:49,0,255:52	+	28	23	5	0.18	Ht


	my($line) = $_;
	chomp($line);
	my @campos = split /\t/,$line;
	my $geneid = $campos[0];
	my $snpid = $campos[1];
	
		my @snpid_info = split /_/,$snpid;
		my $posabs = $snpid_info[1];
	
	my $sample = $campos[2];
	my $type = $campos[3];
	my $chr = $campos[4];
	my $pos = $campos[5];
	my $id = $campos[6];
	my $ref = $campos[7];
	my $alt = $campos[8];
	my $qual = $campos[9];
	my $filter = $campos[10];
	my $info = $campos[11];
	my $format = $campos[12];
	my $fvalues = $campos[13];
	my $strand = $campos[14];
	my $depth = $campos[15];
	my $refcnt = $campos[16];
	my $altcnt = $campos[17];
	my $ratio = $campos[18];
	my $homohetero = $campos[19];
	
	$hs_anysample_genotype{$geneid}{$pos}{$type}="$ref" if ($type eq 'SNP'); # default;
	$hs_evaluated_genotype{$geneid}{$pos}{$type}{$sample}="$chr\t$posabs";

	if ( $type eq 'SNP') {
		if ( $homohetero eq 'Hm' ) {
			if ($ratio < SNP_HOMOZ_MIN_RATIO ) {
				$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}++;
				$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}="$ratio\t$depth\t$qual\t1\t?:?";
				$hs_auxiliar{$geneid}{$pos}{$type}{$sample}=1;
				print STATS $line,"\t1\n";
			} elsif ($depth < SNP_HOMOZ_MIN_DEPTH) {
				$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}++;
				$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}="$ratio\t$depth\t$qual\t2\t?:?";
				$hs_auxiliar{$geneid}{$pos}{$type}{$sample}=2;
				print STATS $line,"\t2\n";
				next;
			} elsif ($qual < SNP_HOMOZ_MIN_QUAL) {
				$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}++;
				$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}="$ratio\t$depth\t$qual\t3\t?:?";
				$hs_auxiliar{$geneid}{$pos}{$type}{$sample}=3;
				print STATS $line,"\t3\n";
				next;
			}
		} elsif ( $homohetero eq 'Ht' ) {
			if ($ratio < SNP_HETER_MIN_RATIO ) {
				$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}++;
				$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}="$ratio\t$depth\t$qual\t4\t?:?";
				$hs_auxiliar{$geneid}{$pos}{$type}{$sample}=4;
				print STATS $line,"\t4\n";
				next;
			} elsif ($depth < SNP_HETER_MIN_DEPTH) {
				$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}++;
				$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}="$ratio\t$depth\t$qual\t5\t?:?";
				$hs_auxiliar{$geneid}{$pos}{$type}{$sample}=5;
				print STATS $line,"\t5\n";
				next;
			} elsif ($qual < SNP_HETER_MIN_QUAL) {
				$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}++;
				$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}="$ratio\t$depth\t$qual\t6\t?:?";
				$hs_auxiliar{$geneid}{$pos}{$type}{$sample}=6;
				print STATS $line,"\t6\n";
				next;
			}
		} 
	}
	elsif ( $type eq 'INDEL' ) {
		if ( $homohetero eq 'Hm' ) {
			if ($ratio < INDEL_HOMOZ_MIN_RATIO ) {
				$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}++;
				$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}="$ratio\t$depth\t$qual\t7\t?:?";
				$hs_auxiliar{$geneid}{$pos}{$type}{$sample}=7;
				print STATS $line,"\t7\n";
				next;
			} elsif ($depth < INDEL_HOMOZ_MIN_DEPTH) {
				$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}++;
				$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}="$ratio\t$depth\t$qual\t8\t?:?";
				$hs_auxiliar{$geneid}{$pos}{$type}{$sample}=8;
				print STATS $line,"\t8\n";
				next;
			} elsif ($qual < INDEL_HOMOZ_MIN_QUAL) {
				$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}++;
				$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}="$ratio\t$depth\t$qual\t9\t?:?";
				$hs_auxiliar{$geneid}{$pos}{$type}{$sample}=9;
				print STATS $line,"\t9\n";
				next;
			}
		} elsif ( $homohetero eq 'Ht' ) {
			if ($ratio < INDEL_HETER_MIN_RATIO ) {
				$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}++;
				$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}="$ratio\t$depth\t$qual\t10\t?:?";
				$hs_auxiliar{$geneid}{$pos}{$type}{$sample}=10;
				print STATS $line,"\t10\n";
				next;
			} elsif ($depth < INDEL_HETER_MIN_DEPTH) {
				$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}++;
				$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}="$ratio\t$depth\t$qual\t11\t?:?";
				$hs_auxiliar{$geneid}{$pos}{$type}{$sample}=11;
				print STATS $line,"\t11\n";
				next;
			} elsif ($qual < INDEL_HETER_MIN_QUAL) {
				$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}++;
				$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}="$ratio\t$depth\t$qual\t12\t?:?";
				$hs_auxiliar{$geneid}{$pos}{$type}{$sample}=12;
				print STATS $line,"\t12\n";
				next;
			}
		} 
	}	
	
	$hs_discarded_genotype_count{$geneid}{$pos}{$type}{'OK'}++;

 }

close (VCFPLUS);

#
# "rescue" of polymorphisms (second chance for polymorphisms to be selected)
#
my $saved = 0;
my $removed = 0;

foreach my $geneid (sort keys %hs_discarded_genotype_count) {
	foreach my $pos (sort keys %{$hs_discarded_genotype_count{$geneid}}) {
		foreach my $type (sort keys %{$hs_discarded_genotype_count{$geneid}{$pos}}) {
			my $ok = 0;
			my $ko = 0;
			if (exists $hs_discarded_genotype_count{$geneid}{$pos}{$type}{'OK'}) {
				$ok = $hs_discarded_genotype_count{$geneid}{$pos}{$type}{'OK'};
			}
			if (exists $hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'}) {
				$ko = $hs_discarded_genotype_count{$geneid}{$pos}{$type}{'KO'};
			}
			
			next if $ko == 0;
			
			my $ratio = $ok / $ko;
			
			if ($ratio < RESCUE_THRESHOLD_RATIO && $ok < RESCUE_THRESHOLD_SAMPLES) {
				foreach my $sample (sort keys %{$hs_discarded_genotype_temp{$geneid}{$pos}{$type}}) {
					print FILTERED "$geneid\t$pos\t$type\t$sample\t$ratio\t$ok\t$ko\t$hs_discarded_genotype_temp{$geneid}{$pos}{$type}{$sample}\t$hs_evaluated_genotype{$geneid}{$pos}{$type}{$sample}\n";
					print FILTANDS "$geneid\t$pos\t$type\t$sample\t$ratio\t$ok\t$ko\tfiltered\t$hs_auxiliar{$geneid}{$pos}{$type}{$sample}\t$hs_evaluated_genotype{$geneid}{$pos}{$type}{$sample}\n";
					$removed++; 
				}
			} else {
				foreach my $sample (sort keys %{$hs_discarded_genotype_temp{$geneid}{$pos}{$type}}) {
					print SAVED "$geneid\t$pos\t$type\t$sample\t$ratio\t$ok\t$ko\t-:-\t$hs_evaluated_genotype{$geneid}{$pos}{$type}{$sample}\n"; 
					print FILTANDS "$geneid\t$pos\t$type\t$sample\t$ratio\t$ok\t$ko\tsaved\t$hs_evaluated_genotype{$geneid}{$pos}{$type}{$sample}\n";
					$saved++;
				}
			}
			
		}
	}
}

print "$removed removed\n";
print "$saved saved\n";

print "End of script $0\n";

close FILTERED;
close SAVED;
close FILTANDS;
close STATS;

exit;
