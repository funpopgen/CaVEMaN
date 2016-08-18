#CaVEMaN : Causal Variant Evidence Mapping using Non-parametric resampling

##Introduction:

TODO
##NAME

CaVEMaN  - Mapping causal variants using bootstrap resampling.
##SYNOPSIS

CaVEMaN [options]
##DESCRIPTION

TO DO
##OPTIONS

	   --help    Display help information.

	   --version Display version information.

	   --bed CHAR
                 Phenotype file [last argument].

	   --vcf CHAR
                 Genotype file.

	   --out, --o CHAR
                 Output file [stdout].

	   --job-number INT
                 Split the analysis into a number of smaller runs which can run in parallel
                 on  a  cluster.  This option specifies which of the sub-analyses should be
                 run.

	   --genes INT
                 This specifies the number of genes to be analysed in each job.

	   --window INT
                 The size in base pairs of the cis window around  the  transcription  start
                 site of the gene [1,000,000].

	   --perm INT(,INT)
                 Calculated permuted p values, one following number indicates the number of
                 permutations, two comma separated numbers gives the number of permutations
                 and the seed.

	   --spear   Runs analysis based on non-parametric Spearman correlation test of associ‐
                 ation.

	   --bcftools
                 Use the bcftools plugin to calculate dosage from the GT field.

	   --noheader
                 Suppress writing of header line.

##FILE FORMATS

###INPUT FILE FORMATS

       Phenotype
              BED file.

	   Genotype
              VCF format. Format should be one DS field unless bcftools flag is  given.  In
              this  case  bcftools  plugin  dosage  will be used to calulate dosage from GT
              field.

###OUTPUT FILE FORMAT

	   Output contains the chromosome, location, reference and  alternate  alleles  of  the
       SNP, followed by the correlation, P value for association and the CaVEMaN weighting.
