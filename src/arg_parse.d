module arg_parse;

import std.array : split;
import core.stdc.stdlib : exit;
import std.array : array;
import std.algorithm : canFind, countUntil, filter, joiner, map;
import std.conv : to, ConvException;
import std.exception : enforce;
import std.process : pipeShell, Redirect, wait;
import std.range : indexed, iota;
import std.stdio : File, writeln, stderr;
import std.string : chomp;

class Opts
{
  import std.getopt;

  //write appropriate string and quit
  bool version_ = false;
  bool spear = false;
  //phenotype and genotype ids are given
  bool noheader = false;
  //number of genotype columns to skip, and phenotype column
  int genes = 1;
  int jobNumber = 0;
  //permutation numbers and seeds
  bool give_seed = false;
  uint[] perms;
  //file names
  string vcf = "";
  string bed = "";
  string output = "";
  string correct = "";
  string cov = "";

  string[] interval;
  string results;

  bool normal = false;
  string best = "";
  size_t window = 1_000_000;
  bool gt = false;
  long loc = 0;
  bool nocheck = false;
  size_t[] genotypeLocations;
  size_t[] phenotypeLocations;
  size_t[] covLocations;

  string weights = "";

  this(string[] args)
  {
    bool noArgs = args.length == 1;
    try
    {
      // dfmt off
      arraySep = ",";
      auto options = getopt(args,
			    "bed", "Phenotype file [last argument].\n", &bed,
			    "vcf", "Genotype file.\n", &vcf,
			    "out|o", "Output file [stdout].\n", &output,
			    "weights", "Specify CaVEMaN weightings.\n", &weights,
			    "job-number", "Split the analysis into a number of smaller runs which can run in parallel on a cluster. This option specifies which of the sub-analyses should be run.\n", &jobNumber,
			    "genes", "This specifies the number of genes to be analysed in each job.\n", &genes,
			    "perm", "Number of bootstrap samplings, with an optional seed. One following number indicates the number of bootstraps, two comma separated numbers gives the number of bootstraps and the seed.\n", &perms,
			    "window", "The size in base pairs of the cis window around the transcription start site of the gene [1,000,000].\n", &window,
			    "spear", "Runs analysis based on non-parametric Spearman correlation test of association.\n", &spear,
			    "correct", "Specify eQTL file to output single genetic signal bed file.\n", &correct,
			    "cov", "Optional covariates matrix if correcting phenotypes.\n", &cov,
			    "normal", "Map single signal bed file phenotypes onto a normal distribution.\n", &normal,
			    "best", "Produce probabilities for most significant association from results file.\n", &best,
			    "interval", "Given a results file, produce the smallest set of SNPs with a given probability of containing the causal variant for each gene. Takes either 1 or 2 comma separated arguments, the name of the results file and optionally the required probability [0.9].\n", &interval,
			    "nocheck", "Do not attempt to match genotype and phenotype IDs.\n", &nocheck,
			    "noheader", "Suppress writing of header line.\n", &noheader,
			    "version", "Display version information.\n", &version_,
			      );
      // dfmt on
      if (options.helpWanted || noArgs)
      {
        defaultGetoptPrinter("NAME
       CaVEMaN - Causal Variant Evidence Mapping using Non-parametric resampling.

SYNOPSIS
       CaVEMaN [options]

DESCRIPTION
       CaVEMaN estimates the probability that a particular variant is the causal variant underlying a genetic association by observing the properties of the association under repeated resampling. The standard mode will run the CaVEMaN analysis on data in standard vcf and bed format, giving the CaVEMaN score for each variant in a cis window of the gene. When there are multiple eQTL for certain genes, the --correct flag will take the list of eQTL and the bed file and create a \"single signal\" expression bed file. This has a phenotype for each eQTL, regressing out all other eQTL, and can now be analysed using the standard pipeline. Once a results file has been created, a summary results file with the most significant association per gene and the probability that the specific variant is causal can be produced with the --best flag. Finally, the --interval flag will reduce a results file to a set of SNPs with probability x that they contain the causal variant.

OPTIONS
", options.options);
        writeln("
FILE FORMATS

   INPUT FILE FORMATS
       Phenotype
              BED file. This file should only contain genes for which there is an eQTL and should not be compressed.

       Genotype
              VCF format. This file should be compressed with bgzip and tabixed.

       eQTL   With the --correct flag a list of eQTL are required, a tab separated file with the following fields: Gene, Chromosome, Location (bp), Reference allele and Alternate allele (all fields necessary to uniquely specify the eQTL).

       Results.
              The --best and --interval flags require only a results file, in the format described in the next section.

   OUTPUT FILE FORMAT
       Standard run and the --interval flag produce an output file with the following fields: Chromosome, Location (bp), Reference and Alternate alleles of the SNP, followed by the Correlation and P value for association and the CaVEMaN weighting. With the --best flag an extra column is added with the probability of being the causal variant. Running with the --correct flag produces a bed file for further analysis.
");
        exit(0);
      }

      if (version_)
        giveHelp(versionString);

      if (best == "" && interval.length == 0)
      {
        if (bed == "" && args.length > 1)
          bed = args[$ - 1];

        if (perms.length == 0)
          perms = [10_000];

        matchIds();
      }
    }
    catch (Exception e)
    {
      writeln("Error with command: ", e.msg);
      exit(0);
    }
  }

  private void matchIds()
  {
    string[] phenotypeIds;
    try
    {
      auto bedFile = File(bed);

      phenotypeIds = bedFile.readln.chomp.split[4 .. $];
    }
    catch (Exception e)
    {
      stderr.writeln("Failed to read phenotype IDs. ", e.msg);
      exit(0);
    }
    if (nocheck)
    {
      genotypeLocations = iota(phenotypeIds.length).array;
      phenotypeLocations = iota(phenotypeIds.length).array;
    }
    else
    {
      string[] genotypeIds;
      try
      {
        auto pipes = pipeShell("zcat " ~ vcf ~ " | grep -v '##' | head -2", Redirect.stdout);
        scope (exit)
          wait(pipes.pid);

        auto line = pipes.stdout.readln.chomp;

        genotypeIds = line.split[9 .. $].to!(string[]);

        auto formatField = pipes.stdout.readln.chomp.split[8].split(':');
        loc = countUntil(formatField, "DS");
        if (loc == -1)
        {
          loc = countUntil(formatField, "GT");
          gt = true;
          if (loc == -1)
          {
            stderr.writeln("DS and GT fields are both missing from the vcf file.");
            exit(0);
          }
        }
      }
      catch (Exception e)
      {
        stderr.writeln("Failed to read genotype IDs. ", e.msg);
        exit(0);
      }

      phenotypeLocations = genotypeIds.map!(a => phenotypeIds.countUntil(a))
        .filter!(a => a != -1).array.to!(size_t[]);

      genotypeLocations = iota(genotypeIds.length).filter!(
          a => phenotypeIds.canFind(genotypeIds[a])).array;

      if (phenotypeIds.indexed(phenotypeLocations)
          .array != genotypeIds.indexed(genotypeLocations).array)
      {
        stderr.writeln("Failed to match IDs. THIS SHOULD NEVER HAPPEN.");
        exit(0);
      }

      if (genotypeLocations.length == 0 || phenotypeLocations.length == 0)
      {
        stderr.writeln("No individuals to analyse.");
        exit(0);
      }

      if (correct != "" && cov != "")
      {

        string[] covIds;

        try
        {
          covIds = File(cov).readln.chomp.split;

          auto temp = genotypeIds.indexed(genotypeLocations).map!(a => covIds.countUntil(a)).array;
          if (temp.canFind(-1))
          {
            auto missing = genotypeIds.indexed(genotypeLocations)
              .filter!(a => !covIds.canFind(a)).joiner(", ");
            stderr.writeln(
                "Individuals in genotype and phenotype file, but not in covariate file. Missing individuals are ",
                missing, ".");
            exit(0);
          }

          covLocations = temp.to!(size_t[]);
        }
        catch (Exception e)
        {
          stderr.writeln("Failed to read covariate IDs. ", e.msg);
          exit(0);
        }
      }
    }
  }

}

static immutable string versionString = "CaVEMaN, Causal Variant Evidence Mapping using Non-parametric Sampling, version 0.9.0";

void giveHelp(immutable string quitString)
{
  writeln(quitString);
  exit(0);
}
