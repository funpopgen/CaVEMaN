module correct;

import arg_parse : Opts;
import core.stdc.stdlib : exit;
import read_data : getDosage, InputException;
import std.algorithm : makeIndex, map, max, reduce;
import std.array : array, split;
import std.conv : to;
import std.exception : enforce;
import std.format : format;
import std.math : approxEqual, sqrt;
import std.process : pipeShell, Redirect, wait;
import std.random : Random, uniform01, unpredictableSeed;
import std.range : enumerate, indexed, iota;
import std.stdio : File, stderr, stdout, write, writef, writefln, writeln;
import std.string : chomp;

extern (C)
{
  void generateRandomValues(double* randomSample, size_t nInd, double sigma, size_t seed);

  struct gsl_block_struct
  {
    size_t size;
    double* data;
  }

  alias gsl_block = gsl_block_struct;

  struct gsl_matrix
  {
    size_t size1;
    size_t size2;
    size_t tda;
    double* data;
    gsl_block* block;
    int owner;
  }

  struct gsl_vector
  {
    size_t size;
    size_t stride;
    double* data;
    gsl_block* block;
    int owner;
  }

  struct gsl_multifit_linear_workspace
  {
    size_t nmax; /* number of observations */
    size_t pmax; /* number of parameters */
    gsl_matrix* A;
    gsl_matrix* Q;
    gsl_matrix* QSI;
    gsl_vector* S;
    gsl_vector* t;
    gsl_vector* xt;
    gsl_vector* D;
  }

  gsl_vector* gsl_vector_alloc(size_t n);
  void gsl_vector_free(gsl_vector* v);

  double gsl_vector_get(const gsl_vector* v, const size_t i);
  void gsl_vector_set(gsl_vector* v, const size_t i, double x);

  gsl_matrix* gsl_matrix_alloc(size_t n1, size_t n2);
  void gsl_matrix_free(gsl_matrix* m);

  void gsl_matrix_set_all(gsl_matrix* m, double x);

  double gsl_matrix_get(const gsl_matrix* m, const size_t i, const size_t j);
  void gsl_matrix_set(gsl_matrix* m, const size_t i, const size_t j, double x);

  gsl_multifit_linear_workspace* gsl_multifit_linear_alloc(const size_t n, const size_t p);
  void gsl_multifit_linear_free(gsl_multifit_linear_workspace* work);

  int gsl_multifit_linear(const gsl_matrix* X, const gsl_vector* y, gsl_vector* c,
      gsl_matrix* cov, double* chisq, gsl_multifit_linear_workspace* work);

  int gsl_multifit_linear_residuals(const gsl_matrix* X, const gsl_vector* y,
      const gsl_vector* c, gsl_vector* r);

  int gsl_fit_linear(const double* x, const size_t xstride, const double* y,
      const size_t ystride, size_t n, double* c0, double* c1, double* cov00,
      double* cov01, double* cov11, double* sumsq);

  double gsl_cdf_gaussian_Pinv(double P, double sigma);

  double gsl_cdf_ugaussian_Pinv(double P);
  double gsl_stats_sd(const double* data, size_t stride, size_t n);
}

void correct(const Opts opts)
{
  auto eqtlList = getEqtl(opts.eqtl);
  if (opts.verbose)
  {
    stderr.writeln("Finished extracting ",
        eqtlList.byKey.map!(a => eqtlList[a].length).reduce!((a, b) => a + b), " eQTL.");
  }

  auto snps = getSnps(eqtlList, opts.vcf, opts.genotypeLocations, opts.loc, opts.gt);
  if (opts.verbose)
  {
    stderr.writeln("Finished extracting ", snps.length, " SNPs.");
  }

  auto cov = getCov(opts);
  if (opts.verbose)
  {
    stderr.writeln("Finished reading covariates.");
  }

  writeBed(opts, eqtlList, snps, cov);
}

string[][string] getEqtl(string inFile)
{
  string[][string] eqtlList;

  try
  {
    auto eqtlFile = File(inFile);
    foreach (line; eqtlFile.byLine)
    {
      auto lineSplit = line.split;
      enforce(lineSplit.length == 5, new InputException("Line doesn't have 5 columns."));
      auto gene = lineSplit[0].to!string;
      auto snp = format("%-(%s_%)", lineSplit[1 .. $]);
      if (auto p = gene in eqtlList)
      {
        eqtlList[gene] ~= snp;
      }
      else
      {
        eqtlList[gene] = [snp];
      }
    }
  }
  catch (Exception e)
  {
    stderr.writeln("Failed to parse eQTL file. ", e.msg);
    exit(1);
  }

  return eqtlList;
}

double[][string] getSnps(string[][string] eqtlList, string vcfFile,
    const size_t[] locations, const long loc, const bool gt)
{
  double[][string] snps;

  foreach (e; eqtlList.byKey)
  {
    foreach (f; eqtlList[e])
    {
      if (!(f in snps))
      {
        string[] snpBreakdown = f.split("_");
        string tabixCommand = "tabix " ~ vcfFile ~ " " ~ snpBreakdown[0] ~ ":"
          ~ snpBreakdown[1] ~ "-" ~ snpBreakdown[1] ~ " | grep -v '#' | cut -f1,2,4,5,10-";

        auto pipes = pipeShell(tabixCommand, Redirect.stdout);
        scope (exit)
          wait(pipes.pid);

        foreach (line; pipes.stdout.byLine)
        {
          auto splitLine = line.split;
          if (format("%-(%s_%)", splitLine[0 .. 4]) == f)
          {
            snps[f] = splitLine[4 .. $].indexed(locations).map!(a => getDosage(a, loc, gt)).array;
          }
        }

        if (!(f in snps))
        {
          stderr.writeln("SNP ", f, " not found in vcf file.");
          exit(1);
        }
      }
    }
  }
  return snps;
}

double[] getCov(const Opts opts)
{
  double[] cov;
  if (opts.cov == "")
  {
    return cov;
  }
  else
  {
    try
    {
      auto covFile = File(opts.cov);
      covFile.readln;
      foreach (line; covFile.byLine)
      {
        auto tempLine = line.split;
        if (tempLine.length < opts.covLocations.length)
        {
          stderr.writeln("Too few individuals in covariates file.");
          exit(1);
        }
        cov ~= tempLine.indexed(opts.covLocations).map!(a => to!double(a)).array;
      }
    }
    catch (Exception e)
    {
      stderr.writeln("Problem reading covariates file. ", e.msg);
      exit(1);
    }
    return cov;
  }
}

void writeBed(const Opts opts, string[][string] eqtlList, double[][string] snps, double[] cov)
{
  File outFile;
  File bedFile;
  try
  {
    if (opts.output == "")
    {
      outFile = stdout;
    }
    else
    {
      outFile = File(opts.output, "w");
    }

    bedFile = File(opts.bed);
  }
  catch (Exception e)
  {
    stderr.writeln("Trouble opening files. ", e.msg);
    exit(1);
  }

  Random rnd;
  if (opts.perms.length != 0)
  {
    rnd.seed(opts.perms[0]);
  }
  else
  {
    rnd.seed(unpredictableSeed);
  }

  immutable auto nInd = opts.phenotypeLocations.length;
  immutable auto baseCov = cov.length / nInd;
  auto maxEqtl = eqtlList.byKey.map!(a => eqtlList[a].length).reduce!(max) + baseCov;

  auto outcome = gsl_vector_alloc(nInd);
  auto workSpace = gsl_multifit_linear_alloc(nInd, maxEqtl);
  auto residuals = gsl_vector_alloc(nInd);

  double chisq;

  scope (exit)
  {
    gsl_vector_free(outcome);
    gsl_multifit_linear_free(workSpace);
    gsl_vector_free(residuals);
  }

  auto headerLineSplit = bedFile.readln.chomp.split;
  outFile.writefln("%-(%s\t%)\t%-(%s\t%)", headerLineSplit[0 .. 4],
      headerLineSplit[4 .. $].indexed(opts.phenotypeLocations));

  foreach (bedLine; bedFile.byLine)
  {
    auto bedLineSplit = bedLine.split.to!(string[]);
    if (bedLineSplit[3] in eqtlList)
    {
      if (eqtlList[bedLineSplit[3]].length == 1 && baseCov == 0 && !opts.simulate)
      {
        outFile.writef("%-(%s\t%)_%s\t", bedLineSplit[0 .. 4], eqtlList[bedLineSplit[3]][0]);
        if (opts.normal)
        {
          auto values = bedLineSplit[4 .. $].indexed(opts.phenotypeLocations)
            .map!(a => to!double(a)).array;
          normalise(values);
          outFile.writefln("%-(%g\t%)", values);
        }
        else
        {
          outFile.writefln("%-(%s\t%)", bedLineSplit[4 .. $].indexed(opts.phenotypeLocations));
        }
      }
      else
      {
        auto nCov = eqtlList[bedLineSplit[3]].length + baseCov;
        auto covariates = gsl_matrix_alloc(nInd, nCov);
        auto coefficients = gsl_vector_alloc(nCov);
        auto corrMat = gsl_matrix_alloc(nCov, nCov);

        scope (exit)
        {
          gsl_matrix_free(covariates);
          gsl_matrix_free(corrMat);
          gsl_vector_free(coefficients);
        }

        foreach (i; 0 .. nInd)
        {
          gsl_matrix_set(covariates, i, 0, 1);
          gsl_vector_set(outcome, i, bedLineSplit[opts.phenotypeLocations[i] + 4].to!double);
        }

        foreach (i, e; cov)
        {
          gsl_matrix_set(covariates, i % nInd, i / nInd + 1, e);
        }

        auto snpKeys = eqtlList[bedLineSplit[3]];

        foreach (e; snpKeys)
        {
          auto j = 1 + baseCov;
          foreach (i, f; enumerate(snpKeys))
          {
            if (f != e)
            {
              foreach (k; 0 .. nInd)
              {
                gsl_matrix_set(covariates, k, j, snps[f][k]);
              }
              j++;
            }
          }

          gsl_multifit_linear(covariates, outcome, coefficients, corrMat, &chisq, workSpace);
          gsl_multifit_linear_residuals(covariates, outcome, coefficients, residuals);

          auto values = gslToArray(residuals);

          if (!opts.simulate)
          {
            if (opts.normal)
            {
              normalise(values);
            }

            outFile.writef("%-(%s\t%)_%s\t", bedLineSplit[0 .. 4], e);
            outFile.writefln("%-(%g\t%)", values);
          }
          else
          {

            outFile.writef("%-(%s\t%)_%s\t", bedLineSplit[0 .. 4], e);

            double c0, effectSize, cov00, cov01, cov11, sigma;

            gsl_fit_linear(snps[e].ptr, 1, values.ptr, 1, values.length, &c0,
                &effectSize, &cov00, &cov01, &cov11, &sigma);

            sigma = sqrt(sigma / (nInd - 2));

            if (opts.verbose)
            {
              stderr.writeln("For gene ", bedLine.split[3], " the effect size estimate is ",
                  effectSize, "; sigma estimate is ", sigma, ".");
            }

            values = iota(nInd).map!(a => gsl_cdf_gaussian_Pinv(uniform01(rnd), sigma)).array;

            foreach (i; 0 .. values.length)
            {
              values[i] += effectSize * snps[e][i];
            }

            gsl_fit_linear(snps[e].ptr, 1, values.ptr, 1, values.length, &c0,
                &effectSize, &cov00, &cov01, &cov11, &sigma);

            if (opts.verbose)
            {
              stderr.writeln("After: For gene ", bedLine.split[3], " the effect size estimate is ",
                  effectSize, "; sigma estimate is ", sqrt(sigma / (nInd - 2)), ".");
            }

            if (opts.normal)
            {
              normalise(values);
            }
            outFile.writefln("%-(%g\t%)", values);

          }
        }
      }

    }
  }
}

void normalise(ref double[] residuals)
{
  size_t[] orderBuffer = new size_t[](residuals.length);

  makeIndex(residuals, orderBuffer);

  size_t count = 0;
  double previous = residuals[orderBuffer[count]];
  residuals[orderBuffer[0]] = count.to!double;

  foreach (e; orderBuffer[1 .. $])
  {
    if (residuals[e] > previous)
    {
      count++;
    }
    previous = residuals[e];
    residuals[e] = count.to!double;
  }

  double[] inverseNormal = iota(count + 1).map!(
      a => gsl_cdf_ugaussian_Pinv(1.0 * (a + 1) / (count + 2)).to!double).array;

  foreach (ref e; residuals)
    e = inverseNormal[e.to!size_t];
}

double[] gslToArray(gsl_vector* vec)
{
  double[] array = new double[](vec.size);
  foreach (i, e; array)
  {
    array[i] = gsl_vector_get(vec, i);
  }
  return array;
}

@system unittest
{

  double[] residuals = [1.0, 4, 4.5, 2, 1, 1, 1];
  double[] resultsFromR = [
    -0.8416212, 0.2533471, 0.8416212, -0.2533471, -0.8416212, -0.8416212, -0.8416212
  ];

  normalise(residuals);

  assert(approxEqual(residuals, resultsFromR));

  residuals = [-1.8770457188144, 0.211096107415374, -0.56072969755841,
    1.11564407213545, 1.99740513952311, -0.0887294570420509,
    -0.610132655912594, -0.378524183553565, -2.41530373819483,
    -0.186576831967272, -0.577841481720977, -0.760309403824646,
    0.500262066216114, -0.0672838126269592, -0.438653775896711,
    -0.119212500084674, 1.37251408213678, 2.24643382422438, -1.04865674605464, 0.649746098602992];

  resultsFromR = [-1.30917171678578, 0.430727299295458, -0.430727299295457,
    0.876142849246841, 1.30917171678578, 0.180012369792705,
    -0.712443032389489, -0.180012369792705, -1.66839119394708,
    -0.0597170997853229, -0.565948821932863, -0.876142849246841,
    0.565948821932863, 0.302980448056207, -0.302980448056207,
    0.0597170997853226, 1.06757052387814, 1.66839119394708, -1.06757052387814, 0.712443032389489];

  normalise(residuals);
  assert(approxEqual(residuals, resultsFromR));

  residuals = [6, 3, 2, 6, 9, 3, 1, 8, 3, 5, 3, 2, 4, 4, 1, 4, 8, 6, 2, 10];

  resultsFromR = [0.2533471031358, -0.524400512708041, -0.841621233572914,
    0.2533471031358, 0.841621233572914, -0.524400512708041,
    -1.2815515655446, 0.524400512708041, -0.524400512708041, 0,
    -0.524400512708041, -0.841621233572914, -0.2533471031358,
    -0.2533471031358, -1.2815515655446, -0.2533471031358, 0.524400512708041,
    0.2533471031358, -0.841621233572914, 1.2815515655446];

  normalise(residuals);
  assert(approxEqual(residuals, resultsFromR));
}
