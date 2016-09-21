module calculation;

import std.exception : enforce;
import std.math : fabs, sqrt;
import arg_parse : Opts;

enum double EPSILON = 0.00000001; //comparison for X>=Y is done X > Y - epsilon 

version (unittest)
{
  import std.math : approxEqual;
}

class VarianceException : Exception
{
  //thrown if variable is constant
  pure nothrow this(string s)
  {
    super(s);
  }
}

pure nothrow extern (C)
{
  //call GSL to calculate P values from T statistics
  double gsl_cdf_tdist_P(double x, double nu);
}

unittest
{
  // Checks GSL gives right p value for t statistic
  assert(approxEqual(gsl_cdf_tdist_P(-1.6, 7), 0.07681585));
}

pure ref double[] rank(ref double[] rankArray)
{
  //ranks array, giving ties mean rank
  import std.algorithm : makeIndex;

  immutable size_t len = rankArray.length;
  auto orderIndex = new size_t[](len);
  makeIndex!("a < b")(rankArray, orderIndex);

  double sumrank = 0.0;
  size_t dupcount = 0;
  double avgrank;

  foreach (i, ref e; orderIndex)
  {
    sumrank += i;
    dupcount++;
    if (i == (len - 1) || rankArray[e] != rankArray[orderIndex[i + 1]])
    {
      avgrank = sumrank / dupcount + 1;
      foreach (ref j; orderIndex[(i - dupcount + 1) .. (i + 1)])
        rankArray[j] = avgrank;
      sumrank = 0;
      dupcount = 0;
    }
  }
  return rankArray;
}

unittest
{
  //Simple test of ranking with ties
  double[] vector = [10, 9, 2, 9, 3];

  assert(rank(vector) == [5, 3.5, 1, 3.5, 2]);
}

pure void transform(ref double[] vector)
{
  //transforms array so mean =0 sum of squares = 1
  int n = 0;
  double mean = 0;
  double M2 = 0;
  double delta;

  foreach (ref e; vector)
  {
    n++;
    delta = e - mean;
    mean += delta / n;
    M2 += delta * (e - mean);
  }

  enforce(M2 != 0, new VarianceException(""));

  M2 = sqrt(M2);

  foreach (ref e; vector)
    e = (e - mean) / M2;
}

unittest
{
  //Checks that transform works on randomly generated vector
  import std.algorithm : reduce;
  import std.random : uniform;

  double[] x = new double[](10);
  foreach (ref e; x)
    e = uniform(0.0, 10.0);

  transform(x);
  auto mean = 0.0.reduce!((a, b) => a + b)(x);

  assert(approxEqual(mean, 0.0));
  assert(approxEqual(0.0.reduce!((a, b) => a + (b - mean) * (b - mean))(x), 1));
}

pure nothrow double[2] correlation(ref double[] vector1, ref double[] vector2)
{
  //calculates correlation, t stat and p value for two arrays
  import std.numeric : dotProduct;

  double[2] results;
  results[0] = dotProduct(vector1, vector2);
  results[1] = gsl_cdf_tdist_P(
      -fabs(results[0] * sqrt((vector1.length - 2) / (1 - results[0] * results[0]))),
      vector1.length - 2) * 2;
  return results;
}

unittest
{
  //Check correlation of phenotype with 3rd row genotype against estimates from R
  double[2] corFromR = [-0.2863051, 0.4225695];

  double[] genotype = [0.115, 2, 0.0964, 1, 1, 1, 0, 1, 0, 0.0563];
  double[] phen = [
    -1.3853088072, -0.785797093643, 1.14540423638, -0.785797093643, 1.03820492508,
    -1.25652676836, -0.787662180447, -2.05355237841, -0.245457234103, 1.14277217712
  ];

  transform(rank(phen));
  transform(rank(genotype));
  double[2] cor = correlation(genotype, phen);

  assert(approxEqual(cor[0], corFromR[0]));
  assert(approxEqual(cor[1], corFromR[1]));
}

size_t[] genPerms(Opts opts, size_t nInd)
{
  import std.array : array;
  import std.algorithm : map;
  import std.range : iota;
  import std.random : rndGen, uniform;

  if (opts.perms.length > 1)
    rndGen.seed(opts.perms[1]);

  return iota(opts.perms[0] * nInd).map!(a => uniform(0, nInd)).array;

}
