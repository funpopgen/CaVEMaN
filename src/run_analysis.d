module run_analysis;

import std.algorithm : map, sum;
import std.array : array;
import std.conv : ConvException, to;
import std.math : fabs;
import std.numeric : dotProduct;
import std.range : chunks, enumerate, indexed, zip;
import std.stdio : File, writeln;

import calculation : correlation, rank, Opts, transform, VarianceException;
import read_data : Phenotype, Genotype, readGenotype, readWeights;

version (unittest)
{
  import std.digest.sha;
  import std.range : put;
  import std.file : exists, remove;
}

class InputException : Exception
{
  pure this(string s)
  {
    super(s);
  }
}

struct CorPlace
{
  double cor = 0;
  size_t[] place;
}

double corCalc(ref double[] genChunk, ref double[] phenChunk, bool ttest)
{
  try
  {
    if (ttest)
      transform(genChunk);
    else
      transform(rank(genChunk));
    return fabs(dotProduct(genChunk, phenChunk));
  }
  catch (VarianceException e)
  {
    return 0;
  }
}

void CaVEMaN(Phenotype phenotype, size_t[] perms, File outFile, Opts opts)
{
  auto genotype = readGenotype(opts, phenotype.chromosome, phenotype.location,
      phenotype.values.length, phenotype.geneName);

  auto snpWeights = getWeights(opts, perms, phenotype, genotype);

  writeResults(opts, snpWeights, phenotype, genotype, outFile);
}

auto getWeights(Opts opts, size_t[] permIndices, Phenotype phenotype, Genotype[] genotypes)
{
  immutable size_t nInd = phenotype.values.length;
  immutable nPerm = opts.perms[0];

  auto perms = phenotype.values.indexed(permIndices).array;

  foreach (ref perm; chunks(perms, nInd))
  {
    try
    {
      if (!opts.spear)
        transform(perm);
      else
        transform(rank(perm));
    }
    catch (VarianceException)
    {
      perm[] = 0;
    }
  }

  CorPlace[][] trackCor = new CorPlace[][](nPerm, 10);

  size_t countLine = 0;

  foreach (ref genotype; genotypes)
  {
    auto permGenotype = genotype.values.indexed(permIndices).array;
    auto simplePerm = zip(chunks(permGenotype, nInd), chunks(perms, nInd)).map!(
        e => corCalc(e[0], e[1], !opts.spear));
    foreach (i, e; simplePerm.enumerate)
    {
      foreach (j; 0 .. 10)
      {
        if (trackCor[i][j].cor < e)
        {
          trackCor[i][(j + 1) .. 10] = trackCor[i][j .. 9].dup;
          trackCor[i][j] = CorPlace(e, [countLine]);
          break;
        }
        else if (trackCor[i][j].cor == e)
        {
          trackCor[i][j].place ~= countLine;
          break;
        }
      }
    }
    countLine++;
  }

  auto snpWeights = new double[](genotypes.length);

  snpWeights[] = 0;

  double[] weights = readWeights(opts);

  foreach (ref e; trackCor)
  {
    foreach (j, ref f; enumerate(e))
    {
      if (f.cor != 0)
      {
        foreach (g; f.place)
        {
          snpWeights[g] += weights[j];
        }
      }
    }
  }

  auto total = snpWeights.sum;

  snpWeights = snpWeights.map!(a => a / total).array;

  return snpWeights;
}

void writeResults(Opts opts, double[] snpWeights, Phenotype phenotype,
    Genotype[] genotypes, File outFile)
{
  import std.array : join;

  string cor;
  immutable size_t nInd = phenotype.values.length;
  size_t countLine = 0;

  if (!opts.spear)
    transform(phenotype.values);
  else
    transform(rank(phenotype.values));

  foreach (ref genotype; genotypes)
  {
    if (!opts.spear)
      transform(genotype.values);
    else
      transform(rank(genotype.values));

    cor = correlation(genotype.values, phenotype.values)[].to!(string[]).join("\t");

    outFile.writeln(phenotype.geneName, genotype.snpId, cor, "\t",
        snpWeights[countLine].to!string);

    countLine++;
  }
}
