module run_analysis;

import calculation : correlation, Opts, transform, VarianceException;
import read_data : Genotype, Phenotype, readGenotype, readWeights;
import std.algorithm : map, sum;
import std.array : array;
import std.conv : to;
import std.math : fabs;
import std.numeric : dotProduct;
import std.range : chunks, enumerate, indexed, zip;
import std.stdio : File, stderr, writefln, writeln;

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

double corCalc(ref double[] genChunk, ref double[] phenChunk)
{
  try
  {
    transform(genChunk);
    return fabs(dotProduct(genChunk, phenChunk));
  }
  catch (VarianceException e)
  {
    return 0;
  }
}

void caveman(Phenotype phenotype, const size_t[] perms, File outFile, const Opts opts)
{
  auto genotype = readGenotype(opts, phenotype.chromosome, phenotype.location, phenotype.geneName);

  if (opts.verbose)
  {
    stderr.writeln("Extracted ", genotype.length, " usable genotypes.");
  }

  auto snpWeights = getWeights(opts, perms, phenotype, genotype);

  writeResults(snpWeights, phenotype, genotype, outFile);
}

auto getWeights(const Opts opts, const size_t[] permIndices,
    Phenotype phenotype, Genotype[] genotypes)
{
  immutable size_t nInd = phenotype.values.length;
  immutable nPerm = opts.perms[0];

  auto perms = phenotype.values.indexed(permIndices).array;

  foreach (ref perm; chunks(perms, nInd))
  {
    try
    {
      transform(perm);
    }
    catch (VarianceException)
    {
      perm[] = 0;
    }
  }

  CorPlace[][] trackCor = new CorPlace[][](nPerm, 10);

  foreach (countLine, ref genotype; enumerate(genotypes))
  {
    auto permGenotype = genotype.values.indexed(permIndices).array;
    auto simplePerm = zip(chunks(permGenotype, nInd), chunks(perms, nInd)).map!(
        e => corCalc(e[0], e[1]));
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

  immutable auto total = snpWeights.sum;

  snpWeights = snpWeights.map!(a => a / total).array;

  return snpWeights;
}

void writeResults(double[] snpWeights, Phenotype phenotype,
    Genotype[] genotypes, File outFile)
{
  transform(phenotype.values);

  foreach (i, ref genotype; enumerate(genotypes))
  {
    transform(genotype.values);
    auto cor = correlation(genotype.values, phenotype.values);

    outFile.writefln("%s%s%s\t%s\t%s", phenotype.geneName, genotype.snpId,
        cor[0], cor[1], snpWeights[i]);

  }
}
