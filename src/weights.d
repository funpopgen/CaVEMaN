module weights;

import arg_parse : Opts;
import core.stdc.stdlib : exit;
import std.algorithm : count, filter, map, sort, sum, uniq;
import std.array : split;
import std.conv : to;
import std.range : iota;
import std.stdio : File, readln, stderr, stdout, writeln;

struct Gene
{
  string[] causal;
  double pVal;
  double causalP;
  double caveman;
  size_t rank;

  this(string[] line)
  {
    causal = line[0].split("_")[($ - 4) .. $];
    caveman = line[7].to!double;
    pVal = line[6].to!double;
    if (causal == line[1 .. 5])
    {
      causalP = pVal;
    }
  }

  void update(double cave, double p)
  {
    caveman = cave;
    pVal = p;
  }
}

void getWeights(const Opts opts)
{
  File inFile;
  File weightFile;
  File rankFile;

  try
  {
    inFile = File(opts.results);
  }
  catch (Exception e)
  {
    stderr.writeln("Failed to read results file. ", e.msg);
    exit(1);
  }

  inFile.readln;

  Gene[] genes;
  double[] allPVals;
  double pVal;
  allPVals.assumeSafeAppend;
  string currentGene = "";

  foreach (line; inFile.byLine)
  {
    auto splitLine = line.split;
    auto gene = splitLine[0].to!string;
    pVal = splitLine[6].to!double;

    if (currentGene != gene)
    {
      if (genes.length > 0)
      {
        genes[$ - 1].rank = allPVals.sort!().uniq.count!(a => a < genes[$ - 1].causalP);
      }
      currentGene = gene;
      genes ~= Gene(splitLine.to!(string[]));
      allPVals.length = 0;
    }
    else if (genes[$ - 1].pVal > pVal)
    {
      genes[$ - 1].update(splitLine[7].to!double, pVal);
    }
    if (genes[$ - 1].causal == splitLine[1 .. 5])
    {
      genes[$ - 1].causalP = pVal;
    }
    allPVals ~= pVal;

  }

  genes[$ - 1].rank = allPVals.sort!().uniq.count!(a => a < genes[$ - 1].causalP);

  genes.sort!((a, b) => a.caveman < b.caveman);

  try
  {
    weightFile = File(opts.weights, "w");
  }
  catch (Exception e)
  {
    stderr.writeln("Failed to open file to write weights. ", e.msg);
    exit(1);
  }

  try
  {
    rankFile = File(opts.rank, "w");
  }
  catch (Exception e)
  {
    stderr.writeln("Failed to open file to write ranks. ", e.msg);
    exit(1);
  }

  double[size_t] rankCounts;

  foreach (ref e; genes)
  {
    if (e.rank in rankCounts)
    {
      rankCounts[e.rank]++;
    }
    else
    {
      rankCounts[e.rank] = 1.0;
    }
  }

  if (opts.verbose)
  {
    foreach (e; rankCounts.keys.sort!())
    {
      stderr.writeln(e, "\t", rankCounts[e]);
    }
  }

  double totalRanks = iota(0, 10).filter!(a => a in rankCounts)
    .map!(a => rankCounts[a]).sum;

  foreach (e; 0 .. 10)
  {
    rankFile.writeln(e, "\t", e in rankCounts ? rankCounts[e] / totalRanks : 0);
  }

  weightFile.writeln("0\t0");

  foreach (e; iota(1, 41, 2))
  {
    weightFile.writeln(genes[e * genes.length / 40].caveman, "\t",
        genes[(e - 1) * genes.length / 40 .. (e + 1) * genes.length / 40].count!(a => a.rank == 0)
        .to!double / (genes.length / 20).to!double);
  }

  weightFile.writeln("1\t1");
}
