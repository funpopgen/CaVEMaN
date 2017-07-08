module weights;

import arg_parse : Opts;
import core.stdc.stdlib : exit;
import std.algorithm : count, sort, uniq;
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
    causal = line[0].split("_")[1 .. $];
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
  allPVals.assumeSafeAppend;
  string currentGene = "";
  
  foreach (line; inFile.byLine)
  {
    auto splitLine = line.split;
    auto gene = splitLine[0].to!string;
    auto pVal = splitLine[6].to!double;

    if (currentGene != gene)
    {
      genes[$ - 1].rank = allPVals.sort!().uniq.count!(a => a < pVal);
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

  genes.sort!((a, b) => a.caveman > b.caveman);
  
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

  size_t[size_t] rankCounts;

  foreach (ref e; genes)
  {
    if (e.rank in rankCounts)
    {
      rankCounts[e.rank]++;
    }
    else
    {
      rankCounts[e.rank] = 1;
    }
  }

  foreach(e; rankCounts.keys.sort!())
  {
    rankFile.writeln(e, "\t", rankCounts[e]);
  }

  foreach (e; iota(1, 21, 3))
  {
    weightFile.writeln(genes[e * genes.length / 20].caveman,
		       "\t",
		       genes[0 .. (e * genes.length / 20 + 1)].count!(a => a.rank == 0).to!double / (e * genes.length / 20 + 1).to!double);
  }
}
