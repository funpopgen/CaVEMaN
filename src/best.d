module best;

import arg_parse : Opts;
import core.stdc.stdlib : exit;
import read_data : InputException;
import std.algorithm : map, max, min, reduce, sort;
import std.array : array, split;
import std.conv : to;
import std.exception : enforce;
import std.math : fabs;
import std.range : zip;
import std.stdio : File, readln, stderr, stdout, writeln;
import std.string : chomp;

extern (C)
{
  void interpolate(double* weights_x, double* weights_y, size_t N,
      double* caveman, double* predictions, size_t Npred);
}

struct Gene
{
  string[] lines;
  double caveman;
  double pVal;
  double cor;

  this(string line, double caveman_, double pval_, double cor_)
  {
    lines = [line];
    caveman = caveman_;
    pVal = pval_;
    cor = cor_;
  }

  void update(string line)
  {
    lines ~= line;
  }
}

void best(const Opts opts)
{
  File inFile;
  File outFile;

  try
  {
    inFile = File(opts.best);
  }
  catch (Exception e)
  {
    stderr.writeln("Failed to read results file. ", e.msg);
    exit(1);
  }

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
  }
  catch (Exception e)
  {
    stderr.writeln("Failed to open results file for writing. ", e.msg);
    exit(1);
  }

  outFile.writeln(inFile.readln.chomp, "\tProbability");

  Gene[string] genes;

  foreach (line; inFile.byLine)
  {
    auto splitLine = line.split;

    if (splitLine.length != 8)
    {
      stderr.writeln(
          "The following line appears to be truncated, this could be because of a previous job failing:");
      stderr.writeln(line);
    }
    {
      auto gene = splitLine[0].to!string;

      if (!(gene in genes))
      {
        genes[gene] = Gene(line.to!string, splitLine[7].to!double, splitLine[6].to!double, fabs(splitLine[5].to!double));
      }
      else
      {
        auto pVal = splitLine[6].to!double;
        if (genes[gene].pVal > pVal)
        {
          genes[gene] = Gene(line.to!string, splitLine[7].to!double, pVal, fabs(splitLine[5].to!double));
        }
        else if (genes[gene].pVal == pVal)
        {
	  auto cor = fabs(splitLine[5].to!double);
	  if (genes[gene].cor < cor)
	  {
	    genes[gene] = Gene(line.to!string, splitLine[7].to!double, pVal, cor);
	  }
	  else if (genes[gene].cor == cor)
	  {
	    genes[gene].update(line.to!string);
	  }
        }
      }
    }

  }

  double[] weights_x;
  double[] weights_y;

  if (opts.weights != "")
  {
    try
    {
      auto weightFile = File(opts.weights);
      foreach (line; weightFile.byLine)
      {
        auto splitLine = line.split.to!(double[]);
        enforce(splitLine.length == 2,
            new InputException("Some rows in weights file do not have 2 columns."));
        weights_x ~= splitLine[0];
        weights_y ~= splitLine[1];
      }
    }
    catch (Exception e)
    {
      stderr.write("Failed to read weights from file. ", e.msg);
      exit(1);
    }
  }
  else
  {
    // Average over all UK10K simulations
    // weights_x = [0.0, 0.0252030825, 0.038378465, 0.0485059625, 0.05789589,
    //   0.0671174, 0.0762736125, 0.085624775, 0.0955709875,
    //   0.106504425, 0.118648825, 0.132308225, 0.147419575, 0.1648905, 0.185910125,
    //   0.210655575, 0.241188625, 0.278071125, 0.323693, 0.379883925, 0.441420325, 1];
    // weights_y = [0.0, 0.178678038379531, 0.182632177373508, 0.190306992609437, 0.203553660270078,
    //   0.214641080312722, 0.23152359295054, 0.250035536602701, 0.267623649801023,
    //   0.302358624609264, 0.32646389994315, 0.356828193832599,
    //   0.405002131590166, 0.45316275764037, 0.511727078891258, 0.575955662924542, 0.660034110289937,
    //   0.760483297796731, 0.843803297328027, 0.930490405117271, 0.982802728823195, 1];
    // Blood
    // weights_x = [0.0, 0.02115411, 0.03142417, 0.03935235, 0.04585157, 0.05300084, 0.0601453, 0.06714237, 0.0748771,
    //   0.08302726, 0.09167466, 0.1021211, 0.1135114, 0.126931, 0.1425238,
    //   0.1639797, 0.1883623, 0.2180209, 0.2571315, 0.3038833, 0.3686424, 1];
    // weights_y = [0.0, 0.192878338278932, 0.21037037037037, 0.216617210682493, 0.228486646884273,
    //   0.248888888888889, 0.260740740740741, 0.264094955489614, 0.267062314540059,
    //   0.302670623145401, 0.318518518518519, 0.373333333333333,
    //   0.403560830860534, 0.449554896142433, 0.49406528189911, 0.58962962962963, 0.644444444444444,
    //   0.692878338278932, 0.807121661721068, 0.891691394658754, 0.968888888888889, 1];
    // Fat
    // weights_x = [0.0, 0.0254682825, 0.03931539, 0.050174775, 0.059942925,
    //   0.0695934, 0.0789265775, 0.08862447, 0.098663675,
    //   0.109898175, 0.122507525, 0.1367777, 0.15123185, 0.16863025, 0.18972085,
    //   0.214642625, 0.245544125, 0.2831702, 0.328188125, 0.384124675, 0.442341125, 1];
    // weights_y = [0.0, 0.167924528301887, 0.169811320754717, 0.184905660377358, 0.197169811320755,
    //   0.216509433962264, 0.242924528301887, 0.24622641509434, 0.261320754716981,
    //   0.304245283018868, 0.32311320754717, 0.367452830188679,
    //   0.398584905660377, 0.463679245283019, 0.514622641509434, 0.577830188679245, 0.650943396226415,
    //   0.768976897689769, 0.834040546911834, 0.936792452830189, 0.985849056603774, 1];
    // LCL
    // weights_x = [0.0, 0.02600062, 0.03972796, 0.0502141, 0.05987642, 0.06942618, 0.07932034, 0.08915016,
    //   0.0997241, 0.1113732, 0.124406, 0.1384526, 0.1548056, 0.173301, 0.195358,
    //   0.2215012, 0.2543894, 0.2928926, 0.340405, 0.3952966, 0.454124, 1];
    // weights_y = [0.0, 0.172151898734177, 0.179822710004221, 0.187420852680456,
    //   0.194596876319122, 0.215189873417722, 0.233755274261603,
    //   0.250738708315745, 0.287040945546644,
    //   0.322648671446647, 0.324472573839662, 0.362447257383966,
    //   0.403123680878008, 0.457154917686788, 0.523206751054852, 0.582032897511598, 0.668776371308017,
    //   0.778387505276488, 0.857745884339384, 0.938370620514985, 0.986919831223629, 1];
    // Skin
    // weights_x = [0.0, 0.02610012, 0.0392731025, 0.0495225875, 0.05899701,
    //   0.06800834, 0.076990945, 0.0859244025, 0.095813275,
    //   0.106393625, 0.118231725, 0.1318257, 0.146937375, 0.1640685, 0.184559925,
    //   0.209414, 0.24004685, 0.27599925, 0.322719125, 0.37739245, 0.43393565, 1];
    // weights_y = [0.0, 0.169871794871795, 0.175213675213675, 0.197115384615385, 0.200854700854701,
    //   0.203525641025641, 0.22008547008547, 0.239978621058258, 0.253739316239316,
    //   0.283119658119658, 0.30715811965812, 0.355235042735043,
    //   0.387286324786325, 0.443910256410256, 0.505077498663816, 0.577457264957265, 0.672542735042735,
    //   0.763354700854701, 0.848824786324786, 0.935363247863248, 0.981303418803419, 1];

    // Excluding blood (outlier)

    weights_x = [0, 0.02584941, 0.03948495, 0.0499885, 0.059601455, 0.06904272, 0.0784325, 0.08803015, 0.098245775,
      0.109414, 0.1219489, 0.13584255, 0.15114695, 0.168945, 0.19031085,
      0.21567395, 0.247001, 0.2844645, 0.330803, 0.38620805, 0.4449283, 1];

    weights_y = [0, 0.169600754479723, 0.179217104228895, 0.185191007703191, 0.197767646596447,
      0.21411727715768, 0.229366451815752, 0.24697374626631, 0.268196824398679,
      0.304983493161453, 0.319710782772713, 0.359107339305359,
      0.40003143665514, 0.453316567117259, 0.514541738720327, 0.587486244301211, 0.661845621757585,
      0.767174972488602, 0.850337997170256, 0.934129853796573, 0.98538195535995, 1];
  }

  weights_y.sort!();

  if (opts.verbose)
  {
    stderr.writeln(genes.length, " genes in the results file.");
  }

  auto caveman = genes.byKey.map!(a => genes[a].caveman).array;

  auto maxCaveman = reduce!(min, max)(caveman);

  if (weights_x[0] > maxCaveman[0])
  {
    weights_x[0] = maxCaveman[0];
  }

  if (weights_x[$ - 1] < maxCaveman[1])
  {
    weights_x[$ - 1] = maxCaveman[1];
  }

  double[] predictions = new double[](caveman.length);

  interpolate(weights_x.ptr, weights_y.ptr, weights_x.length, caveman.ptr,
      predictions.ptr, caveman.length);

  double[double] cavemanLookup;

  foreach (e; zip(caveman, predictions))
  {
    cavemanLookup[e[0]] = e[1];
  }

  foreach (e; genes.keys.sort!())
  {
    auto prediction = cavemanLookup[genes[e].caveman];
    foreach (line; genes[e].lines)
    {
      outFile.writeln(line, "\t", prediction);
    }
  }
}
