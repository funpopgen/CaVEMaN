import std.stdio : File, readln, stdout, stderr;
import std.conv : to;
import std.array : split, array;
import std.algorithm : sort;
import std.string : chomp;

extern (C)
{
  double interpolateSingleValue(double* weights_x, double* weights_y, size_t N, double interval);
}

struct Snp
{
  string line;
  double bootstrap;

  this(string snpLine, double boot)
  {
    line = snpLine;
    bootstrap = boot;
  }
}

void interval(string[] options, string output)
{
  // auto weights_x = [
  //   0.0, 0.03326427, 0.04492822, 0.05476801, 0.06435118, 0.07377145, 0.08316422, 0.09303571, 0.1037076,
  //   0.1155559, 0.128699, 0.1432871, 0.159728, 0.1792016, 0.2024936, 0.2308155,
  //   0.2651624, 0.3064091, 0.3584656, 0.4142872, 0.506331
  // ];

  // auto weights_y = [
  //   0.0, 0.169600754479723, 0.179217104228895, 0.185191007703191, 0.197767646596447,
  //   0.21411727715768, 0.229366451815752, 0.24697374626631, 0.268196824398679,
  //   0.304983493161453, 0.319710782772713, 0.359107339305359,
  //   0.40003143665514, 0.453316567117259, 0.514541738720327, 0.587486244301211, 0.661845621757585,
  //   0.767174972488602, 0.850337997170256, 0.934129853796573, 0.98538195535995
  // ];
  auto interval = cast(double[]) import("intervals");
  auto threshold = options.length > 1 ? options[1].to!double : 0.9;

  threshold = interpolateSingleValue(&interval[140664], &interval[0], 140664, threshold);

  version (unittest)
  {
  }
  else
  {
    stderr.writeln("Equivalent CaVEMaN threshold is ", threshold, ".");
  }

  File resultsFile;
  File outFile;

  try
  {
    resultsFile = File(options[0]);
  }
  catch (Exception e)
  {
    stderr.writeln("Failed to open results file. ", e.msg);
  }

  try
  {
    if (output == "")
    {
      outFile = stdout;
    }
    else
    {
      outFile = File(output, "w");
    }
  }
  catch (Exception e)
  {
    stderr.writeln("Failed to open output file. ", e.msg);
  }

  outFile.write(resultsFile.readln);

  Snp[] results;

  string gene = "";

  foreach (line; resultsFile.byLine)
  {
    auto splitLine = line.split;

    if (splitLine[0].to!string != gene)
    {
      results.sort!("a.bootstrap > b.bootstrap");
      auto count = 0.0;
      foreach (e; results)
      {
        outFile.writeln(e.line);
        count += e.bootstrap;
        // auto prob = interpolateSingleValue(weights_x.ptr, weights_y.ptr,
        //     weights_x.length, e.bootstrap);
        // outFile.writeln(e.line, "\t", prob);
        // count += prob;
        if (count > threshold)
          break;
      }
      gene = splitLine[0].to!string;
      results = [Snp(line.idup, splitLine[$ - 1].to!double)];
    }
    else
    {
      results ~= Snp(line.idup, splitLine[$ - 1].to!double);
    }
  }

  results.sort!("a.bootstrap > b.bootstrap");
  auto count = 0.0;
  foreach (e; results)
  {
    outFile.writeln(e.line);
    count += e.bootstrap;

    // auto prob = interpolateSingleValue(weights_x.ptr, weights_y.ptr,
    //     weights_x.length, e.bootstrap);
    // outFile.writeln(e.line, "\t", prob);
    // count += prob;
    if (count > threshold)
      break;
  }

}
