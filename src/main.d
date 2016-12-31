/* The GPL v3 License

   Copyright (C) 2016 University of Geneva.
   #
   # Author: Andrew Brown <andrew.brown@unige.ch>

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

import std.conv : to;
import std.range : enumerate;
import std.stdio : File, stderr, writeln;

import arg_parse : Opts;
import read_data : Phenotype, readBed, makeOut;
import run_analysis : CaVEMaN;
import calculation : genPerms;
import correct : correct;
import best : best;

version (STATICLINKED)
{
  pragma(msg, "Statically linked");
}
else
{
  pragma(lib, "gsl");
  pragma(lib, "gslcblas");
}

version (unittest)
  void main()
{
  writeln("All unit tests completed successfully.");
}

else
  void main(string[] args)
{
  pragma(msg, "CaVEMaN");

  auto opts = new Opts(args.to!(string[]));

  if (opts.correct != "")
  {
    if (opts.verbose)
    {
      stderr.writeln("Producing bed file corrected for multiple eQTL.");
    }

    correct(opts);
  }
  else if (opts.best != "")
  {
    if (opts.verbose)
    {
      stderr.writeln("Converting CaVEMaN score to causal probability.");
    }

    best(opts);
  }
  else
  {
    if (opts.verbose)
    {
      stderr.writeln("Running CaVEMaN analysis.");
    }

    auto phenotype = readBed(opts);

    if (opts.verbose)
    {
      stderr.writeln("Read ", phenotype.length, " phenotypes.");
    }

    auto permutations = genPerms(opts, phenotype[0].values.length);

    auto outFile = makeOut(opts);

    foreach (ref e; phenotype.enumerate)
    {
      if (opts.verbose)
      {
        stderr.writeln("Analysing gene ", e[1].geneName, " (", e[0] + 1,
            " out of ", phenotype.length, ").");
      }
      CaVEMaN(e[1], permutations, outFile, opts);
    }
  }
}

unittest
{

  import core.stdc.stdlib : exit;
  import std.digest.sha;
  import std.file : exists, remove;
  import std.range : put;
  import std.uuid;

  auto testFile1 = randomUUID.toString;
  auto testFile2 = randomUUID.toString;

  while (testFile1.exists || testFile2.exists)
  {
    testFile1 = randomUUID.toString;
    testFile2 = randomUUID.toString;
  }

  scope (exit)
  {
    if (testFile1.exists)
      testFile1.remove;
    if (testFile2.exists)
      testFile2.remove;
  }

  // ./bin/CaVEMaN --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 100000,4

  auto args = [
    "prog", "--bed", "data/phenotype.bed", "--job-number", "1", "--genes", "10",
    "--vcf", "data/genotype.vcf.gz", "--perm", "100000,4", "--out", testFile1
  ];

  auto opts = new Opts(args.to!(string[]));

  auto phenotype = readBed(opts);

  auto permutations = genPerms(opts, phenotype[0].values.length);

  auto outFile = makeOut(opts);

  foreach (ref e; phenotype)
  {
    CaVEMaN(e, permutations, outFile, opts);
  }

  outFile.close;

  SHA1 hash;
  hash.start;
  put(hash, File(testFile1).byChunk(1024));

  // LDC and DMD deal with precision of floating point numbers slightly differently

  version (LDC)
  {
    assert(toHexString(hash.finish) == "27909CDDEDD68A8DB35C8F1CC8421EB7C56825E2");
  }

  version (DigitalMars)
  {
    assert(toHexString(hash.finish) == "FA9DDAD12C9DC98C447AB0660D574C2E061DF869");
  }
  stderr.writeln("Passed: CaVEMaN.");

  // ./bin/CaVEMaN --best testtemp
  opts.best = testFile1;
  opts.output = testFile2;

  best(opts);

  hash.start;
  put(hash, File(testFile2).byChunk(1024));

  assert(toHexString(hash.finish) == "C530F1AC9E5D57C43F68FA64C7861780C2F742EE");
  stderr.writeln("Passed: extract best.");

  // ./bin/CaVEMaN --correct data/eQTL --bed data/phenotype.bed --vcf data/genotype.vcf.gz

  opts.correct = "data/eQTL";
  opts.output = testFile1;

  correct(opts);

  hash.start;
  put(hash, File(testFile1).byChunk(1024));

  assert(toHexString(hash.finish) == "FDED43C25211773C54FB6F854FFED8D0A0CFEA9C");
  stderr.writeln("Passed: correct phenotypes.");

  opts.normal = true;

  correct(opts);

  hash.start;
  put(hash, File(testFile1).byChunk(1024));

  assert(toHexString(hash.finish) == "BA58F5A4E604A5185270E074CC9BC754DD582C7E");
  stderr.writeln("Passed: correct with normalisation.");

  // ./bin/CaVEMaN --correct data/eQTL --bed data/phenotype.bed --vcf data/genotype.vcf.gz --cov data/covariates

  opts.normal = false;
  opts.cov = "data/covariates";
  opts.covLocations = [1, 0, 2, 3, 4, 5, 6, 7, 9, 8];

  correct(opts);

  hash.start;
  put(hash, File(testFile1).byChunk(1024));

  assert(toHexString(hash.finish) == "798E1AD6FF67BCEE6C96B19E8297F107439E6609");
  stderr.writeln("Passed: correct with covariates.");
}
