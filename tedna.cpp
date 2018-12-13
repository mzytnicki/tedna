/**
Copyright (C) 2013 INRA-URGI
This file is part of TEDNA, a short reads transposable elements assembler
TEDNA is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
GNU Affero General Public License for more details.
See the GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License
along with this program.
**/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include "optionparser.h"
#include "assembler.hpp"

enum  optionIndex {UNKNOWN, INPUT1, INPUT2, INSERT, KMER, OUTPUT, THRESHOLD, PROCESSORS, REPEAT_FREQUENCY, MIN_FREQUENCY, FREQUENCY_DIF, SMALL_GRAPH, BIG_GRAPH, NB_SMALL_GRAPH, MAX_PATHS, EROSION, BUBBLE_SIZE, MIN_LTR, MAX_LTR, MAX_IDENTITY, MIN_OVERLAP, MAX_OVERLAP, SHORT_KMER, INDEL_PEN, MISMATCH_PEN, SIZE_PEN, MAX_PEN, MIN_IDENTITY, MERGE_MAX_NB, MERGE_MAX_NODES, MIN_SCAFFOLD, MAX_SCAFFOLD, SCAFFOLD_MAX_EV, MAX_EVIDENCES, MIN_TE_SIZE, MAX_TE_SIZE, FASTA_INPUT, BYTES_PER_THREAD, MAX_KMERS, MAX_READS, CHECK, HELP, VERSION};
const option::Descriptor usage[] = {
	{UNKNOWN,          0, "" , ""                  , option::Arg::None    , "USAGE: tedna [options]\n\n" "Compulsory options:"},
	{INPUT1,           0, "1", "file1"             , option::Arg::Required, "  -1, --file1  \tFirst FASTQ file."},
	{KMER,             0, "k", "kmer"              , option::Arg::Numeric,  "  -k, --kmer   \tK-mer size."},
	{OUTPUT,           0, "o", "output"            , option::Arg::Required, "  -o, --output \tOutput file."},
	{UNKNOWN,          0, "" ,  ""                 , option::Arg::None    , "\nOther options:"},
	{INPUT2,           0, "2", "file2"             , option::Arg::Required, "  -2, --file2        \tSecond FASTQ file."},
	{INSERT,           0, "i", "insert"            , option::Arg::Numeric,  "  -i, --insert \tInsert size."},
	{THRESHOLD,        0, "t", "threshold"         , option::Arg::Numeric,  "  -t, --threshold    \tK-mer frequency threshold   (default: ad hoc)."},
	{MIN_TE_SIZE,      0, "m", "min-te-size"       , option::Arg::Numeric,  "  -m, --min-te-size  \tMinimum TE size             (default: 500)."},
	{MAX_TE_SIZE,      0, "M", "max-te-size"       , option::Arg::Numeric,  "  -M, --max-te-size  \tMaximum TE size             (default: 30000)."},
	{PROCESSORS,       0, "p", "processors"        , option::Arg::Numeric,  "  -p, --processors   \tNumber of processors        (default: 2)."},
	{HELP,             0, "h", "help"              , option::Arg::None    , "  --help   \tPrint usage and exit."},
	{VERSION,          0, "v", "version"           , option::Arg::None    , "  --version\tPrint version and exit."},
	{UNKNOWN,          0, "" ,  ""                 , option::Arg::None    , "\nObscure options:"},
	{UNKNOWN,          0, "" ,  ""                 , option::Arg::None    , "\n  k-mer frequency:"},
	{REPEAT_FREQUENCY, 0, "" , "repeat-frequency"  , option::Arg::Numeric,  "  --repeat-frequency   \tMinimum number of repetitions      (default: 2)."},
	{MIN_FREQUENCY,    0, "" , "min-frequency"     , option::Arg::Numeric,  "  --min-frequency      \tMinimum k-mer frequency            (default: 3)."},
	{FREQUENCY_DIF,    0, "" , "frequency-dif"     , option::Arg::Numeric,  "  --frequency-dif      \tMaximum k-mer frequency difference (default: 2.5)."},
	{UNKNOWN,          0, "" ,  ""                 , option::Arg::None    , "\n  graph:"},
	{SMALL_GRAPH,      0, "" , "small-graph"       , option::Arg::Numeric,  "  --small-graph        \tMinimum graph size                 (default: 300)."},
	{BIG_GRAPH,        0, "" , "big-graph"         , option::Arg::Numeric,  "  --big-graph          \tMaximum graph size                 (default: 100000)."},
	{NB_SMALL_GRAPH,   0, "" , "small-graph-count" , option::Arg::Numeric,  "  --small-graph-count  \tStop after N small graphs          (default: 10000), 0: never stop."},
	{MAX_PATHS,        0, "" , "max-paths"         , option::Arg::Numeric,  "  --max-paths          \tMaximum # paths                    (default: 100), 0: never stop."},
	{EROSION,          0, "" , "erosion"           , option::Arg::Numeric,  "  --erosion            \tErosion strength                   (default: 100)."},
	{BUBBLE_SIZE,      0, "" , "bubble-size"       , option::Arg::Numeric,  "  --bubble-size        \tSize of the bubbles                (default: 1000)."},
	{UNKNOWN,          0, "" ,  ""                 , option::Arg::None    , "\n  LTR elements:"},                                      
	{MIN_LTR,          0, "" , "min-ltr"           , option::Arg::Numeric,  "  --min-ltr            \tMinimum LTR size                   (default: 50)."},
	{MAX_LTR,          0, "" , "max-ltr"           , option::Arg::Numeric,  "  --max-ltr            \tMaximum LTR size                   (default: 5000)."},
	{UNKNOWN,          0, "" ,  ""                 , option::Arg::None    , "\n  duplicate removal:"},                                       
	{MAX_IDENTITY,     0, "" , "duplicate-id"      , option::Arg::Numeric,  "  --duplicate-id       \tMaximum id. to remove duplicate    (default: 30%)."},
	{UNKNOWN,          0, "" ,  ""                 , option::Arg::None    , "\n  merge:"},                                             
	{MIN_OVERLAP,      0, "" , "min-overlap"       , option::Arg::Numeric,  "  --min-overlap        \tMinimum overlap to merge TEs       (default: 20)."},
	{MAX_OVERLAP,      0, "" , "max-overlap"       , option::Arg::Numeric,  "  --max-overlap        \tMaximum overlap to merge TEs       (default: 500)."},
	{SHORT_KMER,       0, "" , "short-kmer"        , option::Arg::Numeric,  "  --short-kmer         \tSmall k-mer size                   (default: 15)."},
	{INDEL_PEN,        0, "" , "indel-pen"         , option::Arg::Numeric,  "  --indel-pen          \tIndel penalty                      (default: 30)."},
	{MISMATCH_PEN,     0, "" , "mismatch-pen"      , option::Arg::Numeric,  "  --mismatch-pen       \tMismatch penalty                   (default: 10)."},
	{SIZE_PEN,         0, "" , "size-pen"          , option::Arg::Numeric,  "  --size-pen           \tSize penalty                       (default: 1)."},
	{MAX_PEN,          0, "" , "max-pen"           , option::Arg::Numeric , "  --max-pen            \tMaximum penalty                    (default: 200)."},
	{MIN_IDENTITY,     0, "" , "min-id"            , option::Arg::Numeric , "  --min-id             \tMinimum identity                   (default: 80)."},
	{MERGE_MAX_NB,     0, "" , "merge-max-nb"      , option::Arg::Numeric , "  --merge-max-nb       \tMaximum number of repeats          (default: 10000), 0: do not use."},
	{MERGE_MAX_NODES,  0, "" , "merge-max-nodes"   , option::Arg::Numeric , "  --merge-max-nodes    \tMaximum number of neighbor/node    (default: 10), 0: do not use."},
	{UNKNOWN,          0, "" ,  ""                 , option::Arg::None    , "\n  scaffolding:"},                                       
	{MIN_SCAFFOLD,     0, "" , "min-scaffold"      , option::Arg::Numeric , "  --min-scaffold       \tMinimum number of evidences/scaff. (default: 100)."},
	{MAX_SCAFFOLD,     0, "" , "max-scaffold"      , option::Arg::Numeric , "  --max-scaffold       \tMaximum number of evidences/scaff. (default: 10000)."},
	{SCAFFOLD_MAX_EV,  0, "" , "scaffold-max-nb"   , option::Arg::Numeric , "  --scaffold-max-nb    \tMaximum number of neighbor/node    (default: 5), 0: do not use."},
	{UNKNOWN,          0, "" ,  ""                 , option::Arg::None    , "\n  input reading:"},                                             
	{FASTA_INPUT,      0, "" , "fasta"             , option::Arg::None    , "  --fasta-input        \tInput file is in FASTA format      (default: not set)."},
	{BYTES_PER_THREAD, 0, "" , "bytes-per-thread"  , option::Arg::Numeric , "  --bytes-per-thread   \tNumber of bytes read per thread    (default: 10000000)."},
	{MAX_READS,        0, "" , "max-reads"         , option::Arg::Numeric , "  --max-reads          \tMaximum number of reads read       (default: 0), 0: read all."},
	{CHECK,            0, "" , "check"             , option::Arg::Optional, "  --check              \tCheck if a sequence is assembled   (default: none)."},
	{UNKNOWN,          0, "" , ""                  , option::Arg::None    , "\nExample:\n  ./tedna -1 left.fastq -2 right.fastq -k 61 -i 300 -o output.fasta"},
	{0,0,0,0,0,0}
};

int main(int argc, char **argv) {

	argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
	option::Stats  stats(usage, argc, argv);
	option::Option *options = new option::Option[stats.options_max];
	option::Option *buffer  = new option::Option[stats.buffer_max];
	option::Parser parse(usage, argc, argv, options, buffer);

	if (parse.error())
		return 1;

	if (options[HELP] || argc == 0) {
		option::printUsage(std::cout, usage);
		return 0;
	}

	if (options[VERSION] || argc == 0) {
		cout << "Using tedna version " << Globals::VERSION << endl;
		return 0;
	}

	bool unknownOption = false;
	for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next()) {
		cout << "Error: unknown option: " << opt->name << endl;
		unknownOption = true;
	}
	if (unknownOption) {
		return 1;
	}

	if (! options[INPUT1]) {
		cout << "Error: first input FASTQ file is missing." << endl;
		option::printUsage(std::cout, usage);
		return 1;
	}
	ifstream f(options[INPUT1].arg);
    if (f.good()) {
        f.close();
	}
	else {
		cout << "Error: cannot open first input FASTQ file ('" << options[INPUT1].arg << "')." << endl;
		return 1;
	}
	if (options[INPUT2]) {
    ifstream g(options[INPUT2].arg);
    if (g.good()) {
        g.close();
    }
    else {
      cout << "Error: cannot open second input FASTQ file ('" << options[INPUT2].arg << "')." << endl;
      return 1;
    }
	}
	if (! options[OUTPUT]) {
		cout << "Error: output file is missing." << endl;
		option::printUsage(std::cout, usage);
		return 1;
	}
	if ((! options[INSERT]) && (options[INPUT2])) {
		cout << "Error: insert size is missing." << endl;
		option::printUsage(std::cout, usage);
		return 1;
	}

	const char* fileName1      = options[INPUT1].arg;
	const char* fileName2      = (options[INPUT2])? options[INPUT2].arg: nullptr;
	const char *outputFileName = options[OUTPUT].arg;
	int         insertSize     = (options[INSERT])? atoi(options[INSERT].arg): 0;
	int         thresholdPc    = -1;

	if (options[KMER])
		Globals::KMER = strtoul(options[KMER].arg, NULL, 0);
	if (options[MIN_TE_SIZE])
		Globals::MIN_TE_SIZE = atoi(options[MIN_TE_SIZE].arg);
	if (options[THRESHOLD])
		thresholdPc = atoi(options[THRESHOLD].arg);
	if (options[PROCESSORS])
		Globals::NB_THREADS = atoi(options[PROCESSORS].arg);
	if (options[MAX_TE_SIZE])
		Globals::MAX_TE_SIZE = atoi(options[MAX_TE_SIZE].arg);
	if (options[MIN_TE_SIZE])
		Globals::MIN_TE_SIZE = atoi(options[MIN_TE_SIZE].arg);
	if (options[REPEAT_FREQUENCY])
		Globals::NB_REPETITIONS = atoi(options[REPEAT_FREQUENCY].arg);
	if (options[MIN_FREQUENCY])
		Globals::MIN_COUNT = atoi(options[MIN_FREQUENCY].arg);
	if (options[FREQUENCY_DIF])
		Globals::FREQUENCY_DIFFERENCE = atoi(options[FREQUENCY_DIF].arg);
	if (options[MIN_TE_SIZE])
		Globals::MIN_TE_SIZE = atoi(options[MIN_TE_SIZE].arg);
	if (options[SMALL_GRAPH])
		Globals::MIN_NB_NODES = atoi(options[SMALL_GRAPH].arg);
	if (options[BIG_GRAPH])
		Globals::MAX_NB_NODES = atoi(options[BIG_GRAPH].arg);
	if (options[NB_SMALL_GRAPH])
		Globals::NB_SMALL_GRAPHS = strtoul(options[NB_SMALL_GRAPH].arg, NULL, 0);
	if (options[MAX_PATHS])
		Globals::MAX_PATHS = atoi(options[MAX_PATHS].arg);
	if (options[EROSION])
		Globals::EROSION_STRENGTH = atoi(options[EROSION].arg);
	if (options[BUBBLE_SIZE])
		Globals::BUBBLE_SIZE = atoi(options[BUBBLE_SIZE].arg);
	if (options[MIN_LTR])
		Globals::MIN_LTR_SIZE = atoi(options[MIN_LTR].arg);
	if (options[MAX_LTR])
		Globals::MAX_LTR_SIZE = atoi(options[MAX_LTR].arg);
	if (options[MAX_IDENTITY])
		Globals::MAX_IDENTITY = atoi(options[MAX_IDENTITY].arg) / 100.0;
	if (options[MIN_OVERLAP])
		Globals::MIN_MERGE_SIZE = atoi(options[MIN_OVERLAP].arg);
	if (options[MAX_OVERLAP])
		Globals::MAX_MERGE_SIZE = atoi(options[MAX_OVERLAP].arg);
	if (options[SHORT_KMER])
		Globals::SHORT_KMER_SIZE = strtoul(options[SHORT_KMER].arg, NULL, 0);
	if (options[INDEL_PEN])
		Globals::PENALTY_INDEL = atoi(options[INDEL_PEN].arg);
	if (options[MISMATCH_PEN])
		Globals::PENALTY_MISMATCH = atoi(options[MISMATCH_PEN].arg);
	if (options[SIZE_PEN])
		Globals::PENALTY_SIZE = atoi(options[SIZE_PEN].arg);
	if (options[MAX_PEN])
		Globals::MAX_PENALTY = atoi(options[MAX_PEN].arg);
	if (options[MIN_IDENTITY])
		Globals::MIN_IDENTITY = atoi(options[MIN_IDENTITY].arg) / 100.0;
	if (options[MERGE_MAX_NB])
		Globals::MERGE_MAX_NB = atoi(options[MERGE_MAX_NB].arg);
	if (options[MERGE_MAX_NODES])
		Globals::MERGE_MAX_NODES = atoi(options[MERGE_MAX_NODES].arg);
	if (options[MIN_SCAFFOLD])
		Globals::MIN_SCAFFOLD_KMERS = strtoul(options[MIN_SCAFFOLD].arg, NULL, 0);
	if (options[MAX_SCAFFOLD])
		Globals::MAX_SCAFFOLD_COUNTS = strtoul(options[MAX_SCAFFOLD].arg, NULL, 0);
	if (options[SCAFFOLD_MAX_EV])
		Globals::SCAFFOLD_MAX_EV = strtoul(options[SCAFFOLD_MAX_EV].arg, NULL, 0);
	if (options[FASTA_INPUT])
		Globals::FASTA_INPUT = true;
	if (options[BYTES_PER_THREAD])
		Globals::SIZE_THREAD = atol(options[BYTES_PER_THREAD].arg);
	if (options[MAX_READS])
		Globals::NB_READS = strtoul(options[MAX_READS].arg, NULL, 0);
	if (options[CHECK])
		Globals::CHECK = options[CHECK].arg;

	if (Globals::NB_BLOCKS * 32 - 1 < Globals::KMER) {
		cerr << "Cannot use k-mer of size " << Globals::KMER << " with current install. Please re-compile with option 'make k=" << Globals::KMER << "'." << endl;
		return 1;
	}

	Assembler assembler(fileName1, fileName2, outputFileName, insertSize, thresholdPc);
	assembler.assemble();
	assembler.dump();

	delete[] options;
	delete[] buffer;
	cout << "Successfully done." << endl;
}
