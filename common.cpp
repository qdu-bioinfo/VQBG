/*
 *  common.h
 */
#include "common.h"
#include <cstring>
#include <cstdlib>
#include <errno.h>
#include <libgen.h>

// general
bool g_help = false;
bool g_debug = false;
int g_kmer_length = 25;

int decision_read = 0;
float g_ratio = 0.01;
int g_window_length = 1000;
int g_re_gene = 150;
float g_com_kmer_ratio = 0.3;
int avr_kmer_dep = 0;
int first_p = 1;
int second_p = 10;
int dep_search = 0;

// assemble
std::vector<std::string> reads_file;
std::string kmers_file = "";
std::string reads_type = "fastq";
std::string output_filename = "paths.fasta";
