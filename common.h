#ifndef COMMON_H
#define COMMON_H

/*
 *  common.h
 */

#include <vector>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <bits/typesizes.h>
#include <bits/types.h>


// general options
extern int g_kmer_length;
extern bool g_help;
extern bool g_debug;

extern int decision_read;
extern float g_ratio;
extern int g_window_length;
extern int g_re_gene;
extern float g_com_kmer_ratio;
extern int avr_kmer_dep;
extern int first_p;
extern int second_p;
extern int dep_search;

// reads
extern std::vector<std::string> reads_file;
extern std::string reads_type;

extern std::string output_filename;

#define OPT_DEBUG			313

#endif

