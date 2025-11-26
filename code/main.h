#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include "common.h"

#include <getopt.h>

static const char *short_options = "k:i:t:o:a:q";

static struct option long_options[] = {
        // general options
        {"help",no_argument,0,'h'},
        {"kmer_length",required_argument,0,'k'},
        {"reads",required_argument,0,'i'},
        {"trunk_filename",required_argument, 0,'t'},
        {"output_filename",required_argument, 0,'o'},
        {"fasta",required_argument, 0,'a'},
        {"fastq",required_argument, 0,'q'},
        {"debug",no_argument,0,OPT_DEBUG},
        {0,0,0,0} // terminator

};

int parse_options(int argc, char* argv[]);