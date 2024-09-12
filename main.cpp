#include "main.h"
#include "common.h"
#include "sequence_graph.h"
#include "kmer_hash.h"
#include "utility.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cstdio>
#include <cstdarg>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <errno.h>
#include <libgen.h>
#include <ctime>

#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>

using namespace std;


const int MAX_STR = 10340;

void load_data(std::vector<std::string>& data, std::string file, int k) {

    data.reserve(100000);
    std::fstream in;
    in.open(file.c_str(), std::fstream::in);

    if (!in.is_open()) {
        std::cerr << "[error] File " << file << " can't be opened." << std::endl;
        exit(1);
    }

    char c_line[MAX_STR];
    if(reads_type == "fasta"){
        while(!in.eof()) {
            in.getline(c_line, MAX_STR);
            if (c_line[0] == '>') {
                in.getline(c_line, MAX_STR);
                std::string sequence(c_line);
                if (sequence.length() >= g_kmer_length){
                    data.push_back(sequence);
                }
            }
        }
    }
    else if(reads_type == "fastq"){
        if (k == 1){
            while(!in.eof()) {
                in.getline(c_line, MAX_STR);
                if (c_line[0] == '@') {
                    in.getline(c_line, MAX_STR);
                    std::string sequence(c_line);
                    if (sequence.length() >= g_kmer_length){
                        data.push_back(sequence);
                    }
                }
                if (c_line[0] == '+'){
                    in.getline(c_line, MAX_STR);
                }
            }
            cout << "input-left-read: " << data.size() << endl;
        } else if (k == 2){
            while(!in.eof()) {
                in.getline(c_line, MAX_STR);
                if (c_line[0] == '@') {
                    in.getline(c_line, MAX_STR);
                    std::string sequence(c_line);
                    if (sequence.length() >= g_kmer_length){
                        std::string re_sequence = revcomp(sequence);
                        data.push_back(re_sequence);
                    }
                }
                if (c_line[0] == '+'){
                    in.getline(c_line, MAX_STR);
                }
            }
            cout << "input-right-read: " << data.size() << endl;
        }
    }

    in.close();
}

std::string usage() {

    std::stringstream usage_info;
    usage_info
            << std::endl
            << "===============================================================================" << std::endl
            << " Usage: Assemble [--reads/--kmers] <filename>  [opts] " << std::endl
            << "===============================================================================" << std::endl
            << " **Required :" << std::endl
            << " --reads/-i <string>           " << ": the name of the file containing reads" << std::endl;
    //usage_info
    //  << " or " << std::endl
    //  << " --kmers <string>	       " << ": the name of the file containing kmers" << std::endl;
    usage_info
            << std::endl
            << " ** Optional :" << std::endl
            << " --kmer_length/-k <int>        " << ": length of kmer, default: 25." << std::endl
            << " --fasta/-a                    " << ": input reads file is in fasta format." << std::endl
            << " --fastq/-q                    " << ": input reads file is in fastq format." << std::endl
            << " --output_filename/-o <string> " << ": Name of the output file, default: paths.fasta." << std::endl
            << " --help/-h                     " << ": display the help information."<< std::endl
            << std::endl;
    usage_info
            << "================================================================================" << std::endl
            << std::endl;

    return (usage_info.str());
}

int parse_options(int argc, char* argv[]) {

    int option_index = 0 ;
    int next_option;

    while((next_option = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1) {
        switch (next_option){
            case -1:
                break;
            case 'k':
                g_kmer_length = atoi(optarg);
                break;
            case 'i':
                reads_file.push_back(optarg);
                break;
            case 'o':
                output_filename = optarg;
                break;
            case 'a':
                reads_type = "fasta";
                break;
            case 'q':
                reads_type = "fastq";
                break;
            case 'h':
                g_help = true;
                break;
            case OPT_DEBUG:
                g_debug = true;
                break;
            default:
                std::cout << usage();
                exit(1);
        }

    }
    if (g_help) {
        std::cout << usage() ;
        exit(1);
    }

    if (reads_file.size() == 0) {
        std::cerr << "Error : --input option needs an argument!! " << std::endl;
        std::cout << usage() ;
        exit(1);
    }

    if (g_kmer_length > 32) {
        errAbort(const_cast<char *>("Length of kmer can not be excess 32!\n"));
    }


    return 0;
}


int main(int argc, char* argv[]){

    int parse_ret = parse_options(argc, argv);
    if (parse_ret)
        return parse_ret;

    time_t begin = time(NULL);

    std::vector<std::string> data;

    std::cerr << "Begin loading reads ..." << std::endl;
    if(reads_type == "fasta"){
        load_data(data, reads_file.back() ,1);
    } else if (reads_type == "fastq"){
        int k = 1;
        for(const std::string& file : reads_file){
            std::cerr << file << std::endl;
            std::cerr << k << std::endl;
            if (file != "")
                load_data(data, file, k);
            k++;
        }
    }

    cout << data.size() <<" reads have been loaded !" << endl;

    KmerMap kmerMap(g_kmer_length);
    // build the hash
    if (!data.empty()) {
        kmerMap.get_hash(data);
    }else {
        errAbort(const_cast<char *>("Building kmer hash failed."));
    }


    int average;
    kmer_int_type_t seed_kmer = kmerMap.get_seed_kmer(average);
    cout << intval_to_kmer(seed_kmer,g_kmer_length) << "   " << average << endl;

    std::cout << "Begin loading sequence graph ..." << std::endl;
    time_t graph_begin = time(NULL);
    Sequence_graph sequenceGraph;
    sequenceGraph.output_graph(kmerMap,seed_kmer,average,data);
    time_t graph_end = time(NULL);
    std::cout << "Done sequence graph. (elapsed time: " << (graph_end-graph_begin) << " s)" << std::endl;

    time_t end = time(NULL);
    std::cout << " Elapsed time: " << (end-begin) <<" s"<< std::endl;
}