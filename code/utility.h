//
// Created by Tian on 2023/7/14.
//

#ifndef ABEGIN_UTILITY_H
#define ABEGIN_UTILITY_H

#include "utility.h"

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <cstdio>
#include <cstdarg>

typedef unsigned long long kmer_int_type_t;
char int_to_base(int baseval);
int base_to_int(char nucleotide);

std::string revcomp (const std::string& kmer);
kmer_int_type_t kmer_to_intval(const std::string& kmer);
kmer_int_type_t kmer_to_intval(const std::string& kmer,unsigned int kmer_length);
std::string intval_to_kmer (kmer_int_type_t intval, unsigned int kmer_length);
bool contains_non_gatc(const std::string& kmer, unsigned int kmer_length);
kmer_int_type_t revcomp_val(kmer_int_type_t kmer, unsigned int kmer_length);
kmer_int_type_t get_DS_kmer_val(kmer_int_type_t kmer_val, unsigned int kmer_length);

void errAbort(char* format, ...);

#endif //ABEGIN_UTILITY_H
