//
// Created by Tian on 2023/7/14.
//

#include "utility.h"

static char _int_to_base [4] = {'G', 'A', 'T', 'C'};
unsigned char _base_to_int [256] = {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //   0-19
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  20-39
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  40-59
        255, 255, 255, 255, 255,   1, 255,   3, 255, 255, 255,   0, 255, 255, 255, 255, 255, 255, 255, 255, //  60-79
        255, 255, 255, 255,   2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   1, 255,   3, //  80-99
        255, 255, 255,   0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   2, 255, 255, 255, // 100-119
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 120-139
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 140-159
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 160-179
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 180-209
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 200-219
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 220-239
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255                      // 240-255
};

std::string revcomp (const std::string& kmer) {
    std::string revstring;
    for (int i = kmer.size() -1; i >= 0; --i) {
        char c = kmer[i];
        char revchar;
        switch (c) {
            case 'g':
                revchar = 'c';
                break;
            case 'G':
                revchar = 'C';
                break;
            case 'a':
                revchar = 't';
                break;
            case 'A':
                revchar = 'T';
                break;
            case 't':
                revchar = 'a';
                break;
            case 'T':
                revchar = 'A';
                break;
            case 'c':
                revchar = 'g';
                break;
            case 'C':
                revchar = 'G';
                break;
            default:
                revchar = 'N';
        }
        revstring += revchar;
    }
    return (revstring);
}

char int_to_base (int baseval) {
    if (baseval < 0 || baseval > 3) {
        std::cerr << "Error, baseval out of range 0-3" << std::endl;
        exit(1);
    }
    return(_int_to_base[baseval]);
}

int base_to_int (char nucleotide) {
    switch (nucleotide) {
        case 'G':
        case 'g':
            return(0);
        case 'A':
        case 'a':
            return(1);
        case 'T':
        case 't':
            return(2);
        case 'C':
        case 'c':
            return(3);
        default:
            return(-1);
    }
}


//将kmer的每个字符转化成数字并将所有数字合到一起
kmer_int_type_t kmer_to_intval (const std::string& kmer) {
    kmer_int_type_t kmer_val = 0;
    for (unsigned int i = 0; i < kmer.length(); ++i) {
        char c = kmer[i];
        int val = base_to_int(c);
        kmer_val = kmer_val << 2;
        kmer_val |= val;
    }
    return (kmer_val);
}
kmer_int_type_t kmer_to_intval(const std::string& kmer, unsigned int kmer_length) {
    kmer_int_type_t kmer_val = 0;
    for (unsigned int i = 0; i < kmer_length; ++i) {
        char c = kmer[i];
        //将kmer中的每一个字符对应成数字
        int val = base_to_int(c);
        kmer_val = kmer_val << 2;
        kmer_val |= val;
    }
    return (kmer_val);
}

std::string intval_to_kmer (kmer_int_type_t intval, unsigned int kmer_length) {
    // std::string kmer = "";
    std::string kmer(kmer_length, ' ');
    for (unsigned int i = 1; i <= kmer_length; ++i) {
        int base_num = intval & 3ll;//去最后位的字符
        // char base = int_to_base(base_num);
        kmer[kmer_length-i] = _int_to_base[base_num];//取最后位的字符放在字符串的最后一位
        intval = intval >> 2;//右移两位
        // kmer = base + kmer;
    }
    return (kmer);
}

bool contains_non_gatc (const std::string& kmer, unsigned int kmer_length) {
    for (unsigned int i = 0; i < kmer_length; ++i) {
        unsigned char c = kmer[i];
        if (_base_to_int[c] > 3)
            return(true);
    }
    return(false);
}

kmer_int_type_t revcomp_val(kmer_int_type_t kmer, unsigned int kmer_length) {
    kmer_int_type_t rev_kmer = 0;
    kmer = ~kmer;
    for (unsigned int i = 0; i < kmer_length; i++) {
        int base = kmer & 3;
        rev_kmer = rev_kmer << 2;
        rev_kmer += base;
        kmer = kmer >> 2;
    }
    return rev_kmer;
}

kmer_int_type_t get_DS_kmer_val(kmer_int_type_t kmer_val, unsigned int kmer_length) {
    kmer_int_type_t rev_kmer = revcomp_val(kmer_val, kmer_length);
    if (rev_kmer > kmer_val)
        kmer_val = rev_kmer;
    return(kmer_val);
}

void errAbort(char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[Error] ");
    vfprintf(stderr, format, args);
    va_end(args);
    exit(1);
}
