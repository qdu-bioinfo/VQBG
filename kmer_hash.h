#ifndef KMER_HASH_H     // 添加这行
#define KMER_HASH_H
#include "utility.h"
#include "common.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <errno.h>
#include <libgen.h>
#include <ctime>

#include <boost/unordered_map.hpp>

using namespace std;

typedef pair<kmer_int_type_t, size_t> kmer_occurence_pair_t;
typedef pair<size_t, vector<pair<size_t,size_t>>> read_pair;
class KmerMap{
private:

    typedef typename boost::unordered_map<kmer_int_type_t, read_pair> kmer_map_base_t;
    typedef typename boost::unordered_map<kmer_int_type_t, read_pair>::iterator kmer_map_base_iterator_t;

    class kmer_sorter_by_count_desc_t {
    public:
        kmer_sorter_by_count_desc_t(KmerMap& k) : kmer_map_(k) {};

        bool operator() (const kmer_int_type_t& i,
                         const kmer_int_type_t& j) {
            return ((kmer_map_.get_kmer_count(i) > kmer_map_.get_kmer_count(j))
                    || (kmer_map_.get_kmer_count(i) == kmer_map_.get_kmer_count(j) && i > j) );
        }

        bool operator() (const kmer_occurence_pair_t& i,
                         const kmer_occurence_pair_t& j) {
            return ( (i.second > j.second)
                     || (i.second == j.second && i.first > j.first) );
        }

        bool operator() (const kmer_map_base_iterator_t& i,
                         const kmer_map_base_iterator_t& j) {
            size_t count_i = i->second.first;
            size_t count_j = j->second.first;
            return ( (count_i > count_j)
                     || ( count_i == count_j && i->first > j->first));
        }

        bool operator() (const std::string& i, const std::string& j) {
            kmer_int_type_t val_i = kmer_to_intval(i,g_kmer_length);//kmer的字符值转化为对应的二进制整数
            kmer_int_type_t val_j = kmer_to_intval(j,g_kmer_length);
            return ( (kmer_map_.get_kmer_count(val_i) > kmer_map_.get_kmer_count(val_j))
                     || (kmer_map_.get_kmer_count(val_i) == kmer_map_.get_kmer_count(val_j) && val_i > val_j) );
        }

    private:
        KmerMap& kmer_map_;
    };

    kmer_map_base_t kmer_map;
    int kmer_length;
    size_t count;

public:
    KmerMap () {}

    KmerMap (int kmer_length = g_kmer_length) {
        this->kmer_length = kmer_length;
    }


    size_t get_size() {
        return kmer_map.size();
    }

    bool empty() {
        return kmer_map.empty();
    }

    kmer_map_base_iterator_t find_kmer(kmer_int_type_t kmer_val) {
        return kmer_map.find(kmer_val);
    }

    inline size_t get_kmer_count(kmer_int_type_t kmer_val) {

        kmer_map_base_iterator_t it = find_kmer(kmer_val);

        if (it != kmer_map.end())
            return it->second.first;
        else
            return 0;
    }

    inline size_t get_kmer_count(const string& kmer) {

        kmer_int_type_t kmer_val = kmer_to_intval(kmer,kmer_length);
        return (get_kmer_count(kmer_val));
    }

    // check a kmer exist or not
    bool exists(const kmer_int_type_t kmer_val) {

        return (get_kmer_count(kmer_val) > 0);
    }

    bool exists(const std::string& kmer) {

        kmer_int_type_t kmer_val = kmer_to_intval(kmer,kmer_length);
        return (exists(kmer_val));
    }

    std::vector<pair<size_t,size_t>> get_kmer_read(kmer_int_type_t kmer_val){
        std::vector<pair<size_t,size_t>> null_read;
        kmer_map_base_iterator_t it = find_kmer(kmer_val);
        if (it != kmer_map.end())
            return it->second.second;
        else
            return null_read;
    }

    std::vector<pair<size_t,size_t>> get_kmer_read(const string& kmer){
        kmer_int_type_t kmer_val = kmer_to_intval(kmer,kmer_length);
        return (get_kmer_read(kmer_val));
    }

    void get_hash(std::vector<std::string>& data) {

        size_t data_size = data.size();
        std::cerr << "Beginning kmer hash ..." << std::endl;
        time_t beg = time(NULL);

        if (data.empty())
            return;

        for (size_t i = 0; i < data_size; ++i) {

            const std::string& sequence = data[i];
            int n = 0;
            for (size_t j = 0; j <= sequence.length()-kmer_length; ++j) {

                const std::string& kmer = sequence.substr(j, kmer_length);
                //检查kmer是否只包含ATGC
                if (contains_non_gatc(kmer, kmer_length)){
                    n++;
                    continue;
                }
                //获得kmer的二进制数字形式
                kmer_int_type_t kmer_val = kmer_to_intval(kmer, kmer_length);
                kmer_map[kmer_val].first++;
//                if (j == n){
                    pair<size_t,size_t> read_infor;
                    read_infor.first = i;
                    read_infor.second = j;
                    kmer_map[kmer_val].second.push_back(read_infor);
//                }

            }
        }

        time_t end = time(NULL);
        std::cerr << "Kmer hash finished, get " << kmer_map.size() << " kmers! (elapsed time: " << (end-beg) << " s)" << std::endl;

    }

    int get_forward_candidates(kmer_int_type_t seed_kmer, std::vector<kmer_occurence_pair_t>& candidates) {

        candidates.clear();
        kmer_int_type_t forward_prefix = (seed_kmer << (33-kmer_length)*2) >> (32-kmer_length)*2;
        int sum = 0 , candidate_average = 0;
        for (kmer_int_type_t i = 0; i < 4; ++i) {

            kmer_occurence_pair_t candidate;
            candidate.first = forward_prefix | i;
            candidate.second = get_kmer_count(candidate.first);

//            std::cerr << intval_to_kmer(candidate.first,g_kmer_length) << "  : " << candidate.second << std::endl;

            if (candidate.second != 0) {
                candidates.push_back(candidate);
                sum = sum + candidate.second;
            }

        }

        if (candidates.size() != 0){
            candidate_average = sum / candidates.size();
        }
        kmer_sorter_by_count_desc_t sorter(*this);
        sort (candidates.begin(), candidates.end(), sorter);

        return candidate_average;
    }

    int get_reverse_candidates (kmer_int_type_t seed_kmer, std::vector<kmer_occurence_pair_t>& candidates) {

        candidates.clear();
        kmer_int_type_t reverse_suffix = seed_kmer >> 2;
        int sum = 0 , candidate_average = 0;

        for (kmer_int_type_t i = 0; i < 4; ++i) {

            kmer_occurence_pair_t candidate;
            candidate.first = (i << (kmer_length*2-2)) | reverse_suffix;
            candidate.second = get_kmer_count(candidate.first);

            if (candidate.second != 0) {
                candidates.push_back(candidate);
                sum = sum + candidate.second;
            }

        }
        if (candidates.size() != 0){
            candidate_average = sum / candidates.size();
        }


        kmer_sorter_by_count_desc_t sorter(*this);
        sort (candidates.begin(), candidates.end(), sorter);

        return candidate_average;
    }


    kmer_int_type_t get_seed_kmer(int& average){
        kmer_map_base_iterator_t it;
        int max = 0;
        int sum = 0;
        kmer_int_type_t seed_kmer = 0;
        for (it = kmer_map.begin(); it != kmer_map.end(); it++) {
            if (it->second.first > max){
                max = it->second.first;
                seed_kmer = it->first;
            }
            sum = sum + it->second.first;
        }
        average = sum / kmer_map.size();
        return seed_kmer;
    }

    int get_kmer(std::vector<kmer_int_type_t>& kmer_list){
        kmer_map_base_iterator_t it;
        int sum = 0;
        for (it = kmer_map.begin(); it != kmer_map.end(); it++) {
//            if(get_kmer_count(it->first)<=1)
//                continue;
            kmer_list.push_back(it->first);
            sum += it->second.first;
        }
        int average = sum / kmer_map.size();
        time_t beg = time(NULL);
        std::cerr << "Sorting kmers ..." << std::endl;
        kmer_sorter_by_count_desc_t sorter(*this);
        sort(kmer_list.begin(), kmer_list.end(), sorter);
        time_t end = time(NULL);
        std::cerr << "Done sorting. (elapsed time: " << (end-beg) << " s)" << std::endl;
        return average;
    }

    std::vector<kmer_int_type_t> get_all_kmer(){
        std::vector<kmer_int_type_t> kmer_list;
        kmer_map_base_iterator_t it;
        int sum = 0;
        for (it = kmer_map.begin(); it != kmer_map.end(); it++) {
            kmer_list.push_back(it->first);
        }
        return kmer_list;
    }

    void kmer_ca(){
        kmer_map_base_iterator_t it;
        std::ofstream nodefile("6polio_kmer.fasta");
        for (it = kmer_map.begin(); it != kmer_map.end(); it++) {
            nodefile << "kmer : " << intval_to_kmer(it->first,g_kmer_length) <<endl;
            std::vector<kmer_occurence_pair_t> candidates;
            int candidate_average = get_forward_candidates(it->first, candidates);
            for (int i = 0; i < candidates.size(); ++i) {
                nodefile << "cand : " << intval_to_kmer(candidates[i].first,g_kmer_length) <<endl;
            }
        }
        nodefile.close();
    }
    int get_kmer_read_num(){
        kmer_map_base_iterator_t it;
        std::ofstream nodefile("5HIV_kmer_and_reads.fasta");
        int num = 0;
        for (it = kmer_map.begin(); it != kmer_map.end(); it++) {
            nodefile << intval_to_kmer(it->first,g_kmer_length) << " : ";
            for (int i = 0; i < it->second.second.size(); ++i) {
                nodefile << it->second.second[i].first << "(" << it-> second.second[i].second << ")" <<"  ";
            }
            nodefile << "" << endl;
            num = num + it->second.second.size();
        }
        nodefile.close();
        return num;

    }

    void get_sorter(std::vector<kmer_int_type_t>& userKmer){
        kmer_sorter_by_count_desc_t sorter(*this);
        sort(userKmer.begin(), userKmer.end(), sorter);
        std::ofstream outfile("6P_userKmers.txt");
        for (int i = 0; i < userKmer.size(); ++i) {
            outfile << intval_to_kmer(userKmer[i],g_kmer_length) << "     " << get_kmer_count(userKmer[i]) << endl;
        }
        outfile.close();

    }

};
#endif