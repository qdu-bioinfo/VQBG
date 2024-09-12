//
// Created by Tian on 2023/7/14.
//
#ifndef ABEGIN_SEQUENCE_GRAPH_H
#define ABEGIN_SEQUENCE_GRAPH_H
#include "utility.h"
#include "common.h"
#include "kmer_hash.h"
#include <string>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <functional>
#include <unordered_map>
#include <unordered_set>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/set.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>


class Sequence_graph{

private:
    typedef int node_idx_t;

    struct VectorHash {
        template <typename T>
        std::size_t operator()(const std::vector<T>& vec) const {
            std::size_t seed = vec.size();
            for (const auto& value : vec) {
                seed ^= std::hash<T>()(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };


    class Node{

        friend class boost::serialization::access;

        template<class Archive>
        void serialize(Archive & ar) {
            ar & id;
            ar & sequence;
            ar & parents;
            ar & children;
            ar & coverage;
            ar & paired_node;
            ar & re_paired_node;
            ar & node_read;
            ar & c_node_id;
            ar & edge_coverage;
            ar & node_layer;
        }

    public:
        Node(): sequence("") {};

        Node(const std::string& mysequence) {
            sequence = mysequence;
        }

        Node(const Node& node) {
            id = node.id;
            sequence = node.sequence;
            children = node.children;
            parents = node.parents;
            coverage = node.coverage;
            paired_node = node.paired_node;
            re_paired_node = node.re_paired_node;
            node_read = node.node_read;
            c_node_id = node.c_node_id;
            edge_coverage = node.edge_coverage;
            node_layer = node.node_layer;
        }
        void set_id(size_t myid) {
            id = myid;
        }
        size_t get_id() {
            return id;
        }

        void set_sequence(const std::string& myseq) {
            sequence = myseq;
        }

        std::string get_sequence() {
            return sequence;
        }

        int get_coverage(){
            return coverage;
        }

        bool add_child(node_idx_t child) {
            if (child < 0)
                return false;
            if (!children.empty()) {
                for (size_t i = 0; i < children.size(); ++i) {
                    if (children[i] == child ) // if exist already
                        return false;
                }
            }
            this->children.push_back(child);
            return true;
        }

        bool add_parent(node_idx_t parent) {
            if (parent < 0)
                return false;
            if (! parents.empty()) {
                for (size_t i = 0; i < parents.size(); ++i) {
                    if (parents[i] == parent )
                        return false;
                }
            }
            this->parents.push_back(parent);
            return true;
        }

        bool is_child(node_idx_t child) {
            if (child < 0)
                return false;
            std::vector<node_idx_t>::iterator it = children.begin();
            for ( ; it != children.end(); ++it ) {
                if (*it == child)
                    return true;
            }
            return false;
        }

        bool is_parent(node_idx_t parent) {
            if (parent < 0)
                return false;
            std::vector<node_idx_t>::iterator it = parents.begin();
            for ( ; it != parents.end(); ++it ) {
                if (*it == parent)
                    return true;
            }
            return false;
        }

        bool delete_child(node_idx_t child) {
            if (child < 0)
                return false;
            std::vector<node_idx_t>::iterator it = children.begin();
            for ( ; it != children.end(); ++it) {
                if (*it == child)
                    break;
            }
            if (it != children.end()) {
                children.erase(it);
                return true;
            } else {
                return false;
            }
        }

        bool delete_parent(node_idx_t parent) {
            if (parent < 0)
                return false;
            std::vector<node_idx_t>::iterator it = parents.begin();
            for ( ; it != parents.end(); ++it) {
                if (*it == parent)
                    break;
            }
            if (it != parents.end()) {
                parents.erase(it);
                return true;
            } else {
                return false;
            }
        }

        void clear_children() {
            children.clear();
        }

        void clear_parents() {
            parents.clear();
        }

        void clear() {
            sequence.clear();
            children.clear();
            parents.clear();
        }

    public:
        size_t id;
        std::string sequence;
        std::vector<node_idx_t> parents;
        std::vector<node_idx_t> children;
        double coverage;
        std::vector<size_t> node_read;//点所属的reads
        std::vector<size_t> c_node_id;
        int node_layer;
        std::vector<pair<size_t,size_t>> edge_coverage;//first与当前点连接的点，second与这个点连接的边的丰度
        std::unordered_map<std::vector<size_t>,pair<size_t,std::vector<size_t>>,VectorHash> paired_node;
        std::unordered_map<std::vector<size_t>,size_t,VectorHash> re_paired_node;

    };

public:

    Sequence_graph(){
        size_ = 0;
    }

    size_t get_size(){
        return size_;
    }

    bool is_used(kmer_int_type_t intval) {
        return (used_kmers_.find(intval) != used_kmers_.end());
    }

    void add_used_kmer(kmer_int_type_t intval, size_t cov) {
        used_kmers_[intval] = cov;
    }

    void delete_used_kmer(kmer_int_type_t intval){
        used_kmers_.erase(intval);
    }


    size_t get_kmer_count(kmer_int_type_t intval) {
        if (is_used(intval))
            return used_kmers_[intval];
        else
            return 0;
    }

    size_t get_kmer_count(const std::string& kmer) {

        kmer_int_type_t intval = kmer_to_intval(kmer);
        return get_kmer_count(intval);
    }



    std::string forward_extend(KmerMap& kmer_map, kmer_int_type_t kmer_val,int& type) {
        //获得向右扩展后的序列
        kmer_int_type_t intval = kmer_val;
        std::string str = intval_to_kmer(intval, g_kmer_length);
        std::vector<kmer_occurence_pair_t> candidates;
        int m = 0;
        while (true) {
            m++;
            int candidate_average = kmer_map.get_forward_candidates(intval, candidates);
            if (candidates.empty()) {
                type = 0;
                break;
            }
            kmer_int_type_t candidate;
            bool flag = false;
//            if (is_used(candidates[0].first) && m!=1){
//                int base_num = candidates[0].first & 3ll;
//                char base = int_to_base(base_num);
//                str += base;
//                break;
//            }
            for (size_t i = 0; i < candidates.size(); ++i) {
                if (!is_used(candidates[i].first)&& candidates[i].second > candidate_average*0.05) {
                    flag = true;
                    candidate = candidates[i].first;
                    break;
                }
            }
            bool flag2 = true;
            if (m > 1){
                for (size_t i = 0; i < candidates.size(); ++i) {
                    if (is_used(candidates[i].first)) {
                        flag2 = false;
                        candidate = candidates[i].first;
                        break;
                    }
                }
            }

            if (!flag2){
                int base_num = candidate & 3ll;
                char base = int_to_base(base_num);
                str += base;
                type = 1;
                break;
            }
            if (!flag){
                type = 2;
                break;
            }
            size_t cov = kmer_map.get_kmer_count(candidate);
            add_used_kmer(candidate, cov);

            int base_num = candidate & 3ll;
            char base = int_to_base(base_num);
            str += base;
            intval = candidate;
        }
        return str;
    }
    std::string forward_extend(KmerMap& kmer_map, kmer_int_type_t kmer_val) {
        //获得向右扩展后的序列
        kmer_int_type_t intval = kmer_val;
        std::string str = intval_to_kmer(intval, g_kmer_length);
        std::vector<kmer_occurence_pair_t> candidates;
        while (true) {
            int candidate_average = kmer_map.get_forward_candidates(intval, candidates);
            if (candidates.empty()) break;
            kmer_int_type_t candidate;
            bool flag = false;
            for (size_t i = 0; i < candidates.size(); ++i) {
                if (!is_used(candidates[i].first)) {
                    flag = true;
                    candidate = candidates[i].first;
                    break;
                }
            }
            if (!flag){
                break;
            }
            size_t cov = kmer_map.get_kmer_count(candidate);
            add_used_kmer(candidate, cov);

            int base_num = candidate & 3ll;
            char base = int_to_base(base_num);
            str += base;
            intval = candidate;
        }
        return str;
    }
    std::string reverse_extend(KmerMap& kmer_map, kmer_int_type_t kmer_val) {
        //获得向左扩展后的序列
        kmer_int_type_t intval = kmer_val;
        std::string str = intval_to_kmer(intval, g_kmer_length);
        std::vector<kmer_occurence_pair_t> candidates;
        while (true) {
            int candidate_coverage = kmer_map.get_reverse_candidates(intval, candidates);
            if (candidates.empty()) break;
            kmer_int_type_t candidate;
            bool flag = false;
            for (size_t i = 0; i < candidates.size(); ++i) {
                //被使用返回true ,没被使用返回false
                if (!is_used(candidates[i].first)) {
                    flag = true;
                    candidate = candidates[i].first;
                    break;
                }
            }

            if (!flag) break;  // all candidates have been used before
            size_t cov = kmer_map.get_kmer_count(candidate);
            add_used_kmer(candidate, cov);

            int base_num = (candidate >> (g_kmer_length*2-2)) & 3ll;
            char base = int_to_base(base_num);
            str = base + str;

            intval = candidate;
        }
        return str;
    }

    string get_trunk(KmerMap& kmer_map, kmer_int_type_t seed ,int& average) {
        add_used_kmer(seed, 0);
        std::string left;
        std::string right;
        int type;
        left = reverse_extend(kmer_map, seed);//向前扩展
        right = forward_extend(kmer_map, seed,type);//向后扩展
        std::string trunk = left + right.substr(g_kmer_length);
        return trunk;
    }

    std::vector<std::set<string>> get_all_nodes(KmerMap& kmer_map,string& trunk,std::map<string,pair<size_t,size_t>>& span_part){
        std::vector<std::set<string>> nodes;
        nodes.resize(trunk.length()-g_kmer_length+1);
        for (int i = 0; i < trunk.length()-g_kmer_length+1; ++i) {
            string kmer = trunk.substr(i,g_kmer_length);
            nodes[i].insert(kmer);
        }
        for (int i = 0; i < nodes.size(); ++i) {
            for (auto it: nodes[i]) {
                std::vector<kmer_occurence_pair_t> candidates;
                kmer_int_type_t intval = kmer_to_intval(it);

                int candidate_average = kmer_map.get_forward_candidates(intval, candidates);
//                cout << candidate_average << " : " <<  candidate_average*0.02 <<endl;
                if (candidates.size() > 1){
                    for (int j = 0; j < candidates.size(); ++j) {
                        if (is_used(candidates[0].first) && j == 0){
                            continue;
                        }
                        if (candidates[j].second > candidate_average*0.05){
//                            cout << i << endl;
                            int type;
                            string extend_str = forward_extend(kmer_map,intval,type);
//                            cout << i << " : " << extend_str <<endl;
                            if (extend_str.length() > g_kmer_length){
                                string final_kmer = extend_str.substr(extend_str.length()-g_kmer_length);

                                int final_k;
                                for (int k = i+1; k < nodes.size(); ++k) {
                                    if (nodes[k].find(final_kmer)!=nodes[k].end()){
                                        final_k = k;
                                        break;
                                    }
                                }
                                if (final_k - i == extend_str.length()-g_kmer_length){
                                    for (int k = 1; k < extend_str.length()-g_kmer_length; ++k) {
                                        string kmer = extend_str.substr(k,g_kmer_length);
                                        nodes[i+k].insert(kmer);
                                    }
//                                    cout << i << " --- " << final_k <<endl;
                                } else{
                                    if (type == 1 && final_k - i < 50 && final_k - i > 0){
//                                        cout << "dddddddddddd" <<endl;
                                        for (int k = 1; k < extend_str.length()-g_kmer_length-1; ++k) {
                                            string kmer = extend_str.substr(k,g_kmer_length);
                                            nodes[i+k].insert(kmer);
                                            pair<size_t,size_t> two;
                                            two.first = i+k;//当前kmer所在的位置
                                            two.second = i+k+1;//当前kmer应该连接的kmer在哪个位置
                                            span_part[kmer] = two;
//                                            cout << " +1 : "<< two.first << "  " << two.second <<endl;
                                        }
                                        pair<size_t,size_t> two;
                                        two.first = i+extend_str.length()-g_kmer_length-1;
                                        two.second = final_k;
                                        span_part[extend_str.substr(extend_str.length()-g_kmer_length-1,g_kmer_length)]=two;
                                        nodes[i+extend_str.length()-g_kmer_length-1].insert(extend_str.substr(extend_str.length()-g_kmer_length-1,g_kmer_length));
//                                        cout << " +1 : "<< two.first << "  " << two.second <<endl;
//                                        cout << i << " : " << extend_str << " | " << extend_str.length() << " | " << final_k << " == " << final_k - i <<endl;
                                    } else if (type == 0 && i >= nodes.size()-g_kmer_length){
//                                        cout << "0000 : " << extend_str << endl;
                                        for (int k = 1; k < extend_str.length()-g_kmer_length+1; ++k) {
                                            string kmer = extend_str.substr(k,g_kmer_length);
                                            nodes[i+k].insert(kmer);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return nodes;
    }

    int countDifferentCharacters(std::string& str1,std::string& str2) {
        int count = 0;
        // Assuming both strings are of equal length
        for (size_t i = 0; i < str1.length(); ++i) {
            if (str1[i] != str2[i]) {
                ++count;
            }
        }

        return count;
    }
    //处理开始节点
    void not_used_candi(KmerMap& kmer_map,std::vector<kmer_int_type_t>& kmer_list,std::vector<std::set<string>>& nodes){
        cout << "Being not used and not candidate..." << endl;
        int m = 0;
        std::vector<string> first_in_kmer;
//        std::ofstream n_u_c("15first_kmer.txt");
        for (int i = 0; i < kmer_list.size(); ++i) {
            if (used_kmers_.find(kmer_list[i])==used_kmers_.end()){
                std::vector<kmer_occurence_pair_t> candidates;
                int candidate_coverage = kmer_map.get_reverse_candidates(kmer_list[i], candidates);
                if (candidates.empty()){
                    string str1 = *nodes[0].begin();
                    string str2 = intval_to_kmer(kmer_list[i],g_kmer_length);
                    int count = countDifferentCharacters(str1,str2);
//&& kmer_map.get_kmer_count(kmer_list[i]) > 1
                    if (count < 4){
                        first_in_kmer.push_back(intval_to_kmer(kmer_list[i],g_kmer_length));
//                        cout << intval_to_kmer(kmer_list[i],g_kmer_length) << "  " << kmer_map.get_kmer_count(kmer_list[i]) << "  " << count<<endl;
//                        n_u_c << intval_to_kmer(kmer_list[i],g_kmer_length) << "  " << kmer_map.get_kmer_count(kmer_list[i]) << "  " << count<<endl;
                    }
                    m++;
                }
            }
        }
        for (int i = 0; i < first_in_kmer.size(); ++i) {
            std::vector<kmer_occurence_pair_t> candidates;
            kmer_int_type_t intval = kmer_to_intval(first_in_kmer[i]);
            int candidate_average = kmer_map.get_forward_candidates(intval, candidates);
            if (candidates.size() > 0) {
                for (int j = 0; j < candidates.size(); ++j) {
                    if (is_used(candidates[0].first) && j == 0) {
                        continue;
                    }
                    if (candidates[j].second > candidate_average * 0.2 && candidates[j].second > 1) {
                        int type;
                        string extend_str = forward_extend(kmer_map, intval, type);
//                        cout << extend_str <<endl;
                        if (extend_str.length() > g_kmer_length) {
                            string final_kmer = extend_str.substr(extend_str.length() - g_kmer_length);
                            int final_k;
                            for (int k = 0; k < nodes.size(); ++k) {
                                if (nodes[k].find(final_kmer) != nodes[k].end()) {
                                    final_k = k;
                                    break;
                                }
                            }
                            if (final_k == extend_str.length()-g_kmer_length) {
                                add_used_kmer(kmer_to_intval(first_in_kmer[i]),kmer_map.get_kmer_count(first_in_kmer[i]));
                                for (int k = 0; k < extend_str.length() - g_kmer_length; ++k) {
                                    string kmer = extend_str.substr(k, g_kmer_length);
                                    nodes[k].insert(kmer);
                                }
//                                cout << final_k << endl;
                            }
                        }
                    }
                }
            }
        }
//        n_u_c.close();
//        cout<< m << endl;
    }
    void delete_error_kmer(KmerMap& kmerMap,std::vector<std::set<string>>& nodes){
        cout << "Begin delete error kmer..." <<endl;
        int error = 0;
        for (int i = 0; i < nodes.size(); ++i) {
            int max = 0;
            for (auto it: nodes[i]){
                if (kmerMap.get_kmer_count(it) > max){
                    max = kmerMap.get_kmer_count(it);
                }
            }
            for (auto it = nodes[i].begin(); it != nodes[i].end(); ) {
                if (kmerMap.get_kmer_count(*it) < max * 0.01) {
                    it = nodes[i].erase(it);  // 删除元素并返回下一个元素的迭代器
//                    cout << *it << endl;
                    error++;
                } else {
                    ++it;  // 移动到下一个元素
                }
            }
        }
        cout << "error_num : " << error <<endl;
    }

    node_idx_t add_node(Node& node) {
        node.id = size_;
        node_set_.push_back(node);
        return (size_++);
    }
    //将点连接起来
    vector<std::map<string,node_idx_t>> get_bubble_graph(std::vector<std::set<string>>& nodes,std::map<string,pair<size_t,size_t>>& span_part){//将点连成图
        vector<std::map<string,node_idx_t>> node_id_position;
        node_id_position.resize(nodes.size());
        vector<vector<int>> two;
        vector<int> one;
        for (int i = 0; i < nodes.size(); ++i) {
            std::map<string,node_idx_t> str_id;
            one.clear();
            for (auto it: nodes[i]) {
                Node node;
                node.sequence = it;
                node.node_layer = i;
                node_idx_t p = add_node(node);
                one.push_back(p);
                str_id[it] = p;
            }
            two.push_back(one);
            node_id_position[i] = str_id;
        }
        for (int i = 0; i < two.size()-1; ++i) {
            one = two[i];
            for (int j = 0; j < one.size(); ++j) {
                for (int k = 0; k < two[i+1].size(); ++k) {
                    if (node_set_[one[j]].sequence.substr(1) == node_set_[two[i+1][k]].sequence.substr(0,g_kmer_length-1)){
                        node_set_[one[j]].add_child(two[i+1][k]);
//                        node_set_[two[i+1][k]].add_parent(one[j]);
                    }
                }
            }
        }
        for (int i = 0; i < two.size()-1; ++i) {
            for (int j = 0; j < two[i].size(); ++j) {
                if (node_set_[two[i][j]].children.empty()){
                    for (auto it2 : nodes[span_part[node_set_[two[i][j]].sequence].second]) {
                        if (node_set_[two[i][j]].sequence.substr(1)==it2.substr(0,g_kmer_length-1)){
                            for (int k = 0; k < two[span_part[node_set_[two[i][j]].sequence].second].size(); ++k) {
                                if (node_set_[two[span_part[node_set_[two[i][j]].sequence].second][k]].sequence == it2){
                                    node_set_[two[i][j]].add_child(two[span_part[node_set_[two[i][j]].sequence].second][k]);
//                                    node_set_[node_id_position[span_part[it].second][it2]].add_parent(node_id_position[i][it]);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        return node_id_position;
    }

    void set_parents(std::vector<Node>& nodes) {
        cout << "Begin add parents..." <<endl;
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (!nodes[i].parents.empty())
                nodes[i].parents.clear();
        }

        for (size_t i = 0; i < nodes.size(); ++i) {
            std::vector<node_idx_t>::const_iterator it;
            for (it = nodes[i].children.begin(); it != nodes[i].children.end(); ++it) {
                nodes[*it].add_parent(i);
            }
        }
    }
    //合并某些点
    std::vector<Node> unique_connect(std::vector<Node>& sim,vector<std::map<string,node_idx_t>>& node_id_position,std::vector<std::set<string>>& nodes,
                                     std::map<size_t,std::set<size_t>>& layer_nodes){//合并单进单出的节点，且同层其他点没有多分支
        cout << "Begin merge single points ..." <<endl;
        std::vector<Node> sim2;
        std::vector<Node> initial = sim;
        std::vector<size_t> delete_node;
        for (int i = 0; i < sim.size(); ++i) {
//             && nodes[node_position[i]].size() == 1
//            cout << i <<endl;
            int m = 0;
            while (sim[i].children.size() == 1 && sim[sim[i].children[0]].parents.size() == 1){
//                nodes[node_position[i]].size() > 1
                bool flag = true;
                if (m == 0){
                    if (nodes[sim[sim[i].children[0]].node_layer].size() > 1){
                        for (auto it1 : nodes[sim[i].node_layer]){
                            if (it1!=sim[i].sequence){
                                if (initial[node_id_position[sim[i].node_layer][it1]].children.size() > 0){
                                    if (initial[node_id_position[sim[i].node_layer][it1]].children.size() > 1 || initial[initial[node_id_position[sim[i].node_layer][it1]].children[0]].parents.size() > 1){
                                        flag = false;
                                    }
                                }

                            }
                        }
                    }
                    for (auto it2 : nodes[sim[sim[i].children[0]].node_layer]) {

                        if (initial[node_id_position[sim[sim[i].children[0]].node_layer][it2]].children.size() > 0) {
                            if (initial[node_id_position[sim[sim[i].children[0]].node_layer][it2]].children.size() > 1 || initial[initial[node_id_position[sim[sim[i].children[0]].node_layer][it2]].children[0]].parents.size() > 1) {
                                flag = false;
                                break;
                            }
                        }

                    }
                } else{

                    for (auto it2 : nodes[sim[sim[i].children[0]].node_layer]) {

                        if (initial[node_id_position[sim[sim[i].children[0]].node_layer][it2]].children.size() > 0){
                            if (initial[node_id_position[sim[sim[i].children[0]].node_layer][it2]].children.size() > 1 || initial[initial[node_id_position[sim[sim[i].children[0]].node_layer][it2]].children[0]].parents.size() > 1){
                                flag = false;
                                break;
                            }
                        }

                    }
                }
                if (!flag){
                    break;
                }
                m++;
                size_t child_idx = sim[i].children[0];
                sim[i].sequence = sim[i].sequence + sim[sim[i].children[0]].sequence.substr(g_kmer_length-1);
                for (int j = 0; j < sim[sim[i].children[0]].children.size(); ++j) {
                    for (int k = 0; k < sim[sim[sim[i].children[0]].children[j]].parents.size(); ++k) {
                        if (sim[sim[sim[i].children[0]].children[j]].parents[k] == sim[i].children[0]){
                            sim[sim[sim[i].children[0]].children[j]].parents[k] = i;
                        }
                    }
                }
                sim[i].children = sim[sim[i].children[0]].children;
//                for (int j = 0; j < sim[i].children.size(); ++j) {
//                    cout << i << " : " << sim[i].children[j] << endl;
//                }
                sim[child_idx].children.clear();
                sim[child_idx].parents.clear();
                delete_node.push_back(child_idx);

            }
        }
        int new_count = 0;
        for (int i = 0; i < sim.size(); ++i) {
            if (sim[i].children.size()!=0 || sim[i].parents.size()!=0){
                layer_nodes[sim[i].node_layer].insert(new_count);
                sim2.push_back(sim[i]);
                new_count++;
            }
        }
        cout << sim2.size() <<endl;
        return sim2;
    }
    void re_id(std::vector<Node>& graph){//重置点的id
        for (node_idx_t i = 0; i < graph.size(); ++i) {
            for (node_idx_t j = 0; j < graph[i].children.size(); ++j) {
                for (node_idx_t k = 0; k < graph.size(); ++k) {
                    if (graph[k].id == graph[i].children[j] && k != graph[k].id){
                        graph[i].children[j] = k;
                    }
                }
            }
        }
        for (int i = 0; i < graph.size(); ++i) {
            graph[i].id = i;
        }
    }

    //获取序列中kmer丰度的和
    size_t str_sum(KmerMap& kmerMap,std::string& sequence){
        int sum = 0;
        for (size_t j = 0; j < sequence.length()-g_kmer_length+1; j++) {
            std::string kmer = sequence.substr(j, g_kmer_length);
            //获得kmer的二进制数字形式
//            kmer_int_type_t kmer_val = kmer_to_intval(kmer, g_kmer_length);

            int cov = kmerMap.get_kmer_count(kmer);
            sum = sum + cov;
        }

        size_t average = sum / (sequence.length() - g_kmer_length + 1);

        return average;
    }
    //根据kmer丰度获取点的丰度
    void add_kmer_mean_value(KmerMap& kmerMap,std::vector<Node>& nodes){
        cout << "Begin loading add coverage..." << endl;
        for (int i = 0; i < nodes.size(); ++i) {
            if (nodes[i].sequence.length() >= g_kmer_length){
                int cov = str_sum(kmerMap,nodes[i].sequence);
                nodes[i].coverage = cov;
            }
        }
    }
    void correction_ratio(std::vector<Node>& nodes, std::map<size_t,std::set<size_t>>& layer_nodes,int& i){//确保每层点的丰度之和为1
        double sum;
        for (auto it: layer_nodes[nodes[i].node_layer]) {
            sum = sum + nodes[it].coverage;
        }
        if (sum < 1){
            double min_ra = 1;
            int min_i;
            for (auto it: layer_nodes[nodes[i].node_layer]) {
                if (min_ra > nodes[it].coverage){
                    min_ra = nodes[it].coverage;
                    min_i = it;
                }
            }
            nodes[min_i].coverage =  nodes[min_i].coverage - sum + 1;
        } else{
            double max_ra = 0;
            int max_i;
            for (auto it: layer_nodes[nodes[i].node_layer]) {
                if (max_ra < nodes[it].coverage){
                    max_ra = nodes[it].coverage;
                    max_i = it;
                }
            }
            nodes[max_i].coverage =  nodes[max_i].coverage - sum + 1;
        }
    }
    //点的丰度归一化
    void get_coverage_normalization(std::vector<Node>& nodes,std::map<size_t,std::set<size_t>>& layer_nodes){
        for (int i = 0; i < nodes.size(); ++i) {
            int sum = 0;
            if (nodes[i].coverage > 1){
                for (auto it: layer_nodes[nodes[i].node_layer]) {
                    sum = sum + nodes[it].coverage;
                }
                for (auto it: layer_nodes[nodes[i].node_layer]) {
                    nodes[it].coverage = nodes[it].coverage / sum;
                    nodes[it].coverage = std::round(nodes[it].coverage * 100) / 100.0;
                }
                //检查分配比例之和是否为1
                correction_ratio(nodes,layer_nodes,i);
            }
        }
    }
    std::vector<Node> topologicalSort(std::vector<Node>& graph) {
        std::vector<Node> result;
        size_t id_counter = 0;
        std::vector<int> inDegree(graph.size(), 0);
        // 计算每个节点的入度
        for (const Node& node : graph) {
            for (node_idx_t child : node.children) {
                inDegree[child]++;
            }
        }
        std::queue<node_idx_t> q;
        // 将所有入度为0的节点加入队列
        for (int i = 0; i < graph.size() ; i++) {
            if (inDegree[i] == 0) {
                q.push(i);
            }
        }
        while (!q.empty()) {
            node_idx_t current = q.front();//访问队列中第一个元素
//            cout << current << endl;
            q.pop();//队列中第一个元素出队列
            graph[current].set_id(id_counter);
            result.push_back(graph[current]);
            result.back().set_id(id_counter);
            id_counter++;
            // 遍历当前节点的所有子节点
            for (node_idx_t child : graph[current].children) {
                inDegree[child]--;
                if (inDegree[child] == 0) {
                    q.push(child);
                }
            }
        }
        // 检查是否存在环
        if (result.size() != graph.size()) {
            std::cout << "There is a cycle in the graph." << std::endl;
            // 返回一个空的结果列表表示拓扑排序失败
            return std::vector<Node>();
        }

        return result;
    }
    void all_left_paired_info(KmerMap& kmerMap,std::vector<std::string>& data,std::vector<Node>& after_gra,std::vector<Node>& result) {
        cout << "Begin paired information..." << endl;
        for (int i = 0; i < result.size(); ++i) {
            for (int j = 0; j < result[i].sequence.length() - g_kmer_length + 1; ++j) {
                string kmer = result[i].sequence.substr(j, g_kmer_length);
                std::vector <pair<size_t, size_t>> reads = kmerMap.get_kmer_read(kmer);
                for (int k = 0; k < reads.size(); ++k) {
//                    if (read_and_node.find(reads[j].first) == read_and_node.end()) {
//                        continue;
//                    }
                    string read = data[reads[k].first];
                    if (read.length() < result[i].sequence.length()){
                        break;
                    }

                    if (read.find(result[i].sequence) == std::string::npos){
                        continue;
                    }

                    bool flag = true;
                    Node node = result[i];
                    result[i].node_read.push_back(reads[k].first);
                    while (flag) {
                        if (node.children.size() == 0)
                            break;
                        for (int l = 0; l < node.children.size(); l++) {
                            if(read.find(after_gra[node.children[l]].sequence) != std::string::npos){
                                flag = true;
                                read_and_node[reads[k].first].push_back(after_gra[node.children[l]].id);
                                node = result[after_gra[node.children[l]].id];
                                break;
                            }
                            flag = false;
                        }
                    }
                }
            }
        }
    }
    void get_paired_info(KmerMap& kmerMap,std::vector<std::string>& data,std::vector<Node>& after_gra,std::vector<Node>& result){
        cout << "Begin paired information..." << endl;
        for (int i = 0; i < result.size(); ++i) {
            for (int j = 0; j < result[i].sequence.length()-g_kmer_length+1; ++j) {
                string kmer = result[i].sequence.substr(j,g_kmer_length);
                std::vector<pair<size_t,size_t>> reads = kmerMap.get_kmer_read(kmer);
                for (int k = 0; k < reads.size(); ++k) {
                    if (read_and_node.find(reads[j].first) != read_and_node.end()){
//                        cout << "llllllllll" <<endl;
                        continue;
                    }
                    string read = data[reads[k].first].substr(reads[k].second);
                    string node_seq = result[i].sequence.substr(j);
                    if (node_seq.length() > read.length())
                        break;
                    string read_seq = read.substr(0,node_seq.length());

                    if (node_seq == read_seq){//判断当前首个点是否能和reads比对上
                        result[i].node_read.push_back(reads[k].first);
                        read = read.substr(node_seq.length());
                    } else{
                        continue;
                    }
                    Node node = result[i];
                    bool flag = true;
                    while (flag){
                        if (node.children.size() == 0)
                            break;
                        for (int l = 0; l < node.children.size(); l++) {
                            string node_str = after_gra[node.children[l]].sequence.substr(g_kmer_length-1);
                            if (node_str.length() > read.length()){
//                                if (node_str.find(read)!= std::string::npos){
//                                    read_and_node[reads[k].first].push_back(result[node.children[l]].id);
//                                }
                                flag = false;
                                continue;
                            } else if(node_str.length() == read.length()){
                                if (read == node_str){
                                    read_and_node[reads[k].first].push_back(after_gra[node.children[l]].id);
                                    flag = false;
                                    break;
                                }
                                flag = false;
                                continue;
                            }
                            string read_str = read.substr(0,node_str.length());
                            if(read_str == node_str){
                                flag = true;
                                read = read.substr(node_str.length());
                                read_and_node[reads[k].first].push_back(after_gra[node.children[l]].id);
//                                if (i == 3){
//                                    cout << i << " : " << node.children[l] <<endl;
//                                }
                                node = result[after_gra[node.children[l]].id];
                                break;
                            }
                            flag = false;
                        }
                    }
                }
            }
        }
    }
    void add_paired_node(std::vector<Node>& result,std::vector<std::string>& data){
        cout << "Begin add paired node..." << endl;
        for (int i = 0; i < result.size(); ++i) {
            for (int j = 0; j < result[i].node_read.size(); ++j) {
                auto it1 = read_and_node.find(result[i].node_read[j]);
                if(it1 != read_and_node.end()){
                    std::vector<size_t>& paired_reads = it1->second;
                    result[i].paired_node[paired_reads].first++;
                    result[i].paired_node[paired_reads].second.push_back(result[i].node_read[j]);
                }
                auto it2 = read_and_node.find(result[i].node_read[j]+(data.size()/2));
                if(it2 != read_and_node.end()){
                    std::vector<size_t>& re_paired_reads = it2->second;
                    result[i].re_paired_node[re_paired_reads]++;
                }

            }
        }
    }
    //判断其中的包含关系-互相包含，a包含b或b包含a都可以-比对两个vector中的元素是否是包含关系
    int paired_nodes_include(std::vector<size_t> nodes1,std::vector<size_t> nodes2){
        int num , flag,i;
        num = (nodes1.size() > nodes2.size()) ? nodes2.size() : nodes1.size();
        for (i = 0; i < num; ++i) {
            if (nodes1[i] != nodes2[i])
                break;
        }
        if (i == num)
            flag = 1;
        else
            flag = 0;
        return flag;
    }
    //将被包含的配对组合与包含其的组合合并
    std::vector<std::vector<std::vector<size_t>>> merge_duplicate_points(std::vector<Node>& topology_node){
        cout << "Begin merge duplicate points..." << endl;
        std::vector<std::vector<std::vector<size_t>>> keys_to_remove;
        keys_to_remove.resize(100000);
        for (int i = 0; i < topology_node.size(); ++i) {
            for (auto it1 = topology_node[i].paired_node.begin(); it1 != topology_node[i].paired_node.end(); ++it1) {//遍历所有组合两两之间进行比对，看是否两种组合之间存在包含关系
                for (auto it2 = std::next(it1); it2 != topology_node[i].paired_node.end(); ++it2) {
                    if (it1->first.size() != it2->first.size() && it1->first[0] == it2->first[0]){
                        if (paired_nodes_include(it1->first,it2->first) == 1){
                            if (it1->first.size() > it2->first.size()){
                                keys_to_remove[i].push_back(it2->first);
                            } else if (it1->first.size() < it2->first.size()){
                                keys_to_remove[i].push_back(it1->first);
                            }
                        }
                    }
                }
            }
        }
        return keys_to_remove;
    }
    std::vector<std::vector<std::vector<size_t>>> re_merge_duplicate_points(std::vector<Node>& topology_node){
        cout << "Begin re merge duplicate points..." << endl;
        std::vector<std::vector<std::vector<size_t>>> keys_to_remove;
        keys_to_remove.resize(100000);
        for (int i = 0; i < topology_node.size(); ++i) {
            for (auto it1 = topology_node[i].re_paired_node.begin(); it1 != topology_node[i].re_paired_node.end(); ++it1) {//遍历所有组合两两之间进行比对，看是否两种组合之间存在包含关系
                for (auto it2 = std::next(it1); it2 != topology_node[i].re_paired_node.end(); ++it2) {
                    if (it1->first.size() != it2->first.size() && it1->first[0] == it2->first[0]){
                        if (paired_nodes_include(it1->first,it2->first) == 1){
                            if (it1->first.size() > it2->first.size()){//计算合并的丰度
                                topology_node[i].re_paired_node[it1->first] =
                                        topology_node[i].re_paired_node[it1->first]+topology_node[i].re_paired_node[it2->first];
                                keys_to_remove[i].push_back(it2->first);
                            } else if (it1->first.size() < it2->first.size()){
                                topology_node[i].re_paired_node[it2->first] =
                                        topology_node[i].re_paired_node[it2->first]+topology_node[i].re_paired_node[it1->first];
                                keys_to_remove[i].push_back(it1->first);
                            }
                        }
                    }
                }
            }
        }
        return keys_to_remove;
    }
    void delete_paired_nodes(std::vector<Node>& topology_node,std::vector<std::vector<std::vector<size_t>>>& keys_to_remove){
        for (int i = 0; i < keys_to_remove.size(); ++i) {
            std::vector<std::vector<size_t>> keys = keys_to_remove[i];
            for (int j = 0; j < keys.size(); ++j) {
                std::vector<size_t> key = keys[j];
                auto it = topology_node[i].paired_node.find(key);
                if (it != topology_node[i].paired_node.end()) {
                    topology_node[i].paired_node.erase(it);
                }
            }

        }
    }
    void re_delete_paired_nodes(std::vector<Node>& topology_node,std::vector<std::vector<std::vector<size_t>>>& keys_to_remove){
        for (int i = 0; i < keys_to_remove.size(); ++i) {
            std::vector<std::vector<size_t>> keys = keys_to_remove[i];
            for (int j = 0; j < keys.size(); ++j) {
                std::vector<size_t> key = keys[j];
                auto it = topology_node[i].re_paired_node.find(key);
                if (it != topology_node[i].re_paired_node.end()) {
                    topology_node[i].re_paired_node.erase(it);
                }
            }

        }
    }

    //双末端相同位置的点
    void get_paired_position(KmerMap& kmerMap,std::vector<std::string>& data,std::vector<Node>& result,std::map<size_t,size_t>& read_first_node,std::map<size_t,vector<size_t>>& node_pair){
        for (int i = 0; i < result.size(); ++i) {
            string str;
            str = result[i].sequence;
            for (int k = 0; k < str.length() - g_kmer_length + 1; ++k) {
                string kmer = str.substr(k, g_kmer_length);
                std::vector<pair<size_t, size_t>> reads = kmerMap.get_kmer_read(kmer);
                for (int l = 0; l < reads.size(); ++l) {
                    read_first_node[reads[l].first] = i;
                    node_pair[i].push_back(reads[l].first);
                }
            }
        }
    }

    void adjust_layer(std::vector<Node>& after_gra, std::map<size_t,std::set<size_t>>& layer_nodes) {
        for (auto& node : layer_nodes) {
            std::set<size_t>& node_set = node.second; // 获取对应节点的 set 引用
            std::set<size_t> updated_set; // 用于存储更新后的 set

            for (auto it = node_set.begin(); it != node_set.end(); ++it) {
                size_t updated_value = after_gra[*it].id; // 获取更新后的值
                updated_set.insert(updated_value); // 将更新后的值插入到新的 set 中
            }

            node_set = updated_set; // 将更新后的 set 赋值回原始的 set 中
        }
    }

    void dfs(std::vector<Node>& nodeList, std::vector<Node>& after_gra, size_t currentNode, std::vector<size_t>& path,
             std::unordered_set<size_t>& visited, std::unordered_map<std::vector<size_t>, size_t, VectorHash>& allPaths,
             int& maxLength) {
        size_t currentLength = 0;
        for (size_t node : path) {
            currentLength += nodeList[node].sequence.length()-g_kmer_length+1;
        }
        if (currentLength >= maxLength) {
            allPaths[path]++;
            return;
        }
        visited.insert(currentNode);
        path.push_back(currentNode);
        bool hasChildren = false;
        for (size_t neighbor : nodeList[currentNode].children) {
            if (visited.find(after_gra[neighbor].id) == visited.end()) {
                hasChildren = true;
                dfs(nodeList, after_gra,after_gra[neighbor].id,path,visited,allPaths,maxLength);
            }
        }
        if (!hasChildren) {
            allPaths[path]++;
        }
        visited.erase(currentNode);
        path.pop_back();
    }
    std::unordered_map<std::vector<size_t>,size_t,VectorHash> findAllPathsWith50Nodes(std::vector<Node>& nodeList,std::vector<Node>& after_gra,
                                                                                      size_t& startNode,int& maxLength,int& m) {
        std::vector<size_t> path;
        std::unordered_set<size_t> visited;
        std::unordered_map<std::vector<size_t>, size_t, VectorHash> allPaths;

        if (m == 0){
            while (maxLength > 0) {
                dfs(nodeList, after_gra, startNode, path, visited, allPaths, maxLength);
                if (allPaths.size() <= 20) {
                    break;
                }
                allPaths.clear();
                maxLength -= 10;
            }
        } else{//确保同一层开始的节点在同一层结束
            dfs(nodeList, after_gra, startNode, path, visited, allPaths, maxLength);
        }

        return allPaths;
    }

    vector<pair<size_t,size_t>> points_decision(std::vector<node_idx_t>& v1,std::vector<node_idx_t>& v2,std::vector<Node>& after_gra){
        vector<pair<size_t,size_t>> same_path_nodes;
        double cha;
        for (int i = 0; i < v1.size(); ++i) {
            pair<size_t,size_t> same_path_node;
            for (int j = 0; j < v2.size(); ++j) {
                if (after_gra[v1[i]].coverage > after_gra[v2[j]].coverage){
                    cha = after_gra[v1[i]].coverage - after_gra[v2[j]].coverage;
                } else if (after_gra[v1[i]].coverage < after_gra[v2[j]].coverage){
                    cha = after_gra[v2[j]].coverage - after_gra[v1[i]].coverage;
                } else{
                    cha = 0;
                }
                if (cha < 0.015){
//                    cout << after_gra[v1[i]].id << " ~~ " << after_gra[v2[j]].id << " cc " << cha <<endl;
                    same_path_node.first = after_gra[v1[i]].id;
                    same_path_node.second = after_gra[v2[j]].id;
//                    cout << same_path_node.first << " -- " << same_path_node.second <<endl;
                    break;
                }
            }
            if (same_path_node.first!=0 && same_path_node.second!=0){
                same_path_nodes.push_back(same_path_node);
            }
        }
        return same_path_nodes;
    }

    //删除错误的小路径
    void retain_simplified_points(KmerMap& kmerMap,std::vector<Node>& nodes,std::vector<Node>& after_gra,std::unordered_map<std::vector<size_t>,size_t,VectorHash>& allPaths,
                                  std::map<size_t,std::set<size_t>>& layer_nodes,std::vector<std::string>& data){
//        cout << "Begin delete small path ..." <<endl;

        vector<vector<size_t>> retain_points;
        map<pair<size_t,size_t>,size_t> same_roads;
        for (auto path : allPaths) {
            std::vector<size_t> key = path.first;
            vector<pair<size_t,size_t>> sam_road;
            for (int i = 0; i < key.size()-1; ++i) {
                if (nodes[key[i]].children.size() > 1){
                    for (int j = i; j >= 1 ; --j) {
                        if (nodes[key[j]].parents.size() == nodes[key[i]].children.size()) {
//                            if (layer_nodes[after_gra[nodes[key[j]].parents[0]].node_layer].size() == layer_nodes[after_gra[nodes[key[i]].children[0]].node_layer].size()) {//考虑他们所在的层的点数是否相等
//                                cout << key[i] << " : " << key[j] <<endl;
                                sam_road = points_decision(nodes[key[j]].parents, nodes[key[i]].children, after_gra);
                            if (sam_road.size()!=0){
                                break;
                            }

//                            }
                        }
                    }
//                    if (flag){
//                        break;
//                    }
                    for (int i = 0; i < sam_road.size(); ++i) {
//                cout << sam_road[i].first << " -- " << sam_road[i].second <<endl;
                        same_roads[sam_road[i]]++;
                    }
                }
            }
        }
        vector<pair<size_t,size_t>> keep_one;
        if (same_roads.size()!=0){
            // 嵌套循环遍历两次
            for (auto entry1 = same_roads.begin(); entry1!= same_roads.end();entry1++) {
                for (auto entry2 = std::next(entry1); entry2!= same_roads.end();entry2++) {
                    double one_cha;
                    if (entry1->first.first == entry2->first.first){
                        if (abs(nodes[entry1->first.second].coverage -nodes[entry1->first.first].coverage) > abs(nodes[entry2->first.second].coverage-nodes[entry1->first.first].coverage)){
                            keep_one.push_back(entry1->first);
                        } else{
                            keep_one.push_back(entry2->first);
                        }
                    }
                }
            }
            for (int i = 0; i < keep_one.size(); ++i) {
                same_roads.erase(keep_one[i]);
            }
        }
        unordered_map<std::vector<size_t>,size_t,VectorHash> del_points;
        for (auto& path : allPaths) {
            vector<size_t> key_d = path.first;
            string str_d = nodes[key_d[0]].sequence;
            for (int idd = 1; idd < key_d.size(); ++idd) {//获取小路径的序列
                str_d = str_d + nodes[key_d[idd]].sequence.substr(g_kmer_length-1);
            }
            size_t support_num = 0;
            string kmer = str_d.substr(0,g_kmer_length);
            std::vector<pair<size_t, size_t>> reads = kmerMap.get_kmer_read(kmer);
            for (int r = 0; r < reads.size(); ++r) {
                string read = data[reads[r].first].substr(reads[r].second);
                if (read.length() > str_d.length()){
                    if (read.substr(0,str_d.length()) == str_d){
                        support_num++;
                    }
                }
            }
/*            for (int dt = 0; dt < data.size(); ++dt) {
                if (data[dt].find(str_d)!=std::string::npos){
                    support_num++;
                }
                if (support_num > 7){
                    break;
                }
            }*/
            path.second = support_num;
            if (nodes[key_d[key_d.size()-1]].children.empty() || nodes[key_d[0]].parents.empty()){
                if (support_num == 0){//对于头和尾的小路径，删除reads支持数为0的
                    del_points[key_d]++;
                }
            } else{
                if (support_num <= 6){//对于其他中间路径，删除reads支持数小于等于6的
                    del_points[key_d]++;
                }
            }
        }

        if (del_points.size() == allPaths.size()){
            del_points.clear();
            for (auto it : same_roads) {
                for (auto path : allPaths) {//检查每条路径中是否同时出现两个在一条路上的点
                    std::vector <size_t> key = path.first;
                    for (int l = 0; l < key.size(); ++l) {
                        if (it.first.first == key[l]) {
                            bool n_same = false;
                            for (int m = l + 1; m < key.size(); ++m) {
                                if (it.first.second == key[m]) {
                                    n_same = true;
                                    break;
                                }
                            }
                            if (!n_same){
                                del_points[key]++;
                            }
                            break;
                        }
                    }
                    for (int l = 1; l < key.size(); ++l) {
                        if (it.first.second == key[l]) {
                            bool n_same2 = false;
                                for (int m = l-1; m >=0 ; --m) {
                                    if (it.first.first == key[m]) {
                                        n_same2 = true;
                                        break;
                                    }
                                }

                            if (!n_same2){
                                del_points[key]++;
                            }
                            break;
                        }
                    }
                }
            }
        }

        if (del_points.size() < allPaths.size()){
            for (auto itd : del_points) {
                vector<size_t> del_vector = itd.first;
                auto it = allPaths.find(del_vector);
                if (it != allPaths.end()) {
                    allPaths.erase(it);
                } else {
                    std::cout << "Key not found in map." << std::endl;
                }
            }
        }
//        cout << "retain_points : " << retain_points.size() <<endl;
//        return retain_points;
    }

    vector<vector<size_t>> compression_nodes(KmerMap& kmerMap,std::vector<Node>& nodeList, std::vector<Node>& after_gra,std::map<size_t,std::set<size_t>>& layer_nodes,
                                             vector<vector<size_t>>& sim_layer,vector<size_t>& layer_sim_node,std::vector<std::string>& data,int& limit_length){
        cout << "Begin compression nodes..." <<endl;
        vector<size_t> start;
        vector<vector<size_t>> simplified_points;
        for (int i = 0; i < nodeList.size(); ++i) {
            if (nodeList[i].parents.empty()){
                start.push_back(i);
            }
        }
        int m =0;

        int before_sum = 0;

        size_t sim_m = 0;

        while(1){

            m++;
        set<size_t> new_start;
        bool stop_flag = false;
        vector<size_t> sim_layer_node;
        int maxLength = limit_length*0.4;
        bool layer_end_flag = true;
        int min_layer = nodeList[start[0]].node_layer;
        for (int i = 1; i < start.size(); ++i) {
            if (nodeList[start[0]].node_layer != nodeList[start[i]].node_layer){
                if (nodeList[start[i]].node_layer < min_layer){
                    min_layer = nodeList[start[i]].node_layer;
                }
                layer_end_flag = false;
          }
        }
        if(!layer_end_flag){

            for (int sta = 0; sta < start.size(); ++sta) {
                if (nodeList[start[sta]].node_layer != min_layer){
                    start.erase(start.begin()+sta);
                }
            }

        }
        for (int i = 0; i < start.size(); ++i) {

            std::unordered_map<std::vector<size_t>,size_t,VectorHash> allPaths = findAllPathsWith50Nodes(nodeList,after_gra,start[i],maxLength,m);
            before_sum = before_sum + allPaths.size();
            retain_simplified_points(kmerMap,nodeList,after_gra,allPaths,layer_nodes,data);
                for (auto path : allPaths) {
                    std::vector <size_t> key = path.first;

                    vector<size_t> points;
                    string st = "";
                    for (size_t node : key) {
                        points.push_back(node);
                        st = st + nodeList[node].sequence.substr(g_kmer_length-1);

                    }

                    sim_layer_node.push_back(sim_m);
                    layer_sim_node.push_back(m);

                    sim_m++;
                    simplified_points.push_back(points);
                    if (nodeList[key[key.size()-1]].children.empty()){
                        stop_flag = true;
                        continue;
                    }
                    for (int j = 0; j < nodeList[key[key.size()-1]].children.size(); ++j) {
                        new_start.insert(after_gra[nodeList[key[key.size()-1]].children[j]].id);
                    }
                }
//            }

        }
        sim_layer.push_back(sim_layer_node);
        start.clear();
        if (stop_flag){
            break;
        }
        for (auto it : new_start) {//更新start节点信息
            start.push_back(it);
        }
        }
        cout << before_sum << " -- " << simplified_points.size() <<endl;

        return simplified_points;
    }

    bool same_children(size_t& a,size_t& b,vector<Node>& result){
        map<size_t,size_t> c;
        bool flag = false;
        for (int i = 0; i < result[a].children.size(); ++i) {
            c[result[a].children[i]]++;
        }
        for (int i = 0; i < result[b].children.size(); ++i) {
            c[result[b].children[i]]++;
        }
        for (auto it : c){
            if (it.second > 1){
                flag = true;
            }
        }
        return flag;
    }

    void alone_possess_nodes(vector<vector<size_t>>& sim_layer, vector<vector<pair<size_t,size_t>>>& parent_different_node,
                             vector<Node>& sim_graph1,std::vector<Node>& result){//获取压缩图中特有的点，即同一层中这个点只在压缩图中的当前点出现过
        map<size_t,size_t> alone_used;

        for (int i = 0; i < sim_layer.size(); ++i){//遍历每一层的点
            for (int j = 0; j < sim_layer[i].size(); j++) {
                if (alone_used.find(sim_layer[i][j])==alone_used.end()){//记录所有同一层中，最后一个点相同的节点即孩子节点相同的节点
                    vector<size_t> same_end_p;
                    same_end_p.push_back(sim_layer[i][j]);
                    alone_used[sim_layer[i][j]]++;
                    for (int k = 0; k < sim_layer[i].size(); k++) {
                        if (sim_layer[i][j]!=sim_layer[i][k]){
                            if (sim_graph1[sim_layer[i][j]].c_node_id[sim_graph1[sim_layer[i][j]].c_node_id.size()-1] == sim_graph1[sim_layer[i][k]].c_node_id[sim_graph1[sim_layer[i][k]].c_node_id.size()-1]
                            || same_children(sim_graph1[sim_layer[i][j]].c_node_id[sim_graph1[sim_layer[i][j]].c_node_id.size()-1],sim_graph1[sim_layer[i][k]].c_node_id[sim_graph1[sim_layer[i][k]].c_node_id.size()-1],result)){
                                same_end_p.push_back(sim_layer[i][k]);
                                alone_used[sim_layer[i][k]]++;
                            }
                        }
                    }
                    difference_vec(same_end_p,parent_different_node,sim_graph1);//对比两个孩子节点相同的节点的区别
                }
            }
        }

    }

    void difference_vec(vector<size_t>& local_children,vector<vector<pair<size_t,size_t>>>& different_node,vector<Node>& sim_graph1){
        map<size_t,size_t> alone_used;
        for (int i = 0; i < local_children.size(); ++i) {
            for (int k = 0; k < sim_graph1[local_children[i]].c_node_id.size(); ++k) {
                alone_used[sim_graph1[local_children[i]].c_node_id[k]]++;
            }
        }

        for (int i = 0; i < local_children.size(); ++i) {
            bool spe_flag = false;
            for (int k = 0; k < sim_graph1[local_children[i]].c_node_id.size(); ++k) {
                if (alone_used[sim_graph1[local_children[i]].c_node_id[k]] == 1){
                    spe_flag = true;
                    pair<size_t,size_t> spe;
                    spe.first = k;//位置
                    spe.second = sim_graph1[local_children[i]].c_node_id[k];//点的id
                    different_node[local_children[i]].push_back(spe);
                }
            }
        }
    }
    bool areEqualWithOneOrLessErrors(const std::string& s1, const std::string& s2) {
        // 如果长度不同，直接返回false
        if (s1.length() != s2.length()) {
            return false;
        }

        int errorCount = 0;

        // 比较两个字符串的每个字符
        for (size_t i = 0; i < s1.length(); ++i) {
            if (s1[i] != s2[i]) {
                errorCount++;
                if (errorCount > 1) {
                    return false;  // 错误超过一个字符，直接返回false
                }
            }
        }

        // 如果错误字符数小于等于1，则返回true
        return errorCount <= 1;
    }

    bool merge_length_less_hundred(KmerMap& kmerMap,std::vector<Node>& result,std::vector<std::string>& data,vector<size_t>& p,vector<size_t>& c){
        string p_str = result[p[0]].sequence;
        string surplus_str = "";
        for (int i = 1; i < p.size()/2; ++i) {
            p_str += result[p[i]].sequence.substr(g_kmer_length-1);
        }
        for (int i = p.size()/2; i < p.size(); ++i) {
            surplus_str += result[p[i]].sequence.substr(g_kmer_length-1);
        }
        for (int i = 0; i < c.size(); ++i) {
            surplus_str += result[c[i]].sequence.substr(g_kmer_length-1);
        }
        int support = 0;
        for (int i = 0; i < p_str.length()-g_kmer_length+1; ++i) {
            string kmer = p_str.substr(i,g_kmer_length);
            std::vector<pair<size_t, size_t>> reads = kmerMap.get_kmer_read(kmer);
            for (int r = 0; r < reads.size(); ++r) {
                string read = data[reads[r].first].substr(reads[r].second);
                string str1 = p_str.substr(i);
                string str2 = str1 + surplus_str;
                if (str2 == read.substr(0, str2.length())) {
                    support++;
                }
            }
        }
        bool flag;
        cout << "ss" << p[0] << " : " <<support <<endl;
        if (support > 1){
            flag = true;
        } else{
            flag = false;
        }
        return flag;
    }

    void create_sim_graph(KmerMap& kmerMap,std::vector<Node>& result,std::vector<Node>& after_gra,vector<vector<size_t>>& sim_gra,vector<Node>& sim_graph1
                          ,std::map<size_t,std::set<size_t>>& layer_nodes,vector<vector<size_t>>& sim_layer,vector<size_t>& layer_sim_node
                          ,vector<vector<pair<size_t,size_t>>>& parent_different_node,std::vector<std::string>& data){
        cout << "Begin create sim graph..." <<endl;
        sim_graph1.resize(sim_gra.size());

        for (int i = 0; i < sim_gra.size(); ++i) {
            sim_graph1[i].id = i;
            sim_graph1[i].node_layer = layer_sim_node[i];
            sim_graph1[i].c_node_id = sim_gra[i];
            string str = result[sim_gra[i][0]].sequence;
            for (int j = 1; j < sim_gra[i].size(); ++j) {
                str = str + result[sim_gra[i][j]].sequence.substr(g_kmer_length-1);
            }
            sim_graph1[i].sequence = str;
        }

        parent_different_node.resize(sim_graph1.size());
        alone_possess_nodes(sim_layer,parent_different_node,sim_graph1,result);
        vector<vector<pair<size_t,size_t>>> children_different_node;

        //简化图中点的丰度
        for (int i = 0; i < sim_graph1.size(); ++i) {
            double cov = 0;
            int cov_sum = 0;
            double min_cov = 1;
            for (int j = 0; j < sim_graph1[i].c_node_id.size(); ++j) {
                if (result[sim_graph1[i].c_node_id[j]].coverage < min_cov){
                    min_cov = result[sim_graph1[i].c_node_id[j]].coverage;
                }
            }
            for (int j = 0; j < sim_graph1[i].c_node_id.size(); ++j) {
                if (abs(result[sim_graph1[i].c_node_id[j]].coverage-min_cov) < 0.02){
                    cov = cov + result[sim_graph1[i].c_node_id[j]].coverage;
                    cov_sum++;
                }
            }
            sim_graph1[i].coverage = cov/cov_sum;
        }

        for (size_t i = 0; i < sim_gra.size(); ++i) {//给点之间添加关系
//            cout << i <<endl;
            vector<size_t> local_children;
//            }
            double c_sum_c = 0;
            for (size_t j = i + 1; j < sim_gra.size(); ++j) {//确定当前节点所有可能连接的点
                for (int k = 0; k < result[sim_gra[i][sim_gra[i].size() - 1]].children.size(); ++k) {
                    if (after_gra[result[sim_gra[i][sim_gra[i].size() - 1]].children[k]].id == sim_gra[j][0]){
                        local_children.push_back(j);
                        c_sum_c = c_sum_c + sim_graph1[j].coverage;
                        break;
                    }
                }
            }

            if (local_children.size() == 1){
                pair<size_t,size_t> edge;
                edge.first = local_children[0];
                edge.second = 0;
                sim_graph1[i].edge_coverage.push_back(edge);
                sim_graph1[i].add_child(local_children[0]);
                sim_graph1[local_children[0]].add_parent(i);
            } else if(local_children.size() > 1) {
                bool connect = false;
                children_different_node.clear();
                children_different_node.resize(sim_graph1.size());
                difference_vec(local_children, children_different_node, sim_graph1);
                for (int l = 0; l < local_children.size(); ++l) {
                    int support = 0;
                    string str = sim_graph1[i].sequence + sim_graph1[local_children[l]].sequence.substr(g_kmer_length - 1);
                    string kmer = str.substr(0,g_kmer_length);
                    std::vector<pair<size_t, size_t>> reads = kmerMap.get_kmer_read(kmer);
                    for (int r = 0; r < reads.size(); ++r) {
                        string read = data[reads[r].first].substr(reads[r].second);
                        if (read.length() > str.length()){
                            if (read.substr(0,str.length()) == str){
                                support++;
                            }
                        }
                    }
                    if (support > 2) {
                        connect = true;
                        pair <size_t, size_t> edge;
                        edge.first = local_children[l];
                        edge.second = support;
                        sim_graph1[i].edge_coverage.push_back(edge);
                        sim_graph1[i].add_child(local_children[l]);
                        sim_graph1[local_children[l]].add_parent(i);
                    }
                }
                if (!connect) {
                    for (int l = 0; l < local_children.size(); ++l) {
                        int support = 0;
                        string str = sim_graph1[i].sequence + sim_graph1[local_children[l]].sequence.substr(g_kmer_length - 1);
                        for (int d = 0; d < data.size(); ++d) {
                            if (str.length() > data[d].length()) {
                                if (str.find(data[d]) != std::string::npos) {
                                    support++;
                                }
                            } else {
                                if (data[d].find(str) != std::string::npos) {
                                    support++;
                                }
                            }
                        }
                        if (support > 0) {
                            connect = true;
                            pair <size_t, size_t> edge;
                            edge.first = local_children[l];
                            edge.second = support;
                            sim_graph1[i].edge_coverage.push_back(edge);
                            sim_graph1[i].add_child(local_children[l]);
                            sim_graph1[local_children[l]].add_parent(i);
                        }
                    }
                }
                /*if (!connect) {
                    for (int l = 0; l < local_children.size(); ++l) {
                        if (parent_different_node[i].size() > 0 && children_different_node[local_children[l]].size() > 0) {
                            int support = 0;
                            string str1;

                            str1 = result[sim_graph1[local_children[l]].c_node_id[0]].sequence;
                            for (int j = 1; j <= children_different_node[local_children[l]][0].first; ++j) {
                                str1 = str1 + result[sim_graph1[local_children[l]].c_node_id[j]].sequence.substr(
                                        g_kmer_length - 1);
                            }
                            string str2;
                            for (int d = 0; d < parent_different_node[i].size(); ++d) {
                                str2 = "";
                                str2 = result[sim_graph1[i].c_node_id[parent_different_node[i][d].first]].sequence;
                                for (int j = parent_different_node[i][d].first + 1; j < sim_graph1[i].c_node_id.size(); ++j) {
                                    str2 = str2 + result[sim_graph1[i].c_node_id[j]].sequence.substr(g_kmer_length - 1);

                                }
                                if (str2.length() + str1.length() - g_kmer_length + 1 < data[0].length()) {
                                    break;
                                }
                            }
                            string str3 = "";
                            if (parent_different_node[i][parent_different_node[i].size() - 1].first < sim_graph1[i].c_node_id.size() - 1) {//最后一个特殊点开始到最后的序列
                                for (int j = parent_different_node[i][parent_different_node[i].size() - 1].first + 1; j < sim_graph1[i].c_node_id.size(); ++j) {
                                    str3 = result[sim_graph1[i].c_node_id[j]].sequence.substr(g_kmer_length - 1);
                                }
                            }
                            if (str2.length() + sim_graph1[local_children[l]].sequence.substr(g_kmer_length - 1).length() < data[0].length()) {
                                string str6 = str2 + sim_graph1[local_children[l]].sequence.substr(g_kmer_length - 1);
                                for (int u = 0; u < data.size(); ++u) {
                                    if (data[u].find(str6) != std::string::npos) {
                                        support++;
                                    }
                                }
                            } else if (str2.length() + str1.length() - g_kmer_length + 1 < data[0].length()) {
                                string str4 = str2 + str1.substr(g_kmer_length - 1);
                                for (int j = 0; j < str2.length() - str3.length() - 2 * g_kmer_length + 2; ++j) {
                                    string kmer = str2.substr(j, g_kmer_length);
                                    std::vector <pair<size_t, size_t>> reads = kmerMap.get_kmer_read(kmer);
                                    for (int r = 0; r < reads.size(); ++r) {
                                        string read = data[reads[r].first].substr(reads[r].second);
                                        string str5 = str2.substr(j);
                                        if (read.length() > str5.length() && str5 == read.substr(0, str5.length())) {
                                            read = read.substr(str5.length());
                                            if (sim_graph1[local_children[l]].sequence.length() - g_kmer_length + 1 > read.length()) {
                                                if (sim_graph1[local_children[l]].sequence.substr(g_kmer_length - 1,read.length()) == read) {
                                                    support++;
                                                }
                                            } else {
                                                if (sim_graph1[local_children[l]].sequence.substr(g_kmer_length - 1) ==
                                                    read.substr(0, sim_graph1[local_children[l]].sequence.length() -g_kmer_length + 1)) {
                                                    support++;
                                                }
                                            }

                                        }
                                    }
                                }
                            }

                            if (support > 10) {
                                connect = true;
                                pair <size_t, size_t> edge;
                                edge.first = local_children[l];
                                edge.second = support;
                                sim_graph1[i].edge_coverage.push_back(edge);
                                sim_graph1[i].add_child(local_children[l]);
                                sim_graph1[local_children[l]].add_parent(i);
                            }
                        }
                    }

                }*/
                if (!connect) {
                    for (int l = 0; l < local_children.size(); ++l) {
                        pair <size_t, size_t> edge;
                        edge.first = local_children[l];
                        edge.second = 0;
                        sim_graph1[i].edge_coverage.push_back(edge);
                        sim_graph1[i].add_child(local_children[l]);
                        sim_graph1[local_children[l]].add_parent(i);
                    }
                }
            }
        }

//        cout << sim_layer.size() <<endl;
        //给特殊点添加边
        for (int i = 0; i < sim_graph1.size(); ++i) {
            if (sim_graph1[i].children.empty() && !result[sim_graph1[i].c_node_id[sim_graph1[i].c_node_id.size()-1]].children.empty()){
                bool em_child_flag = false;
                vector<int> em_child_i;
                vector<int> em_child_i_pos;
                for (int cn = 0; cn < result[sim_graph1[i].c_node_id[sim_graph1[i].c_node_id.size()-1]].children.size(); ++cn) {
                    for(int j = i+1; j < sim_graph1.size(); j++){//向后寻找没有父节点的节点中是否有可以相连的
                        if (sim_graph1[j].parents.empty()){
                            for (int si = 0; si < sim_graph1[j].c_node_id.size(); ++si) {
                                if (sim_graph1[i].c_node_id[sim_graph1[i].c_node_id.size()-1]
                                == sim_graph1[j].c_node_id[si]){//对于特殊点找到了与之对应的连接处
                                    em_child_i.push_back(j);
                                    em_child_i_pos.push_back(si);
                                    em_child_flag = true;
                                    break;
                                }
                            }
                            if(em_child_flag){
                                break;
                            }
                        }
                    }
                }
                if(em_child_flag){
                    for (int ci = 0; ci < em_child_i.size(); ++ci) {
                        if (front_back_parts_same(sim_graph1[i].c_node_id,sim_graph1[em_child_i[ci]].c_node_id,em_child_i_pos[ci])){
                            sim_graph1[i].c_node_id.erase(sim_graph1[i].c_node_id.begin()+sim_graph1[i].c_node_id.size()-em_child_i_pos[ci]-1
                                    ,sim_graph1[i].c_node_id.end());
                            string spe_str = result[sim_graph1[i].c_node_id[0]].sequence;
                            for (int ct = 1; ct < sim_graph1[i].c_node_id.size(); ++ct) {
                                spe_str += result[sim_graph1[i].c_node_id[ct]].sequence.substr(g_kmer_length-1);
                            }
                            sim_graph1[i].sequence = spe_str;
                            pair <size_t, size_t> edge;
                            edge.first = em_child_i[ci];
                            edge.second = 0;
                            sim_graph1[i].edge_coverage.push_back(edge);
                            sim_graph1[i].add_child(em_child_i[ci]);
                            sim_graph1[em_child_i[ci]].add_parent(i);
                        }
                    }
                }
            }
        }
//        in_degree_zero(sim_layer,sim_graph1,result);//删除没有父节点的，可能是错误的节点
    }

    void in_degree_zero(vector<vector<size_t>>& sim_layer,vector<Node>& sim_graph1,std::vector<Node>& result){//去除无父亲的节点，这些节点是错误的
        for (int i = 5; i < sim_layer.size()-1; ++i) {
            for (int j = 0; j < sim_layer[i].size(); ++j) {
                if (sim_graph1[sim_layer[i][j]].parents.empty() && !sim_graph1[sim_layer[i][j]].children.empty() && !result[sim_graph1[sim_layer[i][j]].c_node_id[0]].parents.empty()){
                    for (int c = 0; c < sim_graph1[sim_layer[i][j]].children.size(); ++c) {
                        for (int k = 0; k < sim_graph1[sim_graph1[sim_layer[i][j]].children[c]].parents.size(); ++k) {
                            if (sim_graph1[sim_graph1[sim_layer[i][j]].children[c]].parents[k] == sim_layer[i][j]){
                                sim_graph1[sim_graph1[sim_layer[i][j]].children[c]].parents.erase(sim_graph1[sim_graph1[sim_layer[i][j]].children[c]].parents.begin() + k);
                            }
                        }

                    }
                    sim_graph1[sim_layer[i][j]].children.clear();
                    sim_graph1[sim_layer[i][j]].edge_coverage.clear();
                }
            }
        }
        for (int i = sim_graph1.size()-1; i >= 0; --i) {
            if (sim_graph1[i].children.empty() && !sim_graph1[i].parents.empty() && !result[sim_graph1[i].c_node_id[sim_graph1[i].c_node_id.size()-1]].children.empty()){
                for (int j = 0; j < sim_graph1[i].parents.size(); ++j) {
                    for (int k = 0; k < sim_graph1[sim_graph1[i].parents[j]].children.size(); ++k) {
                        if (sim_graph1[sim_graph1[i].parents[j]].children[k] == i){
                            sim_graph1[sim_graph1[i].parents[j]].children.erase(sim_graph1[sim_graph1[i].parents[j]].children.begin() + k);
                        }
                    }
                    for (int k = 0; k < sim_graph1[sim_graph1[i].parents[j]].edge_coverage.size(); ++k) {
                        if (sim_graph1[sim_graph1[i].parents[j]].edge_coverage[k].first == i){
                            sim_graph1[sim_graph1[i].parents[j]].edge_coverage.erase(sim_graph1[sim_graph1[i].parents[j]].edge_coverage.begin() + k);
                        }
                    }
                }
                sim_graph1[i].parents.clear();
            }
        }
    }

    bool front_back_parts_same(vector<size_t>& a,vector<size_t>& b,int& position){
        bool same_flag = true;
        int j = 0;
        for (int i = a.size()-position-1; i < a.size(); ++i) {
            if (a[i] != b[j]){
                same_flag = false;
            }
            j++;
        }
        return same_flag;
    }

    bool not_any_unique(vector<Node>& sim_graph1,node_idx_t& i){//i节点的父亲和孩子都不只于他相连
        bool flag = true;

        if(sim_graph1[i].parents.size() == 1){
            flag = false;
        }

        return flag;
    }

    bool consistent_with_abundance(vector<Node>& sim_graph1,size_t& i,node_idx_t& c){//判断根据当前孩子与父节点之间是否可能存在的连接
        bool flag = false;
        if(abs(sim_graph1[i].coverage-sim_graph1[c].coverage) < 0.03){

            flag = true;
        }
        double sum_c = 0;
        for (int j = 0; j < sim_graph1[i].children.size(); ++j) {
            sum_c = sum_c + sim_graph1[sim_graph1[i].children[j]].coverage;
        }
        if (abs(sum_c-sim_graph1[i].coverage) < 0.05){

            flag = true;
        }
        for (int j = 0; j < sim_graph1[i].children.size(); ++j) {
            if (sim_graph1[i].children[j]!=c){
                double two = sim_graph1[sim_graph1[i].children[j]].coverage + sim_graph1[c].coverage;
                if (abs(two-sim_graph1[i].coverage) < 0.03){

                    flag = true;
                    break;
                }
            }
        }
        return flag;
    }

    bool one_node_dif_par(vector<Node>& sim_graph1,size_t& i,node_idx_t& c){
        bool flag = false;
        if (sim_graph1[c].parents.size() > 1){
            for (int k = 0; k < sim_graph1[c].parents.size(); ++k) {
                if (sim_graph1[sim_graph1[c].parents[k]].children.size() > 1) {
                    for (int j = 0; j < sim_graph1[sim_graph1[c].parents[k]].children.size(); ++j) {
                        bool only_flag = false;
                        for (auto pe: sim_graph1[sim_graph1[c].parents[k]].paired_node) {
                            std::vector <size_t> key = pe.first;
//                        cout << i << " : " << key[0] << "  " <<sim_graph1[i].children[j] <<endl;
                            if (key[0] == sim_graph1[sim_graph1[c].parents[k]].children[j]) {
                                only_flag = true;
                            }
                        }
                        if (!only_flag){
                            if (sim_graph1[sim_graph1[c].parents[k]].children[j] == c){
                                flag = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
        return flag;
    }

    void based_on_pairing_information_delete_edge(vector<Node>& sim_graph1){//根据配对信息删除边
        cout << "Being based_on_pairing_information_delete_edge ... " <<endl;
        for (size_t i = 0; i < sim_graph1.size(); ++i) {
            if (sim_graph1[i].children.size() > 1){
                vector<size_t> del;
                for (int j = 0; j < sim_graph1[i].children.size(); ++j) {
                    bool only_flag = false;
                    for (auto pe : sim_graph1[i].paired_node) {
                        std::vector<size_t> key = pe.first;
//                        cout << i << " : " << key[0] << "  " <<sim_graph1[i].children[j] <<endl;
                        if (key[0] == sim_graph1[i].children[j]){
                            only_flag = true;
                        }
                    }
                    if (!only_flag){
//                        cout << i << " : " << sim_graph1[i].children[j] << "---"<<endl;

                        if (not_any_unique(sim_graph1,sim_graph1[i].children[j]) && !consistent_with_abundance(sim_graph1,i,sim_graph1[i].children[j])){//此点不是某个点唯一相连的点
                            del.push_back(sim_graph1[i].children[j]);
//                            cout << i << " : " << sim_graph1[i].children[j] <<endl;
                        }
                    }
                }
                if (del.size() < sim_graph1[i].children.size()){
                    for (int j = 0; j < del.size(); ++j) {
//                        cout << i << " : " << del[j] <<endl;
                        for (int k = 0; k < sim_graph1[i].children.size(); ++k) {
                            if(sim_graph1[i].children[k] == del[j]){
                                sim_graph1[i].children.erase(sim_graph1[i].children.begin()+k);
                            }
                        }
                        for (int k = 0; k < sim_graph1[del[j]].parents.size(); ++k) {
                            if (sim_graph1[del[j]].parents[k] == i){
                                sim_graph1[del[j]].parents.erase(sim_graph1[del[j]].parents.begin()+k);
                            }
                        }
                        for (int k = 0; k < sim_graph1[i].edge_coverage.size(); ++k) {
                            if(sim_graph1[i].edge_coverage[k].first == del[j]){
                                sim_graph1[i].edge_coverage.erase(sim_graph1[i].edge_coverage.begin()+k);
                            }
                        }
                    }
                }

                for (int j = 0; j < sim_graph1[i].edge_coverage.size(); ++j) {//给每条边以评分
                    bool only_flag = false;
                    for (auto pe: sim_graph1[i].paired_node) {
                        std::vector <size_t> key = pe.first;
                        if (key[0] == sim_graph1[i].edge_coverage[j].first) {
                            only_flag = true;
                        }
                    }
                    if (!only_flag) {
                        if (consistent_with_abundance(sim_graph1,i,sim_graph1[i].children[j])){//无reads支持但丰度相近
                            sim_graph1[i].edge_coverage[j].second = 2;
                        } else{//无reads支持但丰度不相近
                            sim_graph1[i].edge_coverage[j].second = 1;
                        }
                    } else{//有reads支持
                        if (consistent_with_abundance(sim_graph1,i,sim_graph1[i].children[j])){//有reads支持但丰度相近
                            sim_graph1[i].edge_coverage[j].second = 4;
                        } else{//有reads支持但丰度不相近
                            sim_graph1[i].edge_coverage[j].second = 3;
                        }
//                        sim_graph1[i].edge_coverage[j].second = 3;
                    }
                }

            } else{
                for (int j = 0; j < sim_graph1[i].edge_coverage.size(); ++j) {
                    sim_graph1[i].edge_coverage[j].second = 3;
                }
            }
        }
        //调整首个节点的丰度
        for (int i = 0; i < sim_graph1.size(); ++i) {
            if(sim_graph1[i].parents.size() == 0 && sim_graph1[i].children.size()!= 0){
                bool cov_flag = false;
                double s_c_c = 0;
                for (int j = 0; j < sim_graph1[i].children.size(); ++j) {
                    if (sim_graph1[sim_graph1[i].children[j]].parents.size() == 1){
                        s_c_c += sim_graph1[sim_graph1[i].children[j]].coverage;
                    }
                }
                sim_graph1[i].coverage = s_c_c;

            }
        }
    }

    //双末端相同位置的点
    void get_sim_paired_position(KmerMap& kmerMap,std::vector<std::string>& data,std::vector<Node>& result,vector<Node>& sim_graph1,std::map<size_t,size_t>& read_first_node
                                 ,vector<vector<pair<size_t,size_t>>>& parent_different_node,std::map<size_t,vector<size_t>>& node_pair){
        for (int i = 0; i < sim_graph1.size(); ++i) {
            string str;
            if (parent_different_node[i].size() > 0){
//                for (int j = 0; j < sim_graph1[i].c_node_id.size(); ++j) {
                for (int j = 0; j < parent_different_node[i].size(); ++j) {
                    str = result[parent_different_node[i][j].second].sequence;
                    for (int k = 0; k < str.length() - g_kmer_length + 1; ++k) {
                        string kmer = str.substr(k, g_kmer_length);
                        std::vector<pair<size_t, size_t>> reads = kmerMap.get_kmer_read(kmer);
                        for (int l = 0; l < reads.size(); ++l) {
                            if (reads[l].second == 0){
                                read_first_node[reads[l].first] = i;
                                node_pair[i].push_back(reads[l].first);
                            }
                        }
                    }
                }
            } else{
                double min_cov = 1;
                for (int j = 0; j < sim_graph1[i].c_node_id.size(); ++j) {
                    if (result[sim_graph1[i].c_node_id[j]].coverage < min_cov){
                        min_cov = result[sim_graph1[i].c_node_id[j]].coverage;
                    }
                }
                for (int j = 0; j < sim_graph1[i].c_node_id.size(); ++j) {
                    if (abs(result[sim_graph1[i].c_node_id[j]].coverage-min_cov) < 0.02){
                        str = result[sim_graph1[i].c_node_id[j]].sequence;
                        for (int k = 0; k < str.length() - g_kmer_length + 1; ++k) {
                            string kmer = str.substr(k, g_kmer_length);
                            std::vector<pair<size_t, size_t>> reads = kmerMap.get_kmer_read(kmer);
                            for (int l = 0; l < reads.size(); ++l) {
                                if (reads[l].second == 0){
                                    read_first_node[reads[l].first] = i;
                                    node_pair[i].push_back(reads[l].first);
                                }
                            }
                        }
                    }
                }
            }

        }
    }

    void get_pair_end(std::vector<std::map<size_t,size_t>>& same_position,std::vector<std::string>& data,vector<Node>& result,std::map<size_t,size_t>& read_first_node,std::map<size_t,vector<size_t>>& node_pair){
        for (int i = 0; i < result.size(); ++i) {
            if(result[i].coverage != 1){
                std::vector<size_t> value = node_pair[i];
                std::map<size_t,size_t> local_pos;
                for (int j = 0; j < value.size(); ++j) {
                    if (value[j] < data.size()/2){
                        if (read_first_node.find(value[j]+data.size()/2)!=read_first_node.end()){
                            if (result[read_first_node[value[j]+data.size()/2]].coverage!=1){
                                local_pos[read_first_node[value[j]+data.size()/2]]++;
                            }
                        }
                    }

                }
                same_position[i] = local_pos;
            }

        }
    }

    size_t same_two(std::vector<size_t>& c1,std::vector<size_t>& c2){
        size_t minSize = std::min(c1.size(), c2.size());
        size_t index = 0;
        while (index < minSize && c1[index] == c2[index]) {
            index++;
        }
        return index;
    }

    void get_used_reads(KmerMap& kmerMap,std::vector<std::string>& user_data,std::vector<std::string>& data){
        for (auto it: used_kmers_) {
            std::vector<pair<size_t,size_t>> reads = kmerMap.get_kmer_read(it.first);
            for (int i = 0; i < reads.size(); ++i) {
                user_data.push_back(data[reads[i].first]);
            }
        }
    }
    std::vector<Node> handle_nodes(std::vector<Node>& after_gra){
        int size = after_gra.size();
        std::vector<Node> after_gra2;
        for (int i = after_gra.size(); i >= 0; --i) {
            if (after_gra[i].children.empty() && i < size - 150){
                for (int j = 0; j < after_gra[i].parents.size(); ++j) {
                    for (int k = 0; k < after_gra[after_gra[i].parents[j]].children.size(); ++k) {
                        if (after_gra[after_gra[i].parents[j]].children[k] == i){
                            after_gra[after_gra[i].parents[j]].children.erase(after_gra[after_gra[i].parents[j]].children.begin()+k);
                            break;
                        }
                    }
                }
//                cout << i << endl;
                after_gra[i].parents.clear();
            }
        }
        for (int i = 0; i < after_gra.size(); ++i) {
            if (after_gra[i].parents.empty() && !after_gra[i].children.empty() && after_gra[i].node_layer!=0){
//                cout << i << endl;
                for (int j = 0; j < after_gra[i].children.size(); ++j) {
                    for (int k = 0; k < after_gra[after_gra[i].children[j]].parents.size(); ++k) {
                        if (after_gra[after_gra[i].children[j]].parents[k] == i){
                            after_gra[after_gra[i].children[j]].parents.erase(after_gra[after_gra[i].children[j]].parents.begin()+k);
                            break;
                        }
                    }
                }
                after_gra[i].children.clear();
            }
        }
        for (int i = 0; i < after_gra.size(); ++i) {
            if (!after_gra[i].parents.empty() || !after_gra[i].children.empty()){
                after_gra2.push_back(after_gra[i]);
            }
        }
        return after_gra2;
    }

    void get_local_pair_nodes(KmerMap& kmerMap,vector<Node>& result,vector<Node>& sim_graph1,vector<vector<pair<size_t,size_t>>>& parent_different_node,std::vector<std::string>& data){
        cout << "Being get left pair nodes ..." <<endl;
        for (int i = 0; i < sim_graph1.size(); ++i) {
            bool long_flag = true;
            for (int j = 0; j < sim_graph1[i].children.size(); ++j) {
                if (sim_graph1[i].sequence.length() + sim_graph1[sim_graph1[i].children[j]].sequence.length() -g_kmer_length+1 > 250){
                    long_flag = false;
                    break;
                }
            }
            if (long_flag){
                std::vector<size_t> paired_reads;
                for (int j = 0; j < sim_graph1[i].sequence.length()-g_kmer_length+1; ++j) {
                    string kmer = sim_graph1[i].sequence.substr(j,g_kmer_length);
                    std::vector<pair<size_t,size_t>> reads = kmerMap.get_kmer_read(kmer);
                    for (int k = 0; k < reads.size(); ++k) {
                        if (read_and_node.find(reads[j].first) != read_and_node.end()){
//                        cout << "llllllllll" <<endl;
                            continue;
                        }
                        string read = data[reads[k].first].substr(reads[k].second);
                        string node_seq = sim_graph1[i].sequence.substr(j);
                        if (node_seq.length() > read.length())
                            break;
                        string read_seq = read.substr(0,node_seq.length());

                        if (node_seq == read_seq){//判断当前首个点是否能和reads比对上
                            sim_graph1[i].node_read.push_back(reads[k].first);
                            read = read.substr(node_seq.length());
                        } else{
                            continue;
                        }
                        Node node = sim_graph1[i];
                        bool flag = true;
                        while (flag){
                            if (node.children.size() == 0)
                                break;
                            for (int l = 0; l < node.children.size(); l++) {
                                string node_str = sim_graph1[node.children[l]].sequence.substr(g_kmer_length-1);
                                if (node_str.length() > read.length()){
//                                if (node_str.find(read)!= std::string::npos){
//                                    read_and_node[reads[k].first].push_back(result[node.children[l]].id);
//                                }
                                    flag = false;
                                    continue;
                                } else if(node_str.length() == read.length()){
                                    if (read == node_str){
//                                        read_and_node[reads[k].first].push_back(sim_graph1[node.children[l]].id);
                                        paired_reads.push_back(node.children[l]);
                                        flag = false;
                                        break;
                                    }
                                    flag = false;
                                    continue;
                                }
                                string read_str = read.substr(0,node_str.length());
                                if(read_str == node_str){
                                    flag = true;
                                    read = read.substr(node_str.length());
//                                    read_and_node[reads[k].first].push_back(sim_graph1[node.children[l]].id);
                                    paired_reads.push_back(node.children[l]);
//                                if (i == 3){
//                                    cout << i << " : " << node.children[l] <<endl;
//                                }
                                    node = sim_graph1[node.children[l]];
                                    break;
                                }
                                flag = false;
                            }
                        }
                        if (paired_reads.size() > 0){
                            sim_graph1[i].paired_node[paired_reads].first++;
                            sim_graph1[i].paired_node[paired_reads].second.push_back(reads[k].first);
                        }
                        paired_reads.clear();
                    }
                }
            } else{
                vector<size_t> local_children;
                for (int j = 0; j < sim_graph1[i].children.size(); ++j) {
                    local_children.push_back(sim_graph1[i].children[j]);
                }
                vector<vector<pair<size_t,size_t>>> children_different_node;
                children_different_node.clear();
                children_different_node.resize(sim_graph1.size());
                difference_vec(local_children,children_different_node,sim_graph1);
                for (int l = 0; l < local_children.size(); ++l) {
                    if (parent_different_node[i].size() > 0 && children_different_node[local_children[l]].size()>0){
                        int support = 0;
                        string str1;

                        str1 = result[sim_graph1[local_children[l]].c_node_id[0]].sequence;
                        for (int j = 1; j <= children_different_node[local_children[l]][0].first; ++j) {
                            str1 = str1 + result[sim_graph1[local_children[l]].c_node_id[j]].sequence.substr(g_kmer_length-1);
                        }

                        string str2;
                        for (int d = 0; d < parent_different_node[i].size() ; ++d) {
                            str2 = "";
                            str2 = result[sim_graph1[i].c_node_id[parent_different_node[i][d].first]].sequence;
                            for (int j = parent_different_node[i][d].first+1; j < sim_graph1[i].c_node_id.size(); ++j) {//截取能包含两端两个特殊节点的最考前的开始位置
                                str2 = str2 + result[sim_graph1[i].c_node_id[j]].sequence.substr(g_kmer_length-1);

                            }
                            if (str2.length() + str1.length() - g_kmer_length+1 < 250){
                                break;
                            }
                        }
                        string str3 = "";
                        if(parent_different_node[i][parent_different_node[i].size()-1].first < sim_graph1[i].c_node_id.size()-1){
                            for (int j = parent_different_node[i][parent_different_node[i].size()-1].first+1; j < sim_graph1[i].c_node_id.size(); ++j) {
                                str3 = result[sim_graph1[i].c_node_id[j]].sequence.substr(g_kmer_length-1);
                            }
                        }
                        if (str2.length() + sim_graph1[local_children[l]].sequence.substr(g_kmer_length-1).length() < 250){
                            string str6 = str2 + sim_graph1[local_children[l]].sequence.substr(g_kmer_length-1);
                            for (int u = 0; u < data.size(); ++u) {
                                if (data[u].find(str6)!=std::string::npos){
                                    vector<size_t> p_n;
                                    p_n.push_back(local_children[l]);
//                                            cout << i << " : " << sim_graph1[i].children[j] <<endl;
                                    sim_graph1[i].paired_node[p_n].first++;
                                    sim_graph1[i].paired_node[p_n].second.push_back(u);
                                }
                            }
                        } else if (str2.length() + str1.length() - g_kmer_length+1 < 250) {
                            string str4 = str2 + str1.substr(g_kmer_length - 1);
                            for (int j = 0; j < str2.length() - str3.length() - g_kmer_length + 2; ++j) {
                                string kmer = str2.substr(j, g_kmer_length);
                                std::vector <pair<size_t, size_t>> reads = kmerMap.get_kmer_read(kmer);
                                for (int r = 0; r < reads.size(); ++r) {
                                    string read = data[reads[r].first].substr(reads[r].second);
                                    string str5 = str2.substr(j);
                                    if (str5 == read.substr(0, str5.length())) {
                                        read = read.substr(str5.length());
                                        if (sim_graph1[local_children[l]].sequence.length() - g_kmer_length + 1 > read.length()) {
//                                            cout << sim_graph1[local_children[l]].sequence.substr(g_kmer_length-1) << endl;
//                                            cout << "-- " << read <<endl;areEqualWithOneOrLessErrors(
                                            if (sim_graph1[local_children[l]].sequence.substr(g_kmer_length - 1,read.length()) == read) {
                                                vector < size_t > p_n;
                                                p_n.push_back(local_children[l]);
//                                            cout << i << " : " << sim_graph1[i].children[j] <<endl;
                                                sim_graph1[i].paired_node[p_n].first++;
                                                sim_graph1[i].paired_node[p_n].second.push_back(reads[r].first);
                                            }
                                        } else {
//                                            cout << sim_graph1[local_children[l]].sequence.substr(g_kmer_length-1) << endl;
//                                            cout << read <<endl;areEqualWithOneOrLessErrors(
                                            if (sim_graph1[local_children[l]].sequence.substr(g_kmer_length - 1) ==
                                                read.substr(0, sim_graph1[local_children[l]].sequence.length()-g_kmer_length + 1)) {
                                                vector < size_t > p_n;
                                                p_n.push_back(local_children[l]);
//                                            cout << i << " : " << sim_graph1[i].children[j] <<endl;
                                                sim_graph1[i].paired_node[p_n].first++;
                                                sim_graph1[i].paired_node[p_n].second.push_back(reads[r].first);
                                            }
                                        }

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }

    bool same_score(std::vector<std::pair<size_t, size_t>>& edge){
        map<size_t,size_t> same_num;
        bool flag;
        for (int i = 0; i < edge.size(); ++i) {
            same_num[edge[i].second]++;
        }
        if (same_num.size() == 1){
            flag = true;
        } else{
            flag = false;
        }
        return flag;
    }

    bool not_p_c_unique(vector<Node>& sim_graph1,size_t & i){//i节点的父亲和孩子都不只于他相连
        bool flag = true;
        for (int j = 0; j < sim_graph1[i].children.size(); ++j) {
            if (sim_graph1[sim_graph1[i].children[j]].parents.size() == 1){
                flag = false;
                break;
            }
        }
        for (int j = 0; j < sim_graph1[i].parents.size(); ++j) {
            if (sim_graph1[sim_graph1[i].parents[j]].children.size() == 1){
                flag = false;
                break;
            }
        }
        return flag;
    }

    std::vector<Node> merge_orders_in_and_out(std::vector<Node>& result,vector<Node>& sim_graph1){
        cout << "Begin merge single points ..." <<endl;
        std::vector<Node> sim2;
        std::vector<Node> initial = sim_graph1;
        for (int i = 0; i < sim_graph1.size(); ++i) {
            int num = 0;
            while (sim_graph1[i].children.size() == 1 && sim_graph1[sim_graph1[i].children[0]].parents.size() == 1){//互为唯一相连的节点
                num++;
                size_t child_idx = sim_graph1[i].children[0];
                sim_graph1[i].sequence = sim_graph1[i].sequence + sim_graph1[sim_graph1[i].children[0]].sequence.substr(g_kmer_length-1);
                for (int c = 0; c < sim_graph1[sim_graph1[i].children[0]].c_node_id.size(); ++c) {
                    sim_graph1[i].c_node_id.push_back(sim_graph1[sim_graph1[i].children[0]].c_node_id[c]);
                }
                for (int j = 0; j < sim_graph1[sim_graph1[i].children[0]].children.size(); ++j) {
                    for (int k = 0; k < sim_graph1[sim_graph1[sim_graph1[i].children[0]].children[j]].parents.size(); ++k) {
                        if (sim_graph1[sim_graph1[sim_graph1[i].children[0]].children[j]].parents[k] == sim_graph1[i].children[0]){
                            sim_graph1[sim_graph1[sim_graph1[i].children[0]].children[j]].parents[k] = i;
                        }
                    }
                }
                sim_graph1[i].children = sim_graph1[sim_graph1[i].children[0]].children;
                sim_graph1[i].coverage += sim_graph1[sim_graph1[i].children[0]].coverage;
                sim_graph1[child_idx].children.clear();
                sim_graph1[child_idx].parents.clear();
            }
            sim_graph1[i].coverage = sim_graph1[i].coverage/num;
        }
        int new_count = 0;
        for (int i = 0; i < sim_graph1.size(); ++i) {
            if (sim_graph1[i].children.size()!=0 || sim_graph1[i].parents.size()!=0 || result[sim_graph1[i].c_node_id[0]].parents.empty()){
                sim2.push_back(sim_graph1[i]);
                new_count++;
            }
        }
//        cout << sim2.size() <<endl;
        return sim2;
    }

    vector<pair<vector<size_t>,double>> find_paths(std::vector<Node>& result,vector<Node>& sim_graph1,std::vector<std::map<size_t,size_t>>& same_position_sim){
            cout << "Begin find paths ... " << endl;
            vector<pair<vector<size_t>,double>> paths;
            int m = 0;
            double sum_cov = 0;
            vector<Node> sim_graph = sim_graph1;
            map<size_t,size_t> use_node;
        int start_num = 0;
        double start_aver = 1;
        for (int i = 0; i < sim_graph1.size(); ++i) {
            if (sim_graph1[i].parents.size() == 0 && sim_graph1[i].children.size()!= 0){
                start_num++;
//                start_aver = start_aver + sim_graph1[i].coverage;
//                cout << i << " : " << sim_graph1[i].coverage << "  " <<endl;
            }
            if (sim_graph1[i].parents.size() == 0 && sim_graph1[i].children.size()!= 0 && sim_graph1[i].coverage > 0 && start_aver > sim_graph1[i].coverage){
                start_aver = sim_graph1[i].coverage;
            }
        }
        double limit_cov;
        if (start_aver > 0.05){
            limit_cov = 0.08;
        } else{
            limit_cov = start_aver;
        }
//        cout << "limit cov : " << limit_cov <<endl;
            for (int i = 0; i < sim_graph1.size(); ++i) {
                double min_start = 1;
                int min_start_i;
                bool start_flag = false;
                if (start_num > 8){
                    for (int j = 0; j < sim_graph1.size(); ++j) {
                        if (sim_graph1[j].parents.size() == 0 && sim_graph1[j].children.size()!= 0 && sim_graph1[j].coverage > limit_cov && min_start > sim_graph1[j].coverage ){
                            min_start = sim_graph1[j].coverage;
                            min_start_i = j;
                            start_flag = true;
                        }
                    }
                } else{
                    for (int j = i; j < sim_graph1.size(); ++j) {
                        if (sim_graph1[j].parents.size() == 0 && sim_graph1[j].children.size()!= 0 && sim_graph1[j].coverage > limit_cov){
                            min_start_i = j;
                            start_flag = true;
                            break;
                        }
                    }
                }

                if (!start_flag){
                    continue;
                }
//                cout << min_start_i << "  start_coverage : " << min_start <<endl;

                while(sim_graph1[min_start_i].coverage > limit_cov){
                    vector<double> path_node_coverage;
                    if (start_num > 8){
                        if (sim_graph1[min_start_i].coverage < 0.02 || use_node[min_start_i] >= 1){
                            break;
                        }
                    } else{
                        if (sim_graph1[min_start_i].coverage < 0.02 && use_node[min_start_i] >= 1){
                            break;
                        }
                    }

                    if (abs(sum_cov - 1) < 0.02 || sum_cov > 1.02){
                        break;
                    }

                    m++;
                    Node node = sim_graph1[min_start_i];
                    vector<size_t> path;
                    path.push_back(node.id);
                    use_node[node.id]++;
                    double path_cov = 10;
                    while (node.children.size()>0){
                        if (node.children.size() == 1){
                            node = sim_graph1[node.children[0]];
                            path.push_back(node.id);
                            use_node[node.id]++;
                            if (node.coverage > 0){
                                path_node_coverage.push_back(node.coverage);
                            }
                        } else{
                            bool p_flag = false;
                            node_idx_t pair_i;
                            int p_i;
                            vector<int> pair_num;
                            for (int s = path.size()-1; s >= 1; --s) {//根据配对节点判断

                                if (sim_graph1[path[s]].parents.size() == node.children.size()) {
                                    for (int c = 0; c < node.children.size(); ++c) {
                                        if (same_position_sim[path[s-1]].find(node.children[c]) != same_position_sim[path[s-1]].end() && sim_graph1[node.children[c]].coverage != 0) {
                                            pair_num.push_back(node.children[c]);
                                        }
                                    }
                                    if (pair_num.size() > 0) {//可以通过配对信息找到唯一与这条路匹配的路径
                                        p_flag = true;
                                        p_i = s;
                                        break;
                                    }
                                }
                            }
                            bool find_flag = false;
                            vector<int> cov_num;
                            double child_min_cov = 1;
                            for (int cmc = 0; cmc < node.children.size(); ++cmc) {
                                if (sim_graph1[node.children[cmc]].coverage < child_min_cov){
                                    child_min_cov = sim_graph1[node.children[cmc]].coverage;
                                }
                            }
                            for (int k = path.size() - 1; k >= 1; --k) {//根据点的丰度判断
                                if (sim_graph1[path[k]].parents.size() == node.children.size()) {
                                    for (int j = 0; j < node.children.size(); ++j) {
                                        if (abs(sim_graph1[path[k-1]].coverage - sim_graph1[node.children[j]].coverage) <
                                                    child_min_cov * 0.2 && sim_graph1[node.children[j]].coverage!=0) {
                                            find_flag = true;
                                            cov_num.push_back(node.children[j]);
                                        }
                                    }
                                    if (cov_num.size() > 0) {
                                        break;
                                    }
                                }
                            }

                            if (find_flag && p_flag){//既符合配对条件也符合丰度条件
                                if (pair_num.size() == 1 && cov_num.size() == 1 && pair_num[0] == cov_num[0]){
                                    node = sim_graph1[pair_num[0]];
                                    path.push_back(node.id);
                                    use_node[node.id]++;
                                    if (node.coverage > 0){
                                        path_node_coverage.push_back(node.coverage);
                                    }
                                } else if (pair_num.size() == 1 && cov_num.size() == 1 && pair_num[0] != cov_num[0]){
                                    node = sim_graph1[pair_num[0]];
                                    path.push_back(node.id);
                                    use_node[node.id]++;
                                    if (node.coverage > 0){
                                        path_node_coverage.push_back(node.coverage);
                                    }
                                } else if (pair_num.size() == 1 && cov_num.size() != 1){
                                    node = sim_graph1[pair_num[0]];
                                    path.push_back(node.id);
                                    use_node[node.id]++;
                                    if (node.coverage > 0){
                                        path_node_coverage.push_back(node.coverage);
                                    }
                                } else if (pair_num.size() != 1 && cov_num.size() == 1){
                                    node = sim_graph1[cov_num[0]];
                                    path.push_back(node.id);
                                    use_node[node.id]++;
                                    if (node.coverage > 0){
                                        path_node_coverage.push_back(node.coverage);
                                    }
                                } else if (pair_num.size() > 1 && cov_num.size() > 1){
                                    double min_cov = 1;
                                    int min_n_i;
                                    for (int pa = 0; pa < pair_num.size(); ++pa) {
                                        if (sim_graph1[pair_num[pa]].coverage < min_cov){
                                            min_cov = sim_graph1[pair_num[pa]].coverage;
                                            min_n_i = pair_num[pa];
                                        }
                                    }
                                    node = sim_graph1[min_n_i];
                                    path.push_back(node.id);
                                    use_node[node.id]++;
                                    if (node.coverage > 0){
                                        path_node_coverage.push_back(node.coverage);
                                    }
                                }

                            } else if (p_flag && !find_flag){

                                if (pair_num.size() == 1){
                                    node = sim_graph1[pair_num[0]];
                                    path.push_back(node.id);
                                    use_node[node.id]++;
                                    if (node.coverage > 0){
                                        path_node_coverage.push_back(node.coverage);
                                    }
                                } else{
                                    double min_cov = 1;
                                    int min_n_i;
                                    for (int pa = 0; pa < pair_num.size(); ++pa) {
                                        if (sim_graph1[pair_num[pa]].coverage < min_cov){
                                            min_cov = sim_graph1[pair_num[pa]].coverage;
                                            min_n_i = pair_num[pa];
                                        }
                                    }
                                    node = sim_graph1[min_n_i];
                                    path.push_back(node.id);
                                    use_node[node.id]++;
                                    if (node.coverage > 0){
                                        path_node_coverage.push_back(node.coverage);
                                    }

                                }

                            } else if (!p_flag && find_flag){
                                if (cov_num.size() == 1){
                                    node = sim_graph1[cov_num[0]];
                                    path.push_back(node.id);
                                    use_node[node.id]++;
                                    if (node.coverage > 0){
                                        path_node_coverage.push_back(node.coverage);
                                    }
                                } else{
                                    double min_cov = 1;
                                    int min_n_i;
                                    for (int pa = 0; pa < cov_num.size(); ++pa) {
                                        if (sim_graph1[cov_num[pa]].coverage < min_cov){
                                            min_cov = sim_graph1[cov_num[pa]].coverage;
                                            min_n_i = cov_num[pa];
                                        }
                                    }
                                    node = sim_graph1[min_n_i];
                                    path.push_back(node.id);
                                    use_node[node.id]++;
                                    if (node.coverage > 0){
                                        path_node_coverage.push_back(node.coverage);
                                    }
                                }
                            } else{

                                    bool final_flag = false;
                                    bool final_flag2 = false;
                                    node_idx_t c2;
                                    if (path_cov!=10){
//                                        cout << "!10  " << node.id <<endl;
                                        for (int j = 0; j < node.children.size(); ++j) {
//                                            cout << abs(path_cov - sim_graph1[node.children[j]].coverage) << "  " << sim_graph1[node.children[j]].id << " : "<< sim_graph1[node.children[j]].coverage<<endl;
                                            if (abs(path_cov - sim_graph1[node.children[j]].coverage) < 0.02 && sim_graph1[node.children[j]].coverage!=0) {

                                                final_flag2 = true;
                                                c2 = node.children[j];

                                            }
                                        }

                                    }
                                    if (final_flag2) {//与路径丰度相近
                                        node = sim_graph[c2];
                                        path.push_back(node.id);
                                        use_node[node.id]++;
                                        if (node.coverage > 0){
                                            path_node_coverage.push_back(node.coverage);
                                        }
//                                        cout << node.id <<endl;
                                    } else {//与路径丰度不相近
                                        for (int k = path.size() - 1; k >= 1; --k) {
                                            node_idx_t c;
                                            for (int j = 0; j < node.children.size(); ++j) {
                                                if (abs(sim_graph1[path[k]].coverage - sim_graph1[node.children[j]].coverage) < 0.02 &&
                                                    sim_graph1[node.children[j]].coverage != 0) {//当前点是否与路径中的某个点丰度相似
                                                    final_flag = true;
                                                    c = node.children[j];
                                                }
                                            }
                                            if (final_flag) {
                                                node = sim_graph[c];
                                                path.push_back(node.id);
                                                use_node[node.id]++;
                                                if (node.coverage > 0) {
                                                    path_node_coverage.push_back(node.coverage);
                                                }
                                                break;
                                            }
                                        }
//                                    cout << "sssssssssss  " << node.id <<endl;
                                        if (!final_flag) {//在当前路径中没有相似丰度的点，看丰度和
                                            double s_sim_cov = 0;
                                            for (int j = 0; j < node.children.size(); ++j) {
                                                s_sim_cov = s_sim_cov + sim_graph1[node.children[j]].coverage;
                                            }
                                            if (abs(s_sim_cov - node.coverage) < 0.03) {//此处可能是个分支
                                                double min_c_c = 1;
                                                int min_c_c_i;
                                                bool min_c_c_flag = false;
                                                for (int j = 0; j < node.children.size(); ++j) {
                                                    if (sim_graph1[node.children[j]].coverage != 0 &&
                                                        sim_graph1[node.children[j]].coverage < min_c_c) {//先选取两个可能的分支中丰度小的那一条
                                                        min_c_c = sim_graph1[node.children[j]].coverage;
                                                        min_c_c_i = j;
                                                        min_c_c_flag = true;
                                                    }
                                                }
                                                if (min_c_c_flag) {//分支中还有非0的分支
                                                    node = sim_graph1[node.children[min_c_c_i]];
                                                    path.push_back(node.id);
                                                    use_node[node.id]++;
                                                    if (node.coverage > 0) {
                                                        path_node_coverage.push_back(node.coverage);
                                                    }
                                                if (node.id < 100){
                                                    path_cov = node.coverage;
//                                                    cout << node.id << "  :  " << path_cov << endl;
                                                }

                                                } else {//可能不能这样
                                                    double max_f = 0;
                                                    int max_f_i;
                                                    for (int c = 0; c < node.children.size(); ++c) {
                                                        if (sim_graph[node.children[c]].coverage > max_f) {
                                                            max_f = sim_graph[node.children[c]].coverage;
                                                            max_f_i = c;
                                                        }
                                                    }
                                                    for (int e = 0; e < sim_graph1[node.id].edge_coverage.size(); ++e) {
                                                        if (sim_graph1[node.id].edge_coverage[e].first == node.children[max_f_i]){
                                                            sim_graph1[node.id].edge_coverage[e].second = 0;
                                                        }
                                                    }
//                                                    cout << " 1dddd " << node.id << endl;
                                                    node = sim_graph1[node.children[max_f_i]];
                                                    path.push_back(node.id);
                                                    use_node[node.id]++;
                                                    if (node.coverage > 0) {
                                                        path_node_coverage.push_back(node.coverage);
                                                    }
                                                }

//                                    cout << "33  " << node.id << " : " << min_cov <<endl;
                                            } else {
                                                double max_f = 0;
                                                int max_f_i;
                                                for (int c = 0; c < node.children.size(); ++c) {
                                                    if (sim_graph[node.children[c]].coverage > max_f) {
                                                        max_f = sim_graph[node.children[c]].coverage;
                                                        max_f_i = c;
                                                    }
                                                }
                                                for (int e = 0; e < sim_graph1[node.id].edge_coverage.size(); ++e) {
                                                    if (sim_graph1[node.id].edge_coverage[e].first == node.children[max_f_i]){
                                                        sim_graph1[node.id].edge_coverage[e].second = 0;
                                                    }
                                                }
//                                                cout << " 2dddd " << node.id << endl;
                                                node = sim_graph[node.children[max_f_i]];

                                                path.push_back(node.id);
                                                use_node[node.id]++;
                                                if (node.coverage > 0) {
                                                    path_node_coverage.push_back(node.coverage);
                                                }
                                            }

                                        }
                                    }
                            }
                        }
                    }
                    std::sort(path_node_coverage.begin(), path_node_coverage.end());

                    if (path_cov == 10){
                        path_cov = 0;
                        double num_n = 0;

                            for (int pa = 3; pa < path_node_coverage.size()*0.5; ++pa) {
                                if (path_node_coverage[pa]>0){
                                    path_cov = path_cov + path_node_coverage[pa];
                                    num_n++;
                                }
                            }

                        path_cov = path_cov/num_n;
                    }
                    pair<vector<size_t>,double> p;
                    p.first = path;
                    p.second = path_cov;
                    paths.push_back(p);

//                    cout << "min_cov : " << path_cov <<endl;
                    sum_cov = sum_cov + path_cov;
                    for (int j = 0; j < path.size(); ++j) {
                        if (abs(sim_graph1[path[j]].coverage - path_cov) < 0.03 || sim_graph1[path[j]].coverage < path_cov){
                            sim_graph1[path[j]].coverage = 0;
                        } else{
                            sim_graph1[path[j]].coverage = sim_graph1[path[j]].coverage - path_cov;
                        }
                    }

                }
            }
//            cout << " not used nodes num : " << sim_graph1.size() - use_node.size() <<endl;
//            cout << " now path : " << paths.size() <<endl;
//            cout << " ------------------------ "<<endl;
            for (int i = 0; i < sim_graph1.size()*0.9; ++i) {//根据剩余没用过的节点再找路
                vector<size_t> path2;
//                if (sim_graph1[i].children.size() > 0){
//                    cout << sim_graph1[i].edge_coverage[0].second <<endl;
//                }
                if (use_node.find(i)==use_node.end() && sim_graph1[i].children.size() == 1 &&
                use_node.find(sim_graph1[i].children[0])==use_node.end() && sim_graph1[sim_graph1[i].children[0]].parents.size() == 1 &&
                (sim_graph1[i].edge_coverage[0].second == 4 || sim_graph1[i].edge_coverage[0].second == 3)){
//                    cout << i <<endl;
                    Node node = sim_graph1[i];
                    path2.push_back(node.id);
                    use_node[node.id]++;
                    double path_cov2 = sim_graph1[i].coverage;
                    int used_count = use_node[node.id];
                    int used_num = 0;
                    while(!node.children.empty()){
//                        cout << node.id << " : " << use_node[node.id] << endl;
                        if (use_node[node.id] == used_count){
                            used_num++;
                        } else{
                            used_count = use_node[node.id];
                            used_num = 0;
                        }
                        if(use_node[node.id] >= 10 || (used_count >= 4 && used_num >= 5)){
//                            cout << " 101010101010 " <<endl;
                            break;
                        }
                        if (node.children.size() == 1){
                            node = sim_graph1[node.children[0]];
                            path2.push_back(node.id);
//                            cout << "22  " << node.id <<endl;
                            use_node[node.id]++;
                        } else{
                            bool is_used_flag = false;
                            for (int c = 0; c < node.children.size(); ++c) {
                                if (use_node.find(node.children[c])==use_node.end()){
                                    is_used_flag = true;
                                    node = sim_graph1[node.children[c]];
                                    path2.push_back(node.id);
//                                    cout << "55  " << node.id <<endl;
                                    use_node[node.id]++;
                                }
                            }
                            if (!is_used_flag){
                                bool p_flag = false;
                                for (int s = path2.size()-1; s >= 1; --s) {
                                    int pair_num = 0;
                                    node_idx_t pair_i;
                                    for (int c = 0; c < node.children.size(); ++c) {
                                        if (same_position_sim[path2[s]].find(node.children[c]) !=
                                            same_position_sim[path2[s]].end()) {
                                            pair_num++;
                                            pair_i = node.children[c];
                                        }
                                    }
                                    if (pair_num == 1) {
                                        node = sim_graph1[pair_i];
                                        p_flag = true;
                                        path2.push_back(node.id);
                                        use_node[node.id]++;
//                                        cout << "66  " << node.id <<endl;
                                        break;
                                    }
                                }
                                if (!p_flag){
                                    double cha = 1;
                                    node_idx_t similar_c;
                                    for (int c = 0; c < node.children.size(); ++c) {
                                        if (abs(sim_graph1[node.children[c]].coverage - path_cov2) < cha){
                                            cha = abs(sim_graph1[node.children[c]].coverage - path_cov2);
                                            similar_c = node.children[c];
                                        }
                                    }
                                    node = sim_graph1[similar_c];
                                    path2.push_back(node.id);
//                                    cout << " 77 " << node.id <<endl;
                                    use_node[node.id]++;
                                }
                            }
                        }

                    }
                    node = sim_graph1[i];
                    used_count = use_node[node.id];
                    used_num = 0;
                    while (!node.parents.empty()){
//                        cout << node.id << " : " << use_node[node.id] << endl;
                        if (use_node[node.id] == used_count){
                            used_num++;
                        } else{
                            used_count = use_node[node.id];
                            used_num = 0;
                        }
                        if(use_node[node.id] > 10 || (used_count >= 4 && used_num >= 5)){
//                            cout << " 22222222222222 " <<endl;
                            break;
                        }
                        if (node.parents.size() == 1){
                            node = sim_graph1[node.parents[0]];
                            path2.insert(path2.begin(),node.id);
                            use_node[node.id]++;
//                            cout << "99  " << node.id <<endl;
                        } else{
                            int not_use=0;
                            int not_p;
                            for (int p = 0; p < node.parents.size(); ++p) {
                                if (use_node.find(node.parents[p])== use_node.end()){
                                    not_use++;
                                    not_p = node.parents[p];
                                }
                            }
                            if(not_use == 1){
                                node = sim_graph1[not_p];
                                path2.insert(path2.begin(),node.id);
//                                cout << "88  " << node.id <<endl;
                                use_node[node.id]++;
                            } else{
                                int same_flag = 0;
                                int same_i;
                                for(int nc = 0; nc < node.parents.size(); nc++){
                                    for(int p2 = 0 ;p2 < path2.size() ;p2++){
                                        if (node.parents.size() == sim_graph1[path2[p2]].children.size()){
                                            if(same_position_sim[node.parents[nc]].find(path2[p2+1])!=same_position_sim[node.parents[nc]].end()){
                                                same_flag++;
                                                same_i = node.parents[nc];
                                            }
                                        }

                                    }
                                }
                                if(same_flag == 1){
                                    node = sim_graph1[same_i];
                                    path2.insert(path2.begin(),node.id);
//                                    cout << "99  " << node.id <<endl;
                                    use_node[node.id]++;
                                } else{
                                    double cha = 1;
                                    node_idx_t similar_p;
                                    for (int p = 0; p < node.parents.size(); ++p) {
                                        if (abs(sim_graph1[node.parents[p]].coverage - path_cov2) < cha){
                                            cha = abs(sim_graph1[node.parents[p]].coverage - path_cov2);
                                            similar_p = node.parents[p];
                                        }
                                    }
                                    node = sim_graph1[similar_p];
                                    path2.insert(path2.begin(),node.id);
//                            cout << "88  " << node.id <<endl;
                                    use_node[node.id]++;
                                }

                            }

                        }
                    }

                    pair<vector<size_t>,double> p;
                    p.first = path2;
                    p.second = path_cov2;
                    paths.push_back(p);

//                    cout << "min_cov : " << path_cov2 <<endl;
                    for (int j = 0; j < path2.size(); ++j) {
                        if (abs(sim_graph1[path2[j]].coverage - path_cov2) < 0.03 || sim_graph1[path2[j]].coverage < path_cov2){
                            sim_graph1[path2[j]].coverage = 0;
                        } else{
                            sim_graph1[path2[j]].coverage = sim_graph1[path2[j]].coverage - path_cov2;
                        }
                    }

                }
            }
            if(sim_graph1.size() > 2*use_node.size()){
//                cout << " --------------more double--------------" <<endl;
                for (int i = 0; i < sim_graph1.size()*0.9; ++i) {//根据剩余没用过的节点再找路
                    vector<size_t> path2;
//                if (sim_graph1[i].children.size() > 0){
//                    cout << sim_graph1[i].edge_coverage[0].second <<endl;
//                }
                    if (use_node.find(i)==use_node.end()){
//                    cout << i <<endl;
                        Node node = sim_graph1[i];
                        path2.push_back(node.id);
                        use_node[node.id]++;
                        double path_cov2 = sim_graph1[i].coverage;
                        int used_count = use_node[node.id];
                        int used_num = 0;
                        while(!node.children.empty()){
//                        cout << node.id << " : " << use_node[node.id] << endl;
                            if (use_node[node.id] == used_count){
                                used_num++;
                            } else{
                                used_count = use_node[node.id];
                                used_num = 0;
                            }
                            if(use_node[node.id] >= 10 || (used_count >= 4 && used_num >= 5)){
//                            cout << " 101010101010 " <<endl;
                                break;
                            }
                            if (node.children.size() == 1){
                                node = sim_graph1[node.children[0]];
                                path2.push_back(node.id);
//                            cout << "22  " << node.id <<endl;
                                use_node[node.id]++;
                            } else{
                                bool is_used_flag = false;
                                for (int c = 0; c < node.children.size(); ++c) {
                                    if (use_node.find(node.children[c])==use_node.end()){
                                        is_used_flag = true;
                                        node = sim_graph1[node.children[c]];
                                        path2.push_back(node.id);
//                                    cout << "55  " << node.id <<endl;
                                        use_node[node.id]++;
                                    }
                                }
                                if (!is_used_flag){
                                    bool p_flag = false;
                                    for (int s = path2.size()-1; s >= 1; --s) {
                                        int pair_num = 0;
                                        node_idx_t pair_i;
                                        for (int c = 0; c < node.children.size(); ++c) {
                                            if (same_position_sim[path2[s]].find(node.children[c]) !=
                                                same_position_sim[path2[s]].end()) {
                                                pair_num++;
                                                pair_i = node.children[c];
                                            }
                                        }
                                        if (pair_num == 1) {
                                            node = sim_graph1[pair_i];
                                            p_flag = true;
                                            path2.push_back(node.id);
                                            use_node[node.id]++;
//                                        cout << "66  " << node.id <<endl;
                                            break;
                                        }
                                    }
                                    if (!p_flag){
                                        double cha = 1;
                                        node_idx_t similar_c;
                                        for (int c = 0; c < node.children.size(); ++c) {
                                            if (abs(sim_graph1[node.children[c]].coverage - path_cov2) < cha){
                                                cha = abs(sim_graph1[node.children[c]].coverage - path_cov2);
                                                similar_c = node.children[c];
                                            }
                                        }
                                        node = sim_graph1[similar_c];
                                        path2.push_back(node.id);
//                                    cout << " 77 " << node.id <<endl;
                                        use_node[node.id]++;
                                    }
                                }
                            }

                        }
                        node = sim_graph1[i];
                        used_count = use_node[node.id];
                        used_num = 0;
                        while (!node.parents.empty()){
//                        cout << node.id << " : " << use_node[node.id] << endl;
                            if (use_node[node.id] == used_count){
                                used_num++;
                            } else{
                                used_count = use_node[node.id];
                                used_num = 0;
                            }
                            if(use_node[node.id] > 10 || (used_count >= 4 && used_num >= 5)){
//                            cout << " 22222222222222 " <<endl;
                                break;
                            }
                            if (node.parents.size() == 1){
                                node = sim_graph1[node.parents[0]];
                                path2.insert(path2.begin(),node.id);
                                use_node[node.id]++;
//                            cout << "99  " << node.id <<endl;
                            } else{
                                int not_use=0;
                                int not_p;
                                for (int p = 0; p < node.parents.size(); ++p) {
                                    if (use_node.find(node.parents[p])== use_node.end()){
                                        not_use++;
                                        not_p = node.parents[p];
                                    }
                                }
                                if(not_use == 1){
                                    node = sim_graph1[not_p];
                                    path2.insert(path2.begin(),node.id);
//                                cout << "88  " << node.id <<endl;
                                    use_node[node.id]++;
                                } else{
                                    int same_flag = 0;
                                    int same_i;
                                    for(int nc = 0; nc < node.parents.size(); nc++){
                                        for(int p2 = 0 ;p2 < path2.size() ;p2++){
                                            if (node.parents.size() == sim_graph1[path2[p2]].children.size()){
                                                if(same_position_sim[node.parents[nc]].find(path2[p2+1])!=same_position_sim[node.parents[nc]].end()){
                                                    same_flag++;
                                                    same_i = node.parents[nc];
                                                }
                                            }

                                        }
                                    }
                                    if(same_flag == 1){
                                        node = sim_graph1[same_i];
                                        path2.insert(path2.begin(),node.id);
//                                    cout << "99  " << node.id <<endl;
                                        use_node[node.id]++;
                                    } else{
                                        double cha = 1;
                                        node_idx_t similar_p;
                                        for (int p = 0; p < node.parents.size(); ++p) {
                                            if (abs(sim_graph1[node.parents[p]].coverage - path_cov2) < cha){
                                                cha = abs(sim_graph1[node.parents[p]].coverage - path_cov2);
                                                similar_p = node.parents[p];
                                            }
                                        }
                                        node = sim_graph1[similar_p];
                                        path2.insert(path2.begin(),node.id);
//                            cout << "88  " << node.id <<endl;
                                        use_node[node.id]++;
                                    }

                                }

                            }
                        }
                        if (path_cov2>0.04){
                            pair<vector<size_t>,double> p;
                            p.first = path2;
                            p.second = path_cov2;
                            paths.push_back(p);
//                            cout << "min_cov : " << path_cov2 <<endl;
                            for (int j = 0; j < path2.size(); ++j) {
                                if (abs(sim_graph1[path2[j]].coverage - path_cov2) < 0.03 || sim_graph1[path2[j]].coverage < path_cov2){
                                    sim_graph1[path2[j]].coverage = 0;
                                } else{
                                    sim_graph1[path2[j]].coverage = sim_graph1[path2[j]].coverage - path_cov2;
                                }
                            }
                        }

                    }
                }
            }

//            cout << sim_graph1.size() << " - " << use_node.size() << " = " << sim_graph1.size()-use_node.size() <<endl;
            used_repeat_paths(use_node,paths);
//            vector<pair<vector<size_t>,double>> c_paths = cut_the_path(paths,sim_graph1);
//
//            cout << c_paths.size() <<endl;
//
//            cout << paths.size() <<endl;
            return paths;
        }
        
    void used_repeat_paths(map<size_t,size_t>& use_node,vector<pair<vector<size_t>,double>>& paths){
        set<int,std::greater<int>> similar;
        for (int i = 0; i < paths.size(); ++i) {
            int used_node_num = 0;
            for (int j = 0; j < paths[i].first.size(); ++j) {
                if (use_node[paths[i].first[j]] > 10){
                    used_node_num++;
                }
            }
//            cout << paths[i].first.size() << " - " << used_node_num <<endl;
            if(used_node_num >= paths[i].first.size()*0.4 || paths[i].first.size() < 10){
                similar.insert(i);
//                cout << i << endl;
            }
        }
        for (auto it : similar) {
            paths.erase(paths.begin()+it);
        }
    }
    
    vector<pair<vector<size_t>,double>> cut_the_path(vector<pair<vector<size_t>,double>>& paths,vector<Node>& sim_graph1){
        vector<pair<vector<size_t>,double>> c_paths;
        cout << "Begin cut paths ..." <<endl;
        for (int i = 0; i < paths.size(); ++i) {
            vector<size_t> cut_element;
            vector<int> cut_position;
            for (int j = 0; j < paths[i].first.size(); ++j) {
                for (int k = 0; k < sim_graph1[paths[i].first[j]].edge_coverage.size(); ++k) {

                    if (sim_graph1[paths[i].first[j]].edge_coverage[k].second == 0 && paths[i].first[j+1] == sim_graph1[paths[i].first[j]].edge_coverage[k].first){
//                        cout <<  paths[i].first[j+1] << " : " << sim_graph1[paths[i].first[j]].edge_coverage[k].first <<endl;
                        cut_position.push_back(paths[i].first[j]);
                        cut_element.push_back(j);
                    }
                }
            }
            if(cut_position.size() > 5){
//                cout << cut_position[cut_position.size()/2] <<endl;
                pair<vector<size_t>,double> p1;
                pair<vector<size_t>,double> p2;
                if (cut_position[cut_position.size()/2] < sim_graph1.size()*0.8){
                    p1.first = std::vector<size_t>(paths[i].first.begin(),paths[i].first.begin()+cut_element[cut_element.size()/2]);
                    p2.first = std::vector<size_t>(paths[i].first.begin()+cut_element[cut_element.size()/2],paths[i].first.end());
                    p1.second = paths[i].second;
                    p2.second = paths[i].second;
                    c_paths.push_back(p1);
                    c_paths.push_back(p2);
                } else{
                    c_paths.push_back(paths[i]);
                }

            } else if (cut_position.size() == 5){
//                cout << cut_position[cut_position.size()-1] <<endl;
                pair<vector<size_t>,double> p1;
                pair<vector<size_t>,double> p2;
                if (cut_position[cut_position.size()-1] < sim_graph1.size()*0.8){
                    p1.first = std::vector<size_t>(paths[i].first.begin(),paths[i].first.begin()+cut_element[cut_element.size()-1]);
                    p2.first = std::vector<size_t>(paths[i].first.begin()+cut_element[cut_element.size()-1],paths[i].first.end());
                    p1.second = paths[i].second;
                    p2.second = paths[i].second;
                    c_paths.push_back(p1);
                    c_paths.push_back(p2);
                } else{
                    c_paths.push_back(paths[i]);
                }
            } else{
                c_paths.push_back(paths[i]);
            }
        }
        return c_paths;
    }

    int countDifferences(std::string &str1, std::string &str2) {
        // 获取两个字符串的长度
        size_t len1 = str1.length();
        size_t len2 = str2.length();

        // 计算短字符串的长度
        size_t minLength = std::min(len1, len2);

        // 统计不同的字符数
        int diffCount = 0;
        for (size_t i = 0; i < minLength; ++i) {
            if (str1[i] != str2[i] && str1[i] != 'N' && str2[i] != 'N') {
                ++diffCount;
            }
        }

        // 返回结果
        return diffCount;
    }

    void similarity_between_roads(vector<pair<string,double>>& paths_str,vector<pair<vector<size_t>,double>>& paths,vector<Node>& sim_graph,vector<vector<size_t>>& sim_layer){
        cout << "Begin similarity between roads ... " <<endl;
        map<int,int,std::greater<int>> similar;
        //路径对齐
        int max_path_l = 0;
        vector<int> add_length;
        add_length.resize(paths_str.size());
        for (int i = 0; i < paths.size(); ++i) {
            int align_seq = 0;
            if (sim_graph[paths[i].first[0]].node_layer > 1){
                if(sim_graph[paths[i].first[0]].node_layer == 2){
                    double min_cha = 1;
                    int min_sl;
                    for (int sl = 0; sl < sim_layer[0].size(); ++sl) {
                        if (abs(sim_graph[sim_layer[0][sl]].coverage - paths[i].second) < min_cha){
                            min_cha = abs(sim_graph[sim_layer[0][sl]].coverage - paths[i].second);
                            min_sl = sim_layer[0][sl];
                        }
                    }
                    paths[i].first.insert(paths[i].first.begin(),min_sl);
                    paths_str[i].first = sim_graph[min_sl].sequence + paths_str[i].first;
                    continue;
                }
                for (int s = 0; s < sim_graph[paths[i].first[0]].node_layer-1; ++s) {
                    align_seq = align_seq + sim_graph[sim_layer[s][0]].sequence.length();
                }
                string in_str = "";
                for (int a = 0; a < align_seq; ++a) {
                    in_str = in_str + 'N';
                }
                paths_str[i].first = in_str + paths_str[i].first;
                add_length[i]=in_str.length();
            }
        }

        for (int i = 0; i < paths_str.size()-1; ++i) {
            for (int j = i+1; j < paths_str.size(); ++j) {
               int num = countDifferences(paths_str[i].first,paths_str[j].first);

                int minLength = std::min(paths_str[i].first.length()-add_length[i], paths_str[j].first.length()-add_length[j]);
//                cout << minLength <<endl;
                if (num < minLength*0.02 && num != 0){
                    similar[i]++;
                    similar[j]++;


                }

            }
//            paths_str[i].second = std::round(paths_str[i].second * 1000.0) / 1000.0;
        }
        std::vector<std::pair<int, int>> vec(similar.begin(), similar.end());
        // 按值从大到小排序
        std::sort(vec.begin(), vec.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            return a.second > b.second;
        });
        if (paths_str.size() >= 10 && vec.size() > 2){
            if (vec[0].second == vec[1].second){
                if (paths_str[vec[0].first].second > paths_str[vec[1].first].second) {
                    paths_str.erase(paths_str.begin()+vec[1].first);
                } else {
                    paths_str.erase(paths_str.begin()+vec[0].first);
                }
            } else{
                paths_str.erase(paths_str.begin()+vec[0].first);
                paths_str.erase(paths_str.begin()+vec[1].first);
            }
        }

        similar.clear();
        for (int i = 0; i < paths_str.size(); ++i) {
            paths_str[i].first.erase(std::remove(paths_str[i].first.begin(), paths_str[i].first.end(), 'N'), paths_str[i].first.end());
        }
        for (int i = 0; i < paths_str.size()-1; ++i) {
            for (int j = i + 1; j < paths_str.size(); ++j) {
                int num = countDifferences(paths_str[i].first, paths_str[j].first);
//                cout << i << "(" << add_length[i] << ")" << " - " << j << " : " << num <<endl;
                int minLength = std::min(paths_str[i].first.length(),
                                         paths_str[j].first.length());
//                cout << minLength << endl;
                if (num < minLength * 0.01 && num != 0) {
                    if(paths_str[i].first.length()==paths_str[j].first.length()){
//                    if (similar.find(i) == similar.end() && similar.find(j) == similar.end()) {
                        if (paths_str[i].second > paths_str[j].second) {
                            similar[j]++;
                        } else {
                            similar[i]++;
                        }
//                    }
                    } else{
                        if (similar.find(i) == similar.end() && similar.find(j) == similar.end()){
                            if (paths_str[i].first.length()-add_length[i] > paths_str[j].first.length()-add_length[j]){
                                similar[j]++;
                            } else{
                                similar[i]++;
                            }
                        }
                    }
                }
            }
        }
        if (paths_str.size() >= 10){
            for (auto it : similar) {
                paths_str.erase(paths_str.begin()+it.first);
            }
        }

    }

    // 计算两个字符串的汉明距离
    int hammingDistance(const std::string& s1, const std::string& s2) {
        int distance = 0;
        for (size_t i = 0; i < s1.size(); ++i) {
            if (s1[i] != s2[i]) {
                ++distance;
            }
        }
        return distance;
    }

// 查找最相似的字符串
    int findMostSimilar(std::string& target, vector<pair<string,double>>& paths_str,vector<int>& paths_num,int& m) {
        int similar_path;
        int minDistance = std::numeric_limits<int>::max();

        for (int i = 0; i < paths_num.size(); i++) {
            string cut_path = paths_str[paths_num[i]].first.substr(m,target.length());
            int distance = hammingDistance(target, cut_path);
            if (distance < minDistance) {
                minDistance = distance;
                similar_path = i;
            }
        }
        return similar_path;
    }

    bool other_path_nuc(vector<pair<string,double>>& paths_str,vector<int>& use_paths,int& p,int& n,char& nucl){
        bool flag = true;
        for (int i = 0; i < use_paths.size(); ++i) {
            char c = paths_str[use_paths[i]].first[n];
            if (use_paths[i]!=p && c == nucl){
                flag = false;
            }
        }
        return flag;
    }

    void output_graph(KmerMap& kmerMap, kmer_int_type_t& seed_kmer,int& average,std::vector<std::string>& data){

        cout << "Begin extend sequence..." << endl;
        string trunk_str;

        trunk_str = get_trunk(kmerMap,seed_kmer,average);
        cout << trunk_str <<endl;
        std::map<string,pair<size_t,size_t>> span_part;
        cout << "Begin get all nodes..." << endl;
        std::vector<std::set<string>> nodes = get_all_nodes(kmerMap,trunk_str,span_part);

        //获取所有kmer
        std::vector<kmer_int_type_t> kmer_all = kmerMap.get_all_kmer();
        //获取首部缺少的分支kmer
        not_used_candi(kmerMap,kmer_all,nodes);
        //删除错误kmer
        delete_error_kmer(kmerMap,nodes);


        //连接成图
        vector<std::map<string,node_idx_t>> node_id_position = get_bubble_graph(nodes,span_part);
        set_parents(node_set_);

        cout << "after node num : " << node_set_.size() <<endl;
        std::vector<Node> initial_gra = node_set_;
        //合并特殊点
        std::map<size_t,std::set<size_t>> layer_nodes;
        std::vector<Node> after_gra = unique_connect(initial_gra,node_id_position,nodes,layer_nodes);
        re_id(after_gra);
        set_parents(after_gra);
        int after_gra_before = after_gra.size();

        //给点添加丰度，并将丰度归一化
        add_kmer_mean_value(kmerMap,after_gra);
        get_coverage_normalization(after_gra,layer_nodes);
        cout << "after before size : " << after_gra.size() <<endl;
        std::vector<Node> after_gra2 = handle_nodes(after_gra);

        if (after_gra2.size()!= after_gra.size()){
            re_id(after_gra2);
            set_parents(after_gra2);
        }
        cout << "after after size : " << after_gra2.size() <<endl;
       //拓扑排序
        std::vector<Node> result = topologicalSort(after_gra2);

        adjust_layer(after_gra2,layer_nodes);

        int limit_length = 0;
        for(int li = 0 ;li < 500; li++){
//            if (limit_length > data[li].length()){
//                limit_length = data[li].length();
//            }
            limit_length = limit_length + data[li].length();
        }
        limit_length = limit_length/500;
        if(limit_length < g_kmer_length*4){
            limit_length = 250;
        }
        cout << "limit_length : " << limit_length <<endl;

        time_t compression_begin = time(NULL);
        vector<vector<size_t>> sim_layer;
        vector<size_t> layer_sim_node;
        vector<vector<size_t>> sim_gra = compression_nodes(kmerMap,result,after_gra2,layer_nodes,sim_layer,layer_sim_node,data,limit_length);
        time_t compression_end = time(NULL);
        cout << "compression_time: " << compression_end - compression_begin << "s" << endl;

        vector<set<size_t>> sim_alone_node;
        vector<vector<pair<size_t,size_t>>> parent_different_node;
        vector<vector<pair<size_t,size_t>>> children_different_node;
        vector<Node> sim_graph1;
        create_sim_graph(kmerMap,result,after_gra2,sim_gra,sim_graph1,layer_nodes,sim_layer,layer_sim_node,parent_different_node,data);

        get_local_pair_nodes(kmerMap,result,sim_graph1,parent_different_node,data);
        std::vector<std::vector<std::vector<size_t>>> sim_keys_to_remove = merge_duplicate_points(sim_graph1);
        delete_paired_nodes(sim_graph1,sim_keys_to_remove);
        std::vector<std::vector<std::vector<size_t>>> sim_re_keys_to_remove = re_merge_duplicate_points(sim_graph1);
        re_delete_paired_nodes(sim_graph1,sim_re_keys_to_remove);
        //删除可能是错误的连边，并给边分数
        based_on_pairing_information_delete_edge(sim_graph1);

        std::map<size_t,size_t> read_first_node_sim;//reads对应的点
        std::map<size_t,vector<size_t>> node_pair_sim;//点对应的reads
        std::vector<std::map<size_t,size_t>> same_position_sim;
        same_position_sim.resize(sim_graph1.size());
        get_sim_paired_position(kmerMap,data,result,sim_graph1,read_first_node_sim,parent_different_node,node_pair_sim);
        get_pair_end(same_position_sim,data,sim_graph1,read_first_node_sim,node_pair_sim);

        vector<pair<vector<size_t>,double>> paths = find_paths(result,sim_graph1,same_position_sim);

        vector<pair<string,double>> paths_str;
        for (int i = 0; i < paths.size(); ++i) {
            string str = "";
            for (int j = 0; j < paths[i].first.size(); ++j) {
                if (j == 0){
                    str = str + sim_graph1[paths[i].first[j]].sequence;
                } else{
                    str = str + sim_graph1[paths[i].first[j]].sequence.substr(g_kmer_length-1);
                }
            }
            pair<string,double> ps;
            ps.first = str;
            ps.second = paths[i].second;
            paths_str.push_back(ps);

        }
        //删除相似度高的contig
        similarity_between_roads(paths_str,paths,sim_graph1,sim_layer);
        cout << "path_num : " << paths_str.size() <<endl;

        cout << "  Strain    Frequency    Length  " <<endl;
        ofstream path_file(output_filename.c_str());

        for (int i = 0; i < paths_str.size(); ++i) {
            path_file << ">" << i << " frequency (" << paths_str[i].second << ")" <<endl;
            path_file << paths_str[i].first <<endl;
            std::cout << std::fixed << std::setprecision(3);
            cout << "  " << i << "         " << paths_str[i].second << "        " << paths_str[i].first.length() << endl;
        }
        path_file.close();

    }

    std::vector<Node> node_set_;
    size_t size_;
    std::map<kmer_int_type_t,size_t> used_kmers_;
    std::map<size_t,size_t> true_nodes;
    std::map<size_t,vector<size_t>> read_and_node;

};
#endif