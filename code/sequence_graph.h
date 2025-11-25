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
            ar & single_node_id;
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
            single_node_id = node.single_node_id;
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
        std::vector<size_t> single_node_id;
        int node_layer;
        std::vector<pair<size_t,double>> edge_coverage;//first与当前点连接的点，second与这个点连接的边的丰度
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

            for (size_t i = 0; i < candidates.size(); ++i) {
                if (!is_used(candidates[i].first) && candidates[i].second > candidate_average*0.02) {
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

    string get_trunk(KmerMap& kmer_map, kmer_int_type_t seed ,int& average,string& st) {
        add_used_kmer(seed, 0);
        std::string left;
        std::string right;
        int type;
        for (int i = 0; i < st.length()-g_kmer_length+1; ++i) {
            kmer_int_type_t intval = kmer_to_intval(st.substr(i,g_kmer_length));
            add_used_kmer(intval,kmer_map.get_kmer_count(intval));
        }
//        left = reverse_extend(kmer_map, seed);//向前扩展
        right = forward_extend(kmer_map, seed,type);//向后扩展

        std::string trunk = st + right;
//        std::string trunk = left + right.substr(g_kmer_length);
        return trunk;
    }

    std::vector<std::set<string>> get_all_nodes(KmerMap& kmer_map,string& trunk,std::map<string,pair<size_t,size_t>>& span_part,int& in_del_num){
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

                                } else{
                                    if (type == 1 && final_k - i < 50 && final_k - i > 0 && extend_str.length()+i < trunk.length()){

                                        for (int k = 1; k < extend_str.length()-g_kmer_length-1; ++k) {
                                            string kmer = extend_str.substr(k,g_kmer_length);
                                            nodes[i+k].insert(kmer);
                                            pair<size_t,size_t> two;
                                            two.first = i+k;//当前kmer所在的位置
                                            two.second = i+k+1;//当前kmer应该连接的kmer在哪个位置
                                            span_part[kmer] = two;

                                        }

                                        pair<size_t,size_t> two;
                                        two.first = i+extend_str.length()-g_kmer_length-1;
                                        two.second = final_k;
                                        span_part[extend_str.substr(extend_str.length()-g_kmer_length-1,g_kmer_length)]=two;
                                        nodes[i+extend_str.length()-g_kmer_length-1].insert(extend_str.substr(extend_str.length()-g_kmer_length-1,g_kmer_length));
                                        in_del_num++;

                                    } else if (type == 0 && i+extend_str.length() > nodes.size()-50 && extend_str.length()-g_kmer_length+i < nodes.size()){

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

        for (int i = 0; i < kmer_list.size(); ++i) {
            if (used_kmers_.find(kmer_list[i])==used_kmers_.end()){
                std::vector<kmer_occurence_pair_t> candidates;
                int candidate_coverage = kmer_map.get_reverse_candidates(kmer_list[i], candidates);
                if (candidates.empty()){
                    string str1 = *nodes[0].begin();
                    string str2 = intval_to_kmer(kmer_list[i],g_kmer_length);
                    int count = countDifferentCharacters(str1,str2);

                    if (count < 4){
                        first_in_kmer.push_back(intval_to_kmer(kmer_list[i],g_kmer_length));

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

    }
    void delete_error_kmer(KmerMap& kmerMap,std::vector<std::set<string>>& nodes){
//        cout << "Begin delete error kmer..." <<endl;
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
//        cout << "Begin add parents..." <<endl;
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
            string kmer = sequence.substr(j, g_kmer_length);
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
//        cout << "Begin loading add coverage..." << endl;
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
                    nodes[it].coverage = std::round(nodes[it].coverage * 1000) / 1000.0;
                }
                //检查分配比例之和是否为1
                correction_ratio(nodes,layer_nodes,i);
            }
        }
    }

    std::vector<Node> handle_nodes(std::vector<Node>& after_gra,std::map<size_t,std::set<size_t>>& layer_nodes){
        int size = after_gra.size();
        std::vector<Node> after_gra2;
        for (int i = after_gra.size(); i >= 0; --i) {
            if (after_gra[i].children.empty() && i < size - 150){

                if (after_gra[i].coverage < 0.5){
//                    cout << "delete single node coverage : " << after_gra[i].coverage <<endl;
                    for (int j = 0; j < after_gra[i].parents.size(); ++j) {
                        for (int k = 0; k < after_gra[after_gra[i].parents[j]].children.size(); ++k) {
                            if (after_gra[after_gra[i].parents[j]].children[k] == i){
                                after_gra[after_gra[i].parents[j]].children.erase(after_gra[after_gra[i].parents[j]].children.begin()+k);
                                break;
                            }
                        }
                    }
                    after_gra[i].parents.clear();
                } else{

                    int dif_init = g_kmer_length;
                    size_t d;
//                    cout << after_gra[i].node_layer <<endl;
                    size_t ls ;
                    bool ls_f = false;
                    for (auto it1 : layer_nodes[after_gra[i].node_layer]) {
                        if (it1!=i && after_gra[it1].children.size() > 0){
                            ls = it1;
                            ls_f = true;
                            break;
                        }
                    }
                    if (ls_f){
                        for (auto it : layer_nodes[after_gra[after_gra[ls].children[0]].node_layer]) {
                            string s1 = after_gra[i].sequence.substr(after_gra[i].sequence.length()-g_kmer_length+1);
                            string s2 = after_gra[it].sequence.substr(0,g_kmer_length-1);
                            int dif = countDifferences(s1,s2);
                            if (dif < dif_init){
                                dif_init = dif;
                                d = it;
                            }
                        }
                        after_gra[i].add_child(d);
                        after_gra[d].add_parent(i);
                    } else{
                        for (int j = 0; j < after_gra[i].parents.size(); ++j) {
                            for (int k = 0; k < after_gra[after_gra[i].parents[j]].children.size(); ++k) {
                                if (after_gra[after_gra[i].parents[j]].children[k] == i){
                                    after_gra[after_gra[i].parents[j]].children.erase(after_gra[after_gra[i].parents[j]].children.begin()+k);
                                    break;
                                }
                            }
                        }
                        after_gra[i].parents.clear();
                    }

                }

            }
        }
        for (int i = 0; i < after_gra.size(); ++i) {
            if (after_gra[i].parents.empty() && !after_gra[i].children.empty() && after_gra[i].node_layer!=0){
                if (after_gra[i].coverage<0.9 || after_gra[i].coverage == 1){
                    for (int j = 0; j < after_gra[i].children.size(); ++j) {
                        for (int k = 0; k < after_gra[after_gra[i].children[j]].parents.size(); ++k) {
                            if (after_gra[after_gra[i].children[j]].parents[k] == i){
                                after_gra[after_gra[i].children[j]].parents.erase(after_gra[after_gra[i].children[j]].parents.begin()+k);
                                break;
                            }
                        }
                    }
                    after_gra[i].children.clear();
                } else{


                    int dif_init = g_kmer_length;
                    size_t d;
                    size_t ls;
                    bool ls_f = false;
                    for (auto it1 : layer_nodes[after_gra[i].node_layer]) {
                        if (it1!=i && after_gra[it1].parents.size() > 0){
                            ls = it1;
                            ls_f = true;
                            break;
                        }
                    }
                    if (ls_f){
                        for (auto it : layer_nodes[after_gra[after_gra[ls].parents[0]].node_layer]) {
                            string s1 = after_gra[it].sequence.substr(1);
                            string s2 = after_gra[i].sequence.substr(0,g_kmer_length-1);
                            int dif = countDifferences(s1,s2);
                            if (dif < dif_init){
                                dif_init = dif;
                                d = it;
                            }
                        }
                        after_gra[d].add_child(i);
                        after_gra[i].add_parent(d);
                    } else{
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

            }
        }
        for (int i = 0; i < after_gra.size(); ++i) {
            if (!after_gra[i].parents.empty() || !after_gra[i].children.empty()){
                after_gra2.push_back(after_gra[i]);
            }
        }
        return after_gra2;
    }


    void adjust_layer(std::vector<Node>& after_gra2, std::map<size_t,std::set<size_t>>& layer_nodes) {
        std::map<size_t,std::set<size_t>> layer_nodes2;
        int lay_num = 0;
        for (int i = 0; i < after_gra2.size(); ++i) {
            layer_nodes2[after_gra2[i].node_layer].insert(after_gra2[i].id);
        }
        layer_nodes.clear();
        for (auto la : layer_nodes2) {
            layer_nodes[lay_num] = la.second;
            for (auto l2 : la.second){
                after_gra2[l2].node_layer = lay_num;
            }
            lay_num++;
        }
    }

    void dfs(std::vector<Node>& after_gra, size_t currentNode, std::vector<size_t>& path,
             std::unordered_set<size_t>& visited, std::unordered_map<std::vector<size_t>, size_t, VectorHash>& allPaths,
             int& maxLength) {
        size_t currentLength = 0;
        for (size_t node : path) {
            currentLength += after_gra[node].sequence.length()-g_kmer_length+1;
        }
        if (currentLength >= maxLength) {
            allPaths[path]++;
            return;
        }
        visited.insert(currentNode);
        path.push_back(currentNode);
        bool hasChildren = false;
        for (size_t neighbor : after_gra[currentNode].children) {
            if (visited.find(after_gra[neighbor].id) == visited.end()) {
                hasChildren = true;
                dfs(after_gra,after_gra[neighbor].id,path,visited,allPaths,maxLength);
            }
        }
        if (!hasChildren) {
            allPaths[path]++;
        }
        visited.erase(currentNode);
        path.pop_back();
    }
    std::unordered_map<std::vector<size_t>,size_t,VectorHash> findAllPathsWith50Nodes(std::vector<Node>& after_gra,
                                                                                      size_t& startNode,int& maxLength,int& m) {
        std::vector<size_t> path;
        std::unordered_set<size_t> visited;
        std::unordered_map<std::vector<size_t>, size_t, VectorHash> allPaths;

        dfs( after_gra, startNode, path, visited, allPaths, maxLength);

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

    void dfs2(std::vector<Node>& after_gra, size_t currentNode, std::vector<size_t>& path,
             std::unordered_set<size_t>& visited, std::unordered_map<std::vector<size_t>, size_t, VectorHash>& allPaths,
              std::set<size_t>& endNodes,int& maxLayer,int& d) {

        // 如果当前节点是结束节点之一，记录路径并返回
        if (endNodes.find(currentNode) != endNodes.end()) {
            allPaths[path]++;
            d = 1;
            return;
        }else if (after_gra[currentNode].node_layer >= maxLayer){
            allPaths[path]++;
            d = 2;
            return;
        }
        visited.insert(currentNode);
        path.push_back(currentNode);
        bool hasChildren = false;
        for (size_t neighbor : after_gra[currentNode].children) {
            if (visited.find(after_gra[neighbor].id) == visited.end()) {
                hasChildren = true;
                dfs2(after_gra,after_gra[neighbor].id,path,visited,allPaths,endNodes,maxLayer,d);
            }
        }
        if (!hasChildren) {
            allPaths[path]++;
        }
        visited.erase(currentNode);
        path.pop_back();
    }

    std::unordered_map<std::vector<size_t>,size_t,VectorHash> findAllPathsWith50Nodes2(std::vector<Node>& after_gra,std::map<size_t,std::set<size_t>>& layer_nodes,
                                                                                      size_t& startNode,std::set<size_t>& endNodes,int maxLayer) {
        std::vector<size_t> path;
        std::unordered_set<size_t> visited;
        std::unordered_map<std::vector<size_t>, size_t, VectorHash> allPaths;
        int d;
        dfs2(after_gra, startNode, path, visited, allPaths, endNodes,maxLayer,d);
//        cout << d <<endl;
        while(d == 2){
            allPaths.clear();
            maxLayer++;
            endNodes = layer_nodes[maxLayer];
            dfs2(after_gra, startNode, path, visited, allPaths, endNodes,maxLayer,d);
        }
        return allPaths;
    }


    //删除错误的小路径
    void delete_error_local_paths(KmerMap& kmerMap,std::vector<Node>& after_gra,vector<pair<vector<size_t>,double>>& start_local_paths
                                  ,vector<pair<vector<size_t>,double>>& after_local_paths ,std::vector<std::string>& data
                                  ,unordered_map<string,int>& used_string,std::unordered_map<std::vector<size_t>,string,VectorHash>& vec_string
                                  ,int& first_level,int& max_level,double& min_c,int& layer,size_t& node_num){
        unordered_map<std::vector<size_t>,size_t,VectorHash> del_points;
        unordered_map<string,int> path_str;
        for (int i = 0; i < start_local_paths.size(); i++) {
            vector <size_t> key_d = start_local_paths[i].first;
            string str_d = after_gra[key_d[0]].sequence;
            for (int idd = 1; idd < key_d.size(); ++idd) {//获取小路径的序列
                str_d = str_d + after_gra[key_d[idd]].sequence.substr(g_kmer_length - 1);
            }
            path_str[str_d]++;
        }
        for (int i = 0; i < start_local_paths.size(); i++) {
            vector<size_t> key_d = start_local_paths[i].first;
            string str_d = after_gra[key_d[0]].sequence;
            for (int idd = 1; idd < key_d.size(); ++idd) {//获取小路径的序列
                str_d = str_d + after_gra[key_d[idd]].sequence.substr(g_kmer_length-1);
            }
            size_t support_num = 0;
            vector<int> sup;
            unordered_map<string,int> one_diff;
            sup.resize(2);
            string kmer = str_d.substr(0,g_kmer_length);
            std::vector<pair<size_t, size_t>> reads = kmerMap.get_kmer_read(kmer);
            for (int r = 0; r < reads.size(); ++r) {
                string read = data[reads[r].first].substr(reads[r].second);
                if (read.length() > str_d.length()){
                    int diff = hammingDistance(read.substr(0,str_d.length()),str_d);
                    if (diff == 0){
                        support_num++;
                        sup[0]++;
                    } else{
                        sup[1]++;
                        if (path_str.find(read.substr(0,str_d.length())) == path_str.end()){
//                            cout << read.substr(0,str_d.length()) <<endl;
                            one_diff[read.substr(0,str_d.length())]++;
                        }

                    }
                }
                else if (read.length() > str_d.length()*0.9){
                    int diff = hammingDistance(read,str_d.substr(0,read.length()));
                    if (diff == 0){
                        support_num++;
//                        cout << i << " ";
                        sup[0]++;
                    }
                }
            }
//            cout << "" <<endl;
            vec_string[key_d] = str_d;
            used_string[str_d]++;
            int max_str_d = 0;
            string max_s;
            bool flag_str = false;
            if (sup[1] > sup[0] && sup[0] > 2){
                for (auto o : one_diff){
                    if (o.second > sup[0]*10 && used_string.find(o.first) == used_string.end() && o.second > max_str_d){
                        max_str_d = o.second;
                        max_s = o.first;
                        flag_str = true;

                    }
                }

            }
            if (flag_str){
                used_string[max_s]++;
                vec_string[key_d] = max_s;
//                cout << sup[1] << "  " << sup[0] <<endl;
            }


                start_local_paths[i].second = support_num;


        }
        std::sort(start_local_paths.begin(), start_local_paths.end(),[](const auto& a, const auto& b) {
            return a.second > b.second; // 按值的大小升序排列
        });
        double sum_path = 0;
        for (int i = 0; i < start_local_paths.size(); i++) {

                sum_path = sum_path + start_local_paths[i].second;

        }

        unordered_map<size_t,size_t> used_local;
        for (int i = 0; i < start_local_paths.size(); ++i) {
            for (int j = 0; j < start_local_paths[i].first.size(); ++j) {
                used_local[start_local_paths[i].first[j]]++;
            }
        }

        if(sum_path > 0){
            if(min_c > 0.01){
                for (int i = start_local_paths.size()-1; i >= 0; --i) {
                    double cov_path = start_local_paths[i].second/sum_path;
// || start_local_paths[i].second < 2  /2
                    if (cov_path < min_c*0.9){
                        bool used_flag = true;
                        for (int d = 0; d < start_local_paths[i].first.size(); ++d) {
                            used_local[start_local_paths[i].first[d]]--;
                            if (used_local[start_local_paths[i].first[d]] == 0 &&
                            (cov_path > min_c*0.5 || after_gra[start_local_paths[i].first[d]].node_layer == max_level
                            || after_gra[start_local_paths[i].first[d]].node_layer == first_level)){
//                              cout << start_local_paths[i].first[d] << " : " << start_local_paths[i].second << "  " << after_gra[start_local_paths[i].first[d]].node_layer << "   " << max_level << endl ;
                                used_flag = false;
                            }
                        }
                        if (used_flag){
                            del_points[start_local_paths[i].first]++;
                        }

                    }

                }
            } else{
                for (int i = start_local_paths.size()-1; i >= 0; --i) {
                    double cov_path = start_local_paths[i].second/sum_path;

                    if (cov_path < 0.05){
                        bool used_flag = true;
                        for (int d = 0; d < start_local_paths[i].first.size(); ++d) {
                            used_local[start_local_paths[i].first[d]]--;
                            if (used_local[start_local_paths[i].first[d]] == 0 && cov_path > 0.01){
                                used_flag = false;
                            }
                        }
                        if (used_flag){
                            del_points[start_local_paths[i].first]++;
                        }

                    }
                }
            }

        }
        for (int i = 0; i < start_local_paths.size(); ++i) {
            if(del_points.find(start_local_paths[i].first) == del_points.end()){
                after_local_paths.push_back(start_local_paths[i]);
            }
        }
        sum_path = 0;
        for (int i = 0; i < after_local_paths.size(); ++i) {
            sum_path = sum_path + after_local_paths[i].second;
        }
        if(sum_path != 0){
            for (int i = 0; i < after_local_paths.size(); ++i) {
                after_local_paths[i].second = after_local_paths[i].second/sum_path;
            }
        }

        if (after_local_paths.size() == 0){
//            cout << layer << " : " << node_num << endl;
            if(start_local_paths.size() < 5){
                after_local_paths = start_local_paths;
            } else{
                for (int i = 0; i < 5; ++i) {
                    after_local_paths.push_back(start_local_paths[i]);
                }
            }
        }

    }


    vector<vector<size_t>> compression_nodes(KmerMap& kmerMap, std::vector<Node>& after_gra,std::map<size_t,std::set<size_t>>& layer_nodes,
                                             vector<vector<size_t>>& sim_layer,vector<size_t>& layer_sim_node,std::vector<std::string>& data,
                                             int& limit_length,std::unordered_map<std::vector<size_t>,string,VectorHash>& vec_string,
                                             vector<double>& node_cov,double& min_c){
        cout << "Begin compression nodes..." <<endl;
        vector<size_t> start;
        vector<vector<size_t>> simplified_points;
        for (int i = 0; i < after_gra.size(); ++i) {
            if (after_gra[i].parents.empty()){
                start.push_back(i);
            }
        }
        int m = 0;

        int before_sum = 0;
        unordered_map<string,int> used_string;
        size_t sim_m = 0;
        while(1){
            m++;
            set<size_t> new_start;
            bool stop_flag = false;
            vector<size_t> sim_layer_node;
            int maxLength = limit_length*0.4;

            string kmer = after_gra[start[0]].sequence.substr(0,g_kmer_length);
            vector<pair<size_t,size_t>> reads_data = kmerMap.get_kmer_read(kmer);
            int local_le_sum = 250;
            int local_num = 0;
            for (int j = 0; j < reads_data.size(); ++j) {
                if (reads_data[j].second == 0){
                    if (data[reads_data[j].first].length() < local_le_sum){
                        local_le_sum = data[reads_data[j].first].length();
                    }

                    local_num++;
                }
            }

            if (local_num != 0){
                maxLength = local_le_sum * 0.4;
            }
            vector<pair<vector<size_t>,double>> start_local_paths;
            for (int i = 0; i < start.size(); ++i) {

                std::unordered_map<std::vector<size_t>,size_t,VectorHash> allPaths = findAllPathsWith50Nodes(after_gra,start[i],maxLength,i);
                for (auto path : allPaths) {
                    std::vector <size_t> key = path.first;
                    pair<vector<size_t>,size_t> one_path;
                    one_path.first = key;
                    one_path.second = path.second;
                    start_local_paths.push_back(one_path);
                    if (after_gra[key[key.size()-1]].children.empty()){
                        stop_flag = true;
                        continue;
                    }

                    for (int j = 0; j < after_gra[key[key.size()-1]].children.size(); ++j) {
                        new_start.insert(after_gra[after_gra[key[key.size()-1]].children[j]].id);
                    }
                }

            }



            bool local_flag = true;
            int local_layer;
            int max_level;
            vector<size_t> local_start;
            if (!stop_flag){

                for(auto st : new_start) {
                    local_start.push_back(st);
                }
                local_layer = after_gra[local_start[0]].node_layer;
                max_level = after_gra[local_start[0]].node_layer;
                for (int j = 1; j < local_start.size(); ++j) {
                    if (after_gra[local_start[0]].node_layer != after_gra[local_start[j]].node_layer){
                        if (after_gra[local_start[j]].node_layer < local_layer){
                            local_layer = after_gra[local_start[j]].node_layer;
                        }
                        if (after_gra[local_start[j]].node_layer > max_level){
                            max_level = after_gra[local_start[j]].node_layer;
                        }
                        local_flag = false;
                    }

                }
                new_start.clear();
            }

            //判断新开始节点是否在同一层中
            if (!local_flag){//路径结束位置不在同一层


                for (int i = 0; i < start.size(); ++i) {//将当前节点作为compressed graph中的一个点，下个开始节点从下一层开始
                    for (int j = 0; j < after_gra[start[i]].children.size(); ++j) {
                        new_start.insert(after_gra[after_gra[start[i]].children[j]].id);
                    }
                }
                bool local_flag2 = true;
                vector<size_t> local_start2;
                for(auto st : new_start) {
                    local_start2.push_back(st);
                }
                size_t min_layer2 = after_gra[local_start2[0]].node_layer;
                size_t max_layer2 = 0;
                for (int j = 0; j < local_start2.size(); ++j) {

                    if (after_gra[local_start2[j]].node_layer < min_layer2
                        && (after_gra[local_start2[j]].parents.size()!=0 || after_gra[local_start2[j]].children.size()!=0)){
                        min_layer2 = after_gra[local_start2[j]].node_layer;
                    }
                    if (after_gra[local_start2[j]].node_layer > max_layer2
                        && (after_gra[local_start2[j]].parents.size()!=0 || after_gra[local_start2[j]].children.size()!=0)){

                        max_layer2 = after_gra[local_start2[j]].node_layer;
                    }


                }

                if (min_layer2 != max_layer2){
                    local_flag2 = false;
                }


                new_start.clear();
                for (auto l : layer_nodes[min_layer2]) {
                    new_start.insert(after_gra[l].id);
                }
                local_start2.clear();
                for(auto st : new_start) {//更新开始节点
                    local_start2.push_back(st);
                }

                std::set<size_t> endNodes;

                endNodes = layer_nodes[max_level+1];
                vector<pair<vector<size_t>,double>> reserve_paths;
                double resever_sum = 0;
                start_local_paths.clear();
                for (int i = 0; i < start.size(); ++i) {

                    std::unordered_map<std::vector<size_t>,size_t,VectorHash> allPaths = findAllPathsWith50Nodes2(after_gra,layer_nodes,start[i],endNodes,max_level+1);
                    for (auto path : allPaths) {
                        std::vector <size_t> key = path.first;
                        pair<vector<size_t>,size_t> one_path;
                        one_path.first = key;
                        one_path.second = path.second;
                        start_local_paths.push_back(one_path);
                    }
                }
                vector<pair<vector<size_t>,double>> after_local_paths;
                delete_error_local_paths(kmerMap,after_gra,start_local_paths,after_local_paths,data,used_string,vec_string,after_gra[start[0]].node_layer,max_level,min_c,m,sim_m);
                for (int i = 0; i < after_local_paths.size(); ++i) {

                    vector <size_t> points;
                    std::vector <size_t> key = after_local_paths[i].first;

                    for (size_t node : key) {
                        points.push_back(node);
                    }
                    sim_layer_node.push_back(sim_m);
                    layer_sim_node.push_back(m);
                    sim_m++;
                    node_cov.push_back(after_local_paths[i].second);
                    simplified_points.push_back(points);
                    if (after_gra[key[key.size()-1]].children.empty()){
                        stop_flag = true;
                        continue;
                    }

                }

                new_start.clear();
//                        cout << "max_level : " << m << "  " << sim_m <<endl;
                int next_level = max_level+1;
                for (auto l : layer_nodes[next_level]) {
                    new_start.insert(after_gra[l].id);
                }


            } else{

                vector<pair<vector<size_t>,double>> after_local_paths;
                delete_error_local_paths(kmerMap,after_gra,start_local_paths,after_local_paths,data,used_string,vec_string,after_gra[start[0]].node_layer,max_level,min_c,m,sim_m);
                vector<pair<vector<size_t>,double>> start_local_paths2 = start_local_paths;
                vector<pair<vector<size_t>,double>> after_local_paths2 = after_local_paths;
                while (after_local_paths.size() > 20){

                    std::set<size_t> endNodes;
                    after_local_paths.clear();
                    max_level = max_level - 1;
                    if (max_level <= after_gra[start[0]].node_layer+1){
                        start_local_paths.clear();
                        start_local_paths = start_local_paths2;
                        after_local_paths = after_local_paths2;
                        break;
                    }
                    endNodes = layer_nodes[max_level];
                    start_local_paths.clear();
                    for (int i = 0; i < start.size(); ++i) {

                        std::unordered_map<std::vector<size_t>,size_t,VectorHash> allPaths = findAllPathsWith50Nodes2(after_gra,layer_nodes,start[i],endNodes,max_level);
                        for (auto path : allPaths) {
                            std::vector <size_t> key = path.first;
                            pair<vector<size_t>,size_t> one_path;
                            one_path.first = key;
                            one_path.second = path.second;
                            start_local_paths.push_back(one_path);

                            if (after_gra[key[key.size()-1]].children.empty()){
                                stop_flag = true;
                                continue;
                            }

                            for (int j = 0; j < after_gra[key[key.size()-1]].children.size(); ++j) {
                                new_start.insert(after_gra[after_gra[key[key.size()-1]].children[j]].id);
                            }
                        }
                    }

                    if (!stop_flag){
                        local_start.clear();
                        for(auto st : new_start) {
                            local_start.push_back(st);
                        }
                        local_layer = after_gra[local_start[0]].node_layer;
                        max_level = after_gra[local_start[0]].node_layer;
                        for (int j = 1; j < local_start.size(); ++j) {
                            if (after_gra[local_start[0]].node_layer != after_gra[local_start[j]].node_layer){
                                if (after_gra[local_start[j]].node_layer < local_layer){
                                    local_layer = after_gra[local_start[j]].node_layer;
                                }
                                if (after_gra[local_start[j]].node_layer > max_level){
                                    max_level = after_gra[local_start[j]].node_layer;
                                }
                                local_flag = false;
                            }

                        }
                        new_start.clear();

                    }
                    if(!local_flag){
                        start_local_paths.clear();
                        start_local_paths = start_local_paths2;
                        after_local_paths = after_local_paths2;

                        break;
                    }

                    delete_error_local_paths(kmerMap,after_gra,start_local_paths,after_local_paths,data,used_string,vec_string,after_gra[start[0]].node_layer,max_level,min_c,m,sim_m);
                    for (int i = 0; i < after_local_paths.size(); ++i) {

                        vector <size_t> points;
                        std::vector <size_t> key = after_local_paths[i].first;

                        for (size_t node : key) {
                            points.push_back(node);
                        }
                        sim_layer_node.push_back(sim_m);
                        layer_sim_node.push_back(m);
                        sim_m++;
                        node_cov.push_back(after_local_paths[i].second);
                        simplified_points.push_back(points);
                        if (after_gra[key[key.size()-1]].children.empty()){
                            stop_flag = true;
                            continue;
                        }

                    }

                    new_start.clear();
//                        cout << "max_level : " << m << "  " << sim_m <<endl;
                    int next_level = max_level;
                    for (auto l : layer_nodes[next_level]) {
                        new_start.insert(after_gra[l].id);
                    }

                }
                for (int i = 0; i < after_local_paths.size(); ++i) {
                    vector<size_t> points;
                    std::vector <size_t> key = after_local_paths[i].first;
                    for (size_t node : key) {
                        points.push_back(node);
                    }

                    sim_layer_node.push_back(sim_m);
                    layer_sim_node.push_back(m);

                    sim_m++;
                    node_cov.push_back(after_local_paths[i].second);
                    simplified_points.push_back(points);
                    if (after_gra[key[key.size()-1]].children.empty()){
                        stop_flag = true;
                        continue;
                    }
                    for (int j = 0; j < after_gra[key[key.size()-1]].children.size(); ++j) {
                        new_start.insert(after_gra[after_gra[key[key.size()-1]].children[j]].id);
                    }
                }

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


    void get_edge_ratio(vector<vector<size_t>>& sim_layer,vector<Node>& sim_graph1){
        for (int i = 0; i < sim_layer.size(); ++i) {

            double level_edge_sum = 0;
            for (int j = 0; j < sim_layer[i].size(); ++j) {
                for (int k = 0; k < sim_graph1[sim_layer[i][j]].edge_coverage.size(); ++k) {
                   level_edge_sum = level_edge_sum + sim_graph1[sim_layer[i][j]].edge_coverage[k].second;
                }
            }
            for (int j = 0; j < sim_layer[i].size(); ++j) {
                for (int k = 0; k < sim_graph1[sim_layer[i][j]].edge_coverage.size(); ++k) {
                    sim_graph1[sim_layer[i][j]].edge_coverage[k].second = sim_graph1[sim_layer[i][j]].edge_coverage[k].second/level_edge_sum;
                }
            }
        }

    }

    void create_sim_graph(KmerMap& kmerMap,std::vector<Node>& after_gra,vector<vector<size_t>>& sim_gra
                          ,vector<Node>& sim_graph1,std::map<size_t,std::set<size_t>>& layer_nodes
                          ,vector<vector<size_t>>& sim_layer,vector<size_t>& layer_sim_node
                          ,vector<vector<pair<size_t,size_t>>>& parent_different_node,std::vector<std::string>& data
                          ,int& min_l,std::unordered_map<std::vector<size_t>,string,VectorHash>& vec_string
                          ,vector<double>& node_cov, vector<vector<int>>& why_connect, double& min_cov){
        cout << "Begin create sim graph..." <<endl;
        sim_graph1.resize(sim_gra.size());
        why_connect.resize(sim_gra.size());

        for (int i = 0; i < sim_gra.size(); ++i) {
            sim_graph1[i].id = i;
            sim_graph1[i].node_layer = layer_sim_node[i];
            sim_graph1[i].single_node_id.push_back(i);
            sim_graph1[i].c_node_id = sim_gra[i];
            string str = after_gra[sim_gra[i][0]].sequence;
            if (vec_string.find(sim_gra[i])!=vec_string.end()){
                sim_graph1[i].sequence = vec_string[sim_gra[i]];
            } else{
                for (int j = 1; j < sim_gra[i].size(); ++j) {
                    str = str + after_gra[sim_gra[i][j]].sequence.substr(g_kmer_length-1);
                }
                sim_graph1[i].sequence = str;
            }

        }


        parent_different_node.resize(sim_graph1.size());
        alone_possess_nodes(sim_layer,parent_different_node,sim_graph1,after_gra);
        vector<vector<pair<size_t,size_t>>> children_different_node;

        //简化图中点的丰度 || sim_graph1[i].node_layer == 1
        for (int i = 0; i < sim_graph1.size(); ++i) {
            double min_alone = 1;

            if (parent_different_node[i].size() != 0
            && after_gra[sim_graph1[i].c_node_id[0]].parents.size() != 0){


                for (int j = 0; j < parent_different_node[i].size(); ++j) {

                    if (after_gra[parent_different_node[i][j].second].coverage < min_alone){
                        min_alone = after_gra[parent_different_node[i][j].second].coverage;
                    }
                }
                sim_graph1[i].coverage = min_alone;
            } else{
                sim_graph1[i].coverage = std::round(node_cov[i]*1000)/1000.0;
            }

        }

        for (int i = 0; i < sim_graph1.size(); ++i) {
            if (parent_different_node[i].size() == 0){

                if(sim_graph1[i].node_layer < sim_layer.size()-1){
                    for (int j = 0; j < sim_layer[sim_graph1[i].node_layer].size(); ++j) {

                        sim_graph1[sim_layer[sim_graph1[i].node_layer][j]].coverage =
                                std::round(node_cov[sim_layer[sim_graph1[i].node_layer][j]]*1000)/1000.0;
                    }
                }else{

                    sim_graph1[i].coverage = std::round(node_cov[i]*1000)/1000.0;
                }
            }
        }


        for (size_t i = 0; i < sim_gra.size(); ++i) {//给点之间添加关系

            vector<size_t> local_children;
            for (size_t j = i + 1; j < sim_gra.size(); ++j) {//确定当前节点所有可能连接的点
                for (int k = 0; k < after_gra[sim_gra[i][sim_gra[i].size() - 1]].children.size(); ++k) {
                    if (after_gra[after_gra[sim_gra[i][sim_gra[i].size() - 1]].children[k]].id == sim_gra[j][0]){
                        local_children.push_back(j);
                        break;
                    }
                }
            }

            if (local_children.size() == 1){
                pair<size_t,double> edge;
                edge.first = local_children[0];
                edge.second = 1;
                sim_graph1[i].edge_coverage.push_back(edge);
                sim_graph1[i].add_child(local_children[0]);
                sim_graph1[local_children[0]].add_parent(i);

                why_connect[i].push_back(1);
            } else if(local_children.size() > 1) {
                bool connect = false;
                children_different_node.clear();
                children_different_node.resize(sim_graph1.size());
                difference_vec(local_children, children_different_node, sim_graph1);
                vector<int> local_children_sup;
                local_children_sup.resize(local_children.size());

                for (int l = 0; l < local_children.size(); ++l) {
                    int support = 0;
                    string str = sim_graph1[i].sequence + sim_graph1[local_children[l]].sequence.substr(g_kmer_length - 1);
                    string kmer = str.substr(0,g_kmer_length);
                    std::vector<pair<size_t, size_t>> reads = kmerMap.get_kmer_read(kmer);
                    for (int r = 0; r < reads.size(); ++r) {
                        string read = data[reads[r].first].substr(reads[r].second);
                        if (read.length() >= str.length()){
                            if (read.substr(0,str.length()) == str){
                                support++;
                            }
                        }
                    }

                    if(support == 0 && str.length() > min_l*0.9 && min_l > g_kmer_length*4){
                        int odd = 1;
                        while(min_l*0.9 <= str.length()){
                            if (odd % 2 == 0){
                                str = str.substr(1);
                            } else{
                                str = str.substr(0,str.length()-1);
                            }

                        }
                        kmer = str.substr(0,g_kmer_length);
                        reads = kmerMap.get_kmer_read(kmer);
                        for (int r = 0; r < reads.size(); ++r) {
                            string read = data[reads[r].first].substr(reads[r].second);
                            if (read.length() >= str.length()){
                                if (read.substr(0,str.length()) == str){
                                    support++;
                                }
                            }
                        }

                    }

                    local_children_sup[l] = support;

                }
                double su_sum = 0;
                for(int su = 0; su < local_children_sup.size(); su++){
                    su_sum = su_sum + local_children_sup[su];
                }
                for (int l = 0; l < local_children_sup.size(); ++l) {
                    if (su_sum != 0){
                        double d_sup = local_children_sup[l]/su_sum;
//
                        if (d_sup >= 0.08 &&((local_children_sup[l] > 1 && su_sum != 1) || (local_children_sup[l] == 1 && su_sum == 1))){
                            connect = true;
                            pair <size_t,double> edge;
                            edge.first = local_children[l];
                            edge.second = local_children_sup[l];
                            sim_graph1[i].edge_coverage.push_back(edge);
                            sim_graph1[i].add_child(local_children[l]);
                            sim_graph1[local_children[l]].add_parent(i);

                            why_connect[i].push_back(2);
                        }
                    }
                }

            }
        }
        int no_child_num = 0;
        for (int i = 0; i < sim_graph1.size(); ++i) {
            if (sim_graph1[i].children.size() == 0){
                no_child_num++;
            }
        }
//        cout << no_child_num << "  " << sim_graph1.size() <<endl;

        for (size_t i = 0; i < sim_gra.size(); ++i) {//给点之间添加关系
            if (sim_graph1[i].children.size() == 0){

                vector < size_t > local_children;

                for (size_t j = i + 1; j < sim_gra.size(); ++j) {//确定当前节点所有可能连接的点
                    for (int k = 0; k < after_gra[sim_gra[i][sim_gra[i].size() - 1]].children.size(); ++k) {
                        if (after_gra[after_gra[sim_gra[i][sim_gra[i].size() - 1]].children[k]].id == sim_gra[j][0]) {
                            local_children.push_back(j);
                            break;
                        }
                    }
                }

                int select_i = 10;
                if (local_children.size() > 1) {
                    double cha_c = 1;

                    for (int j = 0; j < sim_graph1[i].c_node_id.size(); ++j) {

                        if (after_gra[sim_graph1[i].c_node_id[j]].coverage < 1) {
                            for (int l = 0; l < local_children.size(); ++l) {
                                for (int sg = 0; sg < sim_graph1[local_children[l]].c_node_id.size(); ++sg) {
                                    if (after_gra[sim_graph1[local_children[l]].c_node_id[sg]].coverage < 1
                                        && abs(after_gra[sim_graph1[local_children[l]].c_node_id[sg]].coverage
                                               - after_gra[sim_graph1[i].c_node_id[j]].coverage) < cha_c) {
                                        cha_c = abs(after_gra[sim_graph1[local_children[l]].c_node_id[sg]].coverage
                                                    - after_gra[sim_graph1[i].c_node_id[j]].coverage);
                                        select_i = local_children[l];
                                    }
                                }
                            }
                        }
                    }
                }
                if (select_i != 10){
                    pair <size_t,double> edge;
                    edge.first = select_i;
                    edge.second = 0;
                    sim_graph1[i].edge_coverage.push_back(edge);
                    sim_graph1[i].add_child(select_i);
                    sim_graph1[select_i].add_parent(i);

                } else{
                    for (int l = 0; l < local_children.size(); ++l) {
                        pair <size_t,double> edge;
                        edge.first = local_children[l];
                        edge.second = 0;
                        sim_graph1[i].edge_coverage.push_back(edge);
                        sim_graph1[i].add_child(local_children[l]);
                        sim_graph1[local_children[l]].add_parent(i);


                        why_connect[i].push_back(4);
                    }
                }

            }

        }
        for (int i = 1; i < sim_graph1.size(); ++i) {
            if (sim_graph1[i].parents.size() == 0 && after_gra[sim_gra[i][0]].parents.size()!=0){
                vector<size_t> local_parents;
//            }
                double c_sum_c = 0;
                for (int j = i - 1; j > 0; --j) {//确定当前节点所有可能连接的点
                    for (int k = 0; k < after_gra[sim_gra[i][0]].parents.size(); ++k) {
                        if (after_gra[after_gra[sim_gra[i][0]].parents[k]].id == sim_gra[j][sim_gra[j].size()-1]){
                            local_parents.push_back(j);
                            c_sum_c = c_sum_c + sim_graph1[j].coverage;
                            break;
                        }
                    }
                }

                if (local_parents.size() == 0){
                    continue;
                }
                if (local_parents.size() == 1){
                    pair <size_t, size_t> edge;
                    edge.first = i;
                    edge.second = 0;
                    sim_graph1[local_parents[0]].edge_coverage.push_back(edge);

                    sim_graph1[local_parents[0]].add_child(i);
                    sim_graph1[i].add_parent(local_parents[0]);
                } else{
                    double cha_c = 1;
                    int select_i = 10;
                    double cha_min = 1;
                    int cha_l = 10;
                    for (int l = 0; l < local_parents.size(); ++l) {
                        if (abs(sim_graph1[i].coverage - sim_graph1[local_parents[l]].coverage) < cha_min){
                            cha_min = abs(sim_graph1[i].coverage - sim_graph1[local_parents[l]].coverage);
                            cha_l = l;
                        }
                    }


                    if (cha_l != 10 && cha_min < min_cov){

                        pair <size_t, size_t> edge;
                        edge.first = i;
                        edge.second = 0;
                        sim_graph1[local_parents[cha_l]].edge_coverage.push_back(edge);

                        sim_graph1[local_parents[cha_l]].add_child(i);
                        sim_graph1[i].add_parent(local_parents[cha_l]);

                    } else{
                        for (int l = 0; l < local_parents.size(); ++l) {
                            pair <size_t, size_t> edge;
                            edge.first = i;
                            edge.second = 0;
                            sim_graph1[local_parents[l]].edge_coverage.push_back(edge);

                            sim_graph1[local_parents[l]].add_child(i);
                            sim_graph1[i].add_parent(local_parents[l]);

                        }
                    }
                }
            }

        }

        //给特殊点添加边
        for (int i = 0; i < sim_graph1.size(); ++i) {

            if (sim_graph1[i].children.empty() && !after_gra[sim_graph1[i].c_node_id[sim_graph1[i].c_node_id.size()-1]].children.empty()){

                bool em_child_flag = false;
                vector<int> em_child_i;
                vector<int> em_child_i_pos;
                for (int cn = 0; cn < after_gra[sim_graph1[i].c_node_id[sim_graph1[i].c_node_id.size()-1]].children.size(); ++cn) {
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
//                            cout << i << " :::::::::: " << em_child_i_pos[ci] <<endl;
                            if (sim_graph1[i].c_node_id.size() > em_child_i_pos[ci]+1){
                                sim_graph1[i].c_node_id.erase(sim_graph1[i].c_node_id.begin()+sim_graph1[i].c_node_id.size()-em_child_i_pos[ci]-1
                                        ,sim_graph1[i].c_node_id.end());
                                string spe_str = after_gra[sim_graph1[i].c_node_id[0]].sequence;
                                for (int ct = 1; ct < sim_graph1[i].c_node_id.size(); ++ct) {
                                    spe_str += after_gra[sim_graph1[i].c_node_id[ct]].sequence.substr(g_kmer_length-1);
                                }
                                sim_graph1[i].sequence = spe_str;
                                pair <size_t, size_t> edge;
                                edge.first = em_child_i[ci];
                                edge.second = 0;
                                sim_graph1[i].edge_coverage.push_back(edge);
                                sim_graph1[i].add_child(em_child_i[ci]);
                                sim_graph1[em_child_i[ci]].add_parent(i);

//                                cout << i << " -- " << em_child_i[ci] << " 6" <<endl;
                                why_connect[i].push_back(6);
                            }
                        }
                    }
                }
            }
        }
        get_edge_ratio(sim_layer,sim_graph1);
        cout << "Bubble graph nodes : " << sim_graph1.size() << endl;

        add_in_zero(sim_layer,sim_graph1);

    }

    void add_in_zero(vector<vector<size_t>>& sim_layer,vector<Node>& sim_graph1){
        for (int i = 1; i < sim_layer.size()-1; ++i) {
            for (int j = 0; j < sim_layer[i].size(); ++j) {
                if (sim_graph1[sim_layer[i][j]].parents.empty()){
                    if (sim_graph1[sim_layer[i][j]].parents.empty()){
                        double cha_min = 1;
                        int cha_l;
//                        cout << sim_layer[i][j] << " : " << sim_layer[i-1].size() <<endl;
                        for (int k = 0; k < sim_layer[i-1].size(); ++k) {
                            if (abs(sim_graph1[sim_layer[i-1][k]].coverage-sim_graph1[sim_layer[i][j]].coverage) < cha_min){
                                cha_min = abs(sim_graph1[sim_layer[i-1][k]].coverage-sim_graph1[sim_layer[i][j]].coverage);
                                cha_l = sim_layer[i-1][k];
                            }
                        }
                        sim_graph1[cha_l].add_child(sim_layer[i][j]);
                        sim_graph1[sim_layer[i][j]].add_parent(cha_l);
                        pair <size_t, size_t> edge;
                        edge.first = sim_layer[i][j];
                        edge.second = 0;
                        sim_graph1[cha_l].edge_coverage.push_back(edge);
                    }

                }
            }
        }
        for (int i = sim_layer.size()-3; i > 1; --i) {
            for (int j = 0; j < sim_layer[i].size(); ++j) {
                if (sim_graph1[sim_layer[i][j]].children.empty()){
                    double cha_min = 1;
                    int cha_l;
//                    cout << sim_layer[i][j] << " : " << sim_layer[i+1].size() <<endl;
                    for (int k = 0; k < sim_layer[i+1].size(); ++k) {
                        if (abs(sim_graph1[sim_layer[i+1][k]].coverage-sim_graph1[sim_layer[i][j]].coverage) < cha_min){
                            cha_min = abs(sim_graph1[sim_layer[i+1][k]].coverage-sim_graph1[sim_layer[i][j]].coverage);
                            cha_l = sim_layer[i+1][k];
                        }
                    }
                    sim_graph1[cha_l].add_parent(sim_layer[i][j]);
                    sim_graph1[sim_layer[i][j]].add_child(cha_l);
                    pair <size_t, size_t> edge;
                    edge.first = cha_l;
                    edge.second = 0;
                    sim_graph1[sim_layer[i][j]].edge_coverage.push_back(edge);
                }
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

    void determine_nodes_contigs(vector<Node>& sim_graph1,vector<Node>& sim_graph2,vector<vector<size_t>>& sim_layer
                                 ,vector<vector<int>>& why_connect){

        for (int i = 0; i < sim_graph2.size(); ++i) {
            Node node = sim_graph2[i];

            while(sim_graph2[i].children.size() == 1 && sim_graph2[sim_graph2[i].children[0]].parents.size() == 1){
                for (int j = 0; j < sim_graph2[sim_graph2[i].children[0]].children.size(); ++j) {
                    for (int k = 0; k < sim_graph2[sim_graph2[sim_graph2[i].children[0]].children[j]].parents.size(); ++k) {
                        if (sim_graph2[sim_graph2[sim_graph2[i].children[0]].children[j]].parents[k] == sim_graph2[i].children[0]){

                            sim_graph2[sim_graph2[sim_graph2[i].children[0]].children[j]].parents[k] = sim_graph2[i].id;
                        }
                    }
                }
                size_t nc = sim_graph2[i].children[0];
                sim_graph2[i].edge_coverage.clear();
                for (int j = 0; j < sim_graph2[sim_graph2[i].children[0]].children.size(); ++j) {
                    pair<size_t,double> e;
                    e.first = sim_graph2[sim_graph2[i].children[0]].children[j];
                    e.second = 0;
                    sim_graph2[i].edge_coverage.push_back(e);
                }
                sim_graph2[i].sequence = sim_graph2[i].sequence + sim_graph2[sim_graph2[i].children[0]].sequence.substr(g_kmer_length-1);
                sim_graph2[i].children = sim_graph2[sim_graph2[i].children[0]].children;

                why_connect[i].clear();
                for (int j = 0; j < why_connect[sim_graph2[i].children[0]].size(); ++j) {
                    why_connect[i].push_back(why_connect[sim_graph2[i].children[0]][j]);
                }
//                why_connect[i] = why_connect[sim_graph2[i].children[0]];

                sim_graph2[i].single_node_id.push_back(nc);

                sim_graph2[nc].children.clear();
                sim_graph2[nc].parents.clear();
            }

        }
        for (int i = 0; i < sim_graph2.size(); ++i) {
            if (sim_graph2[i].single_node_id.size() > 1){
                double sum = 0;
                for (int j = 0; j < sim_graph2[i].single_node_id.size(); ++j) {
                    sum = sum + sim_graph1[sim_graph2[i].single_node_id[j]].coverage;
                }
                sim_graph2[i].coverage = sum / sim_graph2[i].single_node_id.size();
            }

        }

    }

    vector<pair<vector<size_t>,double>> get_contigs(vector<Node>& sim_graph2
                                                    ,std::vector<std::map<size_t,size_t>>& same_position_sim
                                                    ,map<size_t,size_t>& use_node,vector<vector<size_t>>& node_reason){
        cout << "Begin find paths ... " << endl;

        vector<size_t> node_r;
        node_r.resize(3000);

        vector<pair<vector<size_t>,double>> paths;
        int m = 0;
        double sum_cov = 0;
        vector<Node> sim_graph = sim_graph2;
//        map<size_t,size_t> use_node;
        int start_num = 0;
        double start_aver = 1;
        for (int i = 0; i < sim_graph2.size(); ++i) {
            if (sim_graph2[i].parents.size() == 0 && sim_graph2[i].children.size()!= 0){
                start_num++;
//                start_aver = start_aver + sim_graph1[i].coverage;
//                cout << i << " : " << sim_graph2[i].coverage << "  " <<endl;
            }
            if (sim_graph2[i].parents.size() == 0 && sim_graph2[i].children.size()!= 0 && sim_graph2[i].coverage > 0 && start_aver > sim_graph2[i].coverage){
                start_aver = sim_graph2[i].coverage;
//                cout << sim_graph2[i].coverage <<endl;
            }
        }
        double limit_cov;
        if (start_aver > 0.02){
            limit_cov = 0.02;
        } else{
            limit_cov = start_aver;
        }

        for (int i = 0; i < sim_graph2.size(); ++i) {
            double min_start = 1;
            int min_start_i;
            bool start_flag = false;

                for (int j = 0; j < sim_graph2.size(); ++j) {
                    if (sim_graph2[j].parents.size() == 0 && sim_graph2[j].children.size()!= 0
                    && sim_graph2[j].coverage > limit_cov && min_start > sim_graph2[j].coverage ){//从丰度最小的节点开始
                        min_start = sim_graph2[j].coverage;
                        min_start_i = j;
                        start_flag = true;
                    }
                }

            if (!start_flag){
                continue;
            }


            while(sim_graph2[min_start_i].coverage > limit_cov){
                vector<double> path_node_coverage;

                if (sim_graph2[min_start_i].coverage < 0.01 || use_node[min_start_i] >= 2){
                    break;
                }

                if (abs(sum_cov - 1) < 0.02 || sum_cov > 1.02){
                    break;
                }

                m++;
                Node node = sim_graph2[min_start_i];
                vector<size_t> path;
                node_r.clear();
                node_r.resize(3000);

                path.push_back(node.id);
                use_node[node.id]++;
                double path_cov = 10;
                while (node.children.size()>0){
                    if (node.children.size() == 1){
                        node_r[node.id] = 1;
                        node = sim_graph2[node.children[0]];
                        path.push_back(node.id);


                        use_node[node.id]++;
                        if (node.coverage > 0){
                            path_node_coverage.push_back(node.coverage);
                        }
                    } else{
                        bool p_flag = false;
                        node_idx_t pair_i;
                        int p_i;
                        bool find_flag = false;
                        vector<int> cov_num;
                        double child_min_cov = 1;
                        for (int cmc = 0; cmc < node.children.size(); ++cmc) {
                            if (sim_graph2[node.children[cmc]].coverage < child_min_cov){
                                child_min_cov = sim_graph2[node.children[cmc]].coverage;
                            }
                        }
                        for (int k = path.size() - 1; k >= 1; --k) {//根据点的丰度判断
                            if (sim_graph2[path[k]].parents.size() == node.children.size()) {
                                for (int j = 0; j < node.children.size(); ++j) {
                                    if (abs(sim_graph2[path[k-1]].coverage - sim_graph2[node.children[j]].coverage) <
                                        child_min_cov * 0.02 && sim_graph2[node.children[j]].coverage!=0) {
                                        find_flag = true;
                                        cov_num.push_back(node.children[j]);
                                    }
                                }
                                if (cov_num.size() > 0) {
                                    break;
                                }
                            }
                        }

                        if (find_flag){
                            if (cov_num.size() == 1){
                                node_r[node.id] = 2;

                                node = sim_graph2[cov_num[0]];
                                path.push_back(node.id);

                                use_node[node.id]++;
                                if (node.coverage > 0){
                                    path_node_coverage.push_back(node.coverage);
                                }
                            } else{
                                double min_cov = 1;
                                int min_n_i;
                                for (int pa = 0; pa < cov_num.size(); ++pa) {
                                    if (sim_graph2[cov_num[pa]].coverage < min_cov){
                                        min_cov = sim_graph2[cov_num[pa]].coverage;
                                        min_n_i = cov_num[pa];
                                    }
                                }

                                node_r[node.id] = 3;
                                node = sim_graph2[min_n_i];
                                path.push_back(node.id);



                                use_node[node.id]++;
                                if (node.coverage > 0){
                                    path_node_coverage.push_back(node.coverage);
                                }
                            }
                        } else{
                                    double s_sim_cov = 0;
                                    for (int j = 0; j < node.children.size(); ++j) {
                                        s_sim_cov = s_sim_cov + sim_graph2[node.children[j]].coverage;
                                    }
                                    if (abs(s_sim_cov - node.coverage) < 0.03) {//此处可能是个分支
                                        double min_c_c = 1;
                                        int min_c_c_i;
                                        bool min_c_c_flag = false;
                                        for (int j = 0; j < node.children.size(); ++j) {
                                            if (sim_graph2[node.children[j]].coverage != 0 &&
                                                sim_graph2[node.children[j]].coverage < min_c_c) {//先选取两个可能的分支中丰度小的那一条
                                                min_c_c = sim_graph2[node.children[j]].coverage;
                                                min_c_c_i = j;
                                                min_c_c_flag = true;
                                            }
                                        }
                                        if (min_c_c_flag) {//分支中还有非0的分支
                                            node_r[node.id] = 6;

                                            node = sim_graph2[node.children[min_c_c_i]];
                                            path.push_back(node.id);


                                            use_node[node.id]++;
                                            if (node.coverage > 0) {
                                                path_node_coverage.push_back(node.coverage);
                                                path_cov = node.coverage;
                                            }

                                        } else {//可能不能这样
                                            vector<int> paired_child;
                                            paired_child.resize(node.children.size());
                                            bool p_flag = false;
                                            for (int c = 0; c < node.children.size(); ++c) {
                                                int pair_num = 0;
                                                for (int s = path.size()-1; s >= 1; --s) {
                                                    if (same_position_sim[path[s]].find(node.children[c]) !=
                                                        same_position_sim[path[s]].end()) {
                                                        p_flag = true;
                                                        pair_num++;
                                                    }
                                                }
                                                paired_child[c] = pair_num;
                                            }
                                            int max_p_c=0;
                                            size_t max_p;
                                            for (int p = 0; p < paired_child.size(); ++p) {
                                                if (paired_child[p] > max_p_c){
                                                    max_p_c = paired_child[p];
                                                    max_p = p;
                                                }
                                            }
                                            if (p_flag) {
                                                node_r[node.id] = 11;

                                                node = sim_graph2[node.children[max_p]];
                                                path.push_back(node.id);

                                                use_node[node.id]++;
                                                break;
                                            } else{
                                                node_idx_t c = node.children[0];
                                                double cha = 1;

                                                for (int k = path.size() - 1; k >= 1; --k) {

                                                    for (int j = 0; j < node.children.size(); ++j) {
                                                        if (abs(sim_graph[path[k]].coverage - sim_graph[node.children[j]].coverage) < cha) {//当前点是否与路径中的某个点丰度相似
                                                            cha = abs(sim_graph[path[k]].coverage - sim_graph[node.children[j]].coverage);
                                                            c = node.children[j];
                                                        }
                                                    }
                                                }
                                                    node_r[node.id] = 5;

                                                    node = sim_graph[c];
                                                    path.push_back(node.id);

                                                    use_node[node.id]++;
                                                    if (node.coverage > 0) {
                                                        path_node_coverage.push_back(node.coverage);
                                                    }

                                            }

                                        }
                                    } else {
                                        node_idx_t c = node.children[0];
                                        double cha = 1;

                                        for (int k = path.size() - 1; k >= 1; --k) {

                                            for (int j = 0; j < node.children.size(); ++j) {
                                                if (abs(sim_graph[path[k]].coverage - sim_graph[node.children[j]].coverage) < cha) {//当前点是否与路径中的某个点丰度相似
                                                    cha = abs(sim_graph[path[k]].coverage - sim_graph[node.children[j]].coverage);
                                                    c = node.children[j];
                                                }
                                            }
                                        }
                                        node_r[node.id] = 5;

                                        node = sim_graph[c];
                                        path.push_back(node.id);

                                        use_node[node.id]++;
                                        if (node.coverage > 0) {
                                            path_node_coverage.push_back(node.coverage);
                                        }

                                    }

                        }
                    }
                }
                std::sort(path_node_coverage.begin(), path_node_coverage.end());


                if (path_cov == 10){

                    path_cov = 0;
                    double num_n = 0;

                    for (int j = 0; j < path.size(); ++j) {
                        if(sim_graph[path[j]].children.size() == 1 && sim_graph[path[j]].parents.size() == 1){
                            path_cov = path_cov + sim_graph[path[j]].coverage;
                            num_n++;
                        }
                    }

                    path_cov = path_cov/num_n;
                }

                pair<vector<size_t>,double> p;
                p.first = path;
                p.second = path_cov;
                paths.push_back(p);

                node_reason.push_back(node_r);
                node_r.clear();
                node_r.resize(3000);

                sum_cov = sum_cov + path_cov;
                for (int j = 0; j < path.size(); ++j) {
                    if (abs(sim_graph2[path[j]].coverage - path_cov) < 0.02 || sim_graph2[path[j]].coverage < path_cov){
                        sim_graph2[path[j]].coverage = 0;
                    } else{
                        sim_graph2[path[j]].coverage = sim_graph2[path[j]].coverage - path_cov;
                    }
                }

            }
        }


        for (int i = 0; i < sim_graph2.size(); ++i){//根据剩余没用过的节点再找路

            vector<size_t> path2;
            if (use_node.find(i)==use_node.end() && (sim_graph2[i].children.size() > 0 || sim_graph2[i].parents.size() > 0) &&
                sim_graph2[i].single_node_id.size() > 1){

                Node node = sim_graph2[i];
                path2.push_back(node.id);
                use_node[node.id]++;
                double path_cov2 = sim_graph2[i].coverage;
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

                    if (node.children.size() == 1){
                        node_r[node.id] = 9;

                        node = sim_graph2[node.children[0]];
                        path2.push_back(node.id);

//                            cout << "22  " << node.id <<endl;
                        use_node[node.id]++;
                    } else{
                        bool is_used_flag = false;
                        double cha = 1;
                        node_idx_t similar_c = node.children[0];
                        for (int c = 0; c < node.children.size(); ++c) {
                            if (use_node.find(node.children[c]) == use_node.end()
                            && abs(sim_graph[node.children[c]].coverage - path_cov2) < cha){
                                cha = abs(sim_graph[node.children[c]].coverage - path_cov2);
                                similar_c = node.children[c];
                                is_used_flag = true;
                            }
                        }

                        if (is_used_flag){
                            node_r[node.id] = 10;

                            node = sim_graph2[similar_c];
                            path2.push_back(node.id);
                            use_node[node.id]++;

//                            cout << "55  " << node.id <<endl;
                        } else{

                            vector<int> paired_child;
                            paired_child.resize(node.children.size());
                            bool p_flag = false;
                            for (int c = 0; c < node.children.size(); ++c) {
                                int pair_num = 0;
                                for (int s = path2.size()-1; s >= 1; --s) {
                                    if (same_position_sim[path2[s]].find(node.children[c]) !=
                                        same_position_sim[path2[s]].end()) {
                                        p_flag = true;
                                        pair_num++;
                                    }
                                }
                                paired_child[c] = pair_num;
                            }
                            int max_p_c=0;
                            size_t max_p;
                            for (int p = 0; p < paired_child.size(); ++p) {
                                if (paired_child[p] > max_p_c){
                                    max_p_c = paired_child[p];
                                    max_p = p;
                                }
                            }
                            if (p_flag) {
                                node_r[node.id] = 11;

                                node = sim_graph2[node.children[max_p]];
                                path2.push_back(node.id);

                                use_node[node.id]++;
//                                break;
                            }else{

                                node_idx_t c = node.children[0];
                                double cha = 1;

                                for (int k = path2.size() - 1; k >= 1; --k) {
                                    for (int j = 0; j < node.children.size(); ++j) {
                                        if (sim_graph[path2[k]].parents.size() == 1
                                            && sim_graph[path2[k]].children.size() == 1
                                            && abs(sim_graph[path2[k]].coverage - sim_graph[node.children[j]].coverage) < cha) {//当前点是否与路径中的某个点丰度相似
                                            cha = abs(sim_graph[path2[k]].coverage - sim_graph[node.children[j]].coverage);
                                            c = node.children[j];
                                        }
                                    }
                                }
                                if (cha != 1){
                                    node_r[node.id] = 12;

                                    node = sim_graph[c];
                                    path2.push_back(node.id);

                                    use_node[node.id]++;
                                } else{
                                    double cha = 1;
                                    node_idx_t similar_c;
                                    for (int c = 0; c < node.children.size(); ++c) {
                                        if (abs(sim_graph[node.children[c]].coverage - path_cov2) < cha){
                                            cha = abs(sim_graph[node.children[c]].coverage - path_cov2);
                                            similar_c = node.children[c];
                                        }
                                    }


                                    node_r[node.id] = 12;

                                    node = sim_graph2[similar_c];
                                    path2.push_back(node.id);

                                    use_node[node.id]++;
                                }


                            }
                        }
                    }

                }
                node = sim_graph2[i];
                used_count = use_node[node.id];
                used_num = 0;
                while (!node.parents.empty()){
                    if (use_node[node.id] == used_count){
                        used_num++;
                    } else{
                        used_count = use_node[node.id];
                        used_num = 0;
                    }

                    if (node.parents.size() == 1){
                        node_r[node.id] = 13;
                        node = sim_graph2[node.parents[0]];
                        path2.insert(path2.begin(),node.id);
                        use_node[node.id]++;

                    } else{
                        int not_use=0;
                        int not_p;
                        double cha = 1;

                        for (int p = 0; p < node.parents.size(); ++p) {
                            if (use_node.find(node.parents[p])== use_node.end()
                            && abs(sim_graph[node.parents[p]].coverage - path_cov2) < cha){
                                cha = abs(sim_graph[node.parents[p]].coverage - path_cov2);
                                not_use++;
                                not_p = node.parents[p];
                            }
                        }
                        if(not_use > 0){
                            node_r[node.id] = 14;
                            node = sim_graph2[not_p];
                            path2.insert(path2.begin(),node.id);

                            use_node[node.id]++;
                        } else{

                            node_idx_t c = node.parents[0];
                            double cha = 1;

                            for (int k = 1; k > path2.size(); ++k) {
                                for (int p = 0; p < node.parents.size(); ++p) {
                                    if (sim_graph[path2[k]].parents.size() == 1
                                    && sim_graph[path2[k]].children.size() == 1
                                    && abs(sim_graph[path2[k]].coverage - sim_graph[node.parents[p]].coverage) < cha) {//当前点是否与路径中的某个点丰度相似
                                        cha = abs(sim_graph[path2[k]].coverage - sim_graph[node.parents[p]].coverage);
                                        c = node.parents[p];
                                    }
                                }
                            }
                            if (cha != 1){
                                node_r[node.id] = 12;

                                node = sim_graph[c];
                                path2.insert(path2.begin(),node.id);

                                use_node[node.id]++;
                            } else{
                                double cha = 1;
                                node_idx_t similar_p;
                                for (int p = 0; p < node.parents.size(); ++p) {
                                    if (abs(sim_graph[node.parents[p]].coverage - path_cov2) < cha){
                                        cha = abs(sim_graph[node.parents[p]].coverage - path_cov2);
                                        similar_p = node.parents[p];
                                    }
                                }
                                node_r[node.id] = 15;

                                node = sim_graph2[similar_p];
                                path2.insert(path2.begin(),node.id);
                                use_node[node.id]++;
                            }


                        }
                    }
                }

                pair<vector<size_t>,double> p;
                p.first = path2;
                p.second = path_cov2;
                paths.push_back(p);

                node_reason.push_back(node_r);

                for (int j = 0; j < path2.size(); ++j) {
                    if (abs(sim_graph2[path2[j]].coverage - path_cov2) < 0.01 || sim_graph2[path2[j]].coverage < path_cov2){
                        sim_graph2[path2[j]].coverage = 0;
                    } else{
                        sim_graph2[path2[j]].coverage = sim_graph2[path2[j]].coverage - path_cov2;
                    }
                }

            }
        }

        return paths;
    }


    vector<pair<vector<size_t>,double>> get_contigs2(vector<Node>& sim_graph2,map<size_t,size_t>& use_node){
        cout << "Begin find paths ... " << endl;
        vector<pair<vector<size_t>,double>> paths;
        int m = 0;
        double sum_cov = 0;
        vector<Node> sim_graph = sim_graph2;

        int start_num = 0;
        double start_aver = 1;
        for (int i = 0; i < sim_graph2.size(); ++i) {
            if (sim_graph2[i].parents.size() == 0 && sim_graph2[i].children.size()!= 0){
                start_num++;

            }
            if (sim_graph2[i].parents.size() == 0 && sim_graph2[i].children.size()!= 0 && sim_graph2[i].coverage > 0 && start_aver > sim_graph2[i].coverage){
                start_aver = sim_graph2[i].coverage;

            }
        }
        double limit_cov;
        if (start_aver > 0.05){
            limit_cov = 0.05;
        } else{
            limit_cov = start_aver;
        }
//        cout << "limit cov : " << limit_cov <<endl;
        for (int i = 0; i < sim_graph2.size(); ++i) {
            double min_start = 1;
            int min_start_i;
            bool start_flag = false;

            for (int j = 0; j < sim_graph2.size(); ++j) {
                if (sim_graph2[j].parents.size() == 0 && sim_graph2[j].children.size()!= 0 && sim_graph2[j].coverage > limit_cov && min_start > sim_graph2[j].coverage ){
                    min_start = sim_graph2[j].coverage;
                    min_start_i = j;
                    start_flag = true;
                }
            }

            if (!start_flag){
                continue;
            }


            while(sim_graph2[min_start_i].coverage > limit_cov){
                vector<double> path_node_coverage;

                if (sim_graph2[min_start_i].coverage < 0.01 || use_node[min_start_i] >= 2){
                    break;
                }

                if (abs(sum_cov - 1) < 0.02 || sum_cov > 1.02){
                    break;
                }

                m++;
                Node node = sim_graph2[min_start_i];
                vector<size_t> path;
                path.push_back(node.id);
                use_node[node.id]++;
                double path_cov = 10;
                while (node.children.size()>0){
                    if (node.children.size() == 1){
                        node = sim_graph2[node.children[0]];
                        path.push_back(node.id);
                        use_node[node.id]++;
                        if (node.coverage > 0){
                            path_node_coverage.push_back(node.coverage);
                        }
                    } else{
                        bool p_flag = false;
                        node_idx_t pair_i;
                        int p_i;
                        bool find_flag = false;
                        vector<int> cov_num;
                        double child_min_cov = 1;
                        for (int cmc = 0; cmc < node.children.size(); ++cmc) {
                            if (sim_graph2[node.children[cmc]].coverage < child_min_cov){
                                child_min_cov = sim_graph2[node.children[cmc]].coverage;
                            }
                        }
                        for (int k = path.size() - 1; k >= 1; --k) {//根据点的丰度判断
                            if (sim_graph2[path[k]].parents.size() == node.children.size()) {
                                for (int j = 0; j < node.children.size(); ++j) {
                                    if (abs(sim_graph2[path[k-1]].coverage - sim_graph2[node.children[j]].coverage) <
                                        child_min_cov * 0.2 && sim_graph2[node.children[j]].coverage!=0) {
                                        find_flag = true;
                                        cov_num.push_back(node.children[j]);
                                    }
                                }
                                if (cov_num.size() > 0) {
                                    break;
                                }
                            }
                        }

                        if (find_flag){
                            if (cov_num.size() == 1){
                                node = sim_graph2[cov_num[0]];
                                path.push_back(node.id);
                                use_node[node.id]++;
                                if (node.coverage > 0){
                                    path_node_coverage.push_back(node.coverage);
                                }
                            } else{
                                double min_cov = 1;
                                int min_n_i;
                                for (int pa = 0; pa < cov_num.size(); ++pa) {
                                    if (sim_graph2[cov_num[pa]].coverage < min_cov){
                                        min_cov = sim_graph2[cov_num[pa]].coverage;
                                        min_n_i = cov_num[pa];
                                    }
                                }
                                node = sim_graph2[min_n_i];
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
                                    if (abs(path_cov - sim_graph2[node.children[j]].coverage) < 0.02 && sim_graph2[node.children[j]].coverage!=0) {

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
                                        if (abs(sim_graph2[path[k]].coverage - sim_graph2[node.children[j]].coverage) < 0.02 &&
                                            sim_graph2[node.children[j]].coverage != 0) {//当前点是否与路径中的某个点丰度相似
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
                                if (!final_flag) {//在当前路径中没有相似丰度的点，看丰度和
                                    double s_sim_cov = 0;
                                    for (int j = 0; j < node.children.size(); ++j) {
                                        s_sim_cov = s_sim_cov + sim_graph2[node.children[j]].coverage;
                                    }
                                    if (abs(s_sim_cov - node.coverage) < 0.03) {//此处可能是个分支
                                        double min_c_c = 1;
                                        int min_c_c_i;
                                        bool min_c_c_flag = false;
                                        for (int j = 0; j < node.children.size(); ++j) {
                                            if (sim_graph2[node.children[j]].coverage != 0 &&
                                                sim_graph2[node.children[j]].coverage < min_c_c) {//先选取两个可能的分支中丰度小的那一条
                                                min_c_c = sim_graph2[node.children[j]].coverage;
                                                min_c_c_i = j;
                                                min_c_c_flag = true;
                                            }
                                        }
                                        if (min_c_c_flag) {//分支中还有非0的分支
                                            node = sim_graph2[node.children[min_c_c_i]];
                                            path.push_back(node.id);
                                            use_node[node.id]++;
                                            if (node.coverage > 0) {
                                                path_node_coverage.push_back(node.coverage);
                                            }
                                            path_cov = node.coverage;
                                        } else {//可能不能这样
                                            double max_f = 0;
                                            int max_f_i;
                                            for (int c = 0; c < node.children.size(); ++c) {
                                                if (sim_graph[node.children[c]].coverage > max_f) {
                                                    max_f = sim_graph[node.children[c]].coverage;
                                                    max_f_i = c;
                                                }
                                            }
                                            for (int e = 0; e < sim_graph2[node.id].edge_coverage.size(); ++e) {
                                                if (sim_graph2[node.id].edge_coverage[e].first == node.children[max_f_i]){
                                                    sim_graph2[node.id].edge_coverage[e].second = 0;
                                                }
                                            }
                                            node = sim_graph2[node.children[max_f_i]];
                                            path.push_back(node.id);
                                            use_node[node.id]++;
                                            if (node.coverage > 0) {
                                                path_node_coverage.push_back(node.coverage);
                                            }
                                        }
                                    } else {
                                        double max_f = 0;
                                        int max_f_i;
                                        for (int c = 0; c < node.children.size(); ++c) {
                                            if (sim_graph[node.children[c]].coverage > max_f) {
                                                max_f = sim_graph[node.children[c]].coverage;
                                                max_f_i = c;
                                            }
                                        }
                                        for (int e = 0; e < sim_graph2[node.id].edge_coverage.size(); ++e) {
                                            if (sim_graph2[node.id].edge_coverage[e].first == node.children[max_f_i]){
                                                sim_graph2[node.id].edge_coverage[e].second = 0;
                                            }
                                        }
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

                    if (path_node_coverage.size() <= 4){
                        for (int pa = 0; pa < path_node_coverage.size(); ++pa){
                            if (path_node_coverage[pa] < 1){
                                path_cov = path_cov + path_node_coverage[pa];
                                num_n++;
                            }
                        }

                    } else{
                        for (int pa = 3; pa < path_node_coverage.size()*0.5; ++pa) {
                            if (path_node_coverage[pa]>0){
                                path_cov = path_cov + path_node_coverage[pa];
                                num_n++;
                            }
                        }
                    }

                    path_cov = path_cov/num_n;
                }

                pair<vector<size_t>,double> p;
                p.first = path;
                p.second = path_cov;
                paths.push_back(p);

                sum_cov = sum_cov + path_cov;
                for (int j = 0; j < path.size(); ++j) {
                    if (abs(sim_graph2[path[j]].coverage - path_cov) < 0.03 || sim_graph2[path[j]].coverage < path_cov){
                        sim_graph2[path[j]].coverage = 0;
                    } else{
                        sim_graph2[path[j]].coverage = sim_graph2[path[j]].coverage - path_cov;
                    }
                }

            }
        }
        for (int i = 0; i < sim_graph2.size(); ++i) {//根据剩余没用过的节点再找路

            vector<size_t> path2;
            if (use_node.find(i)==use_node.end() && (sim_graph2[i].children.size() > 0 || sim_graph2[i].parents.size() > 0) &&
                sim_graph2[i].single_node_id.size() > 1){
//                cout << "---------------- first model ----------------  " << i <<endl;
                Node node = sim_graph2[i];
                path2.push_back(node.id);
                use_node[node.id]++;
                double path_cov2 = sim_graph2[i].coverage;
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
                    if(use_node[node.id] >= 10){
// || (used_count >= 4 && used_num >= 5)
                        break;
                    }
                    if (node.children.size() == 1){
                        node = sim_graph2[node.children[0]];
                        path2.push_back(node.id);
//                            cout << "22  " << node.id <<endl;
                        use_node[node.id]++;
                    } else{
                        bool is_used_flag = false;
                        for (int c = 0; c < node.children.size(); ++c) {
                            if (use_node.find(node.children[c])==use_node.end()){
                                is_used_flag = true;
                                node = sim_graph2[node.children[c]];
                                path2.push_back(node.id);
//                                    cout << "55  " << node.id <<endl;
                                use_node[node.id]++;
                            }
                        }
                        if (!is_used_flag){

                            double cha = 1;
                            node_idx_t similar_c;
                            for (int c = 0; c < node.children.size(); ++c) {
                                if (abs(sim_graph2[node.children[c]].coverage - path_cov2) < cha){
                                    cha = abs(sim_graph2[node.children[c]].coverage - path_cov2);
                                    similar_c = node.children[c];
                                }
                            }
                            node = sim_graph2[similar_c];
                            path2.push_back(node.id);
                            use_node[node.id]++;
                        }
                    }

                }
                node = sim_graph2[i];
                used_count = use_node[node.id];
                used_num = 0;
                while (!node.parents.empty()){
                    if (use_node[node.id] == used_count){
                        used_num++;
                    } else{
                        used_count = use_node[node.id];
                        used_num = 0;
                    }
                    if(use_node[node.id] >= 10){
// || (used_count >= 3 && used_num >= 4)
                        break;
                    }
                    if (node.parents.size() == 1){
                        node = sim_graph2[node.parents[0]];
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
                            node = sim_graph2[not_p];
                            path2.insert(path2.begin(),node.id);
//                                cout << "88  " << node.id <<endl;
                            use_node[node.id]++;
                        } else{
                            double cha = 1;
                            node_idx_t similar_p;
                            for (int p = 0; p < node.parents.size(); ++p) {
                                if (abs(sim_graph2[node.parents[p]].coverage - path_cov2) < cha){
                                    cha = abs(sim_graph2[node.parents[p]].coverage - path_cov2);
                                    similar_p = node.parents[p];
                                }
                            }
                            node = sim_graph2[similar_p];
                            path2.insert(path2.begin(),node.id);
                            use_node[node.id]++;
                        }
                    }
                }

                pair<vector<size_t>,double> p;
                p.first = path2;
                p.second = path_cov2;
                paths.push_back(p);

//                    cout << "min_cov : " << path_cov2 <<endl;
                for (int j = 0; j < path2.size(); ++j) {
                    if (abs(sim_graph2[path2[j]].coverage - path_cov2) < 0.03 || sim_graph2[path2[j]].coverage < path_cov2){
                        sim_graph2[path2[j]].coverage = 0;
                    } else{
                        sim_graph2[path2[j]].coverage = sim_graph2[path2[j]].coverage - path_cov2;
                    }
                }

            }
        }
        return paths;
    }


    //获取两个字符串中不同字符的个数
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

    void similarity_between_roads(vector<pair<string,double>>& paths_str,vector<pair<vector<size_t>,double>>& paths
                                  ,vector<Node>& sim_graph,vector<vector<size_t>>& sim_layer,set<int>& delete_p){
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

        if (vec.size() > 2){

            if (vec[0].second == vec[1].second){
                if (paths_str[vec[0].first].second > paths_str[vec[1].first].second) {
                    delete_p.insert(vec[1].first);
                    paths_str.erase(paths_str.begin()+vec[1].first);
                } else {
                    delete_p.insert(vec[0].first);
                    paths_str.erase(paths_str.begin()+vec[0].first);
                }
            } else{
                delete_p.insert(vec[0].first);
                delete_p.insert(vec[1].first);
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
                delete_p.insert(it.first);
                paths_str.erase(paths_str.begin()+it.first);
            }
        }
        vector<pair<vector<size_t>,double>> renew_paths2;
        for (int i = 0; i < paths.size(); ++i) {
            if (delete_p.find(i) == delete_p.end()){
                renew_paths2.push_back(paths[i]);
            }
        }
        paths.clear();
        paths = renew_paths2;
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

    bool user_num(string& trunk){
        map<string,int> u;
        for (int i = 0; i < trunk.length()-g_kmer_length+1; ++i) {
            string kmer = trunk.substr(i,g_kmer_length);
            u[kmer]++;
        }
        bool big_one = true;
        for (auto ur : u){
            if (ur.second > 1){
                big_one = false;
                break;
            }
        }
        return big_one;
    }

    set<size_t> get_gap_node(std::vector<Node>& after_gra){
        set<size_t> gap_node;
        for (int i = 0; i < after_gra.size(); ++i) {
            for (int j = 0; j < after_gra[i].children.size(); ++j) {
                if (after_gra[i].node_layer+1 != after_gra[after_gra[i].children[j]].node_layer){
                    gap_node.insert(i);
                }
            }
        }
        return gap_node;
    }

    void path_coverage2(KmerMap& kmerMap,vector<pair<vector<size_t>,double>>& paths,std::vector<std::string>& data,set<int>& delete_p,
                       map<size_t,size_t>& used_p,vector<Node>& sim_graph1,vector<pair<string,double>>& paths_str
                       ,vector<vector<size_t>>& sim_layer,int& max_layer,map<size_t,vector<size_t>>& lay_nodes){
        map<int,double> decide_nodes;
        vector<set<int>> paths_nodes;
        for (int i = 0; i < lay_nodes[max_layer].size(); ++i) {
            for (int j = 0; j < sim_layer[lay_nodes[max_layer][i]].size(); ++j) {
                if (sim_graph1[sim_layer[lay_nodes[max_layer][i]][j]].children.size() == 1
                && sim_graph1[sim_layer[lay_nodes[max_layer][i]][j]].parents.size() == 1 ){
                   decide_nodes[sim_layer[lay_nodes[max_layer][i]][j]] = sim_graph1[sim_layer[lay_nodes[max_layer][i]][j]].coverage;
                }
            }
        }

        vector<set<double>> group_nodes;
        for (int i = 0; i < paths.size(); ++i) {
//            cout << ">" << i <<endl;
            std::vector<double> data;
            double local_cov = 0;
            for (int j = 0; j < paths[i].first.size(); ++j) {
                if (decide_nodes.find(paths[i].first[j]) != decide_nodes.end()){
                    data.push_back(sim_graph1[paths[i].first[j]].coverage);
//                    cout << paths[i].first[j] << "(" << sim_graph1[paths[i].first[j]].coverage << ") ";
                }
            }
//            cout << "" <<endl;
            set<double> group_node;
            if (data.size() > 0){
                groupNodeDifference(data,local_cov,group_node);
                group_nodes.push_back(group_node);
//            cout << local_cov << endl;
                paths_str[i].second = local_cov;
            }

        }
        vector<double> cov_p;
        double  sum_cov = 0;
        for (int i = 0; i < paths_str.size(); ++i) {
            cov_p.push_back(paths_str[i].second);
            sum_cov = sum_cov + paths_str[i].second;
        }

        if(abs(sum_cov - 1) > 0.2){
            for (int i = 0; i < group_nodes.size(); ++i) {
                double c_cov = 0;
                if (group_nodes[i].size() > 1){
                    for (auto g:group_nodes[i]){
                        if (g != paths_str[i].second){
                            c_cov = g;
//                            cout << c_cov << endl;
                            break;
                        }

                    }
                }
                if (c_cov != 0){
                    sum_cov = sum_cov - paths_str[i].second + c_cov;
                    paths_str[i].second = c_cov;
                }
                if (abs(sum_cov - 1) < 0.2){
                    break;
                }
            }
        }
    }

    void groupNodeDifference(std::vector<double>& data,double& local_cov,set<double>& group_node) {

        std::sort(data.begin(), data.end());
        double threshold = 0.02;
        std::vector<std::vector<double>> groupedData;

        for (double num : data) {
            // 如果没有分组，创建一个新的分组
            if (groupedData.empty()) {
                groupedData.push_back({num});
            } else {
                bool addedToGroup = false;
                // 检查当前数值是否可以添加到现有的分组
                for (auto& group : groupedData) {
                    if (std::abs(group.back() - num) <= threshold) {
                        group.push_back(num);  // 将当前数值加入到该分组
                        addedToGroup = true;
                        break;
                    }
                }

                // 如果没有找到合适的分组，创建一个新的分组
                if (!addedToGroup) {
                    groupedData.push_back({num});
                }
            }
        }

        int mix_i = 0;
        int mix_size = 1;
        for (int i = 0; i < groupedData.size(); ++i) {
            for (int j = 0; j < groupedData[i].size(); ++j) {
                if (groupedData[i][j] < mix_size){
                    mix_size = groupedData[i][j];
                    mix_i = i;
                }
            }
            group_node.insert(groupedData[i][groupedData[i].size()/2]);

        }

        local_cov = groupedData[mix_i][groupedData[mix_i].size()/2];


    }

    // 归一化处理
    vector<double> normalize(const vector<double>& weightedMean) {
        vector<double> normalizedValues;
        double sum = 0;
        for (size_t i = 0; i < weightedMean.size(); ++i) {
            sum = sum + weightedMean[i];
        }
        // 按照线性归一化公式进行归一化
        for (size_t i = 0; i < weightedMean.size(); ++i) {
            double normalized = weightedMean[i] / sum;
            normalizedValues.push_back(normalized);
        }
        return normalizedValues;
    }


    vector<pair<vector<size_t>,double>> adjust_paths(vector<vector<size_t>>& paths,std::vector<Node>& after_gra2, vector<Node>& sim_graph1
            ,vector<vector<size_t>>& sim_layer,vector<double>& paths_cov,double& min_c,vector<double>& initial_cov
            ,vector<vector<size_t>>& new_paths,std::map<size_t,std::set<size_t>>& layer_nodes){
        cout << "Adjust ..." << endl;

        vector<pair<vector<size_t>,double>> renew_paths;
        map<int,int> node_num;

        vector<map<int,int>> paths_nodes;

        for (int i = 0; i < paths.size(); ++i) {

            vector<int> current;
            vector<int> no_current;
//            tf_file << ">" << i << endl;
            map<int,int> path_node;
            for (int j = 0; j < paths[i].size(); ++j) {
                for (int k = 0; k < sim_graph1[paths[i][j]].c_node_id.size(); ++k) {
                    node_num[sim_graph1[paths[i][j]].c_node_id[k]]++;
                }
                path_node[paths[i][j]]++;
            }

            paths_nodes.push_back(path_node);


        }

        for (int i = 0; i < paths.size(); ++i) {//获取路径丰度

            int max = 0;
            double max_i = 0;
            map<double,int> cover_num;
            for (int j = 0; j < paths[i].size(); ++j) {
                for (int k = 0; k < sim_graph1[paths[i][j]].c_node_id.size(); ++k) {

                    if(after_gra2[sim_graph1[paths[i][j]].c_node_id[k]].children.size() == 1
                    && after_gra2[sim_graph1[paths[i][j]].c_node_id[k]].parents.size() == 1){
                        cover_num[after_gra2[sim_graph1[paths[i][j]].c_node_id[k]].coverage]++;
                    }

                }
            }
            int stop_size = cover_num.size()*0.5;
            int count = 0;

            for (auto item:cover_num) {

                if (count >= stop_size)
                    break;
                count++;
                if (item.second > max){
                    max = item.second;
                    max_i = item.first;
                }
            }
            paths_cov.push_back(max_i);

        }

        vector<map<int,int>> change_nodes;
        change_nodes.resize(paths.size());
        //找出错误的点

        for (int i = 0; i < paths.size(); ++i) {

            map<int,int> change_node;
            for (int p = 0; p < paths.size(); ++p) {
                if (i != p){

                    for (int k = 0; k < paths[i].size(); ++k) {
                        bool flag = false;
                        for (int l = 0; l < sim_graph1[paths[i][k]].c_node_id.size(); ++l) {
//
                            if (node_num[sim_graph1[paths[i][k]].c_node_id[l]] > 1
                               && abs(after_gra2[sim_graph1[paths[i][k]].c_node_id[l]].coverage - paths_cov[p]) < min_c
                                && abs(after_gra2[sim_graph1[paths[i][k]].c_node_id[l]].coverage - paths_cov[i]) > min_c
                                && abs(after_gra2[sim_graph1[paths[i][k]].c_node_id[l]].coverage - paths_cov[i]) > abs(after_gra2[sim_graph1[paths[i][k]].c_node_id[l]].coverage - paths_cov[p])){
//                                cout << sim_graph1[paths[i][k]].c_node_id[l] << " ";
                                flag = true;
                            }
                        }
                        if (flag){
                            if ((sim_graph1[paths[i][k]].children.size() == 1 || sim_graph1[paths[i][k]].children.size() == 0)
                            && (sim_graph1[paths[i][k]].parents.size() == 1 || sim_graph1[paths[i][k]].parents.size() == 0)){
//                                file << paths[i][k] << " ";
                                change_node[paths[i][k]]=k;
                            }
                        }
                    }

                }
            }
            change_nodes[i] = change_node;

        }


        vector<vector<vector<size_t>>> change_node_groups;
        change_node_groups.resize(change_nodes.size());
        //按照点的连续关系给错误点分组，找出连续错误的区域
        for (int i = 0; i < change_nodes.size(); ++i) {
            if (change_nodes[i].empty()){
                continue;
            }
            vector<vector<size_t>> change_node_group;
            vector<size_t> group;

            int begin_layer = sim_graph1[change_nodes[i].begin()->first].node_layer;
            group.push_back(change_nodes[i].begin()->first);

            auto change = change_nodes[i].begin();
            ++change;

            for (;change != change_nodes[i].end(); ++change) {
                begin_layer++;
                if (sim_graph1[change->first].node_layer == begin_layer){
                    group.push_back(change->first);
                } else{
                    change_node_group.push_back(group);
                    group.clear();
                    group.push_back(change->first);
                    begin_layer = sim_graph1[change->first].node_layer;
                }

            }
            change_node_group.push_back(group);

            change_node_groups[i] = change_node_group;

        }
//        cout << " 1 =================== " <<endl;
        correction_paths(paths,sim_graph1,sim_layer,min_c, paths_cov,change_nodes,change_node_groups);

//        cout << " 2 =================== " <<endl;
        change_nodes.clear();
        change_nodes.resize(paths.size());
        //找出错误的点
        for (int i = 0; i < paths.size(); ++i) {

            map<int,int> change_node;
            for (int p = 0; p < paths.size(); ++p) {
                if (i != p){

                    for (int k = 0; k < paths[i].size(); ++k) {
                        bool flag = false;
                        for (int l = 0; l < sim_graph1[paths[i][k]].c_node_id.size(); ++l) {
//
                            if (node_num[sim_graph1[paths[i][k]].c_node_id[l]] > 1
                                && abs(after_gra2[sim_graph1[paths[i][k]].c_node_id[l]].coverage - paths_cov[p]) < min_c
                                && abs(after_gra2[sim_graph1[paths[i][k]].c_node_id[l]].coverage - paths_cov[i]) > min_c
                                && abs(after_gra2[sim_graph1[paths[i][k]].c_node_id[l]].coverage - paths_cov[i]) > abs(after_gra2[sim_graph1[paths[i][k]].c_node_id[l]].coverage - paths_cov[p])){

                                flag = true;
                            }
                        }
                        if (flag){
                            if ((sim_graph1[paths[i][k]].children.size() == 1 || sim_graph1[paths[i][k]].children.size() == 0)
                                && (sim_graph1[paths[i][k]].parents.size() == 1 || sim_graph1[paths[i][k]].parents.size() == 0)){

                                change_node[paths[i][k]]=k;
                            }
                        }
                    }

                }
            }
            change_nodes[i] = change_node;

        }

        //相似的路径，记录相似节点个数

        vector<map<int,int>> same_nodes;
        for (int i = 0; i < paths.size(); ++i) {

            map<int,int> same;
            for (int j = 0; j < paths.size(); ++j) {
                vector<size_t> same_node;

                if (i != j){
                    same_node = find_common_nodes(paths[i],paths[j]);
                }
                if (same_node.size() > 0){

                    if (same_node.size() > paths[i].size()*0.5){
                        same[j] = same_node.size();
                    }

                }
            }
            same_nodes.push_back(same);

        }
//        same_file.close();

        set<int> no_save_paths;
        set<int> no_save_paths2;
        for (int i = 0; i < same_nodes.size(); ++i) {

            if (no_save_paths.find(i) == no_save_paths.end() && no_save_paths2.find(i) == no_save_paths2.end()){
//                cout << "--------- " << i <<endl;
                if (same_nodes[i].size() > 1){
                    for (auto s:same_nodes[i]) {
                        if(change_nodes[s.first].size() > paths[i].size() * 0.2){
                            no_save_paths.insert(s.first);
                            no_save_paths2.insert(s.first);
//                            cout << "delete : " << s.first << endl;
                        }
                    }

                } else if (same_nodes[i].size() == 1){

                    auto s = same_nodes[i].begin();
                        if (change_nodes[s->first].size() > change_nodes[i].size() && paths[i].size() >= paths[s->first].size()){
                            no_save_paths.insert(s->first);
                            no_save_paths2.insert(i);
                        } else if (change_nodes[s->first].size() < change_nodes[i].size() && paths[i].size() <= paths[s->first].size()){
                            no_save_paths.insert(i);
                            no_save_paths2.insert(s->first);
                        } else{
                            if (paths[i].size() > paths[s->first].size()){
                                no_save_paths.insert(s->first);
                                no_save_paths2.insert(i);
                            } else{
                                no_save_paths.insert(i);
                                no_save_paths2.insert(s->first);
                            }
                        }
                }

            }

        }


        change_node_groups.clear();
        change_node_groups.resize(change_nodes.size());
        //按照点的连续关系给错误点分组，找出连续错误的区域
        for (int i = 0; i < change_nodes.size(); ++i) {
            if (change_nodes[i].empty()){
                continue;
            }
            vector<vector<size_t>> change_node_group;
            vector<size_t> group;

            int begin_layer = sim_graph1[change_nodes[i].begin()->first].node_layer;
            group.push_back(change_nodes[i].begin()->first);

            auto change = change_nodes[i].begin();
            ++change;

            for (;change != change_nodes[i].end(); ++change) {
                begin_layer++;
                if (sim_graph1[change->first].node_layer == begin_layer){
                    group.push_back(change->first);
                } else{
                    change_node_group.push_back(group);
                    group.clear();
                    group.push_back(change->first);
                    begin_layer = sim_graph1[change->first].node_layer;
                }

            }
            change_node_group.push_back(group);

            change_node_groups[i] = change_node_group;

        }

        if (paths.size() > 10){
            for (int i = 0; i < paths.size(); ++i) {
                if (no_save_paths.find(i) == no_save_paths.end()){
                    new_paths.push_back(paths[i]);
                    pair<vector<size_t>,double> pair_one;
                    pair_one.first = paths[i];
                    pair_one.second = paths_cov[i];
                    renew_paths.push_back(pair_one);
                }
            }

            paths.clear();
            paths = new_paths;
        } else{
            for (int i = 0; i < paths.size(); ++i) {
//                if (no_save_paths.find(i) == no_save_paths.end()){
                    new_paths.push_back(paths[i]);
                    pair<vector<size_t>,double> pair_one;
                    pair_one.first = paths[i];
                    pair_one.second = paths_cov[i];
                    renew_paths.push_back(pair_one);
//                }
            }
        }
        return renew_paths;

    }

    void correction_paths(vector<vector<size_t>>& paths,vector<Node>& sim_graph1,vector<vector<size_t>>& sim_layer,double& min_c
                          ,vector<double>& paths_cov,vector<map<int,int>>& change_nodes,vector<vector<vector<size_t>>>& change_node_groups){
//        ofstream changeFile("15_tihuan.txt");

        for (int i = 0; i < change_node_groups.size(); ++i) {

            for (int j = 0; j < change_node_groups[i].size(); ++j) {
                bool flag = false;
                if (change_node_groups[i][j].size() > 2){

                    if (sim_graph1[change_node_groups[i][j][0]].parents.size() == 1
                        && sim_graph1[sim_graph1[change_node_groups[i][j][0]].parents[0]].children.size() > 1){//将路径中的这部分换成其他孩子的
                        double cha = 1;
                        int cha_index = 0;
                        for (int k = 0; k < sim_graph1[sim_graph1[change_node_groups[i][j][0]].parents[0]].children.size(); ++k) {
                            int parent = sim_graph1[change_node_groups[i][j][0]].parents[0];
                            if (sim_graph1[parent].children[k] != change_node_groups[i][j][0]){
                                if (cha > abs(sim_graph1[sim_graph1[parent].children[k]].coverage - paths_cov[i])){
                                    cha = abs(sim_graph1[sim_graph1[parent].children[k]].coverage - paths_cov[i]);
                                    cha_index = sim_graph1[parent].children[k];
                                }

                            }
                        }
                        if (cha < min_c){
                            flag = true;
                            vector<size_t> time_change_path;
                            Node node = sim_graph1[cha_index];
                            time_change_path.push_back(node.id);

                            while (time_change_path.size() != change_node_groups[i][j].size()){
                                if (node.children.size() == 1){
                                    node = sim_graph1[node.children[0]];
                                    time_change_path.push_back(node.id);
                                } else{
                                    double children_cha = 1;
                                    int children_index;
                                    for (int k = 0; k < node.children.size(); ++k) {
                                        if (children_cha > abs(sim_graph1[node.children[k]].coverage-paths_cov[i])){
                                            children_cha = abs(sim_graph1[node.children[k]].coverage-paths_cov[i]);
                                            children_index = node.children[k];
                                        }
                                    }
                                    node = sim_graph1[children_index];
                                    time_change_path.push_back(node.id);
                                }

                            }

                            for (int k = 0; k < time_change_path.size(); ++k) {
                                paths[i][change_nodes[i][change_node_groups[i][j][k]]] = time_change_path[k];
                            }
                        }
                        else{

                            if (sim_graph1[change_node_groups[i][j][change_node_groups[i][j].size()-1]].children.size() == 1
                                && sim_graph1[sim_graph1[change_node_groups[i][j][change_node_groups[i][j].size()-1]].children[0]].parents.size() > 1){
                                cha = 1;
                                cha_index = 0;
                                for (int k = 0; k < sim_graph1[sim_graph1[change_node_groups[i][j][change_node_groups[i][j].size()-1]].children[0]].parents.size(); ++k) {
                                    int children = sim_graph1[change_node_groups[i][j][change_node_groups[i][j].size()-1]].children[0];
                                    if (sim_graph1[children].parents[k] != change_node_groups[i][j][change_node_groups[i][j].size()-1]){
                                        if (cha > abs(sim_graph1[sim_graph1[children].parents[k]].coverage - paths_cov[i])){
                                            cha = abs(sim_graph1[sim_graph1[children].parents[k]].coverage - paths_cov[i]);
                                            cha_index = sim_graph1[children].parents[k];
                                        }
                                    }
                                }
                                if (cha < min_c){
                                    flag = true;
                                    vector<size_t> time_change_path;
                                    time_change_path.resize(change_node_groups[i][j].size());
                                    int change_index = change_node_groups[i][j].size()-1;
                                    Node node = sim_graph1[cha_index];
                                    time_change_path[change_index] = node.id;
                                    while (1){
                                        change_index--;
                                        if (node.parents.size() == 1){
                                            node = sim_graph1[node.parents[0]];
                                            time_change_path[change_index] = node.id;
                                        } else if (node.parents.size() > 1){
                                            double parents_cha = 1;
                                            int parents_index;
                                            for (int k = 0; k < node.parents.size(); ++k) {
                                                if (parents_cha > abs(sim_graph1[node.parents[k]].coverage-paths_cov[i])){
                                                    parents_cha = abs(sim_graph1[node.parents[k]].coverage-paths_cov[i]);
                                                    parents_index = node.parents[k];
                                                }
                                            }
                                            node = sim_graph1[parents_index];
                                            time_change_path[change_index] = node.id;
                                        }
                                        if (change_index == 0)
                                            break;
                                    }

                                    for (int k = 0; k < time_change_path.size(); ++k) {
//                                        changeFile << "~~"  << time_change_path[k] << "(" << sim_graph1[time_change_path[k]].coverage << ") ";

                                        paths[i][change_nodes[i][change_node_groups[i][j][k]]] = time_change_path[k];
                                    }
//                                    changeFile << "" <<endl;
                                }


                            }

                        }
                    }
                }


            }
        }

    }


    vector<size_t> find_common_nodes(const vector<size_t>& path1, const vector<size_t>& path2) {
        set<size_t> set1(path1.begin(), path1.end()); // 将path1转为集合
        set<size_t> set2(path2.begin(), path2.end()); // 将path2转为集合

        vector<size_t> common_nodes;

        // 使用set_intersection求两个集合的交集
        set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), back_inserter(common_nodes));

        return common_nodes; // 返回交集，即相同的节点
    }


    void output_graph(KmerMap& kmerMap, kmer_int_type_t& seed_kmer,int& average
                      ,std::vector<std::string>& data,string& trunk,int& input_type){

        cout << "Begin extend sequence..." << endl;
        string trunk_str;
        string seed = intval_to_kmer(seed_kmer,g_kmer_length);
        if (trunk.find(seed) != std::string::npos){
            trunk_str = trunk;
        } else{
            trunk_str = revcomp(trunk);
        }
        bool big_one = user_num(trunk);

        if (!big_one){
            size_t index = trunk_str.find(seed);
            string left = trunk_str.substr(0,index);
            trunk_str = get_trunk(kmerMap,seed_kmer,average,left);
        }

        std::map<string,pair<size_t,size_t>> span_part;
        cout << "Begin get all nodes..." << endl;
        int in_del_num = 0;
        std::vector<std::set<string>> nodes = get_all_nodes(kmerMap,trunk_str,span_part,in_del_num);


        //获取所有kmer
        std::vector<kmer_int_type_t> kmer_all = kmerMap.get_all_kmer();
        //获取首部缺少的分支kmer
        not_used_candi(kmerMap,kmer_all,nodes);
        //删除错误kmer
        delete_error_kmer(kmerMap,nodes);

        //连接成图
        vector<std::map<string,node_idx_t>> node_id_position = get_bubble_graph(nodes,span_part);
        set_parents(node_set_);

        cout << "Bubble-structured de Bruijn graph nodes : " << node_set_.size() <<endl;
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

        std::vector<Node> after_gra2 = handle_nodes(after_gra,layer_nodes);

        if (after_gra2.size()!= after_gra.size()){

            re_id(after_gra2);
            set_parents(after_gra2);
        }
        cout << "Simplified Bubble-structured de Bruijn graph nodes : " <<  after_gra2.size() <<endl;

        adjust_layer(after_gra2,layer_nodes);

        int limit_length = 0;
        int min_l = data[0].length();
        for(int li = 0 ;li < data.size(); li++){
            if (min_l > data[li].length()){
                min_l = data[li].length();
            }
            limit_length = limit_length + data[li].length();
        }

        limit_length = limit_length/data.size();
        limit_length = (limit_length+min_l)/2;
        if(limit_length < g_kmer_length*4){
            limit_length = limit_length * 2;
        }


        time_t compression_begin = time(NULL);
        vector<vector<size_t>> sim_layer;
        vector<size_t> layer_sim_node;
        std::unordered_map<std::vector<size_t>,string,VectorHash> vec_string;


        double min_c = 1;
        for (int i = 0; i < after_gra2.size(); ++i) {
            if (after_gra2[i].coverage < min_c){
                min_c = after_gra2[i].coverage;
            }
        }

        vector<double> node_cov;

        vector<vector<size_t>> sim_gra = compression_nodes(kmerMap,after_gra2,layer_nodes,sim_layer
                                                           ,layer_sim_node,data,limit_length,vec_string
                                                           ,node_cov,min_c);

        map<size_t,vector<size_t>> lay_nodes;
        int max_layer = 0;
        for(int i = 0; i < sim_layer.size(); i++){
            lay_nodes[sim_layer[i].size()].push_back(i);
            if (sim_layer[i].size() > max_layer){
                max_layer = sim_layer[i].size();
            }
        }


        time_t compression_end = time(NULL);
        cout << "compression_time: " << compression_end - compression_begin << "s" << endl;

        vector<set<size_t>> sim_alone_node;
        vector<vector<pair<size_t,size_t>>> parent_different_node;
        vector<vector<pair<size_t,size_t>>> children_different_node;
        vector<Node> sim_graph1;
        vector<vector<int>> why_connect;
        create_sim_graph(kmerMap,after_gra2,sim_gra,sim_graph1,layer_nodes,sim_layer,layer_sim_node
                         ,parent_different_node,data,min_l,vec_string,node_cov,why_connect,min_c);


        vector<pair<vector<size_t>,double>> paths;
        vector<Node> sim_graph2;
        map<size_t,size_t> use_node;
        vector<Node> sim_graph3;

        vector<vector<size_t>> node_reason;

        if (input_type == 3){
            cout << "paired_reads" <<endl;
            std::map<size_t,size_t> read_first_node_sim;//reads对应的点
            std::map<size_t,vector<size_t>> node_pair_sim;//点对应的reads
            std::vector<std::map<size_t,size_t>> same_position_sim;
            same_position_sim.resize(sim_graph1.size());
            get_sim_paired_position(kmerMap,data,after_gra,sim_graph1
                                    ,read_first_node_sim,parent_different_node,node_pair_sim);
            get_pair_end(same_position_sim,data,sim_graph1,read_first_node_sim,node_pair_sim);

            sim_graph2 = sim_graph1;
            determine_nodes_contigs(sim_graph1,sim_graph2,sim_layer,why_connect);
            sim_graph3 = sim_graph2;

            paths = get_contigs(sim_graph2,same_position_sim,use_node,node_reason);
        } else if (input_type == 2){
            cout << "single_reads" <<endl;
            sim_graph2 = sim_graph1;
            determine_nodes_contigs(sim_graph1,sim_graph2,sim_layer,why_connect);
            sim_graph3 = sim_graph2;

            paths = get_contigs2(sim_graph2,use_node);
        }


        vector<double> initial_cov;
        vector<vector<size_t>> restore_paths;
        restore_paths.resize(paths.size());
        for (int i = 0; i < paths.size(); ++i) {
            string str = "";
            int max_after2 = 0;
            int max_af_i;
            for (int j = 0; j < paths[i].first.size(); ++j) {

                for (int k = 0; k < sim_graph2[paths[i].first[j]].single_node_id.size(); ++k) {
                    restore_paths[i].push_back(sim_graph2[paths[i].first[j]].single_node_id[k]);
                }
            }
            initial_cov.push_back(paths[i].second);

        }


        vector<vector<size_t>> new_paths;
        vector<double> paths_cov;
        vector<pair<vector<size_t>,double>> renew_paths = adjust_paths(restore_paths,after_gra2,sim_graph1,sim_layer,paths_cov,min_c,initial_cov,new_paths,layer_nodes);

        vector<pair<string,double>> paths_str;
        for (int i = 0; i < restore_paths.size(); ++i) {
            string str = "";

            for (int j = 0; j < restore_paths[i].size(); ++j) {

                if (j == 0){
                    str = str + sim_graph1[restore_paths[i][j]].sequence;
                } else{
                    str = str + sim_graph1[restore_paths[i][j]].sequence.substr(g_kmer_length-1);
                }
            }
            pair<string,double> ps;
            ps.first = str;
            ps.second = renew_paths[i].second;

            paths_str.push_back(ps);
        }

//        删除相似度高的contig
        set<int> delete_p;
        similarity_between_roads(paths_str,renew_paths,sim_graph1,sim_layer,delete_p);
        cout << "path_num : " << paths_str.size() <<endl;


        path_coverage2(kmerMap,renew_paths,data,delete_p,use_node,sim_graph1,paths_str,sim_layer,max_layer,lay_nodes);

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
    std::map<size_t,vector<size_t>> read_and_node;

};
#endif