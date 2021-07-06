// ******************************************************
// Author       : Haoqiang Guo
// Last modified: 2021-04-08 22:42
// Email        : ghaoqi1@lsu.edu
// Filename     : graphdata.h
// Description  : 
// ******************************************************
#ifndef _GRAPHDATA_H_
#define _GRAPHDATA_H_ 
#include<glog/logging.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<map>
#include<string>
#include<vector>
#include<algorithm>

class Graph{
    public:
        Graph(const std::string&);
        std::string get_name();
        void node_sort(const std::vector<std::vector<int>>&);
        //void csr_insert_diag(const std::vector<std::vector<int>>&);
        
        int get_nodes();
        int get_feature_size();
        int hiden = 16;
        void print_data();
    public:
        std::vector<std::vector<int>> data;
        std::vector<std::vector<int>> sorted_data;
        std::vector<std::vector<int>> csr_diag;
        int feature_size;
        int max_num = 0;
        std::string name;
        int nodes;
    private:
        std::string data_path;
        std::map<std::string,int> feature_dim {{"demo1",1},{"cora",1433},{"cora_ml",2879},{"citeseer",3703},{"polblogs",1490},\
                                    {"pubmed",500},{"ppi",50},{"reddit",602},{"yelp",300},{"flickr",500},\
                                    {"amazon",200},{"ogbn",100},{"ppi_large",50}};
        std::map<std::string,int> output_dim {{"demo1",2},{"cora",7},{"cora_ml",7},{"citeseer",6},{"polblogs",2},\
                                    {"pubmed",3},{"ppi",121},{"reddit",41},{"yelp",100},{"flickr",7},\
                                    {"amazon",107},{"ogbn",47},{"ppi_large",121}};
        std::multimap<int,std::map<int,int>,std::greater<int>> pre_sort_data;
        void load_adj(const std::string&);
};

#endif /* _GRAPHDATA_H_ */
