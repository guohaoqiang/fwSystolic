// ******************************************************
// Author       : Haoqiang Guo
// Last modified: 2021-04-12 18:49
// Email        : ghaoqi1@lsu.edu
// Filename     : analysis.h
// Description  : 
// ******************************************************
#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_ 
#include<memory>
#include<iostream>
#include<algorithm>
#include<cmath>
#include<set>
#include<float.h>
#include<glog/logging.h>
#include "acc.h"
#include "graphdata.h"
template<typename T>
class Analysis{
    public:
      Analysis(std::shared_ptr<Acc>,std::shared_ptr<Graph>,int win);  
      void run_w_order_wo_move(int tm,int tk,int tn, int tc);
      void run_baseline_awbgcn_gcnax(int tm,int tk,int tn,int tc);
      void run_wegnn();
      void run_baseline();
      void extract_timing();
      void print();
      void print_final(std::string);
      int row_divider = 1;
    //private:
      std::shared_ptr<Acc> hw;
      std::shared_ptr<Graph> dat;
      int window;
      int height = 0;
      int len = 0;
      std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> tiles \
          = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(); 
      std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> tiles_r \
          = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(); 
      std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> tiles_c \
          = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(); 
      double global_delay = 0;
      
      double global_edp = 0;
      
      double global_energy = 0;
      
      double  l2_l1 = 0;
      double  l1_l2 = 0;
      double onchip_delay = 0;

      double offchip_delay = 0;
      
      double pec_read_delay = 0;  // pec read from L2
      double pec_write_delay = 0; // pec write to L2

      //std::vector<std::vector<int>> tile(2,std::vector<int>()); //only record row and column indices
      
      void to_binary(int n);
      void print_binary(std::vector<int> &bitmap);
      void tab_gen(const std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> &intile, std::vector<int> &d);
      void val();
};

template<typename T>
Analysis<T>::Analysis(std::shared_ptr<Acc> accelerator,\
        std::shared_ptr<Graph> input_data, int win):\
        hw(accelerator),dat(input_data),window(win){
}  

template<typename T>
void Analysis<T>::to_binary(int n){
    int a = n % 2;
    n = n >> 1;
    if(n != 0){
        to_binary(n);
    }
    std::cout<<a;
}

template<typename T>
void Analysis<T>::print_binary(std::vector<int> &bitmap){
    for(auto item:bitmap){
        to_binary(item);
        std::cout<<std::endl;
    }
}

template<typename T>
void Analysis<T>::tab_gen(const std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> &intile, \
        std::vector<int> &d){
    int base = 0;
    for(size_t i=0; i<intile->size(); i++){ //#rows in a tile 
        VLOG(2)<<"#rows = "<<intile->size();
        for(size_t j=0; j<intile->at(i)->size(); j++){ //NNZ in the i-th row
            VLOG(2)<<"#NZ in the "<<j<<"-th row "<<"="<<intile->at(i)->size();
            unsigned int shifts = intile->at(i)->at(j).at(1)%hw->ar;
            //std::cout<<shifts<<std::endl;
            d.at(i) = d.at(i) | (1<<(hw->ar-shifts-1));
            std::cout<<d.at(i)<<std::endl;
        }
    }
    VLOG(2)<<"Print binary map:";
    print_binary(d);
}


template<typename T>
void Analysis<T>::print(){
    for(auto pec:hw->pec_arr){
       if(pec.opt_counts.size()!=pec.opt_delay.size() || \
               pec.opt_counts.size()!=pec.com_delay.size() ||\
               pec.opt_counts.size()!=pec.l1_read.size() ||\
               pec.opt_counts.size()!=pec.l1_write.size())
           LOG(FATAL)<<"Un..."<<std::endl;
       std::cout<<std::endl;
       std::cout<<"-----------------------"<<std::endl;
       //LOG(INFO)<<"total_delay: "<<pec.total_delay.size()<<std::endl;
       pec.print();
    }
}
template<typename T>
void Analysis<T>::val(){
    int base = 0;
    std::vector<int> mark(hw->ar,0);
    std::vector<int> count(hw->ar,0);
    //the inner most vector denotes an non-zero entry <row,col,val>

    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> tile = \
          std::make_shared<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>>(); 
    
    for(size_t i=base; i<dat->data.at(0).size(); i += hw->ar){ // iterate over col-tiles
        for(size_t j=i; j<std::min(i+hw->ar,dat->data.at(0).size()-1); j++){ // iterate over cols in a col-tile
            mark.at(j%hw->ar) = dat->data.at(0).at(j);
            //VLOG(2)<<"mark["<<j%hw->ar<<"]="<<dat->data.at(0).at(j);
        }
        count.assign(hw->ar,0);
        // in adj, #row == #col
        for(size_t k1=0; k1<dat->data.at(0).size()-1; k1 += window){ // iterate over row-tiles
        //for(size_t k1=0; k1<3; k1 += window){ // iterate over row-tiles
            for(size_t k2=k1; k2<std::min(k1+window,dat->data.at(0).size()-1); k2++){
                std::shared_ptr<std::vector<std::vector<int>>> tmp_tile_row = \
                                        std::make_shared<std::vector<std::vector<int>>>();
                // begin fetch one row of a tile
                for(size_t j=i; j<std::min(i+hw->ar,dat->data.at(0).size()-1); j++){ // iterate over cols in a col-tile
                   //std::cout<<"k2 = "<<k2<<" j = "<<j<<std::endl;
                   //VLOG(2)<<"#nz of "<<j<<":"<<dat->data.at(0).at(j+1) - dat->data.at(0).at(j);
                   if(dat->data.at(1).at(mark.at(j%hw->ar)) == k2){ 
                        std::vector<int> tmp_entry(3,0);
                        tmp_entry.at(0) = dat->data.at(1).at(mark.at(j%hw->ar));//row ID
                        tmp_entry.at(1) = j;//col ID
                        tmp_entry.at(2) = dat->data.at(2).at(mark.at(j%hw->ar));//vals
                        tmp_tile_row->push_back(tmp_entry);
                        count.at(j%hw->ar)++;
                        if(count.at(j%hw->ar)<(dat->data.at(0).at(j+1) - dat->data.at(0).at(j)))
                            mark.at(j%hw->ar)++;
                        VLOG(2)<<"k2 = "<<k2<<" j = "<<j<<" r="<<tmp_entry.at(0)<<" c="<<tmp_entry.at(1)<<" v="<<tmp_entry.at(2);
                   }
                } //end fetch one row of a tile
                std::cout<<std::endl;
                tile->push_back(tmp_tile_row);
            }// end fetch one tile 
            //std::cout<<"-----------"<<std::endl; 
            std::vector<int> bits(window,0); 
            //for(auto g:bits)
            //    std::cout<<g<<std::endl;
            tab_gen(tile,bits);
            return;
        }
    }
}
template<typename T>
void Analysis<T>::run_wegnn(){
    
    VLOG(0)<<"Accelerator: ";
    VLOG(0)<<"Total PEs: "<<hw->pes;
    VLOG(0)<<"Aspect Ratio: "<<hw->ar<<"*"<<hw->ac;
    VLOG(0)<<"#PEs in each PE cluster: "<<hw->pec_size;
    VLOG(0)<<"L1 in each PE cluster: "<<hw->l1<<" KB";
    VLOG(0)<<"L2 on chip: "<<hw->l2<<" MB";
    //VLOG(0)<<"Dataset: "<<dat->name<<"\-"<<dat->nodes<<","<<dat->max_num<<","<<dat->sorted_data.at(0).at(1)-dat->sorted_data.at(0).at(0);
    VLOG(0)<<"Dataset: "<<dat->name<<"-"<<dat->nodes<<","<<dat->feature_size<<","<<dat->hiden;

}
template<typename T>
void Analysis<T>::run_baseline(){
    
    VLOG(0)<<"Accelerator: ";
    VLOG(0)<<"Total PEs: "<<hw->pes;
    VLOG(0)<<"Aspect Ratio: "<<hw->ar<<"*"<<hw->ac;
    VLOG(0)<<"#PEs in each PE cluster: "<<hw->pec_size;
    VLOG(0)<<"L1 in each PE cluster: "<<hw->l1<<" KB";
    VLOG(0)<<"L2 on chip: "<<hw->l2<<" MB";
    //VLOG(0)<<"Dataset: "<<dat->name<<"\-"<<dat->nodes<<","<<dat->max_num;
    VLOG(0)<<"Dataset: "<<dat->name<<"-"<<dat->nodes<<","<<dat->feature_size<<","<<dat->hiden;

}
#endif /* _ANALYSIS_H_ */
