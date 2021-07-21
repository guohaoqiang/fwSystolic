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
#include<assert.h>
#include<glog/logging.h>
#include "acc.h"
#include "graphdata.h"
#include "mip.h"

#define TYPE_LENGTH long long
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
      
      void to_binary(TYPE_LENGTH n);
      void print_binary(std::shared_ptr<std::vector<TYPE_LENGTH>> &bitmap);
      void tab_gen(const std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> &intile, std::shared_ptr<std::vector<TYPE_LENGTH>> &d);

      unsigned int find_nz_pos(TYPE_LENGTH n);
      unsigned int count_nnz(TYPE_LENGTH n, unsigned int pos);
      void val();
      void permutate(std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> &intile,\
        std::shared_ptr<std::vector<int>> &vec);
      void print_tile(const std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> &intile);
      void      deconstruct_tab();
      void      deconstruct_self_loop();
      std::shared_ptr<std::vector<unsigned int>> self_loop= std::make_shared<std::vector<unsigned int>>();
      std::shared_ptr<std::vector<std::vector<int>>> tab = std::make_shared<std::vector<std::vector<int>>>();
      long  comp_delays = 0; 

      std::vector<unsigned int> shortcut_steps;
      long debubbling(const std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> &intile, std::shared_ptr<std::vector<TYPE_LENGTH>> &d);
      void val_naive1();// w/o MIP. Verifying the efficiency of MIP.
};
#endif /* _ANALYSIS_H_ */
