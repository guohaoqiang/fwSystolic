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
class Analysis{
    public:
      Analysis(std::shared_ptr<Acc>,std::shared_ptr<Graph>);  
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
      int height = 0;
      int len = 0;
      std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> tiles = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(); 
      std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> tiles_r = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(); 
      std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> tiles_c = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(); 
      double global_delay = 0;
      
      double global_edp = 0;
      double global_edp0 = 0;
      double global_edp1 = 0;
      double global_edp2 = 0;
      double global_edp3 = 0;
      double global_edp4 = 0;
      double global_edp5 = 0;
      double global_edp6 = 0;
      double global_edp7 = 0;
      
      double global_energy = 0;
      double global_energy0 = 0;
      double global_energy1 = 0;
      double global_energy2 = 0;
      double global_energy3 = 0;
      double global_energy4 = 0;
      double global_energy5 = 0;
      double global_energy6 = 0;
      double global_energy7 = 0;
      double global_energy8 = 0;
      
      double  l2_l1 = 0;
      double  l1_l2 = 0;
      double onchip_delay = 0;
      double onchip_delay0 = 0;
      double onchip_delay1 = 0;
      double onchip_delay2 = 0;

      double offchip_delay = 0;
      double offchip_delay1 = 0;
      double offchip_delay0 = 0;
      
      double pec_read_delay = 0;  // pec read from L2
      double pec_write_delay = 0; // pec write to L2


};

#endif /* _ANALYSIS_H_ */
