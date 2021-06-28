// ******************************************************
// Author       : Haoqiang Guo
// Last modified: 2021-04-10 13:21
// Email        : ghaoqi1@lsu.edu
// Filename     : pec.h
// Description  : 
// ******************************************************
#ifndef _PEC_H_
#define _PEC_H_ 
#include<vector>
#include<cmath>
#include<iostream>
#include<glog/logging.h>
class Pec{
    public:
        Pec(const int,const int);
        void set_opt_counts(int);
        void set_opt_delay(int);
        void set_com_delay(int);
        void set_total_delay(double);
        void set_l1_read(int);
        void set_l1_write(int);
        int get_pe_counts();
        int get_l1();
        void print();
    //private:
        int pe_counts; // pes per pe cluster
        int l1; //L1 cache per pe cluste
    //protected:
    public:
        std::vector<int> opt_counts; // multiply counts
        std::vector<int> opt_delay; // computation delay
        std::vector<int> com_delay; // Intra-PEC communication delay
        std::vector<double> total_delay; 
        std::vector<int> l1_read; // pe read L1 
        std::vector<int> l1_write; // pe write L1

        std::vector<unsigned long> l2_read_2pe;
        std::vector<unsigned long> l2_write_from_pe; 
};
#endif /* _PEC_H_ */
