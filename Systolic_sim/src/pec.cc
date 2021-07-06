#include "../include/pec.h"
Pec::Pec(const int pe_c,const int l1_):pe_counts(pe_c),l1(l1_){
}
void Pec::print(){
   LOG(INFO)<<"opt_counts: "<<opt_counts.size()<<std::endl;
   if(opt_counts.size()!=opt_delay.size() || \
           opt_counts.size()!=com_delay.size() ||\
           opt_counts.size()!=l1_read.size() ||\
           opt_counts.size()!=l1_write.size())
       LOG(FATAL)<<"Un..."<<std::endl;
   for(size_t i = 0; i<opt_counts.size(); ++i){
       std::cout<<"opt_counts: ";
       std::cout<<opt_counts.at(i)<<" ";
   }
   std::cout<<std::endl;
   for(size_t i = 0; i<opt_counts.size(); ++i){
       std::cout<<"opt_delay: ";
       std::cout<<opt_delay.at(i)<<" ";
   }
   std::cout<<std::endl;
   for(size_t i = 0; i<opt_counts.size(); ++i){
       std::cout<<"com_delay: ";
       std::cout<<com_delay.at(i)<<" ";
   }
   std::cout<<std::endl;
   for(size_t i = 0; i<opt_counts.size(); ++i){
       std::cout<<"l1_read: ";
       std::cout<<l1_read.at(i)<<" ";
   }
   std::cout<<std::endl;
   for(size_t i = 0; i<opt_counts.size(); ++i){
       std::cout<<"l1_write: ";
       std::cout<<l1_write.at(i)<<" ";
   }
   std::cout<<std::endl;
   LOG(INFO)<<"total delay: "<<std::endl;
   LOG(INFO)<<"total_delay: "<<total_delay.size()<<std::endl;
   for(size_t i = 0; i<opt_counts.size(); ++i){
       std::cout<<"total_delay: ";
       std::cout<<total_delay.at(i)<<" ";
   }
   LOG(INFO)<<"opt_counts: "<<opt_counts.size()<<std::endl;
   std::cout<<std::endl;
}
void Pec::set_opt_counts(int x){
    opt_counts.push_back(x);
}
void Pec::set_opt_delay(int x){
    opt_delay.push_back(x);
}
void Pec::set_com_delay(int x){
    com_delay.push_back(x);
}
void Pec::set_total_delay(double x){
    total_delay.push_back(x);
}
void Pec::set_l1_read(int x){
    l1_read.push_back(x);
}
void Pec::set_l1_write(int x){
    l1_write.push_back(x);
}
/*
std::vector<int> Pec::get_opt_counts(){
    return opt_counts;
}
void Pec::set_opt_delay(int& x){
    opt_delay.push_back(x);
}
void Pec::set_com_delay(int& x){
    com_delay.push_back(x);
}
void Pec::set_total_delay(int& x){
    total_delay.push_back(x);
}
void Pec::set_l1_read(int& x){
    l1_read.push_back(x);
}
void Pec::set_l1_write(int& x){
    l1_write.push_back(x);
}
*/
int Pec::get_pe_counts(){
    return pe_counts;
}
int Pec::get_l1(){
    return l1;
}
