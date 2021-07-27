#include "../include/acc.h"
Acc::Acc(const int l2_, const int l1_, \
        const int pe_t, const int aspr_r, const int aspr_c):\
                        l2(l2_),pes(pe_t),ar(aspr_r),ac(aspr_c){
        //pec_total = (int)pe_t/pec_sz;
        pes = pe_t;
        //l1 = l1_*1024/pec_total; // KB
        //l1_support = (int)l1_*1024*1024/pec_total/4; // word
        //l2_support = (int)l2_*1024*1024/4; // word
        ar = aspr_r;
        ac = aspr_c;
        if(ar<=8)
            bk = ar;
        else
            bk = 32;
        buf_bandwidth = bk * 5.6; // equals bk times 5.6GB/s = 5.6B/ns
        check_pe_total();
        /*
        for(int i=0; i<pec_total; i++){
            pec_arr.emplace_back(pec_sz,l1_);
        }*/
}
void Acc::check_pe_total(){
    if(pes%(ar*ac) != 0 )
        LOG(FATAL)<<"The given #PEs cannot be divided by #PE_SIZE";
}
