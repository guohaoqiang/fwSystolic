#include "../include/analysis.h"

Analysis::Analysis(std::shared_ptr<Acc> accelerator,\
        std::shared_ptr<Graph> input_data, int win):\
        hw(accelerator),dat(input_data),window(win){
        
        for(size_t i=0; i<hw->bk; ++i)
            bk_cts.push_back(0);
}  

void Analysis::to_binary(TYPE_LENGTH n){
    TYPE_LENGTH a = n % 2;
    n = n >> 1;
    if(n != 0){
        to_binary(n);
    }
    std::cout<<a;
}

void Analysis::print_binary(std::shared_ptr<std::vector<TYPE_LENGTH>> &bitmap){
    for(auto item:(*bitmap)){
        to_binary(item);
        std::cout<<std::endl;
    }
}

unsigned int Analysis::find_nz_pos(TYPE_LENGTH n){ // return the first non-zero entry's pos (r->l)
// 0000 0000     return hw->ar
// 0000 0001     return 0
// 0001 0010     return 1
    unsigned int nz_pos = hw->ar;
    TYPE_LENGTH temp = n;
    for(size_t p=0; p<hw->ar; p++){
        //VLOG(2)<<"temp = "<<temp;
        //VLOG(2)<<"temp&1 = "<<(temp & 1);
        if((temp & 1) != 1){
            temp = n>>1;
            //VLOG(2)<<"temp = "<<temp;
        }else{
            nz_pos = p;
            //VLOG(2)<<"nz_pos = "<<nz_pos;
            break;
        }
    }
    return nz_pos;
}

unsigned int Analysis::count_nnz(TYPE_LENGTH n, unsigned int pos){ // #Non-zero entries previous the position "pos" (inclusive)
//0000 1000 pos:2    return 0
//0000 1000 pos:3    return 1
//1001 0100 pos:hw->ar-1    return 3
    if(pos>=hw->ar)
        LOG(FATAL)<<"Un...wrong non-zeron position"<<std::endl;
    else{
        size_t s = 0;
        unsigned int counts = 0;
        TYPE_LENGTH temp = n;
        do{
            if(temp & 1 == 1)
                counts++;
            temp = temp>>1;
            s++;
        }while(s<=pos);
        return counts; // right most bit pos has count 1
    }
}

void Analysis::tab_gen(const std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> &intile, \
        std::shared_ptr<std::vector<TYPE_LENGTH>> &d){
    int base = 0;
/*    
    for(size_t i=0; i<intile->size(); i++){ //#rows in a tile 
        VLOG(2)<<"#rows = "<<intile->size();
        d->push_back(0);
        for(size_t j=0; j<intile->at(i)->size(); j++){ //NNZ in the i-th row
            VLOG(2)<<"#NZ in the "<<j<<"-th row "<<"="<<intile->at(i)->size();
            unsigned int shifts = intile->at(i)->at(j).at(1)%hw->ar;
            //std::cout<<shifts<<std::endl;
            d->at(i) = d->at(i) | (1<<(hw->ar-shifts-1)); //convert non-zero entries in a row to a bit-line
            //std::cout<<d.at(i)<<std::endl;
        }
    }
    VLOG(2)<<"Print binary map:";
    //print_binary(d);
*/  
    for(size_t m = 0; m<intile->size(); m++){
        //unsigned int f_nz_pos = find_nz_pos(d->at(m));
        unsigned int f_nz_pos = hw->ar - (*(intile->at(m)->end()-1)).at(1)%hw->ar;
        VLOG(2)<<"";

        //self_loop->push_back(count_nnz(d->at(m),hw->ar-1));
        self_loop->push_back(intile->at(m)->size());
        
        VLOG(2)<<"";
        std::vector<int> sub_tab;
        for(size_t n = 0; n<intile->size(); n++){
            if(m!=n){
                //unsigned int s_nz_pos = find_nz_pos(d->at(n));
                unsigned int s_nz_pos = hw->ar - (*(intile->at(n)->end()-1)).at(1)%hw->ar;
                //VLOG(2)<<"d.at("<<m<<")="<<d.at(m)<<", d.at("<<n<<")="<<d.at(n);
                /*
                int diff = d->at(m) & d->at(n); 
                //VLOG(2)<<"d.at(m)&d.at(n)="<<(d.at(m)&d.at(n))<<", diff="<<diff;
                unsigned int diff_pos = find_nz_pos(diff);  // the first position that both have non-zero entries.
                */
                size_t cts_m = 0;
                size_t cts_n = 0;
                bool flag = false;
                unsigned int diff_pos = hw->ar;
                for(int i=intile->at(m)->size()-1; i>=0; --i){
                    for(int j=intile->at(n)->size()-1; j>=0; --j){
                        if(intile->at(m)->at(i).at(1)%hw->ar == intile->at(n)->at(j).at(1)%hw->ar){
                            cts_m = intile->at(m)->size() - i;
                            cts_n = intile->at(n)->size() - j;
                            diff_pos = hw->ar - intile->at(m)->at(i).at(1)%hw->ar - 1;
                            flag = true;
                            break;
                        }
                    }
                    if(flag)
                        break;
                }
                
                //return;
                VLOG(2)<<"diff_pos = "<<diff_pos<<"   f_nz_pos = "<<f_nz_pos;
                if(cts_n != 0){
                    //int make_up_0 = count_nnz(d->at(n),hw->ar-1)-count_nnz(d->at(n),diff_pos)\
                            - (count_nnz(d->at(m),hw->ar-1)-count_nnz(d->at(m),diff_pos)+1);
                    //int make_up_0 = intile->at(n)->size() - cts_n\
                            - (intile->at(n)->size() - cts_m + 1)

                    //int make_up = make_up_0>0?make_up_0:0;
                    
                    //sub_tab.push_back(count_nnz(d->at(n),f_nz_pos)+make_up);
                    int opt = intile->at(n)->size() - intile->at(m)->size();
                    sub_tab.push_back(std::max((int)(cts_n-cts_m+1),opt));
                    
                    VLOG(2)<<"m = "<<m<<" n = "<<n<<" , "<<sub_tab.at(sub_tab.size()-1);
                    //tab.at(m).at(n) = count_nnz(d.at(n),f_nz_pos);
                }else if(cts_n == 0){ // no identical non-zero position
                    if(intile->at(m)->size()>=intile->at(n)->size()){
                        sub_tab.push_back(0);
                    }else{
                        sub_tab.push_back(intile->at(n)->size() - intile->at(m)->size());
                    //tab.at(m).at(n) = 1;
                    }
                    VLOG(2)<<"m = "<<m<<" n = "<<n<<" , "<<sub_tab.at(sub_tab.size()-1);
                }else{
                    LOG(FATAL)<<"Un...wrong non-zeron position"<<std::endl;
                }
            }else{
                sub_tab.push_back(-FLT_MAX);
                //tab.at(m).at(n) = 0;
            }
            VLOG(2)<<"sub_tab length: "<<sub_tab.size()<<" sub_tab: "<<sub_tab.at(sub_tab.size()-1);
            //return;
        }
        tab->push_back(sub_tab);
    
    }
}


void Analysis::print(){
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

long Analysis::val(){
    int base = 0;
    std::vector<int> mark(hw->ar,0);
    std::vector<int> count(hw->ar,0);
    //the inner most vector denotes an non-zero entry <row,col,val>

    
    long max_col_tile_delay = 0; 
    for(size_t i=base; i<dat->csr_diag.at(0).size()-1; i += hw->ar){ // iterate over col-tiles
        std::cout<<"col tile: "<<i<<std::endl;
        for(size_t j=i; j<std::min(i+hw->ar,dat->csr_diag.at(0).size()-1); j++){ // iterate over cols in a col-tile
            mark.at(j%hw->ar) = dat->csr_diag.at(0).at(j);
            //VLOG(2)<<"mark["<<j%hw->ar<<"]="<<dat->data.at(0).at(j);
        }
        count.assign(hw->ar,0);
        long col_tile_delay = 0; 
        // in adj, #row == #col
        for(size_t k1=0; k1<dat->csr_diag.at(0).size()-1; k1 += window){ // iterate over row-tiles
            VLOG(2)<<"row tile: "<<k1;
            std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> tile = \
                                    std::make_shared<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>>(); 
        //for(size_t k1=0; k1<3; k1 += window){ // iterate over row-tiles
            for(size_t k2=k1; k2<std::min(k1+window,dat->csr_diag.at(0).size()-1); k2++){
                std::shared_ptr<std::vector<std::vector<int>>> tmp_tile_row = \
                                        std::make_shared<std::vector<std::vector<int>>>();
                // begin fetch one row of a tile
                for(size_t j=i; j<std::min(i+hw->ar,dat->csr_diag.at(0).size()-1); j++){ // iterate over cols in a col-tile
                   //std::cout<<"k2 = "<<k2<<" j = "<<j<<std::endl;
                   //VLOG(2)<<"#nz of "<<j<<":"<<dat->data.at(0).at(j+1) - dat->data.at(0).at(j);
                   if(dat->csr_diag.at(1).at(mark.at(j%hw->ar)) == k2){ 
                        std::vector<int> tmp_entry(3,0);
                        tmp_entry.at(0) = dat->csr_diag.at(1).at(mark.at(j%hw->ar));//row ID
                        tmp_entry.at(1) = j;//col ID
                        tmp_entry.at(2) = dat->csr_diag.at(2).at(mark.at(j%hw->ar));//vals
                        tmp_tile_row->push_back(tmp_entry);
                        count.at(j%hw->ar)++;
                        if(count.at(j%hw->ar)<(dat->csr_diag.at(0).at(j+1) - dat->csr_diag.at(0).at(j)))
                            mark.at(j%hw->ar)++;
                        VLOG(2)<<"k2 = "<<k2<<" j = "<<j<<" r="<<tmp_entry.at(0)<<" c="<<tmp_entry.at(1)<<" v="<<tmp_entry.at(2);
                   }
                } //end fetch one row of a tile
                //std::cout<<std::endl;
                if(tmp_tile_row->size() != 0)
                    tile->push_back(tmp_tile_row);
            }// end fetch one tile 
            //std::cout<<"-----------"<<std::endl; 
            std::shared_ptr<std::vector<TYPE_LENGTH>> bits = std::make_shared<std::vector<TYPE_LENGTH>>(); 
            //for(auto g:bits)
            //    std::cout<<g<<std::endl;
            tab_gen(tile,bits);
            //print_binary(bits);
           
            // LinDA w/ mip 
            std::shared_ptr<std::vector<int>> uni_path = mip(tab,self_loop); 

            permutate(tile,uni_path);    
            //print_tile(tile);
           
            bits->clear(); 
            
            col_tile_delay += debubbling(tile,bits); 

            deconstruct_tab();
            deconstruct_self_loop();
        } // end visiting all row-tiles of a specific col-tile
        comp_delays += col_tile_delay;
        if(col_tile_delay > max_col_tile_delay)
            max_col_tile_delay = col_tile_delay;
    } // end visiting all col-tiles
    std::cout<<"total delays: "<<comp_delays<<std::endl;
    
    long long nnz_total = 0;
    for(size_t i=0; i<hw->bk; ++i){
        nnz_total += bk_cts.at(i);
    }
    float shannon_entropy = 0;
    for(size_t i=0; i<hw->bk; ++i){
        float temp = (float)bk_cts.at(i)/nnz_total;
        shannon_entropy -= temp*log2(temp)/log2(hw->bk);
    }
    std::cout<<"conflicts: "<<conflicts<<std::endl;
    std::cout<<"shannon entropy: "<<shannon_entropy<<std::endl;
    
    return max_col_tile_delay; 
}
void Analysis::print_tile(const std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> &intile){
    for(auto row:(*intile)){
        for(auto entry:(*row)){
            std::cout<<"r = "<<entry.at(0)<<"c = "<<entry.at(1)<<"val = "<<entry.at(2);
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}


void Analysis::permutate(std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> &intile,\
        std::shared_ptr<std::vector<int>> &vec){
    assert(intile->size() == vec->size());
            
    for( size_t vv = 0; vv < intile->size() - 1; ++vv )
    {
            if (vec->at(vv) == vv)
            {
                continue;
            }
            size_t oo;
            for(oo = vv + 1; oo <vec->size(); ++oo)
            {
                if (vec->at(oo) == vv)
                {
                    break;
                }
            }
            std::iter_swap(intile->begin()+vv, intile->begin()+vec->at(vv) );
            std::iter_swap(vec->begin()+vv, vec->begin()+oo);
    } 
}



long Analysis::debubbling(const std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> &intile, \
        std::shared_ptr<std::vector<TYPE_LENGTH>> &d){
    std::cout<<"start debubbling ... "<<std::endl;
    int base = 0;
    long t = 0;
  /*  
    for(size_t i=0; i<intile->size(); i++){ //#rows in a tile 
//        VLOG(2)<<"#rows = "<<intile->size();
        d->push_back(0);
        for(size_t j=0; j<intile->at(i)->size(); j++){ //NNZ in the i-th row
            VLOG(2)<<"#NZ in the "<<i<<"-th row "<<"="<<intile->at(i)->size();
            unsigned int shifts = intile->at(i)->at(j).at(1)%hw->ar;
            //std::cout<<shifts<<std::endl;
            d->at(i) = d->at(i) | (1<<(hw->ar-shifts-1)); //convert non-zero entries in a row to a bit-line
            //std::cout<<d.at(i)<<std::endl;
        }
    }
    VLOG(2)<<"Print binary map:";
    print_binary(d);
   */ 
    std::shared_ptr<vector<std::shared_ptr<std::vector<std::vector<int>>>>> data_flow_in = \
                            std::make_shared<vector<std::shared_ptr<std::vector<std::vector<int>>>>>();
    for(size_t l=0;l<hw->ar;++l){
        data_flow_in->push_back(std::make_shared<std::vector<std::vector<int>>>());
    }
    std::vector<int> sentry(hw->ar,0); // in the same order with "data_flow_in". Top->Down: 0,1,2,...hw->ar-1
    for(size_t m = 0; m<intile->size(); m++){
        VLOG(2)<<"m = "<<m; 
        unsigned int first_y = intile->at(m)->at(0).at(1)%hw->ar;
        /*
        for(size_t s=0;s<hw->ar;++s){ //"first" here refers to from top to down. reverse to the order of bit representation d
            if((d->at(m)>>hw->ar-1-s) & 1){ //identical 1
                first_y = s;
                break;
            }
        }
        */
        unsigned int first_x = data_flow_in->at(first_y)->size();
        VLOG(2)<<"first_x = "<<first_x<<"   first_y = "<<first_y;

        VLOG(2)<<"intile.size = "<<intile->size()<<"   (m).size() = "<<intile->at(m)->size();
        unsigned int max_y = (unsigned)intile->at(m)->at(0).at(1)%hw->ar;
        unsigned int max_x = data_flow_in->at((unsigned)intile->at(m)->at(0).at(1)%hw->ar)->size();
        unsigned int max_pos = 0;
        unsigned int nz_indx = 0;
        bool flag = false;
        for(size_t s=1; s<intile->at(m)->size(); ++s){
            nz_indx = (unsigned)intile->at(m)->at(s).at(1)%hw->ar; //col ID
            VLOG(2)<<"nz_idx="<<nz_indx<<" size="<<data_flow_in->at(nz_indx)->size()<<" max_x="<<max_x;
            if(data_flow_in->at(nz_indx)->size()>max_x){ //!!! IMPORTANT: no identical sign
                max_x = data_flow_in->at(nz_indx)->size(); //available position index
                max_y = nz_indx;
                max_pos = s;
            }
            /*
            VLOG(2)<<"nz_idx="<<nz_indx<<" size="<<data_flow_in->at(nz_indx)->size()<<" max_x="<<max_x;
            if(data_flow_in->at(nz_indx)->size() == max_x && !flag){ 
                max_x = data_flow_in->at(nz_indx)->size(); //available position index
                max_y = nz_indx;
                max_pos = s;
                flag = true;
                //VLOG(2)<<"s = "<<s;
            }
            */
        }
        VLOG(2)<<"max_x = "<<max_x<<"   max_y = "<<max_y;
   
        std::vector<int> bubble(3,-1); //Placeholder. {row_ID,col_ID,Val} = {-1,-1,-1}.
        if((max_y-first_y) > (max_x-first_x)){ //>45 degree
            //top->down;
            data_flow_in->at(first_y)->push_back(intile->at(m)->at(0));//push back the first non-zero entry
            sentry.at(first_y)++;
            VLOG(2)<<" sentry["<<first_y<<"]="<<sentry.at(first_y);
            for(size_t idx=1; idx<intile->at(m)->size();++idx){ //iterate over the remaining non-zero entries in the same row.
                size_t i_y = intile->at(m)->at(idx).at(1)%hw->ar; //col ID
                VLOG(2)<<intile->at(m)->at(idx).at(0)%hw->ar; //col ID
                VLOG(2)<<" sentry["<<i_y<<"]="<<sentry.at(i_y)<<"<?"<<sentry.at(first_y)+idx-1;
                VLOG(2)<<" sentry["<<first_y<<"]="<<sentry.at(first_y)<<" idx = "<<idx;
                while(sentry.at(i_y)<(int)(sentry.at(first_y)+idx-1)){
                    data_flow_in->at(i_y)->push_back(bubble);
                    sentry.at(i_y)++;
                    //VLOG(2)<<" sentry["<<i_y<<"]="<<sentry.at(i_y)<<" sentry["<<first_y<<"]="<<sentry.at(first_y)<<" idx = "<<idx;
                }
                data_flow_in->at(i_y)->push_back(intile->at(m)->at(idx));
                sentry.at(i_y)++;
                VLOG(2)<<" sentry["<<i_y<<"]="<<sentry.at(i_y);
            }
            VLOG(2)<<"";
        }else if((max_y-first_y) <= (max_x-first_x)){ //<=45 degree
            data_flow_in->at(max_y)->push_back(intile->at(m)->at(max_pos)); //max_pos in line 232,238
            sentry.at(max_y)++;
            VLOG(2)<<" sentry["<<max_y<<"]="<<sentry.at(max_y);
            //elements above max_y: down->top
            for(size_t idx=0; idx<max_pos; ++idx){ //iterate over the remaining non-zero entries in the same row.
                VLOG(2)<<"";
                size_t i_y = intile->at(m)->at(idx).at(1)%hw->ar; //col ID
                VLOG(2)<<" sentry["<<i_y<<"]="<<sentry.at(i_y);
                VLOG(2)<<" sentry["<<max_y<<"]="<<sentry.at(max_y)<<" max_pos = "<<max_pos<<" idx = "<<idx;
                while(sentry.at(i_y)<std::max((int)(sentry.at(max_y)-(max_pos-idx)-1),0)){
                    data_flow_in->at(i_y)->push_back(bubble);
                    sentry.at(i_y)++;
                }
                data_flow_in->at(i_y)->push_back(intile->at(m)->at(idx));
                sentry.at(i_y)++;
                VLOG(2)<<" sentry["<<i_y<<"]="<<sentry.at(i_y);
            }

            //elements below max_y: top->down
            VLOG(2)<<"";
            for(size_t idx=0; idx<intile->at(m)->size()-max_pos-1;++idx){ //iterate over the remaining non-zero entries in the same row.
                size_t i_y = intile->at(m)->at(idx+max_pos+1).at(1)%hw->ar; //col ID
                while(sentry.at(i_y)<(sentry.at(max_y)+idx)){
                    data_flow_in->at(i_y)->push_back(bubble);
                    sentry.at(i_y)++;
                }
                data_flow_in->at(i_y)->push_back(intile->at(m)->at(idx+max_pos+1));
                sentry.at(i_y)++;
                VLOG(2)<<" sentry["<<i_y<<"]="<<sentry.at(i_y);
            }
            VLOG(2)<<"";
        }
        for(auto p:(*data_flow_in)){
            if(p->size()>t)
                t = p->size();
        }
        VLOG(2)<<"t: "<<t;
    
    }
    std::cout<<"end debubbling ... "<<std::endl;
    //linda_mem_analysis(data_flow_in);
    c2sr_mem_analysis(data_flow_in);
    return t;
}

void Analysis::linda_mem_analysis(const std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> &dataflow_in){
    long long length = 0;
    for(auto p:(*dataflow_in)){
        if(p->size()>length)
            length = p->size();
    }
    long long counts = -1;
    std::set<int> bk_sets;
    for(size_t i=0; i<length; ++i){ 
        int nnz = 0;
        for(size_t j=0; j<dataflow_in->size(); ++j){ 
            if(i < dataflow_in->at(j)->size() && dataflow_in->at(j)->at(i).at(0) != -1){
                counts++;
                nnz++;
                bk_sets.insert(counts%hw->bk);
                bk_cts.at(counts%hw->bk)++;
            }else{
                continue;
            }
        }
        conflicts += (nnz-bk_sets.size());
        bk_sets.clear();
    }
    /*
    long long nnz_total = 0;
    for(size_t i=0; i<hw->ar; ++i){
        nnz_total += bk_cts.at(i);
    }
    float shannon_entropy = 0;
    for(size_t i=0; i<hw->ar; ++i){
        float temp = (float)bk_cts.at(i)/nnz_total;
        shannon_entropy += temp*log2(temp);
    }
    std::cout<<"conflicts: "<<conflicts<<std::endl;
    std::cout<<"shannon entropy: "<<shannon_entropy<<std::endl;
    */
}

void Analysis::c2sr_mem_analysis(const std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> &dataflow_in){
    long long length = 0;
    for(auto p:(*dataflow_in)){
        if(p->size()>length)
            length = p->size();
    }
    
    std::set<int> bk_sets;
    for(size_t i=0; i<length; ++i){ 
        int nnz = 0;
        for(size_t j=0; j<dataflow_in->size(); ++j){ 
            if(i < dataflow_in->at(j)->size() && dataflow_in->at(j)->at(i).at(0) != -1){
                nnz++;
                bk_sets.insert(dataflow_in->at(j)->at(i).at(0)%hw->bk);
                bk_cts.at(dataflow_in->at(j)->at(i).at(0)%hw->bk)++;
            }
        }
        conflicts += (nnz - bk_sets.size());
        bk_sets.clear();
    }
}
void Analysis::val_naive1(){
    int base = 0;
    std::vector<int> mark(hw->ar,0);
    std::vector<int> count(hw->ar,0);
    //the inner most vector denotes an non-zero entry <row,col,val>

    
    for(size_t i=base; i<dat->csr_diag.at(0).size()-1; i += hw->ar){ // iterate over col-tiles
        std::cout<<"col tile: "<<i<<std::endl;
        for(size_t j=i; j<std::min(i+hw->ar,dat->csr_diag.at(0).size()-1); j++){ // iterate over cols in a col-tile
            mark.at(j%hw->ar) = dat->csr_diag.at(0).at(j);
            //VLOG(2)<<"mark["<<j%hw->ar<<"]="<<dat->data.at(0).at(j);
        }
        count.assign(hw->ar,0);
        // in adj, #row == #col
        std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> tile = \
                                std::make_shared<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>>(); 
    //for(size_t k1=0; k1<3; k1 += window){ // iterate over row-tiles
        for(size_t k2=0; k2<dat->csr_diag.at(0).size()-1; k2++){
            std::shared_ptr<std::vector<std::vector<int>>> tmp_tile_row = \
                                    std::make_shared<std::vector<std::vector<int>>>();
            // begin fetch one row of a tile
            for(size_t j=i; j<std::min(i+hw->ar,dat->csr_diag.at(0).size()-1); j++){ // iterate over cols in a col-tile
               //std::cout<<"k2 = "<<k2<<" j = "<<j<<std::endl;
               //VLOG(2)<<"#nz of "<<j<<":"<<dat->data.at(0).at(j+1) - dat->data.at(0).at(j);
               if(dat->csr_diag.at(1).at(mark.at(j%hw->ar)) == k2){ 
                    std::vector<int> tmp_entry(3,0);
                    tmp_entry.at(0) = dat->csr_diag.at(1).at(mark.at(j%hw->ar));//row ID
                    tmp_entry.at(1) = j;//col ID
                    tmp_entry.at(2) = dat->csr_diag.at(2).at(mark.at(j%hw->ar));//vals
                    tmp_tile_row->push_back(tmp_entry);
                    count.at(j%hw->ar)++;
                    if(count.at(j%hw->ar)<(dat->csr_diag.at(0).at(j+1) - dat->csr_diag.at(0).at(j)))
                        mark.at(j%hw->ar)++;
                    VLOG(2)<<"k2 = "<<k2<<" j = "<<j<<" r="<<tmp_entry.at(0)<<" c="<<tmp_entry.at(1)<<" v="<<tmp_entry.at(2);
               }
            } //end fetch one row of a tile
            //std::cout<<std::endl;
            if(tmp_tile_row->size() != 0)
                tile->push_back(tmp_tile_row);
        }// end fetch one tile 
        //std::cout<<"-----------"<<std::endl; 
        std::shared_ptr<std::vector<TYPE_LENGTH>> bits = std::make_shared<std::vector<TYPE_LENGTH>>(); 
        //for(auto g:bits)
        //    std::cout<<g<<std::endl;
        //print_binary(bits);
       
        // baseline w/o mip
        comp_delays += debubbling(tile,bits); 
        //comp_delays += self_loop->at(0);
        deconstruct_tab();
        deconstruct_self_loop();
    } // end visiting all col-tiles
    std::cout<<"total delays: "<<comp_delays<<std::endl;
}

void Analysis::deconstruct_tab(){
    for(auto i=0;i<tab->size();i++)
        tab->at(i).clear();
    tab->clear();
}

void Analysis::deconstruct_self_loop(){
    self_loop->clear();
}

void Analysis::run_wegnn(){
    
    VLOG(0)<<"Accelerator: ";
    VLOG(0)<<"Total PEs: "<<hw->pes;
    VLOG(0)<<"Aspect Ratio: "<<hw->ar<<"*"<<hw->ac;
    VLOG(0)<<"#PEs in each PE cluster: "<<hw->pec_size;
    VLOG(0)<<"L1 in each PE cluster: "<<hw->l1<<" KB";
    VLOG(0)<<"L2 on chip: "<<hw->l2<<" MB";
    //VLOG(0)<<"Dataset: "<<dat->name<<"\-"<<dat->nodes<<","<<dat->max_num<<","<<dat->sorted_data.at(0).at(1)-dat->sorted_data.at(0).at(0);
    VLOG(0)<<"Dataset: "<<dat->name<<"-"<<dat->nodes<<","<<dat->feature_size<<","<<dat->hiden;

}

void Analysis::run_baseline(){
    
    VLOG(0)<<"Accelerator: ";
    VLOG(0)<<"Total PEs: "<<hw->pes;
    VLOG(0)<<"Aspect Ratio: "<<hw->ar<<"*"<<hw->ac;
    VLOG(0)<<"#PEs in each PE cluster: "<<hw->pec_size;
    VLOG(0)<<"L1 in each PE cluster: "<<hw->l1<<" KB";
    VLOG(0)<<"L2 on chip: "<<hw->l2<<" MB";
    //VLOG(0)<<"Dataset: "<<dat->name<<"\-"<<dat->nodes<<","<<dat->max_num;
    VLOG(0)<<"Dataset: "<<dat->name<<"-"<<dat->nodes<<","<<dat->feature_size<<","<<dat->hiden;

}
void Analysis::val_naive(){
    int base = 0;
    std::vector<int> mark(hw->ar,0);
    std::vector<int> count(hw->ar,0);
    //the inner most vector denotes an non-zero entry <row,col,val>

    
    for(size_t i=base; i<dat->csr_diag.at(0).size()-1; i += hw->ar){ // iterate over col-tiles
        std::cout<<"col tile: "<<i<<std::endl;
        for(size_t j=i; j<std::min(i+hw->ar,dat->csr_diag.at(0).size()-1); j++){ // iterate over cols in a col-tile
            mark.at(j%hw->ar) = dat->csr_diag.at(0).at(j);
            //VLOG(2)<<"mark["<<j%hw->ar<<"]="<<dat->data.at(0).at(j);
        }
        count.assign(hw->ar,0);
        // in adj, #row == #col
        std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>> tile = \
                                std::make_shared<std::vector<std::shared_ptr<std::vector<std::vector<int>>>>>(); 
    //for(size_t k1=0; k1<3; k1 += window){ // iterate over row-tiles
        for(size_t k2=0; k2<dat->csr_diag.at(0).size()-1; k2++){
            std::shared_ptr<std::vector<std::vector<int>>> tmp_tile_row = \
                                    std::make_shared<std::vector<std::vector<int>>>();
            // begin fetch one row of a tile
            for(size_t j=i; j<std::min(i+hw->ar,dat->csr_diag.at(0).size()-1); j++){ // iterate over cols in a col-tile
               //std::cout<<"k2 = "<<k2<<" j = "<<j<<std::endl;
               //VLOG(2)<<"#nz of "<<j<<":"<<dat->data.at(0).at(j+1) - dat->data.at(0).at(j);
               if(dat->csr_diag.at(1).at(mark.at(j%hw->ar)) == k2){ 
                    std::vector<int> tmp_entry(3,0);
                    tmp_entry.at(0) = dat->csr_diag.at(1).at(mark.at(j%hw->ar));//row ID
                    tmp_entry.at(1) = j;//col ID
                    tmp_entry.at(2) = dat->csr_diag.at(2).at(mark.at(j%hw->ar));//vals
                    tmp_tile_row->push_back(tmp_entry);
                    count.at(j%hw->ar)++;
                    if(count.at(j%hw->ar)<(dat->csr_diag.at(0).at(j+1) - dat->csr_diag.at(0).at(j)))
                        mark.at(j%hw->ar)++;
                    VLOG(2)<<"k2 = "<<k2<<" j = "<<j<<" r="<<tmp_entry.at(0)<<" c="<<tmp_entry.at(1)<<" v="<<tmp_entry.at(2);
               }
            } //end fetch one row of a tile
            //std::cout<<std::endl;
            if(tmp_tile_row->size() != 0)
                tile->push_back(tmp_tile_row);
        }// end fetch one tile 
        //std::cout<<"-----------"<<std::endl; 
        std::shared_ptr<std::vector<TYPE_LENGTH>> bits = std::make_shared<std::vector<TYPE_LENGTH>>(); 
        //for(auto g:bits)
        //    std::cout<<g<<std::endl;
        //print_binary(bits);
       
        // baseline w/o mip
        comp_delays += tile->size(); 
        //comp_delays += self_loop->at(0);
        deconstruct_tab();
        deconstruct_self_loop();
    } // end visiting all col-tiles
    std::cout<<"total delays: "<<comp_delays + hw->ar<<std::endl;
}
void Analysis::timing(){
    //Stage 1
    //channel width = 8 Byte, so a load transaction read 4 entries, each of which is 2 Bytes.
    //row-major
    float t_stage_load_Col_tile_w = (dat->feature_size*ceil(hw->ac/4)*8)/hw->dram_bandwidth; //in the units of ns 
    float t_stage_1 = t_stage_load_Col_tile_w;
    
    
    //Stage 2
    //channel width = 4 Bytes. An operand is 2 Bytes.
    float t_load_from_buf_x_one_cycle = ceil(hw->ar/hw->bk/2)*32*hw->bk/hw->buf_bandwidth;
    float t_load_from_buf_w_one_cycle = ceil(hw->ac/hw->bk/2)*32*hw->bk/hw->buf_bandwidth;

    float t_stage_2 = (dat->feature_size + hw->ar + hw->ac) * \
            (t_load_from_buf_x_one_cycle > hw->mac_latency ? t_load_from_buf_x_one_cycle : hw->mac_latency); 
    
    //Stage 3
    float t_stage_3 = ceil(hw->act_latency/hw->mac_latency)*hw->ar;

    //Stage 4
    float t_stage_4 = val();

    float max1 = (t_stage_1 > t_stage_2)? t_stage_1 : t_stage_2;
    float max2 = (t_stage_3 > t_stage_4)? t_stage_3 : t_stage_4;
    float max_stage = (max1 > max2) ? max1 : max2;
    
    float stage_1_stall = abs(max_stage - t_stage_1);
    float stage_2_stall = abs(max_stage - t_stage_2);
    float stage_3_stall = abs(max_stage - t_stage_3);
    float stage_4_stall = abs(max_stage - t_stage_4);


    float t_sub = (ceil(dat->hiden/hw->ac) + 3)*max_stage; // in us
    float t_sub_stall = (stage_1_stall + stage_2_stall +stage_3_stall + stage_4_stall)*\
                        ceil(dat->hiden/hw->ac);

    float t_total = ceil(dat->nodes/hw->ar)*t_sub;
    float t_total_stall = ceil(dat->nodes/hw->ar)*t_sub_stall;


    std::cout<<"Hardware: "<<hw->pes<<" PEs/MACs |"<<" Aspect Ratio: "<<hw->ar<<"X"<<hw->ac<<std::endl;

    std::cout<<"Stage1: "<<t_stage_1<<"  ns  "<<std::endl;
    std::cout<<"Stage2: "<<t_stage_2<<"  ns  "<<std::endl;
    std::cout<<"Stage3: "<<t_stage_3<<"  ns  "<<std::endl;
    std::cout<<"Stage4: "<<t_stage_4<<"  ns  "<<std::endl;
    std::cout<<"t_sub:  "<<t_sub<<"      ns  "<<std::endl;
    
    
    std::cout<<"Total runtime:  "<<t_total/1000000<<"  ms  "<<std::endl;
    std::cout<<"Total sync:  "<<t_total_stall/1000000<<"  ms  "<<std::endl;

    //Only counts dynamic energy
    float e_dram = hw->dram_power * (ceil(dat->nodes/hw->ar)*hw->ar*dat->feature_size*2/hw->dram_bandwidth \
            + t_stage_1 * ceil(dat->hiden/hw->ac) * ceil(dat->nodes/hw->ar)*0.000000001 \
            + ceil(dat->hiden/hw->ac) * ceil(dat->nodes/hw->ar) * dat->nodes * hw->ar * 2 * 2 /hw->dram_bandwidth); 

    // Mu & Sigma
    float e_activation = 2*hw->ac * hw->activation_power * ceil(dat->hiden/hw->ac) * ceil(dat->nodes/hw->ar) * t_stage_2;
    float e_dense = 2*hw->pes * hw->mac_power * ceil(dat->hiden/hw->ac) * ceil(dat->nodes/hw->ar) * t_stage_3;
    float e_sparse = 2*hw->pes * hw->pe_power * ceil(dat->hiden/hw->ac) * ceil(dat->nodes/hw->ar) * t_stage_4;

    float e_total = e_dram + e_activation + e_dense + e_sparse;

    std::cout<<"DRAM Energy: "<<e_dram/1000<<"  nj  "<<std::endl;
    std::cout<<"Dense Energy: "<<e_dense/1000<<"  nj  "<<std::endl;
    std::cout<<"Activation Energy: "<<e_activation/1000<<"  nj  "<<std::endl;
    std::cout<<"Sparse Energy: "<<e_sparse/1000<<"  nj  "<<std::endl;
    std::cout<<"Total Energy:  "<<e_total/1000<<"      nj  "<<std::endl;
    
    std::cout<<"Theoretical DRAM Bandwidth: 94 GB/s"<<std::endl;
    // the second to last item: col ID and Vals
    float real_dram_bandwidth = (dat->nodes*dat->feature_size*2 \
            + ceil(dat->nodes/hw->ar)*dat->feature_size*dat->hiden*2 \
            + (dat->csr_diag.at(0).size()+dat->csr_diag.at(1).size()+dat->csr_diag.at(2).size())*2 \
            + ceil(dat->nodes/hw->ar)*dat->nodes * dat->hiden * 2 * 2)/t_total;
    std::cout<<"Achieved DRAM Bandwidth: "<<real_dram_bandwidth<<" GB/s"<<std::endl;

    bool baseline = false;
    int coefficient = 0;
    if(baseline)
        coefficient = 1;
    else
        coefficient = 2;
    std::cout<<"Theoretical on-chip buffer Bandwidth: "<<hw->buf_bandwidth<<" GB/s"<<std::endl;
    //last two items: col ID and Vals, //read & write
    float real_buffer_bandwidth = (ceil(dat->hiden/hw->ac)*dat->nodes*dat->feature_size*2 \
            + ceil(dat->nodes/hw->ar)*dat->feature_size*dat->hiden*2 \
            + ceil(dat->hiden/hw->ac)*dat->csr_diag.at(1).size()*2*coefficient \
            + ceil(dat->nodes/hw->ar)*dat->nodes * dat->hiden * 2 * 2)/t_total;  
    std::cout<<"Achieved on-chip buffer Bandwidth: "<<real_buffer_bandwidth<<" GB/s"<<std::endl;
}
