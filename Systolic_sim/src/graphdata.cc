#include "../include/graphdata.h"

void Graph::node_sort(const std::vector<std::vector<int>>& in){
    int in_counts = 0;
    for(size_t it=1;it<in.at(0).size();++it){
        //LOG(INFO) << "it = "<<it<<std::endl;
        int key = in.at(0).at(it) - in.at(0).at(it-1);
        
        //initiallizing col ID and val(NNZ of a row)  of the diagnal entry
        std::map<int,int> temp_col_and_val;
        temp_col_and_val[it-1] = key+1;

        //push back rowID
        auto pt = pre_sort_data.insert({key,temp_col_and_val});
        for(int j=0;j<key;++j){
            (*pt).second.insert({in.at(1).at(in_counts+j),in.at(2).at(in_counts+j)}); // {col ID,Val}
        }
        in_counts += key;
        //LOG(INFO) << "in_counts = "<<in_counts<<std::endl;
    }
    //copy to sorted_data
    sorted_data.at(0).push_back(0);
    int base = 0;
    for(auto it:pre_sort_data){
        base += it.first+1;
        sorted_data.at(0).push_back(base);
        //std::sort(it.second.at(0).begin(),it.second.at(0).end());
        auto ele = it.second.begin();
        while(ele != it.second.end()){
            sorted_data.at(1).push_back((*ele).first); //col ID
            sorted_data.at(2).push_back((*ele).second); // Val
            ele++;
        }
    }
    LOG(INFO) << "node_sort over.";
}

Graph::Graph(const std::string& name1):data(3,std::vector<int>()),\
                                    sorted_data(3,std::vector<int>()),\
                                    csr_diag(3,std::vector<int>()){
    
    this->data_path = name1.substr(0,name1.find_last_of("/"));                                    
    this->name = name1.substr(name1.find_last_of("/")+1,-1);                                    
    //std::cout<<this->name<<std::endl;
    //this->name = name;
    //this->data = {{},{}};
    //this->data.resize(2,std::vector<int>());
    
    this->feature_size = feature_dim[this->name];
    //this->feature_size = 16;
    //this->hiden = output_dim[this->name];
    this->hiden = 16;
    
    
    this->load_adj(this->name);    
    this->nodes = data[0].size()-1;
    node_sort(this->data);
    
    csr_diag.at(0).push_back(0);
    int base = 0;
    int new_base = 0;
    int num = 0;
    for(size_t i=1;i<this->data.at(0).size();++i){
        //int expected_sum = this->data.at(0).at(i)-this->data.at(0).at(i-1)+1;
        //VLOG(2)<<"i = "<<i;
        //VLOG(2)<<"base = "<<base;
        //VLOG(2)<<"# = "<<this->data.at(0).at(i)-this->data.at(0).at(i-1);
        std::map<int,int> temp;    
        //diagnal entries
        temp[i-1] = this->data.at(0).at(i)-this->data.at(0).at(i-1);
        //csr_diag.at(1).push_back(i-1);
        num = this->data.at(0).at(i)-this->data.at(0).at(i-1);
        if(num > this->max_num)
            this->max_num = num;
        //std::cout<<"@74 expected num: "<<expected_sum<<std::endl;
        //std::cout<<"@75 result num: "<<temp.size()<<std::endl;
        for(int j=base;j<base+this->data.at(0).at(i)-this->data.at(0).at(i-1);++j){
            temp[this->data.at(1).at(j)] = this->data.at(2).at(j);
            //csr_diag.at(1).push_back(this->data.at(1).at(j));
          //  VLOG(2)<<"this->data.at(1).at("<<j<<") = "<<this->data.at(1).at(j);
        }
        //std::cout<<"@81 expected num: "<<expected_sum<<std::endl;
        //std::cout<<"@82 result num: "<<temp.size()<<std::endl;
        for(auto tmp_entry:temp){
            csr_diag.at(1).push_back(tmp_entry.first);
            csr_diag.at(2).push_back(tmp_entry.second);
            //csr_diag.at(1).push_back(this->data.at(1).at(j));
          //  VLOG(2)<<"this->data.at(1).at("<<j<<") = "<<this->data.at(1).at(j);
        }
        //std::sort(csr_diag.at(1).begin()+base+i-1,csr_diag.at(1).begin()+base+i-1+this->data.at(0).at(i)-this->data.at(0).at(i-1)+1);
        csr_diag.at(0).push_back(new_base+temp.size());
        new_base += temp.size();
        base += this->data.at(0).at(i)-this->data.at(0).at(i-1);
        //std::cout<<"expected num: "<<expected_sum<<std::endl;
        //std::cout<<"result num: "<<temp.size()<<std::endl;
        //int result_sum = temp.size();
        //assert(expected_sum == result_sum);
    }
    this->max_num += 1;
//    print_data();
    LOG(INFO) << "graph construct over.";
    /*
    std::cout<<"NNZ in data: "<<*(this->data.at(0).end()-1)<<std::endl;
    std::cout<<"nodes: "<<this->data.at(0).size()-1<<std::endl;
    std::cout<<"NNZ in data: "<<this->data.at(1).size()<<std::endl;
    std::cout<<"NNZ in csr_diag: "<<*(this->csr_diag.at(0).end()-1)<<std::endl;
    std::cout<<"NNZ in csr_diag: "<<this->csr_diag.at(1).size()<<std::endl;
    std::cout<<"NNZ in csr_diag: "<<this->csr_diag.at(2).size()<<std::endl;
    */
}

void Graph::load_adj(const std::string& name0){
    std::fstream fin;
    fin.open(this->data_path+"/"+"n_"+name0+".csv",std::ios::in);
    //std::cout<<this->data_path<<std::endl;
    //std::cout<<name0<<std::endl;
    //std::cout<<this->data_path+"\/"+"n_"+name0+".csv"<<std::endl;
    std::string line, word;
    int i = 0;
    while(i<3){
        std::getline(fin,line);
        std::stringstream ss(line);
        while(std::getline(ss,word,',')){
            data.at(i).push_back(std::stoi(word));        
        }
        i++;
    }
    fin.close();
}
std::string Graph::get_name(){
    return this->name;
}
int Graph::get_nodes(){
    return nodes;
}
int Graph::get_feature_size(){
    return feature_size;
}
void Graph::print_data(){
    LOG(INFO) << "print start.";
    std::cout<<"The first 5 elements of rowptr: ";
    for(auto it=data.at(0).begin(); it<data.at(0).begin()+5 ; it++)
        std::cout<<(*it)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"The last 5 elements of rowptr: ";
    for(auto it=data[0].end()-5; it!=data[0].end() ; it++)
        std::cout<<(*it)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"The first 5 elements of indies: ";
    for(auto it=data[1].begin(); it<data[1].begin()+5 ; it++)
        std::cout<<(*it)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"The last 5 elements of indies: ";
    for(auto it=data[1].end()-5; it!=data[1].end() ; it++)
        std::cout<<(*it)<<" ";
    std::cout<<std::endl;

    std::cout<<"The first 5 elements of vals: ";
    for(auto it=data[2].begin(); it<data[2].begin()+5 ; it++)
        std::cout<<(*it)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"The last 5 elements of vals: ";
    for(auto it=data[2].end()-5; it!=data[2].end() ; it++)
        std::cout<<(*it)<<" ";
    std::cout<<std::endl;
   /* 
    std::cout<<"pre_sort_data: "<<std::endl;
    for(auto it:pre_sort_data){
        std::cout<<it.first<<": ";
        for(auto loc:it.second)
            std::cout<<loc<<" ";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    */
    std::cout<<"sorted_data: "<<std::endl;
    for(auto it:sorted_data){
        for(auto loc:it)
            std::cout<<loc<<" ";
        std::cout<<std::endl;
    }
    
    std::cout<<"csr_diag: "<<std::endl;
    for(auto it:csr_diag){
        for(auto loc:it)
            std::cout<<loc<<" ";
        std::cout<<std::endl;
    }
    
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Graph name: "<<get_name()<<std::endl;
    std::cout<<"The number of nodes: "<< get_nodes()<<"   Rowptr: "<<data.at(0).size()<<"   Pointer: "<<data.at(1).size()<<std::endl;
    std::cout<<"The size of a node feature: "<<get_feature_size()<<std::endl;
}
