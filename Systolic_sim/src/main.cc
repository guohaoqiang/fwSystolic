#include "../include/main.h"

int main(int argc, char *argv[])
{
  optparse::OptionParser parser = optparse::OptionParser();
  parser.add_option("-n", "--pe_counts").type("int").set_default("16").help("Total number of PEs. n == r * c");
  parser.add_option("-r", "--ar").type("int").set_default("1").help("Aspect ratio R.");
  parser.add_option("-c", "--ac").type("int").set_default("1").help("Aspect ratio C.");
  parser.add_option("-d", "--data_name").type("string").set_default("../data/test_data").help("The name of datasets.");
  parser.add_option("-l1", "--l1_size").type("int").set_default("4").help("The size of L1 Cache for all PE cluster (in MB).");
  parser.add_option("-l2", "--l2_size").type("int").set_default("8").help("The size of L2 (in MB).");
  parser.add_option("-t", "--type").type("string").set_default("b").help("(b)gcnax OR (w)wegnn");
  
  optparse::Values& option = parser.parse_args(argc,argv);
  std::shared_ptr<Graph> adj_csr = std::make_shared<Graph>((std::string)option.get("data_name"));  
  adj_csr->print_data();
  
  //std::shared_ptr<Acc> acc = std::make_shared<Acc>((int)option.get("l2_size"),\
              (int)option.get("l1_size"),(int)option.get("pe_counts"),\
              (int)option.get("ar"),\
              (int)option.get("ac"));
  
  //Analysis res(acc,adj_csr);

  /*
  std::string my_kernel((std::string)option.get("type"));
  if(!my_kernel.compare("w")){
    res.run_wegnn();
  //res.run_w_order_wo_move(3,1,1,1);
  //res.run_w_order_wo_move(1,85,1,1);
  }else{
    res.run_baseline();
  //  res.run_baseline_awbgcn_gcnax(1,2,1,1);
  }
  */
 


  return 0;

}
