#!/bin/bash

GLOG_vmodule=analysis=0 ./sim -d cora -c 64 > cora_64_log.txt 
GLOG_vmodule=analysis=0 ./sim -d cora -c 256 > cora_256_log.txt 


#nohup GLOG_vmodule=analysis=0 ./sim -d cora_ml -c 4 > cora_ml_4_log.txt &
#nohup GLOG_vmodule=analysis=0 ./sim -d cora_ml -c 16 > cora_ml_16_log.txt &
#nohup GLOG_vmodule=analysis=0 ./sim -d cora_ml -c 64 > cora_ml_64_log.txt &
#nohup GLOG_vmodule=analysis=0 ./sim -d cora_ml -c 256 > cora_ml_256_log.txt &

GLOG_vmodule=analysis=0 ./sim -d citeseer -c 4 > citeseer_4_log.txt 
GLOG_vmodule=analysis=0 ./sim -d citeseer -c 16 > citeseer_16_log.txt 
GLOG_vmodule=analysis=0 ./sim -d citeseer -c 64 > citeseer_64_log.txt 
GLOG_vmodule=analysis=0 ./sim -d citeseer -c 256 > citeseer_256_log.txt 

#nohup GLOG_vmodule=analysis=0 ./sim -d yelp -c 4 > yelp_4_log.txt &
#nohup GLOG_vmodule=analysis=0 ./sim -d yelp -c 16 > yelp_16_log.txt &
#nohup GLOG_vmodule=analysis=0 ./sim -d yelp -c 64 > yelp_64_log.txt &
#nohup GLOG_vmodule=analysis=0 ./sim -d yelp -c 256 > yelp_256_log.txt &

GLOG_vmodule=analysis=0 ./sim -d pubmed -c 4 > pubmed_4_log.txt 
GLOG_vmodule=analysis=0 ./sim -d pubmed -c 16 > pubmed_16_log.txt 
GLOG_vmodule=analysis=0 ./sim -d pubmed -c 64 > pubmed_64_log.txt 
GLOG_vmodule=analysis=0 ./sim -d pubmed -c 256 > pubmed_256_log.txt 

GLOG_vmodule=analysis=0 ./sim -d ppi -c 4 > ppi_4_log.txt 
GLOG_vmodule=analysis=0 ./sim -d ppi -c 16 > ppi_16_log.txt 
GLOG_vmodule=analysis=0 ./sim -d ppi -c 64 > ppi_64_log.txt 
GLOG_vmodule=analysis=0 ./sim -d ppi -c 256 > ppi_256_log.txt 
