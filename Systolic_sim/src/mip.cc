#include "../include/mip.h"

//string itos(int i) {stringstream s; s << i; return s.str(); }
// Given an integer-feasible solution 'sol', find the smallest
// sub-tour.  Result is returned in 'tour', and length is
// returned in 'tourlenP'.

void
findsubtour(int      n,
            double** sol,
            int*     tourlenP,
            int*     tour)
{
  VLOG(2)<<"start tracking the tour...";
  bool* seen = new bool[n];
  int bestind, bestlen;
  int i, node, len, start;

  for (i = 0; i < n; i++)
    seen[i] = false;

  start = 0;
  bestlen = n+1;
  bestind = -1;
  node = 0;
  while (start < n) {
    for (node = 0; node < n; node++)
      if (!seen[node])
        break;
    if (node == n)
      break;
    for (len = 0; len < n; len++) {
      tour[start+len] = node;
      seen[node] = true;
      for (i = 0; i < n; i++) {
        if (sol[node][i] > 0.5 && !seen[i]) {
          node = i;
          break;
        }
      }
      if (i == n) {
        len++;
        if (len < bestlen) {
          bestlen = len;
          bestind = start;
        }
        start += len;
        break;
      }
    }
  }

  for (i = 0; i < bestlen; i++)
    tour[i] = tour[bestind+i];
  *tourlenP = bestlen;

  delete[] seen;
}


std::shared_ptr<std::vector<int>> mip(shared_ptr<vector<vector<int>>> &table, \
        shared_ptr<vector<unsigned int>> &self_loops){
  
  VLOG(2)<<"start computing mip...";
  int n = 2*table->size();
  int i;
  //long mini_delay = 0;
  shared_ptr<vector<int>> uni_tour = make_shared<vector<int>>();

  GRBEnv *env = NULL;
  GRBVar **vars = NULL;

  vars = new GRBVar*[n];
  for (i = 0; i < n; i++)
    vars[i] = new GRBVar[n];
  try {
    int j;

    env = new GRBEnv();
    GRBModel model = GRBModel(*env);

    // Must set LazyConstraints parameter when using lazy constraints

    model.set(GRB_IntParam_LazyConstraints, 1);

    // Create binary decision variables
/*
    double table[8][8] = {{100,100,100,100,-100,0,0,1},\
                          {100,100,100,100,2,-100,1,1},\
                          {100,100,100,100,2,1,-100,1},\
                          {100,100,100,100,0,1,1,-100},\
                          {-100,2,2,0,100,100,100,100},\
                          {0,-100,1,1,100,100,100,100},\
                          {0,1,-100,1,100,100,100,100},\
                          {1,1,1,-100,100,100,100,100}};
*/
    for (i = 0; i < n; i++) {
      for (j = 0; j <=i; j++) {
        //vars[i][j] = model.addVar(0.0, 1.0, distance(x, y, i, j),
        //                          GRB_BINARY, "x_"+itos(i)+"_"+itos(j));
        if(i>=table->size() && j<table->size()){
            vars[i][j] = model.addVar(0.0, 1.0, (float)table->at(i%table->size()).at(j),
                                  GRB_BINARY, "x_"+to_string(i)+"_"+to_string(j));
            VLOG(2)<<"table["<<i<<"]["<<j<<"]="<<table->at(i%table->size()).at(j);
        }else{
            vars[i][j] = model.addVar(0.0, 1.0, M,
                                  GRB_BINARY, "x_"+to_string(i)+"_"+to_string(j));
            VLOG(2)<<"table["<<i<<"]["<<j<<"]="<<"M";
        }
        vars[j][i] = vars[i][j];
      }
    }

    // Degree-2 constraints

    for (i = 0; i < n; i++) {
      GRBLinExpr expr = 0;
      for (j = 0; j < n; j++)
        expr += vars[i][j];
      model.addConstr(expr == 2, "deg2_"+to_string(i));
    }

    // Forbid edge from node back to itself

    for (i = 0; i < n; i++)
      vars[i][i].set(GRB_DoubleAttr_UB, 0);

    // Set callback function

    subtourelim cb = subtourelim(vars, n);
    model.setCallback(&cb);

    // Optimize model

    model.optimize();

    // Extract solution

    if (model.get(GRB_IntAttr_SolCount) > 0) {
      double **sol = new double*[n];
      for (i = 0; i < n; i++)
        sol[i] = model.get(GRB_DoubleAttr_X, vars[i], n);

      int* tour = new int[n];
      int len;

      findsubtour(n, sol, &len, tour);
      std::cout<<"len = "<<len<<" n = "<<n<<std::endl;
      assert(len == n);

      //cout << "Tour: ";
      for (i = 0; i < len; i++){
        //cout << tour[i] << " ";
        if(tour[i]<table->size())
            uni_tour->push_back(tour[i]);
      }  
      //cout << endl;

      for (i = 0; i < n; i++)
        delete[] sol[i];
      delete[] sol;
      delete[] tour;

      //mini_delay = comp_delay(self_loops,uni_tour,table);
    }

  } catch (GRBException e) {
    cout << "Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Error during optimization" << endl;
  }

  for (i = 0; i < n; i++)
    delete[] vars[i];
  delete[] vars;
  delete env;
  //VLOG(2)<<"mini_delay = "<<mini_delay;
  return uni_tour;
}

long comp_delay(shared_ptr<vector<unsigned int>> &self_loop,\
        shared_ptr<vector<int>> &uni_tour,\
        shared_ptr<vector<vector<int>>> &table){
    VLOG(2)<<"start computing delay...";
    for(size_t u=0;u<uni_tour->size();++u){
        VLOG(2)<<"u: "<<uni_tour->at(u)<<" self_loop:"<<self_loop->at(u);
    }
    long min_mark = FLT_MAX;
    unsigned int i;
    unsigned int start_label;
    int inx_next;
    int inx_prev;
    for(i=0;i<uni_tour->size();++i){
        inx_next = uni_tour->at(i);
        inx_prev = uni_tour->at((i+uni_tour->size()-1)%uni_tour->size());
        if((self_loop->at(inx_next) - table->at(inx_prev).at(inx_next))<min_mark){
            start_label = i;
            min_mark = (self_loop->at(inx_next) - table->at(inx_prev).at(inx_next));
        }
    }
    
    long min_delay = self_loop->at(uni_tour->at(start_label));
    VLOG(2)<<"start: "<<uni_tour->at(start_label);
    
    for(size_t j=1;j<uni_tour->size();++j){
        inx_next = uni_tour->at((start_label+j)%uni_tour->size());
        inx_prev = uni_tour->at((start_label+j-1)%uni_tour->size());
        min_delay += table->at(inx_prev).at(inx_next);
        VLOG(2)<<"Touring:"<<inx_next<<"    dist:"<<min_delay;
    }
    return min_delay;
}








