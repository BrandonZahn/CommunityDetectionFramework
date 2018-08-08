#include "global.h"
#include "NSGA2.h"

#include "random.h"

#include <vector>
#include <numeric>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

extern nsga2::random_gen rgen; // global common random generator

extern int nvertices;
extern std::vector< std::vector<double> >*  mat;
extern std::vector< std::vector<double> >*  sim;
extern int numEdges;
extern double totalWeight;
extern int popsize;
extern double pmut_comm;
void test_problem1 (std::vector<double>* obj, std::vector<nsga2::community*>* comm) {
    (*obj)[0] = 0;
    (*obj)[1] = 0;

    for(int i =0; i< (*comm).size(); ++i)
    {
        //cout << i << endl;
        (*obj)[0] -=(*comm)[i]->modularity;
        (*obj)[1] -=(*comm)[i]->w_modularity;

    }
    //(*obj)[0] = round((*obj)[0]*10000.0)/10000.0;
    //(*obj)[1] = round((*obj)[1]*10000.0)/10000.0;

    return;
}

void test_population(nsga2::population& pop) {
    cout << "Evaluating population!" << endl;
    pop.evaluate();
}

/*int main(int argc, char *argv[]) {
    
    if (argc<3) {
        cout << "Usage " << argv[0] << " random_seed input_file" << endl;
        return 1;
    }

    
    int seed = atoi(argv[1]);
    rgen.set_seed(seed);

    
    ifstream cin(argv[2]);*/

  int main()
  {
    
  ifstream cin("input.dat");
	string inputData;
	cin >> inputData;
	int seed = stoi(inputData);
	rgen.set_seed(seed);
  cout << seed << endl;
	
    cin >> nvertices; cout << "nvertices: " << nvertices << endl;

    mat = new std::vector< std::vector<double> >(nvertices, std::vector<double>(nvertices));
    sim = new std::vector< std::vector<double> >(nvertices, std::vector<double>(nvertices));
    int i;
    int j;

    std::vector<nsga2::vertex*> vertices;

    numEdges = 0;
    for(i=0; i<nvertices; ++i)
    {
        for(j=0; j<nvertices; ++j)
        {
            cin >> (*mat)[i][j];
            numEdges += (*mat)[i][j];
        }

    }
    totalWeight = 0;
    for(i=0; i<nvertices; ++i)
    {
        for(j=0; j<nvertices; ++j)
        {
            cin >> (*sim)[i][j];
            totalWeight += (*sim)[i][j];
        }

    }
    numEdges = numEdges /2;
    totalWeight = totalWeight/2.0;
    for(i=0; i<nvertices; ++i)
    {
        nsga2::vertex* temp = new nsga2::vertex(i);
        vertices.push_back(temp);
    }
    for(i=0; i<nvertices; ++i)
    {
        vertices[i]->update(mat, mat);
    }


    nsga2::individual* ind = new nsga2::individual;


   
    // ind->initialize(vertices);
    // cout << "Individual has been initialized with: " << endl;
    // for(i =0; i<nvertices; ++i)
    // {
    //     cout << ind->vertexSet[i] << " " ;
    // }
    // cout << endl;
    // cout  << "i.e. it has the following communities: " <<endl;
    // for(i=0; i<ind->comm.size(); ++i)
    // {
    //     cout << "Community : " << ind->comm[i]->number << ", with modularity " << ind->comm[i]->modularity << " Containing: " <<endl;
    //     for(int j=0; j<ind->comm[i]->vertices.size(); ++j)
    //     {
    //         cout << "Vertex : " << ind->comm[i]->vertices[j]->number <<  "Degree: " << ind->comm[i]->vertices[j]->degree << endl;
    //     }
    // }

    // ind->comm[13]->modularity_update(vertices[0]);

    //     cout << "Individual has been initialized with: " << endl;
    // for(i =0; i<nvertices; ++i)
    // {
    //     cout << ind->vertexSet[i] << " " ;
    // }
    // cout << endl;
    // cout  << "i.e. it has the following communities: " <<endl;
    // for(i=0; i<ind->comm.size(); ++i)
    // {
    //     cout << "Community : " << ind->comm[i]->number << ", with modularity " << ind->comm[i]->modularity << " Containing: " <<endl;
    //     for(int j=0; j<ind->comm[i]->vertices.size(); ++j)
    //     {
    //         cout << "Vertex : " << ind->comm[i]->vertices[j]->number <<  "Degree: " << ind->comm[i]->vertices[j]->degree << endl;
    //     }
    // }

    // cout<< numEdges << endl;

    
    cout << "Enter the problem relevant and algorithm relevant parameters ... " << endl;
    cout << "Enter the population size (a multiple of 4) : " << endl;
    cin >> popsize;
    if (popsize<4 || (popsize%4)!= 0) {
        cout << "population size read is : " << popsize << endl;
        cout << "Wrong population size entered, hence exiting" << endl;
        return 1;
    }

    int ngen;
    cout << "Enter the number of generations : " << endl;
    cin >> ngen;
    if (ngen<1) {
        cout << "number of generations read is : " << ngen << endl;
        cout << "Wrong nuber of generations entered, hence exiting" << endl;
        return 1;
    }

    int nobj;
    cout << "Enter the number of objectives : " << endl;
    cin >> nobj;
    if (nobj<1) {
        cout << "number of objectives entered is : " << nobj << endl;
        cout << "Wrong number of objectives entered, hence exiting" << endl;
        return 1;
    }

    int ncon;
    cout << "Enter the number of constraints : " << endl;
    cin >> ncon;
    if (ncon<0) {
        cout << "number of constraints entered is : " << ncon << endl;
        cout << "Wrong number of constraints enetered, hence exiting" << endl;
        return 1;
    }

 

     cout << "Input data successfully entered, now performing initialization" << endl;
  

    cin >> pmut_comm;

    double pcross_comm;

    cin >> pcross_comm;

     nsga2::NSGA2 nsga2;
    
    nsga2.set_nvertices(nvertices);

    double eta_c;
    double eta_m;
    cin >> eta_c;
    cin >> eta_m;

    nsga2.set_nobj(nobj);
    nsga2.set_ncon(ncon);
    nsga2.set_popsize(popsize);
    nsga2.set_ngen(ngen);
    nsga2.set_pcross_comm(pcross_comm);
    nsga2.set_pmut_comm(pmut_comm);
    
    nsga2.set_eta_c(eta_c);
    nsga2.set_eta_m(eta_m);
   

    nsga2.set_function(&test_problem1);
    nsga2.set_popfunction(&test_population);

    try
    {
        nsga2.initialize(vertices);
    }
    catch(exception &e)
    {
        cerr<<e.what() <<endl;
    }
     nsga2.evolve();
    



    cout << "Didnt crash " << endl;
    return 0;
}
