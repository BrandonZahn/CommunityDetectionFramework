#include "NSGA2.h"
//#include "rand.h"

#include <iostream>
#include <numeric>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>

using namespace nsga2;
using namespace std;



struct modCompare{
    bool operator () (community* a, community* b) { return (a->modularity > b->modularity); }
}modCompare;

struct w_modCompare
{
    bool operator () (community* a, community* b) { return (a->w_modularity > b->w_modularity); }
}w_modCompare;


NSGA2::NSGA2() :
    nreal(-1),
    nbin(-1),
    nobj(-1),
    ncon(-1),
    popsize(-1),
    ngen(-1),
    nreport(1),
    pcross_real(-1),
    pcross_bin(-1),
    pmut_real(-1),
    pmut_bin(-1),
    eta_c(-1),
    eta_m(-1),
    epsilon_c(EPS),
    nbits(0),
    limits_realvar(0),
    limits_binvar(0),
    function(0),
    popFunction(0),
    reportFunction(0),
    // choice(0),
    // obj1(0),
    // obj2(0),
    // obj3(0),
    // angle1(0),
    // angle2(0),
    backupFilename("nsga2_backup_pop.data"),
    nbinmut(0),
    nrealmut(0),
    ncommmut(0),
    nbincross(0),
    nrealcross(0),
    ncommcross(0),
    bitlength(0),
    parent_pop(0),
    child_pop(0),
    mixed_pop(0),
    crowd_obj(true) {
}

NSGA2::~NSGA2() {
    if (parent_pop) {
        delete parent_pop;
        parent_pop = 0;
    }
    if (child_pop) {
        delete child_pop;
        child_pop = 0;
    }
    if (mixed_pop) {
        delete mixed_pop;
        mixed_pop = 0;
    }
}


void NSGA2::initialize(vector<vertex*> verts) throw (nsga2exception) {

    cout << "Initializing NSGA-II v0.2.1\n"
         << "Checking configuration" << endl;

    if (nobj < 1)
        throw nsga2exception("Invalid number of objective functions");
    if (ncon < 0)
        throw nsga2exception("Invalid number of constraints");
    if (popsize<4 || (popsize%4)!= 0)
        throw nsga2exception("Invalid size of population");
    if (pcross_comm<0.0 || pcross_comm>1.0)
        throw nsga2exception("Invalid probability of comm crossover");
    if (pmut_comm<0.0 || pmut_comm>1.0)
        throw nsga2exception("Invalid probability of comm mutation");
    
    if (eta_c<=0)
        throw nsga2exception("Invalid distribution index for crossover");
    if (eta_m<=0)
        throw nsga2exception("Invalid distribution index for mutation");
    if (ngen < 1)
        throw nsga2exception("Invalid number of generations");
      if (function == 0)
        throw nsga2exception("Evaluation function not defined");

    init_streams();
    cout << "Streams init successfully" << endl;
    report_parameters(fpt5);
    cout << "Report parameters successful" <<endl;
    nbinmut    = 0;
    nrealmut   = 0;
    nbincross  = 0;
    nrealcross = 0;
   

    parent_pop = new population(popsize,
                                nreal,
                                nbin,
                                ncon,
                                nbits,
                                limits_realvar,
                                limits_binvar,
                                nobj,
                                pmut_real,
                                pmut_bin,
                                eta_m,
				epsilon_c,
                                function);
    cout << "Parent created successfully" << endl;
    child_pop  = new population(popsize,
                                nreal,
                                nbin,
                                ncon,
                                nbits,
                                limits_realvar,
                                limits_binvar,
                                nobj,
                                pmut_real,
                                pmut_bin,
                                eta_m,
				epsilon_c,
                                function);
    mixed_pop  = new population(popsize*2,
                                nreal,
                                nbin,
                                ncon,
                                nbits,
                                limits_realvar,
                                limits_binvar,
                                nobj,
                                pmut_real,
                                pmut_bin,
                                eta_m,
				epsilon_c,
                                function);

    if (popFunction) {
	parent_pop->set_popfunction(popFunction);
	child_pop->set_popfunction(popFunction);
	mixed_pop->set_popfunction(popFunction);
    }

    parent_pop->crowd_obj = crowd_obj;
    child_pop->crowd_obj = crowd_obj;
    mixed_pop->crowd_obj = crowd_obj;

    //randomize();

    bool fromBackup = false;
    if (!fromBackup) {
        parent_pop->initialize(verts);
        cout << "Initialization done, now performing first generation" << endl;

        //parent_pop->decode();
        parent_pop->evaluate();
        parent_pop->fast_nds();
        parent_pop->crowding_distance_all();
        cout << "Crowding distance successful" << endl;

        t = 1;
    } else {
        cout << "Initialization made from backup file" << endl;
    }

    custom_report(*parent_pop);

    report_pop(*parent_pop,fpt1);
    fpt4 << "# gen = " << t << '\n';
    report_pop(*parent_pop,fpt4);

    fpt1.flush();
    fpt4.flush();
    fpt5.flush();

}

void NSGA2::init_streams() {
    fpt1.open("nsga2_initial_pop.out" , ios::out | ios::trunc);
    fpt2.open("nsga2_final_pop.out"   , ios::out | ios::trunc);
    fpt3.open("nsga2_best_pop.out"    , ios::out | ios::trunc);
    fpt4.open("nsga2_all_pop.out"     , ios::out | ios::trunc);
    fpt5.open("nsga2_params.out"      , ios::out | ios::trunc);

    fpt1.setf(ios::scientific);
    fpt2.setf(ios::scientific);
    fpt3.setf(ios::scientific);
    fpt4.setf(ios::scientific); fpt4.precision(16);
    fpt5.setf(ios::scientific);

    fpt1 << "# This file contains the data of initial population\n";
    fpt2 << "# This file contains the data of final population\n";
    fpt3 << "# This file contains the data of final feasible population (if found)\n";
    fpt4 << "# This file contains the data of all generations\n";
    fpt5 << "# This file contains information about inputs as read by the program\n";

    fpt1 << "# of objectives = "    << nobj
         << ", # of constraints = " << ncon
         << ", # of real_var = "    << nreal
         << ", # of bits of bin_var = " << bitlength
         << ", constr_violation, rank, crowding_distance\n";
    fpt2 << "# of objectives = "    << nobj
         << ", # of constraints = " << ncon
         << ", # of real_var = "    << nreal
         << ", # of bits of bin_var = " << bitlength
         << ", constr_violation, rank, crowding_distance\n";
    fpt3 << "# of objectives = "    << nobj
         << ", # of constraints = " << ncon
         << ", # of real_var = "    << nreal
         << ", # of bits of bin_var = " << bitlength
         << ", constr_violation, rank, crowding_distance\n";
    fpt4 << "# of objectives = "    << nobj
         << ", # of constraints = " << ncon
         << ", # of real_var = "    << nreal
         << ", # of bits of bin_var = " << bitlength
         << ", constr_violation, rank, crowding_distance\n";
}

void NSGA2::report_parameters(std::ostream& os) const {
    os << "Population size = " << popsize
       << "\nNumber of generations = " << ngen
       << "\nNumber of objective functions = " << nobj
       << "\nNumber of constraints = " << ncon
       << "\nNumber of real variables = " << nreal;

    if (nreal != 0) {
        for (int i = 0; i<nreal; ++i) {
            os << "\nLower limit of real variable " << (i+1)
               << " = " << limits_realvar[i].first;
            os << "\nUpper limit of real variable " << (i+1)
               << " = " << limits_realvar[i].second;
        }
        os << "\nProbability of crossover of real variable = " << pcross_real;
        os << "\nProbability of mutation of real variable = " << pmut_real;
        os << "\nDistribution index for crossover = " << eta_c;
        os << "\nDistribution index for mutation = " << eta_m;
    }

    os << "\nNumber of binary variables = " << nbin;
    if (nbin != 0) {
        for (int i = 0; i<nbin; ++i) {
            os << "\nNumber of bits for binary variable " << (i+1)
               << " = " << nbits[i];
            os << "\nLower limit of real variable " << (i+1)
               << " = " << limits_binvar[i].first;
            os << "\nUpper limit of real variable " << (i+1)
               << " = " << limits_binvar[i].second;
        }
        os << "Probability of crossover of binary variable = " << pcross_bin;
        os << "Probability of mutation of binary variable = " << pmut_bin;
    }
    os << "\nSeed for random number generator = " << rgen.get_seed() << endl;
}

void NSGA2::report_pop(const population& pop, std::ostream& os) const {
    pop.report(os);
}

void NSGA2::report_pop_obj(const population& pop, std::ostream& os) const {
    pop.report_obj(os);
}

void NSGA2::save_backup() const {
 //    cout << "Saving backup: ";
 //    if (backupFilename == "") {
	// cout << "No backup file set" << endl;
	// return;
 //    }


 //    char tempfilename[L_tmpnam];
 //    char* res = tmpnam(tempfilename);
 //    if (!res) {
	// perror("Could not create temporary file!");
	// return;
 //    }
 //    cout << tempfilename << endl;

 //    ofstream ofs(tempfilename, ios::binary);

 //    ofs.write(reinterpret_cast<const char*>(&t), sizeof(int));
 //    ofs.write(reinterpret_cast<const char*>(&nbinmut), sizeof(int));
 //    ofs.write(reinterpret_cast<const char*>(&nrealmut), sizeof(int));
 //    ofs.write(reinterpret_cast<const char*>(&nbincross), sizeof(int));
 //    ofs.write(reinterpret_cast<const char*>(&nrealcross), sizeof(int));

 //    parent_pop->dump(ofs);

 //    ofs.flush();
 //    ofs.close();

 //    int result = rename(tempfilename, backupFilename.c_str());
 //    if (result)
 //        perror("Could not save backup!");

}

bool NSGA2::load_backup() {
    cout << "Loading backup: ";
    
    if (backupFilename == "") {
	cout << "No backup file set" << endl;
	return false;
    }

    ifstream ifs(backupFilename.c_str(), ios::in | ios::binary);

    if (!ifs.good()) {
        cout << "Could not load backup file" << endl;
        return false;
    }

    ifs.read(reinterpret_cast<char*>(&t),sizeof(int));
    ifs.read(reinterpret_cast<char*>(&nbinmut),sizeof(int));
    ifs.read(reinterpret_cast<char*>(&nrealmut),sizeof(int));
    ifs.read(reinterpret_cast<char*>(&nbincross),sizeof(int));
    ifs.read(reinterpret_cast<char*>(&nrealcross),sizeof(int));

    parent_pop->load(ifs);

    ifs.close();

    return true;
}

// void NSGA2::report_feasible(const population& pop, std::ostream& os) const {
//     pop.report_feasible(os);
// }

void NSGA2::selection(population& oldpop, population& newpop)
    throw (nsga2::nsga2exception) {
    const int N = oldpop.size();
    if (newpop.size() != N)
        throw nsga2::nsga2exception("Selection error: new and old pops don't have the same size");

    std::vector<int> a1(N), a2(N);
    for (int i = 0; i < N; ++i) {
        a1[i] = a2[i] = i;
    }

    int rand;
    for (int i = 0; i < N; ++i) { // this could be done with a shuffle
        rand = rgen.integer(i,N-1);
        std::swap(a1[rand],a1[i]);
        rand = rgen.integer(i,N-1);
        std::swap(a2[rand],a2[i]);
    }

    for (int i = 0; i < N; i+=4) {

        individual& p11 = tournament(oldpop.ind[a1[i  ]], oldpop.ind[a1[i+1]]);
        individual& p12 = tournament(oldpop.ind[a1[i+2]], oldpop.ind[a1[i+3]]);
       
        crossover(p11,p12,newpop.ind[i  ],newpop.ind[i+1]);
       
        individual& p21 = tournament(oldpop.ind[a2[i  ]], oldpop.ind[a2[i+1]]);
        individual& p22 = tournament(oldpop.ind[a2[i+2]], oldpop.ind[a2[i+3]]);
        crossover(p21,p22,newpop.ind[i+2],newpop.ind[i+3]);
    }

}

individual& NSGA2::tournament(individual& ind1, individual& ind2) const {
    int flag = ind1.check_dominance(ind2);
    if (flag == 1) // ind1 dominates ind2
        return ind1;
    else if (flag == -1) // ind2 dominates ind1
        return ind2;
    else if (ind1.crowd_dist > ind2.crowd_dist)
        return ind1;
    else if (ind2.crowd_dist > ind1.crowd_dist)
        return ind2;
    else if (rgen.realu() <= 0.5)
        return ind1;
    else
        return ind2;
}

void NSGA2::crossover(const individual& parent1, const individual& parent2,
                      individual& child1, individual& child2) {


    commcross(parent1,parent2,child1,child2);

    child1.evaluated = false;
    child2.evaluated = false;

}

std::vector<community*> NSGA2::big_list_make(std::vector<community*> first_list, std::vector<community*> second_list)
{
    int i=0;
    int j =0;
    std::vector<community*> out;
    while(i<first_list.size() || j<second_list.size())
    {
        
        if(i<first_list.size())
        {
            out.push_back(first_list[i]);
            ++i;
        }
        if( j<second_list.size())
        {
            out.push_back(second_list[j]);
            ++j;
        }
    }
    return out;
}
// The function to crossover communities
// Each community with all of its member nodes is considered the 'attribute'
void NSGA2::commcross(const individual& parent1, const individual& parent2,
                      individual& child1, individual& child2)
{

    int numSet1 = 0;
    int numSet2 = 0;


    //We make copies of the community lists so we can make modifications without
    //Affecting the parents. One for normal one for weighted
    std::vector<community*> parent1_commCopy; 
    std::vector<community*> parent1_w_commCopy;
    for(int i =0; i<parent1.comm.size(); ++i)
    {
        parent1_commCopy.push_back(parent1.comm[i]->clone());
        parent1_w_commCopy.push_back(parent1.comm[i]->clone());
    }
    std::sort(parent1_commCopy.begin(), parent1_commCopy.end(), modCompare);
    std::sort(parent1_w_commCopy.begin(), parent1_w_commCopy.end(), w_modCompare);
    std::vector<community*> parent2_commCopy;
    std::vector<community*> parent2_w_commCopy;
    for(int i =0; i<parent2.comm.size(); ++i)
    {
        parent2_commCopy.push_back(parent2.comm[i]->clone());
        parent2_w_commCopy.push_back(parent2.comm[i]->clone());
    }
    std::sort(parent2_commCopy.begin(), parent2_commCopy.end(), modCompare);
    std::sort(parent2_w_commCopy.begin(), parent2_w_commCopy.end(), w_modCompare);

    std::vector<community*> biglist;
    std::vector<community*> biglist2;

    // Merge the two lists of communities in an alternating pattern
    // Randomly choose who goes first
    if(rgen.realu() < 0.5)
    {
        biglist = big_list_make(parent1_commCopy, parent2_w_commCopy);
    }
    else
    {
        biglist = big_list_make(parent2_w_commCopy, parent1_commCopy);
    }

    if(rgen.realu() < 0.5)
    {
        biglist2 = big_list_make(parent2_commCopy, parent1_w_commCopy);
    }
    else
    {
        biglist2 = big_list_make(parent1_w_commCopy, parent2_commCopy);
    }


    int commNumber = 1;
    int i=0;
    int community_added = 0;

    while(numSet1 < nvertices){ // While there are vertices not belonging to any community in the first child
        community_added = 0;
               for(int q = 0; q<biglist[i]->vertices.size(); ++q) // For each vertex in the ith community
                {
                    vertex* c_vertex = biglist[i]->vertices[q];

                    if(child1.vertexSet[c_vertex->number] == 0){ // 0 is flag for no community yet
                      
                        child1.vertexSet[c_vertex->number] = commNumber; // Make it in this comunity
                        if(community_added == 0) // If we are creating a new community in doing so
                        { // Add it to the list of communities
                            biglist[i]->number = commNumber;
                            child1.comm.push_back(biglist[i]);
                            community_added = 1;
                        }                        
                        ++numSet1;
                    }
                    else
                    { //Otherwise this vertex is already in another community
                      // Take it out and update the modularity of the community
                         biglist[i]->modularity_update(c_vertex); 
                         q--;
                    }
                   
                }

                ++i;
                if(community_added != 0)
                    ++commNumber;   
        }
    i=0;
    community_added = 0;
    //Same as above for second child
        while(numSet2 < nvertices){
                community_added = 0;                
                for(int q = 0; q<biglist2[i]->vertices.size(); ++q)
                {
                   vertex* c_vertex = biglist2[i]->vertices[q];
                   
                    if(child2.vertexSet[c_vertex->number] == 0){ // 0 is flag for no community yet
                       
                        child2.vertexSet[c_vertex->number] = commNumber;
                        if(community_added == 0)
                        {
                            biglist2[i]->number = commNumber;
                            child2.comm.push_back(biglist2[i]);
                            community_added = 1;
                        }                        
                        ++numSet2;
                    }
                    else
                    {
                         biglist2[i]->modularity_update(c_vertex);
                         q--;

                    }
                   
                }          
                ++i;
                if(community_added!=0)
                    ++commNumber;     
    }   
}




void NSGA2::realcross(const individual& parent1, const individual& parent2,
                      individual& child1, individual& child2) 
{

    int i;
    double rand;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;
    if (rgen.realu() <= pcross_real) {
        nrealcross++;
	for (i=0; i<nreal; i++) {
	    //if (rgen.realu()<=0.5 ) { this is evil, in my opinion

	    if (fabs(parent1.xreal[i]-parent2.xreal[i]) > EPS) {
		
		if (parent1.xreal[i] < parent2.xreal[i]) {
		    y1 = parent1.xreal[i];
		    y2 = parent2.xreal[i];
		} else {
		    y1 = parent2.xreal[i];
		    y2 = parent1.xreal[i];
		}
		
		yl = limits_realvar[i].first;
		yu = limits_realvar[i].second;
		
		rand = rgen.realu();
		beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
		alpha = 2.0 - pow(beta,-(eta_c+1.0));
		if (rand <= (1.0/alpha)) { // This is a contracting crossover
		    betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
		} else {                   // This is an expanding crossover
		    betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
		}
		c1 = 0.5*((y1+y2)-betaq*(y2-y1));
		
		beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
		alpha = 2.0 - pow(beta,-(eta_c+1.0));
		if (rand <= (1.0/alpha)) { // This is a contracting crossover
		    betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
		} else {                   // This is an expanding crossover
		    betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
		}
		c2 = 0.5*((y1+y2)+betaq*(y2-y1));
		
		c1 = min(max(c1,yl),yu);
		c2 = min(max(c2,yl),yu);
		
		if (rgen.realu()<=0.5) {
		    child1.xreal[i] = c2;
		    child2.xreal[i] = c1;
		} else {
		    child1.xreal[i] = c1;
		    child2.xreal[i] = c2;
		}
	    } else {
		child1.xreal[i] = parent1.xreal[i];
		child2.xreal[i] = parent2.xreal[i];
	    }
	    // } else {
	    //     // This is evil, in my opinion
	    //     child1.xreal[i] = parent1.xreal[i];
	    //     child2.xreal[i] = parent2.xreal[i];
	    // }
	}
    } else {
        for (i=0; i<nreal; i++) {
            child1.xreal[i] = parent1.xreal[i];
            child2.xreal[i] = parent2.xreal[i];
        }
    }
}

void NSGA2::bincross(const individual& parent1, const individual& parent2,
                      individual& child1, individual& child2)
{


    int i, j;
    double rand;
    int temp, site1, site2;
    for (i=0; i<nbin; i++) {
        rand = rgen.realu();
        if (rand <= pcross_bin) {
            nbincross++;
            site1 = rgen.integer(0,nbits[i]-1);
            site2 = rgen.integer(0,nbits[i]-1);
            if (site1 > site2) {
                temp = site1;
                site1 = site2;
                site2 = temp;
            }
            for (j=0; j<site1; j++) {
                child1.gene[i][j] = parent1.gene[i][j];
                child2.gene[i][j] = parent2.gene[i][j];
            }
            for (j=site1; j<site2; j++) {
                child1.gene[i][j] = parent2.gene[i][j];
                child2.gene[i][j] = parent1.gene[i][j];
            }
            for (j=site2; j<nbits[i]; j++) {
                child1.gene[i][j] = parent1.gene[i][j];
                child2.gene[i][j] = parent2.gene[i][j];
            }
        } else {
            for (j=0; j<nbits[i]; j++) {
                child1.gene[i][j] = parent1.gene[i][j];
                child2.gene[i][j] = parent2.gene[i][j];
            }
        }
    }
}

struct sort_n {
    const population& pop;
    sort_n(const population& population) : pop(population) {};
    bool operator() (int i, int j) {
        const individual& ind1 = pop.ind[i];
        const individual& ind2 = pop.ind[j];
        if (ind1.rank < ind2.rank)
            return true;
        else if (ind1.rank == ind2.rank &&
                 ind1.crowd_dist > ind2.crowd_dist)
            return true;
        return false;
    };
};

void printme(const individual& ind) {
    cout << ind << endl;
}

void NSGA2::custom_report(population& pop) {
    if (reportFunction)
	(*reportFunction)(pop);
}

void NSGA2::advance() {

    cout << "Advancing to generation " << t+1 << endl;  

    std::pair<int,int> res;

    // create next population Qt
    selection(*parent_pop,*child_pop);
    cout << "Selection Successful" << endl;
    res = child_pop->mutate();
    cout << "Mutation Successful" << endl;
    child_pop->generation = t+1;
    //child_pop->decode();
    child_pop->evaluate();
    // for(int i=0; i<child_pop->ind.size(); ++i){
    //     printme(child_pop->ind[i]);
    // }

    // mutation book-keeping



    // create population Rt = Pt U Qt
    mixed_pop->merge(*parent_pop,*child_pop);
    mixed_pop->generation = t+1;

    // fpt4 << "#Mixed\n";
    // report_pop(*mixed_pop, fpt4);

    mixed_pop->fast_nds();
    //mixed_pop->crowding_distance_all();

    // fpt4 << "#Mixed nfs\n";
    // report_pop(*mixed_pop, fpt4);

    // Pt+1 = empty
    parent_pop->ind.clear();
    child_pop->reset();

    int i = 0;
    // until |Pt+1| + |Fi| <= N, i.e. until parent population is filled
    while (parent_pop->size() + mixed_pop->front[i].size() < popsize) {
        std::vector<int>& Fi = mixed_pop->front[i];
        mixed_pop->crowding_distance(i);           // calculate crowding in Fi
        for (int j = 0; j < Fi.size(); ++j)        // Pt+1 = Pt+1 U Fi
            parent_pop->ind.push_back(mixed_pop->ind[Fi[j]]);
        i += 1;
    }

    mixed_pop->crowding_distance(i);           // calculate crowding in Fi
    std::sort(mixed_pop->front[i].begin(),
              mixed_pop->front[i].end(),
              sort_n(*mixed_pop) );// sort remaining front using <n

    const int extra = popsize - parent_pop->size();
    for (int j = 0; j < extra; ++j) // Pt+1 = Pt+1 U Fi[1:N-|Pt+1|]
        parent_pop->ind.push_back(mixed_pop->ind[mixed_pop->front[i][j]]);

    t += 1;

    // if (popFunction) {
    //   (*popFunction)(*parent_pop);
    // }

    parent_pop->generation = t;
    //custom_report(*parent_pop);

 //    if (t%nreport == 0) {
	// fpt4 << "# gen = " << t << '\n';
	// report_pop(*parent_pop,fpt4);
	// fpt4.flush();
 //    }

    // save a backup
    save_backup();
}

void NSGA2::evolve() {
  while (t < ngen)
    advance();
  report_pop(*parent_pop,fpt2);
  report_pop_obj(*parent_pop, fpt3);
}

