#include "global.h"
//#include "rand.h"

#include <cmath>
#include <iostream>
#include <algorithm>

using namespace nsga2;
using namespace std;

extern random_gen rgen; // global common random generator
int nvertices;
std::vector< std::vector<double> >*  mat;
std::vector< std::vector<double> >*  sim;
int numEdges;
double totalWeight;
int popsize;
int nobj;
double pmut_comm;
vertex::vertex()
{

}

vertex::vertex(int num, Matrix* adj, Matrix* sim)
{
    number = num;


}

vertex::vertex(int num)
{
    number = num;
    degree =0;
    w_degree = 0;

    
}
    

void vertex::update(std::vector<std::vector<double> >* adj, std::vector<std::vector<double> >* sim)
{
    int i, j;
    double temp_degree =0;
    double temp_w_degree=0;
    vertex* temp; 
    //Below simply loops through and adds up the degree and weighted degree
    for(i=0; i<adj->size(); ++i)
    {
        temp_degree += (*adj)[number][i];
        temp_w_degree+= (*sim)[number][i];
    } 
    degree = temp_degree;
    w_degree = temp_w_degree;
}

community::community(int _number)
{
    number = _number;
    size =0;
}


void community::modularity_calc()
{
    double internal =0;
    double total = 0;
    double w_total =0;
    double w_internal=0;
    double i_current;
    double w_i_current;
    internalDegrees.clear();
    w_internalDegrees.clear();
    for(int i =0; i<vertices.size(); ++i)
    {
        i_current = 0;
        w_i_current = 0;

        for(int j=0; j<vertices.size(); ++j)
        {
            i_current += (*mat)[vertices[i]->number][vertices[j]->number];
            w_i_current += (*sim)[vertices[i]->number][vertices[j]->number];
        }
        internalDegrees.push_back(i_current);
        w_internalDegrees.push_back(w_i_current);

        internal += i_current;
        w_internal += w_i_current;
        total+=vertices[i]->degree;
        w_total += vertices[i]->w_degree;

    }
    totalDegree = total;

    internalDegree = internal;
    w_internalDegree = w_internal;
    w_totalDegree = w_total;
    modularity = (internal/2/numEdges - (total/(2*numEdges))*(total/(2*numEdges))); // Formula for Laylas thesis
    w_modularity = (w_internal/2/totalWeight - (w_total/(2*totalWeight))*(w_total/(2*totalWeight))); // Formula for Laylas thesis
}

void community::modularity_update(vertex* v)
{

    for(int i=0; i<internalDegrees.size(); ++i)
    {

        if(vertices[i] == v)
        { // If this is the vertex we are removing
            modularity += (2*totalDegree*v->degree)/(4*numEdges*numEdges) - ((v->degree)/(2*numEdges))*((v->degree)/(2*numEdges)) - internalDegrees[i]/numEdges; // Adjust according to laylas formulas
            w_modularity += (2*w_totalDegree*v->w_degree)/(4*totalWeight*totalWeight) - ((v->w_degree)/(2*totalWeight))*((v->w_degree)/(2*totalWeight)) - w_internalDegrees[i]/totalWeight;
            totalDegree -= v->degree;
            w_totalDegree -= v->w_degree;
            internalDegrees.erase(internalDegrees.begin() + i);
            w_internalDegrees.erase(w_internalDegrees.begin() + i);
            vertices.erase(vertices.begin() + i); 
        }
        else
        { // Otherwise just adjust the degree of this vertex
            w_internalDegrees[i] -=  (*sim)[v->number][vertices[i]->number];
            internalDegrees[i] -= (*mat)[v->number][vertices[i]->number];
        }
    }

}

community* community::clone()
{
    community* new_comm = new community(number);
    new_comm->vertices = vertices;
    
    new_comm->modularity = modularity;
    new_comm->w_modularity = w_modularity;
    new_comm->totalDegree = totalDegree;
    new_comm-> w_totalDegree = w_totalDegree;
    new_comm-> internalDegree = internalDegree;
    new_comm-> w_internalDegree = w_internalDegree;
    new_comm-> size = size;
    new_comm->number = 0;
    new_comm-> internalDegrees = internalDegrees;
    new_comm-> w_internalDegrees = w_internalDegrees;

    return new_comm;
}



individual::individual() throw () :
    rank(0),
    constr_violation(0),
    xreal(0),
    gene(0),
    xbin(0),
    obj(0),
    constr(0),
    crowd_dist(0),
    config(0)
   { for(int i=0; i<nvertices; ++i)
    {
        vertexSet.push_back(0);
    }
}

void individual::reset()
{
    vertexSet.clear();
    comm.clear();



    for(int i=0; i<nvertices; ++i)
    {
        vertexSet.push_back(0);
    }
}


individual::individual(const individual_config& c
                       /*const unsigned int nreal,
                       const unsigned int nbin,
                       const unsigned int ncon,
                       const std::vector<int>& nbits,
                       const unsigned int nobj*/) throw (nsga2::nsga2exception) :
    rank(0),
    constr_violation(0),
    xreal(0),
    gene(0),
    xbin(0),
    obj(0),
    constr(0),
    crowd_dist(0),
    evaluated(false),
    config(&c) {

    for(int i=0; i<nvertices; ++i)
    {
        vertexSet.push_back(0);
    }
    obj.resize(config->nobj,0);
    constr.resize(config->ncon,0);
}

individual::~individual() {
}


void individual::initialize(std::vector<vertex*> verts, int num_comm) throw (nsga2::nsga2exception) {
  // Apply an intensive local search to increase the quality of initialization
  int t;
  double deltaQ;

  while (t < vertexSet.size())
  {
    t = 0;
    vector<int> tempVertexSet;
    for (int i = 0; i < verts.size(); ++i) 
    {    
      tempVertexSet.push_back(rgen.integer(1, num_comm));			
    }
    
    for(int i=1; i<verts.size()+1; ++i)
    {
      community* current_comm = NULL;
      for(int j=0; j<verts.size(); ++j)
      {
        
        if(tempVertexSet[j] == i)
        {
          if(current_comm == NULL)
          {
            current_comm = new community(i);
          }
          current_comm->grow();
          vertex* temp = verts[j];
          current_comm->vertices.push_back(temp);
        }
      }
      if(current_comm !=NULL)
      {
        comm.push_back(current_comm);
      }
    }

    for(int i =0; i<comm.size(); ++i)
    {
      comm[i]->modularity_calc();
    }
    
    for (int i = 0; i < verts.size(); i++)
    {
      int RandVer = 0; // j -> algorithm 3
      vector<int> neighbours; // M
      for (int j = 0; j < verts.size(); j++)
      {
        if ((*mat)[i][j] == 1 && tempVertexSet[j] != tempVertexSet[i])
        {
          neighbours.push_back(j);
        }
      }
      deltaQ = 0.0;
      while (deltaQ < 0.0)
      {
        double QInit = 0.0;
        for(int j =0; j< comm.size(); ++j)
        {
          QInit -=comm[j].modularity;
        }
        
        for(int j =0; j< comm.size(); ++j)
        {
          delete comm[j];
        }
        comm.clear();
        
        RandVer = neighbours[rgen.integer(1, (int) neighbours.size())];
        tempVertexSet[i] = tempVertexSet[RandVer];
        
        for(int it=1; it<verts.size()+1; ++it)
        {
          community* current_comm = NULL;
          for(int j=0; j<verts.size(); ++j)
          {
            
            if(vertexSet[j] == it)
            {
              if(current_comm == NULL)
              {
                current_comm = new community(it);
              }
              current_comm->grow();
              vertex* temp = verts[j];
              current_comm->vertices.push_back(temp);
            }
          }
          if(current_comm !=NULL)
          {
            comm.push_back(current_comm);
          }
        }

        for(int it =0; it<comm.size(); ++it)
        {
          comm[it]->modularity_calc();
        }
        
        double Q = 0.0;
        for(int j =0; j< comm.size(); ++j)
        {
          Q -=comm[j].modularity;
        }
        
        deltaQ = (-1.0 * Q) - (-1.0 * QInit);
      }
      if (deltaQ > 0.0)
      {
        vertexSet[i] = tempVertexSet[RandVer];
      }
      else 
      {
        t++;
      }
    }
  }
}

void individual::decode() {
    int sum;
    for (int i = 0; i < config->nbin; ++i) {
        sum = 0;
        for (int j = 0; j < config->nbits[i]; ++j) {
            sum += (1 << (config->nbits[i]-1-j));  // TODO: check
        }

        xbin[i] = config->limits_binvar[i].first +
            (double)sum*( config->limits_binvar[i].second - config->limits_binvar[i].first) / (double)((1 << (config->nbits[i]))-1); // TODO: check
    }
}

void individual::evaluate() {

    // workaround to respect the signature of test_problem and its (int**)


    (*config->function) (&obj, &comm);

    evaluated = true;
}

// returns:  1 if this < b (this dominates b),
//          -1 if this > b (this is dominated by b),
//           0 if they are nondominated
int individual::check_dominance(const individual& b) const {

        int flag1 = 0, // to check if this has a smaller objective
            flag2 = 0; // to check if b    has a smaller objective

        for (int i=0; i<config->nobj; ++i) {
    	    if (config->nobj > 1) { // Normal multi objective comparison

        		if (obj[i] < b.obj[i]) {
        		    flag1 = 1;
        		} else if (obj[i] > b.obj[i]) {
        		    flag2 = 1;
        		}
    	    } else { // mono objective comparison with an epsilon
        		if (obj[i] < b.obj[i] && fabs(obj[i]-b.obj[i]) > config->epsilon_c) {
        		    flag1 = 1;
        		} else if (obj[i] > b.obj[i] && fabs(obj[i]-b.obj[i]) > config->epsilon_c) {
        		    flag2 = 1;
        		}
    	    }
        }

        if (flag1==1 && flag2==0) {
            // there is at least one smaller objective for this and none for b

            return 1;

        } else if (flag1==0 && flag2==1) {
            // there is at least one smaller objective for b and none for this

            return -1;

        } else {
            // no smaller objective or both have one smaller

            return 0;
        }

    
}

// returns num_mut_real, num_mut_bin
std::pair<int,int> individual::mutate() {
    std::pair<int,int> num_mut = std::make_pair(0,0);
    



    comm_mutate();

    return num_mut;
}

int individual::comm_mutate(){
    vertex* other_vertex;
    vertex* chosen_vertex;
    community* other_comm = NULL;

    for(int i=0; i<comm.size(); ++i)
    {

        int other_community =0;
        if(rgen.realu() <= pmut_comm)
        {

            int this_vertex = rgen.integer(0,comm[i]->vertices.size()-1);
            chosen_vertex = comm[i]->vertices[this_vertex];
        
            
            for(int j=0; j<nvertices; ++j)
            {
                int other_number = rgen.integer(0, nvertices-1);
                

                if((*mat)[chosen_vertex->number][other_number] > 0){

                    other_community = vertexSet[other_number];
                    break;
                }
              
            }
            if(other_community !=0)
            {
                
              
                for(int j=0; j< comm.size(); ++j)
                {

                    if(other_community == comm[j]->number)
                    {
                        other_comm = comm[j];
                        other_comm->vertices.push_back(chosen_vertex);

                        comm[i]->vertices.erase(comm[i]->vertices.begin() + this_vertex);
                        comm[i]->internalDegrees.erase(comm[i]->internalDegrees.begin() + this_vertex);
                        comm[i]->w_internalDegrees.erase(comm[i]->w_internalDegrees.begin() + this_vertex);
            
                    
                        comm[i]->modularity_calc();
                        vertexSet[this_vertex] = other_community;
                        other_comm->modularity_calc();
                        if(comm[i]->vertices.size() == 0)
                        {
                            comm.erase(comm.begin()+i);
                            --i;
                        }
                        break;
                    }
                }


            }
        }
    }

}


int individual::real_mutate() {
    int j;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;
    int num_mut = 0;
    for (j=0; j<config->nreal; j++) {
        if (rgen.realu() <= config->pmut_real) {
            y = xreal[j];
            yl = config->limits_realvar[j].first;
            yu = config->limits_realvar[j].second;
            delta1 = (y-yl)/(yu-yl);
            delta2 = (yu-y)/(yu-yl);
            rnd = rgen.realu();
            mut_pow = 1.0/(config->eta_m+1.0);
            if (rnd <= 0.5) {
                xy = 1.0-delta1;
                val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(config->eta_m+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            } else {
                xy = 1.0-delta2;
                val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(config->eta_m+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }
            y = y + deltaq*(yu-yl);
            if (y<yl)
                y = yl;
            if (y>yu)
                y = yu;
            xreal[j] = y;
            num_mut+=1;
        }
    }
    return num_mut;
}

int individual::bin_mutate() {
    int j, k;
    double prob;
    int num_mut = 0;
    for (j=0; j<config->nbin; j++) {
        for (k=0; k<config->nbits[j]; k++) {
            prob = rgen.realu();
            if (prob <=config->pmut_bin) {
                if (gene[j][k] == 0) {
                    gene[j][k] = 1;
                } else {
                    gene[j][k] = 0;
                }
                num_mut+=1;
            }
        }
    }
    return num_mut;
}


std::ostream& nsga2::operator<< (std::ostream& os, const individual& ind) {

    os << "{Individual rank=" << ind.rank
       << "\nconstr_violation=" << ind.constr_violation;

    os << "\nxreal=[";
    std::vector<double>::const_iterator it;
    for (it = ind.xreal.begin(); it != ind.xreal.end(); ++it) {
        os << *it;
        if (it+1 != ind.xreal.end())
            os << ",";
    }

    os << "]\ngene=";
    std::vector< std::vector<int> >::const_iterator it1;
    for (it1 = ind.gene.begin(); it1 != ind.gene.end(); ++it1) {
        const std::vector<int>& tmp = *it1;
        std::vector<int>::const_iterator it2;
        if (it1 != ind.gene.begin())
            os << "     "; // tab space
        for (it2 = tmp.begin(); it2 != tmp.end(); ++it2) {
            os << *it2;
        }
        //       gene=
        os << '\n';
    }

    os << "xbin=";
    for (it = ind.xbin.begin(); it != ind.xbin.end(); ++it) {
        os << *it;
        if (it+1 != ind.xbin.end())
            os << ",";
    }

    os << "\nobj=";
    for (it = ind.obj.begin(); it != ind.obj.end(); ++it) {
        os << *it;
        if (it+1 != ind.obj.end())
            os << ",";
    }

    os << "\nconstr=";
    for (it = ind.constr.begin(); it != ind.constr.end(); ++it) {
        os << *it;
        if (it+1 != ind.constr.end())
            os << ",";
    }

    os << "\ncrowd_dist=" << ind.crowd_dist;

    os << " }";

    return os;
}

population::population(const int size,
                       const int nreal,
                       const int nbin,
                       const int ncon,
                       const std::vector<int>& nbits,
                       const std::vector< std::pair<double,double> >& limreal,
                       const std::vector< std::pair<double,double> >& limbin,
                       const int nobj,
                       const double pmut_real,
                       const double pmut_bin,
                       const double eta_m,
		       const double epsilon_c,
                       const individual_config::funcType func)
	  throw (nsga2::nsga2exception) :
	  crowd_obj(true),
	  ind_config(),
	  eval_pop_function(NULL) {

    generation = 1;

    ind_config.nobj           = nobj;
    ind_config.ncon           = ncon;

    ind_config.eta_m          = eta_m;
    ind_config.function       = func;
    ind_config.epsilon_c      = epsilon_c;

    for (int i = 0; i < size; ++i) {
        ind.push_back(individual(ind_config));
    }

}




population::~population() {
}

void population::initialize(std::vector<vertex*> verts) throw (nsga2::nsga2exception) {
    std::vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {
        int num_comm = rgen.integer(1, verts.size());
        it->initialize(verts, num_comm);
    }
}

void population::decode() {
    std::vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {
        it->decode();
    }
}

void population::evaluate() {
    normal_evaluate();
}

void population::reset()
{
    for(int i =0; i<ind.size(); ++i)
    {
        ind[i].reset();
    }
}

void population::custom_evaluate() {
    if (eval_pop_function != NULL)
	(*eval_pop_function)(*this);
    else 
	normal_evaluate_openmp();
}

void population::normal_evaluate() {
    for (int i=0; i<popsize; ++i)
    {

        ind[i].evaluate();
    }
}

void population::normal_evaluate_openmp() {
#ifdef USE_OPENMP
#pragma omp parallel for
    for (int i = 0; i < ind.size(); ++i) {
        ind[i].evaluate();
    }
#else
    normal_evaluate();
#endif
}

void population::set_popfunction(individual_config::popFuncType f) {
    eval_pop_function = f;
}

void population::fast_nds() 
{
    front.resize(1);
    front[0].clear();
    //std::vector< std::vector<int> >  F(1);
#pragma omp parallel for
    for (int i = 0; i < ind.size(); ++i) {
	
        std::vector<int> dom;
        int dcount = 0;
	
        individual& p = ind[i];
        // p.dcounter  = 0;
        // p.dominated.clear();
	
        for (int j = 0; j < ind.size(); ++j) {
	    
            individual& q = ind[j];
	    
            int compare = p.check_dominance(q);
            if (compare == 1) { // p dominates q
                //p.dominated.push_back(j);
                dom.push_back(j);
            } else if (compare == -1) { // q dominates p
                //p.dcounter += 1;
                dcount += 1;
            }
        }
	
#pragma omp critical
        {
            p.dcounter  = dcount;
            p.dominated.clear();
            p.dominated = dom;
	    
	    
            if (p.dcounter == 0) {
                p.rank = 1;
                front[0].push_back(i);
            }
        }
	
    }
    
    // using OpenMP can have different orders in the front[0]
    // so let's sort it so that the algorithm is deterministic
    // given a seed
    std::sort(front[0].begin(), front[0].end());    

    int fi = 1;
    while (front[fi-1].size() > 0) {

        std::vector<int>& fronti = front[fi-1];
        std::vector<int> Q;
        for (int i = 0; i < fronti.size(); ++i) {

            individual& p = ind[fronti[i]];

            for (int j = 0; j < p.dominated.size() ; ++j) {


                individual& q = ind[p.dominated[j]];
                q.dcounter -= 1;

                if (q.dcounter == 0) {
                    q.rank = fi+1;
                    Q.push_back(p.dominated[j]);
                }
            }
        }


        fi += 1;
        front.push_back(Q);
    }

}



struct comparator_obj {
    comparator_obj(const population& population, int index) :
        pop(population), m(index) {};
    const population& pop;
    int m;
    bool operator() (int i, int j) {
        return pop.ind[i].obj[m] < pop.ind[j].obj[m];
    };
};

void population::crowding_distance_all() {
    for (int i = 0; i < front.size(); ++i)
        crowding_distance(i);
}

void population::crowding_distance(int fronti) {

    std::vector<int> F = front[fronti];
    if (F.size() == 0 ) return;

    const int l = F.size();

    for (int i = 0; i < l; ++i)
        ind[F[i]].crowd_dist = 0;

    // for (int m = 0; m < ind_config.nobj; ++m) {
    //     std::sort(F.begin(), F.end(), comparator_obj(*this,m));
    //     ind[F[0]].crowd_dist = INF;
    // }

    const int limit = crowd_obj?ind_config.nobj:ind_config.nvertices;
    for (int m = 0; m < limit; ++m) {

        std::sort(F.begin(), F.end(), comparator_obj(*this,m));

        // in the paper dist=INF for the first and last, in the code
        // this is only done to the first one or to the two first when size=2
        // ind[F[0]].crowd_dist = INF;
        // if (l == 2)
        //      ind[F[0]].crowd_dist = ind[F[1]].crowd_dist = INF;
        ind[F[0]].crowd_dist = INF;
        if (l > 1)
            ind[F[l-1]].crowd_dist = INF;
	//cout << "min " << ind[F[0]].xreal[0];
	//cout << "\tmax " << ind[F[l-1]].xreal[0] << endl;

        for (int i = 1; i < l-1; ++i) {
            if (ind[F[i]].crowd_dist != INF) {
                if (crowd_obj && ind[F[l-1]].obj[m] != ind[F[0]].obj[m]) {
                // ind[F[l-1]].xreal[m] != ind[F[0]].xreal[m])

                    ind[F[i]].crowd_dist +=
                        (ind[F[i+1]].obj[m] - ind[F[i-1]].obj[m]) // crowd over obj
                        / (ind[F[l-1]].obj[m] - ind[F[0]].obj[m]);

                 } 
                //  else if (!crowd_obj && ind[F[l-1]].xreal[m] != ind[F[0]].xreal[m]) {
                //     ind[F[i]].crowd_dist +=
                //         (ind[F[i+1]].xreal[m] - ind[F[i-1]].xreal[m]) // crowd over vars
                //         / (ind[F[l-1]].xreal[m] - ind[F[0]].xreal[m]);
                // }
            }
        }
    }

    for (int i=0; i < l; ++i) { // this is deduced from code, not mentioned in paper
        if (ind[F[i]].crowd_dist != INF)
            ind[F[i]].crowd_dist /= ind_config.nobj;
    }
}

void population::merge(const population& pop1, const population& pop2)
    throw (nsga2::nsga2exception) {

    if (size() < pop1.size() + pop2.size())
        throw nsga2::nsga2exception("Merge: target population not big enough");

    std::copy(pop1.ind.begin(), pop1.ind.end(), ind.begin());
    std::copy(pop2.ind.begin(), pop2.ind.end(), ind.begin() + pop1.size());

}

void population::report(std::ostream& os) const {

    std::vector<individual>::const_iterator it;
    os<< "POP SIZE: " << ind.size() << endl;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {

        for (int j = 0; j < ind_config.nobj; ++j)
            os << it->obj[j] << '\t';

        for (int j = 0; j < nvertices; ++j)
            os << it->vertexSet[j] << ' ';
       os<< '\n';
        // os << it->constr_violation << '\t'
        //    << it->rank << '\t'
        //    << it->crowd_dist << '\n';

    }
}
void population::report_obj(std::ostream& os) const {

    std::vector<individual>::const_iterator it;
    
    for (it  = ind.begin();
         it != ind.end();
         ++it) {
        if(it->rank ==1){
                    for (int j = 0; j < ind_config.nobj; ++j)
            os << it->obj[j] << '\t';
             os<< '\n';
        }



  
        // os << it->constr_violation << '\t'
        //    << it->rank << '\t'
        //    << it->crowd_dist << '\n';

    }
}

void population::dump(std::ostream& os) const {

    std::vector<individual>::const_iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {

        if (ind_config.nobj > 0)
            os.write(reinterpret_cast<const char*>(&(it->obj[0])),
                     sizeof(double)*ind_config.nobj);

        if (ind_config.ncon > 0)
            os.write(reinterpret_cast<const char*>(&(it->constr[0])),
                     sizeof(double)*ind_config.ncon);

        if (ind_config.nreal > 0)
            os.write(reinterpret_cast<const char*>(&(it->xreal[0])),
                     sizeof(double)*ind_config.nreal);

        for (int j = 0; j < ind_config.nbin; ++j)
            os.write(reinterpret_cast<const char*>(&(it->gene[j][0])),
                     sizeof(int)*ind_config.nbits[j]);

        os.write(reinterpret_cast<const char*>(&(it->constr_violation)),
                 sizeof(double));
        os.write(reinterpret_cast<const char*>(&(it->rank)),
                 sizeof(int));
        os.write(reinterpret_cast<const char*>(&(it->crowd_dist)),
                 sizeof(double));

    }
}

void population::load(std::istream& os) {

    std::vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {

        if (ind_config.nobj > 0)
            os.read(reinterpret_cast<char*>(&(it->obj[0])),
                     sizeof(double)*ind_config.nobj);

        if (ind_config.ncon > 0)
            os.read(reinterpret_cast<char*>(&(it->constr[0])),
                     sizeof(double)*ind_config.ncon);

        if (ind_config.nreal > 0)
            os.read(reinterpret_cast<char*>(&(it->xreal[0])),
                     sizeof(double)*ind_config.nreal);

        for (int j = 0; j < ind_config.nbin; ++j)
            os.read(reinterpret_cast<char*>(&(it->gene[j][0])),
                     sizeof(int)*ind_config.nbits[j]);

        os.read(reinterpret_cast<char*>(&(it->constr_violation)),
                 sizeof(double));
        os.read(reinterpret_cast<char*>(&(it->rank)),
                 sizeof(int));
        os.read(reinterpret_cast<char*>(&(it->crowd_dist)),
                 sizeof(double));

    }
}

std::pair<int,int> population::mutate() {
    std::pair<int,int>
        num_mut = std::make_pair(0,0),
        tmp = std::make_pair(0,0);
    std::vector<individual>::iterator it;
    for (it  = ind.begin();
         it != ind.end();
         ++it) {
        tmp = it->mutate();
        num_mut.first  += tmp.first;
        num_mut.second += tmp.second;
    }
    return num_mut;
}

std::ostream& nsga2::operator<< (std::ostream& os, const population& pop) {
    os << "Population: {\n";
    std::vector<individual>::const_iterator it;
    for (it = pop.ind.begin(); it != pop.ind.end(); ++it) {
        os << *it;
    }
    os << '}';
    return os;
}
