#ifndef MVM_H
#define MVM_H

#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>

using std::vector;
using namespace std;

class System_site{

	public:
		int vneighbors*;
		double num_neighbor_states*;
		int total_neighbor_states;
		int state;
		double own_wealth;
		double mutualistic_wealth;
		double lattice_wealth;
		double total_wealth;
};

class System::public Lattice{

	private:
	int num_states;
	double Aij**;
	
	public:
	Dynamics(int dim, int L, int Ns, double aij**): Lattice(dim,L), num_states(Ns), Aij(aij){
	
		for (int i = 0; i < num_sites; ++i) {
    			list_sites[i].frac_states = new int[Ns];
    		}
	};
	
	void set_states();
	void initialize_num_neighbor_states();
	void initialize_wealths();
};

void System::set_states(){ //Particular case! Maybe another subclass for this particular cases of set_states;

	int site;
	
	for(int i=0;i<linear_size;i++){
		list_sites[i].state=gsl_rng_uniform_int(r,num_states);
		for(int j=1;j<linear_size;j++){
			site=j*linear_size+i;
			list_sites[site].state=-1;
		}
	}
}

void System::initialize_num_neighbor_states(){

	int neighbor,state;
	vector <int> counted_states;
	
	for(int i=0;i<num_sites;i++){
		if(list_sites[i].state>0){
			for(int j=0;j<num_states;j++) list_sites[i].num_neighbor_states[j]=0;
			counted_states.clear();
			for(int j=0;j<num_neighbors;j++){
				neighbor=list_sites[i].vneighbor[j];
				state=list_sites[neighbor].state;
				if(state>0){
					list_sites[i].num_neighbor_states[state]++;
					if(counted_states.size()==0) counted_states.push_back(state);
					else{
						for(k=0;k<counted_states.size();k++){
					
							if(state!=counted_states[i]) counted_states.push_back(state);
						}
					}
				}
			}
			list_sites[i].total_neighbor_states=counted_states.size();
		}
	}
}

void System::initialize_wealths(){ //After set_states!

	int state;
	int num_states;
	int tot_states;
	
	
	for(int i=0;i<num_sites;i++){
		list_sites[i].own_wealth=1;
		list_sites[i].lattice_wealth=0; //Ordered case.
		list_sites[i].mutualistic_wealth=0;
		state=list_sites[i].state;
		if(state>0){
			tot_states=list_sites[i].total_neighbors_states;
			for(int j=0;j<num_states;j++){
				num_states=list_sites[i].num_neighbor_states[j];
				list_sites[i].mutualistic_wealth=Aij[state][j]*(num_states/double(tot_states));
			}
		}
		list_sites[i].lattice_wealth=list_sites[i].own_wealth+list_sites[i].lattice_wealth+list_sites[i].mutualistic_wealth;
	}
}


class Lattice{
	
	protected:
		const int num_sites;
		const int dimension;
		const int linear_size;
		const int num_neighbors;
	public:
		System_site list_sites*;
		void set_neighbors();
		
		Lattice(int dim, int L): dimension(dim),linear_size(L),num_sites(pow(linear_size,dim)),num_neighbors(2*dim) 
		{
			list_sites =  new System_site [num_sites];
			for (int i = 0; i < num_sites; ++i) {
    				list_sites[i].vneighbors = new int[num_neighbors];
    			}
		};
};

void Lattice::set_neighbors(){

	//Network 1D
	if(dimension==1){
			
		for(int i=1;i<num_sites-1;i++){
		
			list_sites[i].vneighbors[0]=i-1;
			list_sites[i].vneighbors[1]=i+1;
		}		
			
		list_sites[0].vneighbors[0]=num_sites-1;
		list_sites[0].vneighbors[1]=1;
	
		list_sites[num_sites-1].vneighbors[0]=num_sites-2;
		list_sites[num_sites-1].vneighbors[1]=0;
	}
	
	//Network 2D		
	else if(dimension==2){
	
		int site;
			
		for(int isite=0;isite<linear_size;isite++){
			for(int jsite=0;jsite<linear_size;jsite++){
						
				site=isite*linear_size+jsite;
	
				list_sites[site].vneighbors[0]=((isite+1+linear_size)%linear_size)*linear_size+jsite;
				list_sites[site].vneighbors[1]=((isite-1+linear_size)%linear_size)*linear_size+jsite;
				list_sites[site].vneighbors[2]=isite*linear_size+(jsite+1+linear_size)%linear_size;
				list_sites[site].vneighbors[3]=isite*linear_size+(jsite-1+linear_size)%linear_size;
				
				//list_sites[site].row=isite;
				//list_sites[site].column=jsite;
		
			}
		}
	}
}







//---------------------- Lattice class-------------
class Lattice{

		
		Lattice(int dim, int L, int Ns, vector<double>& vfreq_specie): dimension(dim),linear_size(L),num_species(Ns) 
		{
			num_sites=pow(linear_size,dimension);
			num_neighbors=2*dimension;
			list_sites = new Site[num_sites];
			//Site site_aux;
			//list_sites.clear();
			//for(int i=0;i<num_sites;i++) list_sites.push_back(site_aux); //preguntar a Jose mejor forma de hacerlo!
		};
		~Lattice();
		void set_neighbors();
		int set_specie(gsl_rng *r);
		int get_neighbor0(int index);
};

Lattice::~Lattice(){ cout<<"destroyed"<<endl; }
int Lattice::get_sites(){ return num_sites; }
int Lattice::get_specie(int index){ return list_sites[index].specie;	 }
int Lattice::get_neighbor0(int index){ return list_sites[index].vneighbors[0]; }

int Lattice::set_specie(gsl_rng *r){

	double sum;
	double sum_ant;
	double rand_aux;
	
	//for(int j=0;j<num_species;j++) //cout<<freq_specie[j]<<endl;
	
	for(int i=0;i<num_sites;i++){
		sum=0;
		sum_ant=0;
		rand_aux=gsl_rng_uniform(r);
		for(int j=0;j<num_species;j++){
			sum+=freq_specie[j];
			if((rand_aux>=sum_ant)&&(rand_aux<sum)) list_sites[i].specie=j;
			sum_ant=sum;
		}
	}
	
	//for(int i=0;i<num_sites;i++) //cout<<i<<" "<<list_sites[i].specie<<endl;
}

void Lattice::set_neighbors(){
	
	for(int i=0;i<num_sites;i++) list_sites[i].vneighbors.clear();
	//Network 1D
	if(dimension==1){
			
		for(int i=1;i<num_sites-1;i++){
		
			list_sites[i].vneighbors.push_back(i-1);
			list_sites[i].vneighbors.push_back(i+1);
		}		
			
		list_sites[0].vneighbors.push_back(num_sites-1);
		list_sites[0].vneighbors.push_back(1);
	
		list_sites[num_sites-1].vneighbors.push_back(num_sites-2);
		list_sites[num_sites-1].vneighbors.push_back(0);
	}
	
	//Network 2D		
	else if(dimension==2){
	
		int site;
			
		for(int isite=0;isite<linear_size;isite++){
			for(int jsite=0;jsite<linear_size;jsite++){
						
				site=isite*linear_size+jsite;
	
				list_sites[site].vneighbors.push_back(((isite+1+linear_size)%linear_size)*linear_size+jsite);
				list_sites[site].vneighbors.push_back(((isite-1+linear_size)%linear_size)*linear_size+jsite);
				list_sites[site].vneighbors.push_back(isite*linear_size+(jsite+1+linear_size)%linear_size);
				list_sites[site].vneighbors.push_back(isite*linear_size+(jsite-1+linear_size)%linear_size);
				
				//list_sites[site].row=isite;
				//list_sites[site].column=jsite;
		
			}
		}
	}	
}

//------------------Class Dynamics----

class Dynamics: public Lattice{

	private:
		vector<int> vinterfases;
		vector<int> vinterfase_positions;
		int check_interfase(int index);
		void add_interfase(int index);
		void remove_interfase(int index);
	public:
		Dynamics(int dim, int L, int Ns,vector<double>& vfreq_specie):Lattice(dim,L,Ns,vfreq_specie){};
		int num_interfases;
		int num_vinterfases;
		void get_interfases();
		int get_die(gsl_rng *r);
		int get_reproduce(gsl_rng *r, int site_die);
		void reproduce(int sdie,int srepr);
		void clean_variables();
};
void Dynamics::clean_variables(){

	vinterfases.clear();
	vinterfase_positions.clear();
}
void Dynamics::get_interfases(){
	
	vinterfases.clear();
	vinterfase_positions.clear();
	vinterfase_positions.assign(num_sites,0);
	
	//for(int i=0;i<num_sites;i++) //cout<<i<<" "<<get_specie(i)<<" "<<check_interfase(i)<<endl;
	
	for(int i=0;i<num_sites;i++){
		if(check_interfase(i)==1) add_interfase(i);
		else vinterfase_positions[i]=-1;
		
		//for(int j=0;j<vinterfases.size();j++) //cout<<vinterfases[j]<<" ";
		////cout<<endl;	
	}
	
	num_vinterfases=vinterfases.size();
	num_interfases=round(num_vinterfases/2.);
	if((num_interfases%2!=0)&&(num_interfases!=1)) num_interfases++;

	//for(int i=0;i<num_interfases;i++) //cout<<vinterfases[i]<<" ";
	////cout<<endl;	
}

int Dynamics::check_interfase(int index){

	int specie1,specie2;
	int neighbor_aux;
	
	for(int j=0;j<num_neighbors;j++){
		neighbor_aux=list_sites[index].vneighbors[j];
		specie1=list_sites[index].specie;
		specie2=list_sites[neighbor_aux].specie;
		if(specie1!=specie2) return 1;
	}
	return 0;
}
void Dynamics::add_interfase(int index){

	vinterfase_positions[index]=vinterfases.size();
	vinterfases.push_back(index);
}
void Dynamics::remove_interfase(int site){

	int last_index=vinterfases.size()-1;
	int last_element=vinterfases[last_index];
	int current_index=vinterfase_positions[site];
	
	vinterfases[current_index]=last_element;
	vinterfase_positions[last_element]=vinterfase_positions[site];
	vinterfase_positions[site]=-1;
	
	vinterfases.pop_back();
}

int Dynamics::get_die(gsl_rng *r){ return vinterfases[gsl_rng_uniform_int(r, num_vinterfases)]; }
int Dynamics::get_reproduce(gsl_rng *r, int site_die){ return list_sites[site_die].vneighbors[gsl_rng_uniform_int(r, num_neighbors)]; }
void Dynamics::reproduce(int sdie,int srepr){

	/*for(int i=0;i<num_sites;i++) //cout<<i<<"("<<check_interfase(i)<<") ";
	//cout<<endl;
	//cout<<"  ";
	for(int i=0;i<num_sites;i++) //cout<<get_specie(i)<<"    ";
	//cout<<endl;
	for(int i=0;i<num_interfases;i++) //cout<<vinterfases[i]<<" ("<<vinterfase_positions[vinterfases[i]]<<") ";
	//cout<<endl;*/
	
	int neighbor;
	if(list_sites[sdie].specie!=list_sites[srepr].specie){
	
		list_sites[sdie].specie=list_sites[srepr].specie;
		if(check_interfase(sdie)==0) remove_interfase(sdie);
		
		for(int i=0;i<num_neighbors;i++){
			neighbor=list_sites[sdie].vneighbors[i];
			if(vinterfase_positions[neighbor]>=0){
			
				if(check_interfase(neighbor)==0) remove_interfase(neighbor);
			}
			else{
			
				if(check_interfase(neighbor)==1) add_interfase(neighbor);
			}
		}
	}
	num_vinterfases=vinterfases.size();
	num_interfases=round(num_vinterfases/2.);
	if((num_interfases%2!=0)&&(num_interfases!=1)) num_interfases++;
	
	/*//cout<<"Selected >>> "<<sdie<<" "<<srepr<<endl;
	for(int i=0;i<num_sites;i++) //cout<<i<<"("<<check_interfase(i)<<") ";
	//cout<<endl;
	//cout<<"  ";
	for(int i=0;i<num_sites;i++) //cout<<get_specie(i)<<"    ";
	//cout<<endl;
	for(int i=0;i<num_interfases;i++) //cout<<vinterfases[i]<<" ("<<vinterfase_positions[vinterfases[i]]<<") ";
	//cout<<endl;
	//cout<<endl;
	//cout<<endl;
	//cout<<endl;*/
}

//-------------------- Class simulation

class Simulation{

	protected:
		const int tfin;
		const int trials;
		const int refresh_time;
		const int rand_seed;
		
		Dynamics dynamics;
		
		int time_increment;
		time_t itime, ftime;
		fstream FILE;
		const char* CURRENT_FILE_NAME;
		
		Simulation(Simulation&);
		void operator=(Simulation&);
		
	public:
		const int num_points;
		Simulation(Dynamics dynamics_in, int T, int itrials, int irefresh_time,int iseed, char const* FILE_NAME): 
		dynamics(dynamics_in),tfin(T),trials(itrials),refresh_time(irefresh_time),rand_seed(iseed),num_points(int(log10(double(tfin))*2e2)+2),CURRENT_FILE_NAME(FILE_NAME){
		
			//num_points=int(log10(double(tfin))*2e2)+1;
			
			//cout<<num_points;
			
			//FILE.open(FILE_NAME, ios::out | ios::trunc);
		}; 
		
		//General simulation
		void sequential_trial(gsl_rng *r);
		void promediate();
		
		//Determined at the concret simulation
		virtual void initialice_time_variables(gsl_rng *r){};
		virtual void initialice_time_measures() {};
		virtual void take_initial_measures() {};
		virtual void take_time_measures(int ipoint) {};
		virtual bool conditions_run(double t){ return false; };
		virtual double get_time_increment() { return 1.; };
		virtual void sequential_time_step(gsl_rng *r){};
		virtual void initiate_counter(){};
		virtual void initialice_global_variables() {};
		virtual void initialice_global_measures() {};
		virtual void write_measures(int trial) {};
};

void Simulation::promediate(){ //Maybe here the random generator!?

	//// ------ GENERADOR DE NUMEROS ALEATORIOS ------//
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	//if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
	gsl_rng_default_seed=rand_seed; //pongo semilla a mano
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	////----------------------------------------------//
	
	initialice_global_variables();
	//cout<<" variable global"<<endl;
	initialice_global_measures();
	//cout<<" meas glob"<<endl;
	itime=time(NULL);
	
	for(int itrial=0;itrial<trials;itrial++){
	
		sequential_trial(r);
		//cout<<" seq trial"<<endl;
		write_measures(itrial);
	}
	
	gsl_rng_free (r);
}

void Simulation::sequential_trial(gsl_rng *r){

	double t=1;
	int ipoint=0;
	
	initialice_time_variables(r);
	//cout<<" variables(0) "<<endl;
	initialice_time_measures();
	//cout<<" measures(0) "<<endl;
	take_initial_measures();
	//cout<<" measures take(0) "<<endl;
	while(conditions_run(t)){
	
		t+=get_time_increment();
		//cout<<" time_incr "<<endl;
		ipoint=int(log10(t)*2e2)+1;
		//cout<<"Original: "<<t<<" "<<int(log10(t)*2e2)+1<<" "<<log10(t)<<" "<<ipoint<<endl;
		sequential_time_step(r);
		//cout<<" sequential "<<endl;
		take_time_measures(ipoint);
		initiate_counter();
		//cout<<t<<" ";
		//for(int i=0;i<dynamics.num_sites;i++) cout<<dynamics.get_specie(i)<<" ";
		//cout<<"-->"<<dynamics.num_interfases<<endl;
	}
	
	dynamics.clean_variables();
}

//----------

class Simulation_collisions: public Simulation{

	private:
	 	int num_collisions;
		int num_anh;
		int num_coal;
		int num_triplets_coal;
		int num_triplets_anh;
		vector<double> p_anh;
		vector<double> p_coal;
		vector<double> prom_num_interfases;
		vector<double> p_triplets_anh;
		vector<double> p_triplets_coal;
	protected:
		vector<double> points;
	public:
		Simulation_collisions(Dynamics dynamics_in, int T, int itrials, int irefresh_time,int iseed, char const* FILE_NAME)
		: Simulation(dynamics_in,T,itrials,irefresh_time,iseed,FILE_NAME){};
		
		virtual void sequential_time_step(gsl_rng *r);
		virtual void initialice_time_variables(gsl_rng *r);
		virtual void initialice_time_measures();
		//void take_initial_measures() {}; //In this concrete case there aren't initial measures!
		virtual void take_time_measures(int ipoint);
		virtual bool conditions_run(double t);
		virtual double get_time_increment();
		virtual void get_collisions(int index_die);
		virtual void initiate_counter();
		virtual void get_triplets();
		virtual void initialice_global_variables();
		virtual void initialice_global_measures();
		virtual void write_measures(int itrial);
};
	
void Simulation_collisions::initialice_global_variables(){ 

	dynamics.set_neighbors();
}

void Simulation_collisions::initialice_global_measures(){ 

	p_anh.assign(num_points,0);
	p_coal.assign(num_points,0);
	p_triplets_anh.assign(num_points,0);
	p_triplets_coal.assign(num_points,0);
	points.assign(num_points,0);
	prom_num_interfases.assign(num_points,0);
}

void Simulation_collisions::initialice_time_variables(gsl_rng *r){
 
	dynamics.set_specie(r);
	dynamics.get_interfases();
}

void Simulation_collisions::initialice_time_measures(){

	num_coal=0;
	num_anh=0;
	num_collisions=0;
	num_triplets_coal=0;
	num_triplets_anh=0;
	
}

void Simulation_collisions::initiate_counter(){

	num_coal=0;
	num_anh=0;
	num_collisions=0;
	num_triplets_coal=0;
	num_triplets_anh=0;
}

bool Simulation_collisions::conditions_run(double t){

	if((t<tfin)&&(dynamics.num_vinterfases>0)) return true;
	else return false;
}

double Simulation_collisions::get_time_increment(){ return 1./double(dynamics.num_vinterfases); }
void Simulation_collisions::sequential_time_step(gsl_rng *r){
	
	int index_die,index_repr;

	get_triplets();
	index_die=dynamics.get_die(r);
	get_collisions(index_die);
	index_repr=dynamics.get_reproduce(r,index_die);
	dynamics.reproduce(index_die,index_repr);
}

void Simulation_collisions::get_collisions(int index_die){
	
	int specie1,specie2,specie3;
	int neighbor1,neighbor2;
	int comp1,comp2,comp3;
	
	neighbor1=dynamics.list_sites[index_die].vneighbors[0];
	neighbor2=dynamics.list_sites[index_die].vneighbors[1];
	
	specie1=dynamics.list_sites[neighbor1].specie;
	specie2=dynamics.list_sites[index_die].specie;
	specie3=dynamics.list_sites[neighbor2].specie;
	
	comp1=specie1-specie2;
	comp2=specie2-specie3;
	comp3=specie1-specie3;
	
	if((comp1!=0)&&(comp2!=0)){
	
		num_collisions++;
		if(comp3!=0) num_coal++;
		else num_anh++;
	}
}

void Simulation_collisions::get_triplets(){

	int specie1,specie2,specie3;
	int neighbor1,neighbor2;
	int comp1,comp2,comp3;
	
	for(int i=0;i<dynamics.num_vinterfases;i++){
	
		neighbor1=dynamics.list_sites[i].vneighbors[0];
		neighbor2=dynamics.list_sites[i].vneighbors[1];
	
		specie1=dynamics.list_sites[neighbor1].specie;
		specie2=dynamics.list_sites[i].specie;
		specie3=dynamics.list_sites[neighbor2].specie;
	
		comp1=specie1-specie2;
		comp2=specie2-specie3;
		comp3=specie1-specie3;
	
		if((comp1!=0)&&(comp2!=0)){
		
			if(comp3!=0) num_triplets_coal++;
			else num_triplets_anh++;
		}	
	}
}

void Simulation_collisions::take_time_measures(int ipoint){

	p_triplets_anh[ipoint]+=num_triplets_anh;///double(dynamics.num_vinterfases);
	p_triplets_coal[ipoint]+=num_triplets_coal;///double(dynamics.num_vinterfases);
	p_anh[ipoint]+=num_anh;
	p_coal[ipoint]+=num_coal;
	points[ipoint]++;
	prom_num_interfases[ipoint]+=dynamics.num_vinterfases;
}

void Simulation_collisions::write_measures(int itrial){

	int tlogaux,tlog;
	
	ftime=time(NULL);
	if(difftime(ftime, itime)>refresh_time){
		
		FILE.open(CURRENT_FILE_NAME, ios::out | ios::trunc);
		tlogaux=0.;
		for(int i=0;i<num_points;i++){
			tlog=int(pow(10,((i-1)/(2e2))));
			if((points[i]>1)&&(tlog!=tlogaux)){ 
				FILE<<tlog<<" ";
				FILE<<p_anh[i]/points[i]<<" ";
				FILE<<p_coal[i]/points[i]<<" ";
				if((p_anh[i]/points[i]+p_coal[i]/points[i])!=0) FILE<<(p_anh[i]+p_coal[i])/points[i]<<" ";
				else FILE<<1<<" ";
				FILE<<p_triplets_anh[i]/points[i]<<" ";
				FILE<<p_triplets_coal[i]/points[i]<<" ";
				if((p_triplets_anh[i]/points[i]+p_triplets_coal[i]/points[i])!=0) FILE<<(p_triplets_anh[i]+p_triplets_coal[i])/points[i]<<" ";
				else FILE<<1<<" ";
				FILE<<prom_num_interfases[i]/points[i]<<" ";
				FILE<<itrial<<endl;
			}
			tlogaux=tlog;
		}
		FILE.close();
		itime=time(NULL);
	}
}



#endif
