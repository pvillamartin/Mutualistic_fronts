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
		int *vneighbors;								//Vector of neighbors of the site
		int row,column;								//Fila y columna del 
		int *num_neighbor_states;					//Vector of states of the neighbors of the site
		int num_active_neighbors;
		int state;									//State of the site
		double own_wealth;							//Basice wealth
		double mutualistic_wealth;					//Wealth dued to mutualism
		double lattice_wealth;						//Quenqued disorder would affect it. Ordered case would have this =0
		double total_wealth;							//Total wealth of the site sum of the other contributions
};

// -------------Class System

class System{

	private:
	
	//Lattice variables
	const int dimension;								//Dimension of the system
	const int linear_size;							//Linear size of the system
	const int num_sites;								//Number of sites in the system
	const int neighbors_width;						//Radius of nn: n=1->8nn
	const int num_neighbors;							//Number of neighbors of each site

	//Dynamics variables
	const int num_states;									//Number of possible states of each site
	const double p_swap;
	double **Aij;									//Mutualism matrix, each aij determines de strength of mutualism between ij
	vector <int> list_empty_front_sites;				//List of empty sites near the front that can be occupied
	vector <int> list_active_front_sites;					//List of active sites
	int *position_empty_list;						//Position in the vector of empty sites
	int *position_active_list;						//Position in the vector of empty sites
	double own_wealth_value;
	int initial_width;

	//Private funcions for the dynamics of the system
	void get_num_neighbor_states(int index);			//Obtain the number of neighbors of each state
	void add_empty_site(int index);					//Add empty site to list_empty_front_sites
	void add_active_site(int index);					//Add empty site to list_empty_front_sites
	void remove_empty_site(int site);				//Remove empty site from list_empty_front_sites
	void remove_active_site(int site);				//Remove empty site from list_empty_front_sites
	bool check_toinactive_front(int site);
	bool check_toactive_front(int site);
	
	double mututalistic_term(int site, int site_state,int neighbor_state);
	void add_num_neighbor_states(int index, int istate);
	void remove_num_neighbor_states(int index, int istate);

	public:
	//Lattice variables
	System_site *list_sites;							//List of the sites of the system
	
	//Dynamic variables
	int num_empty_sites;								//Number of empty sites close to an active site
	int num_active_sites;	
	
	//Contructor of the class
	System(int dim, int L, int nw, int Ns, double **aij, double pswap, int initial_w): 
	dimension(dim),linear_size(L),neighbors_width(nw),num_sites(pow(linear_size,dim)),num_neighbors((2*nw+1)*(2*nw+1)-1), num_states(Ns), p_swap(pswap), initial_width(initial_w)
	{
		own_wealth_value=0.1;
		//initial_width=nw;
		
		position_empty_list = new int[num_sites];
		position_active_list = new int[num_sites];
		list_sites =  new System_site [num_sites];
		for (int i = 0; i < num_sites; ++i) {
    				list_sites[i].vneighbors = new int[num_neighbors];
    				list_sites[i].num_neighbor_states = new int[num_states];
   		}
   		
   		
   		Aij =  new double *[num_states];
		for (int i = 0; i < num_states; ++i) Aij[i] = new double[num_states];
		for(int i=0;i<num_states;i++){
			for(int j=0;j<num_states;j++) Aij[i][j]=aij[i][j];
		}	
	};

	//Lattice functions
	void set_neighbors();							//Set the neighbors of each site (at the moment just square)
	int get_dimension();
	int get_linear_size();
	
	//System functions
	void set_states(gsl_rng *r);						//Set initial random states
	void initialize_list_empty_sites();				//Initialize empty sites close to active sites
	void initialize_list_active_sites();
	void initialize_num_neighbor_states();			//Initialize the number of neighbors of each state
	void initialize_wealths();						//Initialize the wealths of sites
	void get_mutualistic_wealth(int istate,int index);			//Get the mutualistic wealth depending of the neighbors that sites has at that moment.
	
	//Dynamics fucntions
	int get_empty_site(gsl_rng *r);						//Get the empty site that is going to be occupied
	int get_reproduce_site(gsl_rng *r,int empty_site);	//Get the neighbor of that site that is going to reproduce
	void reproduce(int empty_site, int reproduce_site);	//Reproduce the site and actualice the variables
	void update_neighbor_reproduction(int neighbor,int neighbor_state, int new_state);
	void update_neighbor_swap(int site, int site_state, int old_neighbor_state, int new_neighbor_state);
	void swap(gsl_rng *r);											//Swap to neighbors
	void clean_variables();								//Clean memory for the next trial
	
};

int System::get_dimension() { return dimension; }
int System::get_linear_size() { return linear_size; }

void System::set_neighbors(){

	//Lattice 1D
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
	
	//lattice 2D		
	else if(dimension==2){
	
		int site,ineighbor;
	
		for(int isite=0;isite<linear_size;isite++){
			for(int jsite=0;jsite<linear_size;jsite++){
		
				site=isite*linear_size+jsite;
				ineighbor=0;
				for(int ksite=1;ksite<=neighbors_width;ksite++){
			
					//recorro fila derecha (isite-ksite->site+ksite, jsite+ksite)
					for(int lsite=isite-ksite;lsite<=isite+ksite;lsite++){
						list_sites[site].vneighbors[ineighbor]=((lsite+linear_size)%linear_size)*linear_size+(jsite+ksite+linear_size)%linear_size;
						ineighbor++;
					}
					//recorro fila izquierda (isite-ksite->site+ksite, jsite-ksite)
					for(int lsite=isite-ksite;lsite<=isite+ksite;lsite++){
						list_sites[site].vneighbors[ineighbor]=((lsite+linear_size)%linear_size)*linear_size+(jsite-ksite+linear_size)%linear_size;
						ineighbor++;
					}
					//recorro fila abajo (isite+ksite, jsite-ksite->jsite+ksite)
					for(int lsite=jsite-ksite+1;lsite<=jsite+ksite-1;lsite++){
						list_sites[site].vneighbors[ineighbor]=((isite+ksite+linear_size)%linear_size)*linear_size+(lsite+linear_size)%linear_size;
						ineighbor++;
					}
					//recorro fila arriba (isite-ksite, jsite-ksite->jsite+ksite)
					for(int lsite=jsite-ksite+1;lsite<=jsite+ksite-1;lsite++){
						
						list_sites[site].vneighbors[ineighbor]=((isite-ksite+linear_size)%linear_size)*linear_size+(lsite+linear_size)%linear_size;
						ineighbor++;
					}
				}
				
				list_sites[site].row=isite;
				list_sites[site].column=jsite;
			}
		}	
	}
}

void System::set_states(gsl_rng *r){

	int site;
	
	for(int i=0;i<linear_size;i++){
		for(int j=0;j<initial_width;j++){
			site=j*linear_size+i;
			list_sites[site].state=gsl_rng_uniform_int(r,num_states);
		}
		for(int k=initial_width;k<linear_size;k++){
			site=k*linear_size+i;
			list_sites[site].state=-1;
		}
	}
}

void System::initialize_list_empty_sites(){

	int site;
	list_empty_front_sites.clear();
	
	for(int i=0;i<num_sites;i++) position_empty_list[i]=-1;
	
	for(int i=0;i<linear_size;i++){
		site=initial_width*linear_size+i;
		position_empty_list[site]=list_empty_front_sites.size();
		list_empty_front_sites.push_back(site);
	}
	num_empty_sites=list_empty_front_sites.size();
}

void System::initialize_list_active_sites(){

	int site;
	list_active_front_sites.clear();
	
	for(int i=0;i<num_sites;i++) position_active_list[i]=-1;
	
	for(int i=0;i<linear_size;i++){
		site=(initial_width-1)*linear_size+i;
		//cout<<site<<" "<<list_active_front_sites.size();
		position_active_list[site]=list_active_front_sites.size();
		list_active_front_sites.push_back(site);
	}
	num_active_sites=list_active_front_sites.size();
}

void System::initialize_num_neighbor_states(){

	int site,neighbor;
	int state;
	
	for(int i=0;i<num_sites;i++){
		for(int j=0;j<num_states;j++) list_sites[i].num_neighbor_states[j]=0;
		list_sites[i].num_active_neighbors=0;
	}
	
	for(int i=0;i<linear_size;i++){
		for(int k=0;k<initial_width;k++){
			site=k*linear_size+i;
			state=list_sites[site].state;
			for(int j=0;j<num_neighbors;j++){
				neighbor=list_sites[site].vneighbors[j];
				list_sites[neighbor].num_neighbor_states[state]++;
				list_sites[neighbor].num_active_neighbors++;
			}
		}
	}
}

void System::initialize_wealths(){ //After set_states!

	int state;

	for(int i=0;i<num_sites;i++){
		list_sites[i].own_wealth=own_wealth_value;
		list_sites[i].lattice_wealth=0; //Ordered case 0;
		state=list_sites[i].state;
		if(state>=0){
			get_mutualistic_wealth(state,i);
			list_sites[i].total_wealth=list_sites[i].own_wealth+list_sites[i].lattice_wealth+list_sites[i].mutualistic_wealth;
		}
		else{
			list_sites[i].mutualistic_wealth=0;
			list_sites[i].total_wealth=0;
		}
	}
}

void System::get_mutualistic_wealth(int istate,int index){

	int neighbor,num_state,num_total;

	list_sites[index].mutualistic_wealth=0;
	for(int j=0;j<num_states;j++){
		num_state=list_sites[index].num_neighbor_states[j];
		num_total=list_sites[index].num_active_neighbors;
		list_sites[index].mutualistic_wealth+=Aij[istate][j]*(num_state/double(num_total));
	}
}

int System::get_empty_site(gsl_rng *r){

	int rand_num=gsl_rng_uniform_int(r,num_empty_sites);
	return list_empty_front_sites[rand_num];
}

int System::get_reproduce_site(gsl_rng *r,int empty_site){

	int neighbor;
	double sum_neighbors_weights;
	double rand_num=gsl_rng_uniform(r);
	double sum,sum_ant;
	
	sum_neighbors_weights=0;
	for (int i=0;i<num_neighbors;i++){
	
		neighbor=list_sites[empty_site].vneighbors[i];
		sum_neighbors_weights+=list_sites[neighbor].total_wealth;
	}
	sum_ant=0;
	sum=0;
	for(int i=0;i<num_neighbors;i++){
		neighbor=list_sites[empty_site].vneighbors[i];
		sum+=list_sites[neighbor].total_wealth/sum_neighbors_weights;
		if((rand_num>sum_ant)&&(rand_num<sum)) return neighbor;
		sum_ant=sum;
	}
}

void System::reproduce(int empty_site, int reproduce_site){

	int neighbor,neighbor_state;
	int new_state=list_sites[reproduce_site].state;

	//Change the state of the site
	list_sites[empty_site].state=new_state;
	
	//Remove it from the list of empty sites
	remove_empty_site(empty_site);
	num_empty_sites=list_empty_front_sites.size();
	
	//Changes frac_states and mutualistic_wealth of the neighbors of empty site!
	for(int i=0;i<num_neighbors;i++){
		
		neighbor=list_sites[empty_site].vneighbors[i];
		neighbor_state=list_sites[neighbor].state;
		if(neighbor_state>=0){
		
			update_neighbor_reproduction(neighbor,neighbor_state,new_state);
			if(check_toinactive_front(neighbor)) remove_active_site(neighbor);
		}
		else{
			
			add_num_neighbor_states(neighbor,new_state);
			list_sites[neighbor].num_active_neighbors++;
			if(position_empty_list[neighbor]<0){
			
				add_empty_site(neighbor); //If it is an empty site it is added to the list of empty sites near the front!
			}
		}
	}
	
	if(check_toactive_front(empty_site)) add_active_site(empty_site);
	get_mutualistic_wealth(new_state,empty_site);
	list_sites[empty_site].total_wealth=list_sites[empty_site].own_wealth+list_sites[empty_site].lattice_wealth+list_sites[empty_site].mutualistic_wealth;
	
	//cout<<num_active_sites<<endl;
}

bool System::check_toinactive_front(int site){

	if(list_sites[site].num_active_neighbors==num_neighbors) return true;
	else return false;
}

bool System::check_toactive_front(int site){

	if(list_sites[site].num_active_neighbors!=num_neighbors) return true;
	else return false;
}

void System::update_neighbor_reproduction(int neighbor,int neighbor_state, int new_state){

	double sustract_amount,sum_amount;

	sustract_amount=mututalistic_term(neighbor,neighbor_state,new_state);
	add_num_neighbor_states(neighbor,new_state);
	list_sites[neighbor].num_active_neighbors++;
	sum_amount=mututalistic_term(neighbor,neighbor_state,new_state);
	list_sites[neighbor].mutualistic_wealth+=sum_amount-sustract_amount;
	list_sites[neighbor].total_wealth+=sum_amount-sustract_amount;
}
double System::mututalistic_term(int site, int site_state,int neighbor_state){

	int num_state=list_sites[site].num_neighbor_states[neighbor_state];
	int num_total=list_sites[site].num_active_neighbors;
	return Aij[site_state][neighbor_state]*(num_state/double(num_total));
}

void System::add_num_neighbor_states(int index, int istate) {list_sites[index].num_neighbor_states[istate]++;}
void System::remove_num_neighbor_states(int index, int istate) {list_sites[index].num_neighbor_states[istate]--;}

void System::swap(gsl_rng *r){

	int empty_site,neighbor;
	int site1,site2;
	int state1,state2;
	int irow,icolumn;
	int neighbor_state;
	
	if(gsl_rng_uniform(r)<p_swap){
	
		//cout<<num_active_sites<<endl;
		site1=list_active_front_sites[gsl_rng_uniform_int(r,num_active_sites)];
		site2=list_active_front_sites[gsl_rng_uniform_int(r,num_active_sites)];

		if(state1!=state2){
		 	
		 	//Change the states of the sites
		 	list_sites[site1].state=state2;
			list_sites[site2].state=state1;
			
			//Change neighbors of site1
			for(int i=0;i<num_neighbors;i++){
		
				neighbor=list_sites[site1].vneighbors[i];
				neighbor_state=list_sites[neighbor].state;
				if(neighbor_state>=0) update_neighbor_swap(neighbor,neighbor_state,state1,state2);
				else{
					add_num_neighbor_states(neighbor,state2);
					remove_num_neighbor_states(neighbor,state1);
				}
			}
			//Change neighbors of site2
			for(int i=0;i<num_neighbors;i++){
		
				neighbor=list_sites[site2].vneighbors[i];
				neighbor_state=list_sites[neighbor].state;
				if(neighbor_state>=0) update_neighbor_swap(neighbor,neighbor_state,state2,state1);
				else{
					add_num_neighbor_states(neighbor,state1);
					remove_num_neighbor_states(neighbor,state2);
				}
			}
		}
	}
}

void System::update_neighbor_swap(int site, int site_state, int old_neighbor_state, int new_neighbor_state){

	double sustract_amount=0,sum_amount=0;
	
	sustract_amount+=mututalistic_term(site,site_state,new_neighbor_state);
	sustract_amount+=mututalistic_term(site,site_state,old_neighbor_state);
	add_num_neighbor_states(site,new_neighbor_state);
	remove_num_neighbor_states(site,old_neighbor_state);
	sum_amount+=mututalistic_term(site,site_state,new_neighbor_state);
	sum_amount+=mututalistic_term(site,site_state,old_neighbor_state);
	list_sites[site].mutualistic_wealth+=sum_amount-sustract_amount;
	list_sites[site].total_wealth+=sum_amount-sustract_amount;
}

void System::add_empty_site(int index){

	position_empty_list[index]=list_empty_front_sites.size();
	list_empty_front_sites.push_back(index);
}

void System::add_active_site(int index){

	position_active_list[index]=list_active_front_sites.size();
	list_active_front_sites.push_back(index);
}

void System::remove_empty_site(int site){

	int last_index=list_empty_front_sites.size()-1;
	int last_element=list_empty_front_sites[last_index];
	int current_index=position_empty_list[site];
	
	list_empty_front_sites[current_index]=last_element;
	position_empty_list[last_element]=position_empty_list[site];
	position_empty_list[site]=-1;
	
	list_empty_front_sites.pop_back();
}

void System::remove_active_site(int site){

	int last_index=list_active_front_sites.size()-1;
	int last_element=list_active_front_sites[last_index];
	int current_index=position_active_list[site];
	
	list_active_front_sites[current_index]=last_element;
	position_active_list[last_element]=position_active_list[site];
	position_active_list[site]=-1;
	
	list_active_front_sites.pop_back();
	
	num_active_sites=list_active_front_sites.size();
}

void System::clean_variables(){

	list_empty_front_sites.clear();
	list_active_front_sites.clear();
}

//-------------------- Class simulation

class Simulation{

	protected:
		const int tfin;
		const int trials;
		const int refresh_time;
		const int rand_seed;
		
		System system;
		int *points;
		
		int time_increment;
		time_t itime, ftime;
		fstream FILE;
		const char* CURRENT_FILE_NAME;
		
		Simulation(Simulation&);
		void operator=(Simulation&);
		
	public:
		const int num_points;
		Simulation(System system_in, int T, int itrials, int irefresh_time,int iseed, char const* FILE_NAME): 
		system(system_in),tfin(T),trials(itrials),refresh_time(irefresh_time),rand_seed(iseed),num_points(int(log10(double(tfin))*2e2)+2),CURRENT_FILE_NAME(FILE_NAME){
		}; 
		
		//General simulation
		void sequential_trial(gsl_rng *r);
		void simulate();
		
		//Determined at the concret simulation
		virtual void initialice_time_variables(gsl_rng *r){};
		virtual void initialice_time_measures() {};
		virtual void take_initial_measures() {};
		virtual void take_time_measures(int ipoint) {};
		virtual bool conditions_run(double t){ return false; };
		virtual double get_time_increment() { return 1.; };
		virtual void sequential_time_step(gsl_rng *r){};
		virtual void initiate_counter(){};
		virtual void initialice_global_variables();
		virtual void initialice_global_measures() {};
		virtual void write(int trial);
		virtual void write_measures() {};
};

void Simulation::initialice_global_variables(){ 

	system.set_neighbors();
}

void Simulation::simulate(){

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
	initialice_global_measures();
	itime=time(NULL);
	
	for(int itrial=0;itrial<trials;itrial++){
	
		sequential_trial(r);
		write(itrial);
	}
	
	gsl_rng_free (r);
}

void Simulation::sequential_trial(gsl_rng *r){

	double t=1;
	int ipoint=0;

	initialice_time_variables(r);
	initialice_time_measures();
	take_initial_measures();
	while(conditions_run(t)){
	
		t+=get_time_increment();
		ipoint=int(log10(t)*2e2)+1;
		sequential_time_step(r);
		take_time_measures(ipoint);
		initiate_counter();
	}
	system.clean_variables();
}

void Simulation::write(int itrial){

	int tlogaux,tlog;
	
	ftime=time(NULL);
	if(difftime(ftime, itime)>refresh_time){
		
		FILE.open(CURRENT_FILE_NAME, ios::out | ios::trunc);
		tlogaux=0.;
		for(int i=0;i<num_points;i++){
			tlog=int(pow(10,((i-1)/(2e2))));
			if((points[i]>1)&&(tlog!=tlogaux)){ 
				FILE<<tlog<<" ";
				write_measures();
				FILE<<itrial<<endl;
			}
			tlogaux=tlog;
		}
		FILE.close();
		itime=time(NULL);
	}
}

//-------Visualization trial

class Visualize_evolution: public Simulation{

	public:
		virtual void initialice_time_variables(gsl_rng *r);
		virtual bool conditions_run(double t);
		virtual double get_time_increment();
		virtual void sequential_time_step(gsl_rng *r);	
		virtual void write(int itrial);
		
		Visualize_evolution(System system_in, int T, int itrials, int irefresh_time,int iseed, char const* FILE_NAME)
		: Simulation(system_in,T,itrials,irefresh_time,iseed,FILE_NAME){};
	
};

void Visualize_evolution::initialice_time_variables(gsl_rng *r){

	system.set_states(r);
	system.initialize_list_empty_sites();
	system.initialize_list_active_sites();
	system.initialize_num_neighbor_states();
	system.initialize_wealths();
	
}
void Visualize_evolution::sequential_time_step(gsl_rng *r){

	int empty_site,reproduce_site;
	
	empty_site=system.get_empty_site(r);
	reproduce_site=system.get_reproduce_site(r,empty_site);
	system.reproduce(empty_site,reproduce_site);
	system.swap(r);
	
}

bool Visualize_evolution::conditions_run(double t){

	if((t<tfin)&&(system.num_empty_sites>0)) return true;
	else return false;
}

double Visualize_evolution::get_time_increment(){ return 1./double(system.num_empty_sites); }

void Visualize_evolution::write(int itrial){ //Solo vale para 2!!!

	int site;
	int linear_size=system.get_linear_size();
	int dimension=system.get_dimension();
	FILE.open(CURRENT_FILE_NAME, ios::out | ios::trunc);
	
	if(dimension==2){
		for(int isite=0;isite<linear_size;isite++){
			for(int jsite=0;jsite<linear_size;jsite++){
						
				site=isite*linear_size+jsite;
				FILE<<isite<<" "<<jsite<<" "<<system.list_sites[site].state<<endl;
			}
		}
	}
	
	FILE.close();
}

#endif
