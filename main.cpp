#include "routingengine.h"
using namespace std;

vector<TEST_DESCRIPTOR> test_vector;
int current_i;
int node_number, link_number;
double total_link_energy, total_node_energy;
int dumplinktime = 0;
vector<int> mapping;
STRING_VECTOR node_name;
FLOAT_MATRIX Topology, Capacity, Consumption, TM, Latency;
FLOAT_MATRIX Link_load;
BOOL_MATRIX Link_state;
PARAMETER_SET engine_parameters;
SIMULATOR_OUTPUT simulator_output;

string rocket_fuel_topo_directory = ROCKET_FUEL_TOPO_DIR;
string rocket_fuel_traffic_directory = ROCKET_FUEL_TM_DIR;

void Load(vector<TEST_DESCRIPTOR>& test_vector, int index);
void DumpOutput(int test_descriptor_index);
TEST_DESCRIPTOR string2TEST_DESCRIPTOR(string str);
//void LoadParameter(int parameter);
//void LoadMethod(int method);
void LoadTopology(int topo);
void read_topology(int topo, string filename);
void initialize();
void LoadTM(int topo, int descriptor);
void read_tm(string filename, int scale=1);
void rocketfuel_tm_read(int topo, int scale);
void Abilene_tm_read(int index);
void GEANT_tm_read(int index);
void CERNET_tm_read(int index);
void DumpMatrix(ofstream& fout, FLOAT_MATRIX& a);
void DumpMatrix(ofstream& fout, BOOL_MATRIX& a);
void DumpMatrix(ofstream& fout, INT_MATRIX& a);
void DumpVector(ofstream& fout, vector<double>& a);

#define REAL_RUN
int main(int argc, char* argv[]){
#ifdef REAL_RUN	
	ifstream test_descriptor_list(argv[1]);
#else	
	ifstream test_descriptor_list("..//..//test_list.txt");
#endif
	while(!test_descriptor_list.eof()){
		string test_descriptor;
		getline(test_descriptor_list,test_descriptor);
		test_vector.push_back(string2TEST_DESCRIPTOR(test_descriptor));
	}
	for(int i=0; i<test_vector.size(); i++){
		current_i = i;
		Load(test_vector, i);
		bool sim_result = RunSimulation(Topology,Capacity,Consumption,TM,Latency,engine_parameters,simulator_output);
		if(sim_result == false){
			fail_report();
			log_report("Simulation fais!");
			continue;
		}
		else{
			DumpOutput(i);
		}
	}
	//system("pause");
}



void DumpOutput(int test_descriptor_index){
/*
typedef struct{
	FLOAT_MATRIX _util;					
	FLOAT_MATRIX _sp_delay;
	FLOAT_MATRIX _cg_delay;
	FLOAT_MATRIX _sp_hop;
	FLOAT_MATRIX _cg_hop;
	INT_MATRIX _ls;
	BOOL_MATRIX _path;
	CDF_DELAY _cdf_delay;
	CDF_DELAY _cdf_hop;
	struct UTILITY_INFO _util_info;
	double _avg_time;
	double _lnk_energy;
}SIMULATOR_OUTPUT;
*/
	//different value in different files
	string output_directory = SIMULATION_OUTPUT_DIRECTORY;
	double total_energy = total_link_energy + total_node_energy;
	ofstream fout((output_directory + "\\basic\\" + test_vector[test_descriptor_index]._str+".txt").c_str());
	fout<<"test descriptor: "<<test_vector[test_descriptor_index]._str<<endl;
	fout<<"node number: "<<node_number<<endl;
	fout<<"link number: "<<link_number<<endl;
	fout<<"total energy: "<<total_energy<<endl;
	fout<<"total link energy: "<<total_link_energy<<endl;
	fout<<"total node energy: "<<total_node_energy<<endl;
	fout<<"used link energy: "<<simulator_output._lnk_energy<<endl;
	fout<<"energy saving ratio: "<<1-(simulator_output._lnk_energy + total_node_energy)/total_energy<<endl;
	fout<<"link energy saving ratio: "<<1-simulator_output._lnk_energy/total_link_energy;
	fout.close();

	cout<<"energy saving ratio: "<<1-(simulator_output._lnk_energy + total_node_energy)/total_energy<<endl;

	fout.open((output_directory + "\\basic_info\\" + test_vector[test_descriptor_index]._str+".txt").c_str());
	//fout<<test_vector[test_descriptor_index]._str<<" ";
	fout<<0<<" ";
	fout<<node_number<<" ";
	fout<<link_number<<" ";
	fout<<total_energy<<" ";
	fout<<total_link_energy<<" ";
	fout<<total_node_energy<<" ";
	fout<<simulator_output._lnk_energy<<" ";
	fout<<1-(simulator_output._lnk_energy + total_node_energy)/total_energy<<" ";
	fout<<1-simulator_output._lnk_energy/total_link_energy;
	fout.close();

	
	fout.open((output_directory + "\\util\\" + test_vector[test_descriptor_index]._str+".txt").c_str());
	DumpMatrix(fout,simulator_output._util);
	fout.close();

	fout.open((output_directory + "\\sp_delay\\" + test_vector[test_descriptor_index]._str+".txt").c_str());
	DumpMatrix(fout,simulator_output._sp_delay);
	fout.close();

	fout.open((output_directory + "\\cg_delay\\" + test_vector[test_descriptor_index]._str+".txt").c_str());
	DumpMatrix(fout,simulator_output._cg_delay);
	fout.close();

	fout.open((output_directory + "\\sp_hop\\" + test_vector[test_descriptor_index]._str+".txt").c_str());
	DumpMatrix(fout,simulator_output._sp_hop);
	fout.close();

	fout.open((output_directory + "\\cg_hop\\" + test_vector[test_descriptor_index]._str+".txt").c_str());
	DumpMatrix(fout,simulator_output._cg_hop);
	fout.close();

	fout.open((output_directory + "\\link_state\\" + test_vector[test_descriptor_index]._str+".txt").c_str());
	DumpMatrix(fout,simulator_output._ls);
	fout.close();

	fout.open((output_directory + "\\path\\" + test_vector[test_descriptor_index]._str+".txt").c_str());
	DumpMatrix(fout,simulator_output._path);
	fout.close();

	fout.open((output_directory + "\\cdf_delay\\" + test_vector[test_descriptor_index]._str+".txt").c_str());
	DumpVector(fout,simulator_output._cdf_delay._xaxis);
	DumpVector(fout,simulator_output._cdf_delay._cdf);
	DumpVector(fout,simulator_output._cdf_delay._cdf_ospf);
	fout.close();

	fout.open((output_directory + "\\cdf_hop\\" + test_vector[test_descriptor_index]._str+".txt").c_str());
	DumpVector(fout,simulator_output._cdf_hop._xaxis);
	DumpVector(fout,simulator_output._cdf_hop._cdf);
	DumpVector(fout,simulator_output._cdf_hop._cdf_ospf);
	fout.close();

	fout.open((output_directory + "\\info\\" + test_vector[test_descriptor_index]._str+".txt").c_str());
	fout<<simulator_output._util_info.avg_utility<<" "<<simulator_output._util_info.min_utility<<" ";
	fout<<simulator_output._util_info.max_utility<<" ";
	fout<<simulator_output._util_info_ospf.avg_utility<<" "<<simulator_output._util_info_ospf.min_utility<<" ";
	fout<<simulator_output._util_info_ospf.max_utility<<" ";
	fout<<1-(simulator_output._lnk_energy + total_node_energy)/total_energy<<" ";
	fout<<simulator_output._avg_time<<" ";
	fout<<simulator_output._average_stretch<<" ";
	fout<<simulator_output._average_stretch_hop<<" ";

	if(simulator_output._util_info.max_utility>1)log_report("PALS overflow!");
	if(simulator_output._util_info_ospf.max_utility>1)log_report("ospf overflow!");
	
	fout.close();
}

void DumpVector(ofstream& fout, vector<double>& a){
		for(int j=0; j<a.size(); j++){
			fout<<a[j]<<" ";
		}
		fout<<endl;
}

void DumpMatrix(ofstream& fout, FLOAT_MATRIX& a){
	for(int i=0; i<a.size(); i++){
		for(int j=0; j<a.size(); j++){
			fout<<a[i][j]<<" ";
		}
		fout<<endl;
	}
}
void DumpMatrix(ofstream& fout, BOOL_MATRIX& a){
	for(int i=0; i<a.size(); i++){
		for(int j=0; j<a.size(); j++){
			fout<<a[i][j]<<" ";
		}
		fout<<endl;
	}
}
void DumpMatrix(ofstream& fout, INT_MATRIX& a){
	for(int i=0; i<a.size(); i++){
		for(int j=0; j<a.size(); j++){
			fout<<a[i][j]<<" ";
		}
		fout<<endl;
	}
}

//R7018_30_MCRT_SILO_55_10
//G_342_MCRT_SILO_55_10
TEST_DESCRIPTOR string2TEST_DESCRIPTOR(string str){
	//PALS_A_16_0
	TEST_DESCRIPTOR re;
	re._str = str;

	int index = str.find('_');
	string topo = str.substr(0,index);
	if(topo == "A") re._t = ABILENE;
	else if(topo == "G") re._t = GEANT;
	else if(topo == "C") re._t = CERNET;
	else{re._t = atoi(topo.substr(1,topo.length()-1).c_str());}// rocketfuel topology
	
	str = str.substr(index+1,str.length()-index-1);
	index = str.find('_');
	string descriptor = str.substr(0, index);
	re._d = atoi(descriptor.c_str());
	
	str = str.substr(index+1,str.length()-index-1);
	index = str.find('_');
	string core_graph = str.substr(0, index);
	if(core_graph == "MCRT") re._c = MCRT;
	else if(core_graph == "CG") re._c = CG;

	str = str.substr(index+1,str.length()-index-1);
	index = str.find('_');
	string scheduling = str.substr(0, index);
	if(scheduling == "SILO") re._s = SILO;
	else if(scheduling == "LISO") re._s = LISO;
	else if(scheduling == "RS") re._s = RS;

	str = str.substr(index+1,str.length()-index-1);
	index = str.find('_');
	string upper_bound = str.substr(0, index);
	re._u = atoi(upper_bound.c_str());

	str = str.substr(index+1,str.length()-index-1);
	index = str.find('_');
	string stability_shift = str.substr(0, index);
	re._f = atoi(stability_shift.c_str());


	if(index != -1){
		str = str.substr(index+1,str.length()-index-1);
		index = str.find('_');
		string probabilistic_onoff = str.substr(0, index);
		if(probabilistic_onoff == "PROB")
			re._p = PROB;
		else if(probabilistic_onoff == "NOPROB")
			re._p = NO_PROB;
	}
	else{
		re._p = NO_PROB;// default no probabilistic
	}

	return re;
}

//R7018_30_MCRT_SILO_55_10
//G_342_MCRT_SILO_55_10
void Load(vector<TEST_DESCRIPTOR>& test_vector, int index){
	TEST_DESCRIPTOR current, last;
	current = test_vector[index];
	if(index!=0) last = test_vector[index-1];
	if(index==0 || current._t != last._t){
		LoadTopology(current._t);
		LoadTM(current._t,current._d);
	}
	if(index==0 || current._d != last._d)
		LoadTM(current._t,current._d);

	engine_parameters._cg = current._c;
	engine_parameters._sm = current._s;
	engine_parameters._ub = current._u;
	if(current._t != ABILENE && current._t != GEANT && current._t != CERNET){
		//for rocket fuel topologies, if MLU of input is larger than 50, then upper bound raise accordingly
		if(current._d>50)engine_parameters._ub = current._d;
#ifdef TH_CHANGE
		engine_parameters._mlu = current._d;
#endif
	}
	else{
		engine_parameters._mlu = -1; //inform that it's not Rocketfuel
	}
	

	engine_parameters._ss = current._f;
	engine_parameters._prb = current._p;

	dumplinktime = 0;
}

void LoadTopology(int topo){
	string path;
	if(topo == ABILENE) path = ABILENE_TOPOLOGY_PATH;
	else if(topo == GEANT) path = GEANT_TOPOLOGY_PATH;
	else if(topo == CERNET) path = CERNET_TOPOLOGY_PATH;
	else{
		strstream buf;buf<<topo;string as_num; buf>>as_num;
		path = rocket_fuel_topo_directory + "topo" + as_num + ".txt";
		//path = rocket_fuel_path[find(mapping.begin(), mapping.end(), topo)-mapping.begin()];
	}
	read_topology(topo, path);
}

void read_topology(int topo, string filename)
{
			ifstream fin(filename.c_str());
			
			fin>>node_number;
			initialize();
			for(int i=0; i<node_number; i++) {
				int index; fin>>index; mapping.push_back(index);
				string name; getline(fin, name);
                node_name.push_back(name);
			}
			
			//remap_with_name();

			fin>>link_number;
			for(int i=0; i<link_number; i++){
				int src, des;
				double weight, consumption, bandwidth;
				fin>>src>>des>>weight>>consumption>>bandwidth;
				src = find(mapping.begin(), mapping.end(), src) - mapping.begin();
				des = find(mapping.begin(), mapping.end(), des) - mapping.begin();
				Topology[src][des] = weight; Consumption[src][des] = consumption;
				Capacity[src][des] = bandwidth;				
			}

			

				for(int i=0; i<link_number; i++){
					int src, des;
					double latency;
					fin>>src>>des>>latency;
					src = find(mapping.begin(), mapping.end(), src) - mapping.begin();
					des = find(mapping.begin(), mapping.end(), des) - mapping.begin();
					Latency[src][des] = latency;
				}
				fin>>total_link_energy>>total_node_energy;
				for(int i=0; i<node_number; i++){
					int node;
					double node_power;
					fin>>node>>node_power;
					//Latency[src][des] = latency;
				}
			
			//node name: node_name vector
			//node number: node_number
			//link number: link_number
			//latency: Topology matrix
			//capacity: Capacity matrix
			//link energy: Consumption matrix
			//Latency: link latency
			fin.close();
}

void initialize(){
	Topology.clear(); Capacity.clear(); Consumption.clear(); 
	mapping.clear(); node_name.clear();
	vector<double> empty_vector(node_number,-1);
	FLOAT_MATRIX initializer(node_number, empty_vector);
	Topology = Capacity = Consumption = Latency = initializer;
}

void LoadTM(int topo, int descriptor){
	if(topo == ABILENE) Abilene_tm_read(descriptor);
	else if(topo == GEANT) GEANT_tm_read(descriptor);
	else if(topo == CERNET) CERNET_tm_read(descriptor);
	else{
		rocketfuel_tm_read(topo, descriptor);
	}
}

void rocketfuel_tm_read(int topo, int scale){
	string path = rocket_fuel_traffic_directory;
	strstream buf; buf<<topo; string rtopo; buf>>rtopo;
	path += "gravity_model_traffic_" + rtopo + "_1.txt";
	read_tm(path, scale);
}

void read_tm(string filename, int scale){
	double sc = scale;
	if(scale == 100) sc = 99;//to avoid actually overflow under OSPF

			ifstream fin(filename.c_str());			
			fin>>node_number;
			TM.clear();
			vector<double> tmp(node_number,0);
			TM.insert(TM.begin(),node_number,tmp);

			for(int i=0; i<node_number; i++) {
				int index; fin>>index; 
				string name; getline(fin, name);
			}
			double traffic;
			for(int i=0; i<node_number; i++)
				for(int j=0; j<node_number; j++){
					fin>>traffic;
					if(i==j)continue;
					TM[i][j] = traffic*sc;//Kbps
				}
			fin.close();
}

void Abilene_tm_read(int index){
	string path = ABILENE_TM_PATH_TEMPLATE;
	strstream buf; buf<<index; string test_index; buf>>test_index;
	path += test_index + ".txt";
	read_tm(path);
}

void GEANT_tm_read(int index){
	string path = GEANT_TM_PATH_TEMPLATE;
	strstream buf; buf<<index; string test_index; buf>>test_index;
	path += test_index + ".txt";
	read_tm(path);
}

void CERNET_tm_read(int index){
	string path = CERNET_TM_PATH_TEMPLATE;
	strstream buf; buf<<index; string test_index; buf>>test_index;
	path += test_index + ".txt";
	read_tm(path);
}
