#include <vector>
#include <string>
#include <io.h> 
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <map>
#include <algorithm>
#include <time.h>
#include <strstream>

//#include <winbase.h>

using namespace std;


//R7018_30_MCRT_SILO_55_10
//G_342_MCRT_SILO_55_10
typedef struct{
	int _t; //topology
	int _d; //descriptor, MLU for Rocketfuel, TM index for Abilene/GEANT/CERNET
	int _c; //core graph, MCRT/2CG
	int _s;//scheduling method: SILO/LISO/RS
	int _u;//upper bound utilization
	int _f;//stability shift
	int _p;//probabilistic onoff
	string _str;
} TEST_DESCRIPTOR;

typedef struct{
	double _min;
	double _max;
	int _slot_number;
	vector<double> _cdf;
	vector<double> _cdf_ospf;
	vector<double> _xaxis;
}CDF_DELAY;

typedef struct{
	int _cg; //core graph: MCRT or CG
	int _sm; //scheduling method: SILO/LISO/RS
	int _ub; //upper bound
	int _ss; //stability shift
	int _prb;//probabilistic switching paths, default no
	double _mlu;//current TM mlu
}PARAMETER_SET;

#define STRING_VECTOR vector<string>

typedef vector<map<pair<int,int>,vector<int>>>  LINK_MAP; 
typedef map<pair<int,int>,double>				LINK_BUDGET;
typedef vector<LINK_BUDGET>						TH_MAP;
typedef vector<vector<double>>					FLOAT_MATRIX;
typedef vector<vector<int>>						INT_MATRIX;
typedef vector<vector<bool>>					BOOL_MATRIX;
typedef vector<double>							VEC_FLOAT;

struct UTILITY_INFO{
	double max_utility;
	double min_utility;
	double avg_utility;
	int max_u,max_v;			//MUL link index
};

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
	struct UTILITY_INFO _util_info_ospf;
	double _avg_time;
	double _lnk_energy;
	double _average_stretch;
	double _average_stretch_hop;
}SIMULATOR_OUTPUT;

#define NO_PR		1000			//no probability

#define ABILENE (-1)
#define GEANT (-2)
#define CERNET (-3)
#define MCRT 0
#define CG 1
#define SILO 0
#define LISO 1
#define RS 2
#define NO_PROB 0
#define PROB 1

#define UNCHANGED (-100)

#define ROUND 10

#define DELAY_SLOT_NUMBER 100

//current folder: simulation\\simulator\\simulator\\
//Abilene topo file: simulation\\Abilene\\Abilene_topo.txt

#define ABILENE_TOPOLOGY_PATH "..\\..\\Abilene\\Abilene_topo.txt"
#define GEANT_TOPOLOGY_PATH "..\\..\\GEANT\\GEANT_topo.txt"
#define CERNET_TOPOLOGY_PATH "..\\..\\CERNET_TM\\CERNET_topo.txt"
#define ROCKET_FUEL_TOPO_DIR "..\\..\\topo\\"
#define ROCKET_FUEL_TM_DIR "..\\..\\traffic\\"
#define ROUTING_DIR "..\\..\\simulator_output\\routing\\"
#define CORELINKS_DIR "..\\..\\simulator_output\\corelinks\\"

#define ABILENE_TM_PATH_TEMPLATE "..\\..\\Abilene\\Abilene_tm_" //Abilene_tm_0.txt Abilene_tm_1.txt
#define GEANT_TM_PATH_TEMPLATE "..\\..\\GEANT\\geant_tm_" //geant_tm_1.txt geant_tm_2.txt
#define CERNET_TM_PATH_TEMPLATE "..\\..\\CERNET_TM\\C_"//C_0.txt

#define SIMULATION_OUTPUT_DIRECTORY "..\\..\\simulator_output\\"



#define LOG_FILE "..\\..\\log\\log_file.txt"
#define LOG_MST_UTIL_FILE "..\\..\\log\\utilization_mst.txt"
#define LOG_ORDINARY_UTIL_FILE "..\\..\\log\\utilization_ordinary.txt"
#define FAIL_LIST "..\\..\\log\\fail_list.txt"
#define LOG_PREV_PATH_FILE "..\\..\\log\\prev_path.txt"
void log_report(string out);
void fail_report();

#define LINK_NOT_EXIST -1
#define LINK_ON 1
#define LINK_OFF 0

//#define CORE_GRAPH_ENERGY_AWARE
//#define DOUBLE_ROOT_TREE

#define PROB_FACTOR 1

//#define TM_SP_TH_PERCENTAGE 10
//#define TM_SP_TRESHOLD

#define TH_CHANGE //when mlu>50, then threshold is set as the mlu

#define OVERLAP_AWARE

bool RunSimulation(FLOAT_MATRIX& Topology, FLOAT_MATRIX& Capacity, 
				   FLOAT_MATRIX& Power, FLOAT_MATRIX& TM, FLOAT_MATRIX& Latency,
				   PARAMETER_SET& engine_parameters,SIMULATOR_OUTPUT& output);