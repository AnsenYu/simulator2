#include "mcrt2mrsg.h"

using namespace std;


bool Generate_SPT(FLOAT_MATRIX& Topology, INT_MATRIX& SPT_previous);
void Dijkstra(FLOAT_MATRIX &Topology, int t, vector<int> &previous);
void MRCT_MST(FLOAT_MATRIX &Topology,FLOAT_MATRIX &Power, INT_MATRIX &SPT_previous, BOOL_MATRIX &MST);
int MRCT2(FLOAT_MATRIX &Tpology, FLOAT_MATRIX &Power, INT_MATRIX &SPT_previous, BOOL_MATRIX &MST);
void Link_Pass(LINK_MAP &MST_map, LINK_MAP &SPT_map, INT_MATRIX &MST_pass);
void MST_Statistics(vector<vector<int>> &MST_all_previous,vector<map<pair<int,int>,vector<int>>> &MST_map);
void SP_ST_Overlap(INT_MATRIX &SPT_all_previous, INT_MATRIX &MST_all_previous, BOOL_MATRIX &bOverlap);
void Path_Delay(FLOAT_MATRIX &Tpology, INT_MATRIX &SPT_previous, INT_MATRIX &MST_previous, 
				FLOAT_MATRIX &SPT_DELAY, FLOAT_MATRIX &MST_DELAY,double& max, double& min, bool hop);
void Prev_Utilization(FLOAT_MATRIX &TM, FLOAT_MATRIX &Capacity, BOOL_MATRIX &prev_path, 
					  INT_MATRIX &SPT_all_previous, INT_MATRIX &MST_all_previous,
					  FLOAT_MATRIX &Utilization, UTILITY_INFO &utility_info);
void Budget_Adjust(FLOAT_MATRIX &Utilization, FLOAT_MATRIX &TM, FLOAT_MATRIX &Capacity,
				   LINK_MAP &MST_map, LINK_MAP &SPT_map,INT_MATRIX &MST_pass,
				   BOOL_MATRIX &MST, BOOL_MATRIX &prev_path,
				   TH_MAP &dFlow_th);
void Scheduling(INT_MATRIX &SPT_all_previous, INT_MATRIX &MST_all_previous,
						 BOOL_MATRIX &MST, LINK_MAP &SPT_MAP, LINK_MAP &MST_MAP,
						 FLOAT_MATRIX &Capacity, FLOAT_MATRIX &TM,
						 BOOL_MATRIX &prev_path, BOOL_MATRIX &bOverlap, 
						 TH_MAP &dFlow_th);
void GetLinkStateAndTotalEnergy(FLOAT_MATRIX& Topology, BOOL_MATRIX& MST, 
	              FLOAT_MATRIX& Utilization, INT_MATRIX& link_state, FLOAT_MATRIX& Consumption,
				  double& total_energy);
void GetPacketDelayHop(FLOAT_MATRIX& SPT_DELAY, FLOAT_MATRIX& MST_DELAY,
	                BOOL_MATRIX& prev_path, FLOAT_MATRIX& TM, 
					CDF_DELAY& weighted_delay, double& avg_stretch);
void GetPacketDelay(FLOAT_MATRIX& SPT_DELAY, FLOAT_MATRIX& MST_DELAY,
	                BOOL_MATRIX& prev_path, FLOAT_MATRIX& TM, 
					CDF_DELAY& weighted_delay, double& avg_stretch);
void Install_engine_parameters(PARAMETER_SET& engine_parameters);
void DumpLinkUtilization(FLOAT_MATRIX& Topology, FLOAT_MATRIX& Util, BOOL_MATRIX& MST, bool MST_TRUE);
void DumpPrev_path(FLOAT_MATRIX& Topology, BOOL_MATRIX& prev_path);
void DumpRouting(FLOAT_MATRIX Topology, INT_MATRIX& spt_routing,INT_MATRIX& mst_routing);
void DumpCoreLinks(BOOL_MATRIX MST);

//default setting
double UTILITY_LB = 0.45;
double UTILITY_UB = 0.55;
int CORE_GRAPH = MCRT;
int SCHEDULING = SILO;
int PROBABILISTIC = NO_PROB;
double TM_sp_threshold;

bool RunSimulation(FLOAT_MATRIX& Topology, FLOAT_MATRIX& Capacity, 
				   FLOAT_MATRIX& Power, FLOAT_MATRIX& TM, FLOAT_MATRIX& Latency,
				   PARAMETER_SET& engine_parameters, SIMULATOR_OUTPUT& output)
{
	Install_engine_parameters(engine_parameters);

	INT_MATRIX SPT_previous;
	bool connected = Generate_SPT(Topology,SPT_previous);//shortest path first routing
	if(connected == false) return false;

	BOOL_MATRIX MST;
	switch(CORE_GRAPH){
		case MCRT:
		MRCT_MST(Topology,Power,SPT_previous,MST);//generate core graph
		break;
		case CG:
		int root = MRCT2(Topology,Power,SPT_previous,MST);   //NewMethod: T-ReG(Tree-based Reliable Graph)
#ifndef DOUBLE_ROOT_TREE
		MRSG_from_MCRT(Topology, MST, root,engine_parameters._mlu);
#endif
		break;
	}

	FLOAT_MATRIX subTpology = Topology;
	int N = Topology.size();
	for(int k=0;k!=N;k++)
		for(int j=0;j!=N;j++){
		    if(MST[k][j]==false && subTpology[k][j]!=-1)
				subTpology[k][j]=-1;
		}
	INT_MATRIX MST_previous;
    Generate_SPT(subTpology,MST_previous);//shortest path routing on MCRT

	DumpRouting(Topology, SPT_previous, MST_previous);
	DumpCoreLinks(MST);

	vector<double> TM_vector;
	for(int k=0;k!=N;k++)
		for(int j=0;j!=N;j++){
		    if(Topology[k][j]!=-1)
				TM_vector.push_back(TM[k][j]);
	}
	sort(TM_vector.begin(),TM_vector.end());

//	TM_sp_threshold = TM_vector.at(TM_vector.size()-TM_vector.size()/TM_SP_TH_PERCENTAGE);


	LINK_MAP MST_map;			//ʵ����ֵ�����㷨��ӳ��
	LINK_MAP SPT_map;
	INT_MATRIX MST_pass;		//ÿ�������Ͼ������flow����
	BOOL_MATRIX bOverlap;
	FLOAT_MATRIX SPT_DELAY, MST_DELAY, SPT_HOP, MST_HOP;//<=======output
	MST_Statistics(MST_previous,MST_map);
	MST_Statistics(SPT_previous,SPT_map);
	Link_Pass(MST_map,SPT_map,MST_pass);
	SP_ST_Overlap(SPT_previous,MST_previous,bOverlap);			//ȷ����ЩOD pair��SP��ST����ȫ�غϵ�
	
	double delay_max, delay_min, delay_max_hop, delay_min_hop;
	Path_Delay(Latency,SPT_previous,MST_previous,SPT_DELAY,MST_DELAY,delay_max,delay_min,false);
	Path_Delay(Latency,SPT_previous,MST_previous,SPT_HOP,MST_HOP,delay_max_hop,delay_min_hop,true);

	FLOAT_MATRIX Utilization;
	vector<bool>  tmp_bool(N,true);
	BOOL_MATRIX prev_path(N,tmp_bool);

	UTILITY_INFO es_utility, es_utility_ospf;
	TH_MAP dFlow_th;
	for(int k=0;k!=N;k++){
		map<pair<int,int>,double> tmp_th; //��ÿ��MST link�趨һ����ֵ
		for(int u=0;u!=N;u++)
			for(int v=0;v!=N;v++){
				if(u!=v && MST[u].at(v)){
					pair<int,int> index(u,v);
					bool dirST = MST_map[k].find(index)!=MST_map[k].end();
					if(!dirST) continue;
					tmp_th[index] = 0;
				}
			}
		dFlow_th.push_back(tmp_th);
	}

	clock_t start, finish; long total_time = 0;double average_time;
	for(int k=0;k!=ROUND;k++)
	{
		Utilization.clear();
		
		//��һ��round�ս���ʱ��״̬
		/*void Prev_Utilization(FLOAT_MATRIX &TM, FLOAT_MATRIX &Capacity, BOOL_MATRIX &prev_path, 
					  INT_MATRIX &SPT_all_previous, INT_MATRIX &MST_all_previous,
					  FLOAT_MATRIX &Utilization, UTILITY_INFO &utility_info)*/
		Prev_Utilization(TM,Capacity,prev_path,SPT_previous,MST_previous,Utilization,es_utility);
		DumpLinkUtilization(Topology,Utilization,MST,true);
		DumpLinkUtilization(Topology,Utilization,MST,false);
		DumpPrev_path(Topology,prev_path);

		//cout<<Utilization[17][25]<<endl;

		if(k==0)es_utility_ospf = es_utility;//record the utility under OSPF

		start = clock();
		//Budget����
		Budget_Adjust(Utilization,TM,Capacity,
			MST_map,SPT_map,MST_pass,MST,prev_path,dFlow_th);

		/*void SILO(INT_MATRIX &SPT_all_previous, INT_MATRIX &MST_all_previous,
						 BOOL_MATRIX &MST, LINK_MAP &SPT_MAP, LINK_MAP &MST_MAP,
						 FLOAT_MATRIX &Capacity, FLOAT_MATRIX &TM,
						 BOOL_MATRIX &prev_path, BOOL_MATRIX &bOverlap, 
						 TH_MAP &dFlow_th)*/
        Scheduling(SPT_previous,MST_previous,MST,SPT_map,MST_map,
            Capacity,TM,prev_path,bOverlap,dFlow_th);

		finish = clock();
		total_time += finish - start;
		//system("pause");
	}
	Utilization.clear();
	Prev_Utilization(TM,Capacity,
		             prev_path,  //<=================output
		             SPT_previous,MST_previous,
					 Utilization,//<=================output
					 es_utility  //<=================output
					 );
	average_time = ((double)total_time)/CLOCKS_PER_SEC/ROUND/N; //<=======output: average for each node each time
	cout<<average_time<<endl;
	double total_link_energy = 0;//<=================output 
	vector<int>  tmp_int(N,LINK_OFF);
	INT_MATRIX link_state(N,tmp_int);//<=================output
	GetLinkStateAndTotalEnergy(Topology, MST, Utilization, link_state,Power,total_link_energy);

	CDF_DELAY weighted_delay_ms, weighted_delay_hop;//<=================output
	double weighted_stretch, weighted_stretch_hop; //<===============output
	weighted_delay_ms._max = delay_max; weighted_delay_ms._min = delay_min; 
	weighted_delay_ms._slot_number = DELAY_SLOT_NUMBER;
	weighted_delay_hop._max = delay_max_hop; weighted_delay_hop._min = delay_min_hop; 
	weighted_delay_hop._slot_number = (int)(delay_max_hop - delay_min_hop);

	GetPacketDelay(SPT_DELAY, MST_DELAY, prev_path, TM, weighted_delay_ms, weighted_stretch);
	GetPacketDelayHop(SPT_HOP, MST_HOP, prev_path, TM, weighted_delay_hop, weighted_stretch_hop);


	//dump output varialbes
	output._avg_time = average_time;
	output._cdf_delay = weighted_delay_ms;
	output._cdf_hop = weighted_delay_hop;
	output._cg_delay = MST_DELAY;
	output._cg_hop = MST_HOP;
	output._lnk_energy = total_link_energy;
	output._ls = link_state;
	output._path = prev_path;
	output._sp_delay = SPT_DELAY;
	output._sp_hop = SPT_HOP;
	output._util = Utilization;
	output._util_info = es_utility;
	output._util_info_ospf = es_utility_ospf;
	output._average_stretch = weighted_stretch;
	output._average_stretch_hop = weighted_stretch_hop;
	return true;
}

void Install_engine_parameters(PARAMETER_SET& engine_parameters){
	UTILITY_UB = ((double)engine_parameters._ub)/100;
	UTILITY_LB = ((double)engine_parameters._ub - engine_parameters._ss)/100;
	CORE_GRAPH = engine_parameters._cg;
	SCHEDULING = engine_parameters._sm;
	PROBABILISTIC = engine_parameters._prb;
}

void GetPacketDelay(FLOAT_MATRIX& SPT_DELAY, FLOAT_MATRIX& MST_DELAY,
	                BOOL_MATRIX& prev_path, FLOAT_MATRIX& TM, 
					CDF_DELAY& weighted_delay, double& avg_stretch)
{
// raw data is kept in output file
	int N = SPT_DELAY.size();
	double max = weighted_delay._max+1;
	double min = weighted_delay._min;
	int slot_number = weighted_delay._slot_number;
	double slot_width = (max - min)/slot_number;
	vector<double> slot_traffic(slot_number,0);
	vector<double> ospf_slot_traffic(slot_number,0);
	double total_traffic = 0;
	double stretch_product = 0;

	for(int i=0; i!=N; i++){
		for(int j=0; j!=N; j++){
			if(i==j)continue;
			total_traffic += TM[i][j];
			ospf_slot_traffic[floor((SPT_DELAY[i][j] - min)/slot_width)]+=TM[i][j];
			if(prev_path[i][j] == true){
				slot_traffic[floor((SPT_DELAY[i][j] - min)/slot_width)]+=TM[i][j];
				stretch_product += TM[i][j];
			}
			else{
				slot_traffic[floor((MST_DELAY[i][j] - min)/slot_width)]+=TM[i][j];
				stretch_product += MST_DELAY[i][j]/SPT_DELAY[i][j] * TM[i][j];
			}
		}
	}

	avg_stretch = stretch_product / total_traffic;

	for(int i=1; i!=slot_number; i++){
		slot_traffic[i]+=slot_traffic[i-1]; //cumulative sum
		ospf_slot_traffic[i]+=ospf_slot_traffic[i-1]; //cumulative sum
	}
	weighted_delay._cdf = slot_traffic;
	weighted_delay._cdf_ospf = ospf_slot_traffic;

	weighted_delay._xaxis.push_back(min); 
	for(int i=1; i!=slot_number; i++){
		weighted_delay._xaxis.push_back(min+slot_width*i); 
	}
}

void GetPacketDelayHop(FLOAT_MATRIX& SPT_DELAY, FLOAT_MATRIX& MST_DELAY,
	                BOOL_MATRIX& prev_path, FLOAT_MATRIX& TM, 
					CDF_DELAY& weighted_delay, double& avg_stretch)
{	
	int N = SPT_DELAY.size();
	int max = weighted_delay._max;
	int min = weighted_delay._min;
	int slot_number = max - min + 1;
	weighted_delay._slot_number = slot_number;
	vector<double> slot_traffic(slot_number,0);
	vector<double> ospf_slot_traffic(slot_number,0);
	double total_traffic = 0;
	double stretch_product = 0;

	for(int i=0; i!=N; i++){
		for(int j=0; j!=N; j++){
			if(i==j)continue;
			total_traffic += TM[i][j];
			ospf_slot_traffic[SPT_DELAY[i][j] - min]+=TM[i][j];
			if(prev_path[i][j] == true){
				slot_traffic[SPT_DELAY[i][j] - min]+=TM[i][j];
				stretch_product += TM[i][j];
			}
			else{
				slot_traffic[MST_DELAY[i][j] - min]+=TM[i][j];
				stretch_product += MST_DELAY[i][j]/SPT_DELAY[i][j] * TM[i][j];
			}
		}
	}

	avg_stretch = stretch_product / total_traffic;

	for(int i=1; i!=slot_number; i++){
		slot_traffic[i]+=slot_traffic[i-1];
		ospf_slot_traffic[i]+=ospf_slot_traffic[i-1]; //cumulative sum
	}
	weighted_delay._cdf = slot_traffic;
	weighted_delay._cdf_ospf = ospf_slot_traffic;

	weighted_delay._xaxis.push_back(min); 
	for(int i=1; i!=slot_number; i++){
		weighted_delay._xaxis.push_back(min+i); 
	}
}

void GetLinkStateAndTotalEnergy(FLOAT_MATRIX& Topology, BOOL_MATRIX& MST, 
	              FLOAT_MATRIX& Utilization, INT_MATRIX& link_state, FLOAT_MATRIX& Consumption,
				  double& total_energy)
{
	//link_state initially off
	total_energy = 0;
	int N = Topology.size();
	for(size_t i=0;i!=N;i++){
		link_state[i][i] = LINK_NOT_EXIST;
		for(size_t j=i+1;j!=N;j++){
			if(Topology[i][j] == -1){
				link_state[i][j] = link_state[j][i] = LINK_NOT_EXIST;
				continue;
			}
			if(MST[i][j]){ 
				link_state[i][j] = link_state[j][i] = LINK_ON;
				total_energy += Consumption[i][j];
			}
			else if(Utilization[i][j]!=0 ){ 
				link_state[i][j] = link_state[j][i] = LINK_ON;
				total_energy += Consumption[i][j];
			}
		}
	}
}


//�������н���SPT����ÿ����㰴�д����all_previous��
bool Generate_SPT(FLOAT_MATRIX& Topology, INT_MATRIX& SPT_previous)
{
	int N=Topology.size();
	SPT_previous.clear();
	for(int k=0;k!=N;k++)
	{
		vector<int> previous(N,-1);
		Dijkstra(Topology,k,previous);

		if(find(previous.begin(),previous.end(),-1)!=previous.end()) //not connected
		{	
			log_report("Generate_SPT fail!");
			return false;
		}

		SPT_previous.push_back(previous);
	}
	return true;
}

//�������������е�t���ڵ�����·����
//previous:��¼ÿ������ǰһ��·����Ԥ�ȷ���ã�
void Dijkstra(FLOAT_MATRIX &Topology, int t, vector<int> &previous)
{
	size_t node_num=Topology.size();
	previous[t]=t;
	vector<double> SP_len(node_num,-1);			//���·������
	SP_len[t]=0;
	vector<bool> visited(node_num,false);
	visited[t]=true;
	bool finished=false;
	for(size_t i=0;i!=previous.size();i++)
		if(Topology.at(t).at(i)!=-1)
		{
			SP_len[i]=Topology.at(t).at(i);
			previous[i]=t;
		}
		int count=1;
		while(count!=node_num)
		{
			double  min_latency=-1;	//INF
			int  min_index=-1, min_pre_index=t;
			//�ҵ����·��
			    for(size_t i=0; i!=visited.size();i++)
					if(!visited[i] && SP_len[i]!=-1 && (min_latency==-1 || SP_len[i]<min_latency))
					{
						min_latency=SP_len[i];
						min_index=i;
					}
				if(min_index == -1) break;
				visited[min_index]=true;
				for(size_t k=0;k!=SP_len.size();k++)
				{
					if(Topology.at(min_index).at(k)!=-1 && (SP_len[k]==-1 || SP_len[k]>SP_len[min_index]+Topology.at(min_index).at(k)))
					{	
						previous[k]=min_index;
						SP_len[k]=SP_len[min_index]+Topology.at(min_index).at(k);
					}
				}
				count++;
		}
}

void MRCT_MST(FLOAT_MATRIX &Topology,FLOAT_MATRIX &Power, INT_MATRIX &SPT_previous, BOOL_MATRIX &MST)
{
	//2 approximation��ѡ��median-rooted SPT
	double min_SP_cost = -1;
	int min_index = 0;

	int N = Topology.size();

	//��ÿ������SPT���б�����ѡ��median
	for(int k=0;k!=N;k++)
	{
		vector<int> previous = SPT_previous[k];
		double cost = 0;
		for(int t=0;t!=N;t++)
		{
			if(t==k)continue;
			int v = t, u = v;
			while(u!=k)
			{
				v = u;
				u = previous[u];
#ifdef CORE_GRAPH_ENERGY_AWARE
				cost += Power[u].at(v);
#else
				cost += Topology[u].at(v);
#endif
			}
		}

		if(min_SP_cost==-1 || (cost<min_SP_cost))
		{
			min_SP_cost = cost;
			min_index = k;
		}
	}

	//����Spanning Tree
	vector<bool> tmp_bool(N,false);
	MST.insert(MST.begin(),N,tmp_bool);

	vector<int> previous = SPT_previous[min_index];
	for(int t=0;t!=N;t++)
	{
		if(t==min_index)continue;
		int v = t, u = v;
		while(u!=min_index)
		{
			v = u;
			u = previous[u];
			MST[u].at(v) = true;
			MST[v].at(u) = true;
		}
	}
}

int MRCT2(FLOAT_MATRIX &Tpology, FLOAT_MATRIX &Power, INT_MATRIX &SPT_previous, BOOL_MATRIX &MST)
{
	//2 approximation��ѡ��median-rooted SPT
	double min_SP_cost = -1;
	int min_index = 0;
	double sec_min_SP_cost = -1;
	int sec_min_index = 0;

	int N = Tpology.size();

	//��ÿ������SPT���б�����ѡ��median
	for(int k=0;k!=N;k++)
	{
		vector<int> previous = SPT_previous[k];
		double cost = 0;
		for(int t=0;t!=N;t++)
		{
			if(t==k)continue;
			int v = t, u = v;
			while(u!=k)
			{
				v = u;
				u = previous[u];
#ifdef CORE_GRAPH_ENERGY_AWARE
				cost += Power[u].at(v);
#else
				cost += Tpology[u].at(v);
#endif
			}
		}

		if(min_SP_cost==-1 || (cost<sec_min_SP_cost))
		{
			sec_min_SP_cost = min_SP_cost;
			sec_min_index = min_index;
		}
		if(min_SP_cost==-1 || (cost<min_SP_cost))
		{
			sec_min_SP_cost = min_SP_cost;
			sec_min_index = min_index;
			min_SP_cost = cost;
			min_index = k;
		}
	}

	//����Spanning Tree
	vector<bool> tmp_bool(N,false);
	MST.insert(MST.begin(),N,tmp_bool);

	vector<int> previous = SPT_previous[min_index];
	for(int t=0;t!=N;t++)
	{
		if(t==min_index)continue;
		int v = t, u = v;
		while(u!=min_index)
		{
			v = u;
			u = previous[u];
			MST[u].at(v) = true;
			MST[v].at(u) = true;
		}
	}
	
#ifdef DOUBLE_ROOT_TREE
	previous = SPT_previous[sec_min_index];
	for(int t=0;t!=N;t++)
	{
		if(t==sec_min_index)continue;
		int v = t, u = v;
		while(u!=sec_min_index)
		{
			v = u;
			u = previous[u];
			MST[u].at(v) = true;
			MST[v].at(u) = true;
		}
	}
#endif
	return min_index;
}

//����MST�����ˣ�������ÿ��������ÿ��link�Ͼ�����flowӳ��
//MST_map:ÿ�������һ��map��ÿ��map��������link��ӳ������һ��vector<int>������
void MST_Statistics(vector<vector<int>> &MST_all_previous,vector<map<pair<int,int>,vector<int>>> &MST_map)
{
	int N = MST_all_previous.size();

	for(int k=0;k!=N;k++)
	{
		vector<int> MST_previous = MST_all_previous[k];
		map<pair<int,int>,vector<int>> k_map;

		for(int t=0;t!=N;t++)
		{
			if(t==k)continue;
			int v = t, u = v;

			while(u!=k)
			{
				v = u;
				u = MST_previous[u];

				//k->tҪ����(u,v)��ÿ����·�ǵ���ģ����Իῼ������
				pair<int,int> tmp_pair(u,v);

				map<pair<int,int>,vector<int>>::iterator it = k_map.find(tmp_pair);
				if(it==k_map.end())				//��link��δ���ֹ�
				{
					vector<int> vec_tmp;
					vec_tmp.push_back(t);
					k_map[tmp_pair] = vec_tmp;
				}
				else
				{
					k_map[tmp_pair].push_back(t);
				}

			}
		}
		MST_map.push_back(k_map);
	}
}

//����MST_Statisticsͳ��ÿ����·���������Ͼ�����flow����
//MST_pass:��ÿ����·��ÿ�����򾭹��Ľ����
void Link_Pass(LINK_MAP &MST_map, LINK_MAP &SPT_map, INT_MATRIX &MST_pass)
{
	int N = MST_map.size();
	vector<int> tmp_int(N,0);
	MST_pass.insert(MST_pass.begin(),N,tmp_int);

	for(int k=0;k!=N;k++)
		for(map<pair<int,int>,vector<int>>::iterator it = MST_map[k].begin(); it!=MST_map[k].end(); it++){
#ifdef OVERLAP_AWARE
			pair<int, int>link;link.first = it->first.first; link.second = it->first.second;
			map<pair<int,int>,vector<int>>::iterator it2 = SPT_map[k].find(link);
			if(it2==SPT_map[k].end())
#endif
				MST_pass[it->first.first].at(it->first.second) ++;
		}
}

//ͳ��SP·����ST·���Ƿ���ȫ�غ�
void SP_ST_Overlap(INT_MATRIX &SPT_all_previous, INT_MATRIX &MST_all_previous, BOOL_MATRIX &bOverlap)
{
	int N = SPT_all_previous.size();
	vector<bool> tmp_bool(N,false);
	bOverlap.insert(bOverlap.begin(),N,tmp_bool);

	for(int s=0;s!=N;s++)
	{
		vector<int> previous = SPT_all_previous[s];
		vector<int> MST_previous  = MST_all_previous[s];
		for(int d=0;d!=N;d++)
		{
			if(s==d)continue;
			bool test = true;			//true��ʾ��ȫ�غ�
			int v1=d, u1=v1;
			int v2=d, u2=v2;

			while(u1!=s && u2!=s)
			{
				v1 = u1, u1 = previous[u1];
				v2 = u2, u2 = MST_previous[u2];
				//���SP��ST��ȫ�غϣ����ǵ�˳��Ҳһ������ȫ��ͬ��
				if(u1!=u2 || v1!=v2)			//����ȫ��ͬ��˵��SP��ST����ȫ�غ�
				{
					test = false;
					break;
				}
			}
			if(test && (u1!=s || u2!=s))		//˵��·�����Ȳ�ͬ����ʱҲ�ǲ���ȫ�غϵ�
				test = false;

			bOverlap[s].at(d) = test;
		}
	}
}

//ȷ����������֮���delay
void Path_Delay(FLOAT_MATRIX &Tpology, INT_MATRIX &SPT_previous, INT_MATRIX &MST_previous, 
				FLOAT_MATRIX &SPT_DELAY, FLOAT_MATRIX &MST_DELAY,double& max, double& min, bool hop)
{
	int N = Tpology.size();
	vector<double> tmp(N,0);
	SPT_DELAY.clear();
	MST_DELAY.clear();
	SPT_DELAY.insert(SPT_DELAY.begin(),N,tmp);
	MST_DELAY.insert(MST_DELAY.begin(),N,tmp);

	max = 0; min = 999999999;
	for(int s=0;s!=N;s++)
	{
		vector<int> previous = SPT_previous[s];
		vector<int> previous2 = MST_previous[s];
		for(int d=0;d!=N;d++)
		{
			if(s==d)continue;
			int v = d, u = v;
			while(u!=s)
			{
				v = u;
				u = previous[u];
				if(hop)
					SPT_DELAY[s].at(d)++;
				else
					SPT_DELAY[s].at(d) += Tpology[u].at(v);
			}
			max = max>SPT_DELAY[s].at(d)?max:SPT_DELAY[s].at(d);
			min = min<SPT_DELAY[s].at(d)?min:SPT_DELAY[s].at(d);
			v = d, u = v;
			while(u!=s)
			{
				v = u;
				u = previous2[u];
				if(hop)
					MST_DELAY[s].at(d)++;
				else
					MST_DELAY[s].at(d) += Tpology[u].at(v);
			}
			max = max>MST_DELAY[s].at(d)?max:MST_DELAY[s].at(d);
			min = min<MST_DELAY[s].at(d)?min:MST_DELAY[s].at(d);
		}
	}

}

void Find_Utilization(FLOAT_MATRIX &Utilization, FLOAT_MATRIX &Capacity, UTILITY_INFO &utility_info)
{
	int node_num = Utilization.size();
	utility_info.max_utility=-1,utility_info.min_utility=0, utility_info.avg_utility=0;
	int count=0;
	utility_info.max_u = 0, utility_info.max_v = 0;
	for(size_t i=0;i!=node_num;i++)
		for(size_t j=0;j!=node_num;j++)
			if(Capacity[i].at(j)!=-1)
			{
				Utilization[i].at(j)/=Capacity[i].at(j);
				utility_info.avg_utility+=Utilization[i].at(j);
				count++;
				if(utility_info.max_utility==-1 || (Utilization[i].at(j)>utility_info.max_utility))
				{
					utility_info.max_utility=Utilization[i].at(j);
					utility_info.max_u = i;
					utility_info.max_v = j;
				}

				if(Utilization[i].at(j)<utility_info.min_utility)
					utility_info.min_utility=Utilization[i].at(j);

			}
			else
				Utilization[i].at(j)=-1;

	utility_info.avg_utility/=count;
}

void Prev_Utilization(FLOAT_MATRIX &TM, FLOAT_MATRIX &Capacity, BOOL_MATRIX &prev_path, 
					  INT_MATRIX &SPT_all_previous, INT_MATRIX &MST_all_previous,
					  FLOAT_MATRIX &Utilization, UTILITY_INFO &utility_info)
{
	int N = TM.size();

	//��ʼ��
	vector<double> tmp(N,0);
	Utilization.insert(Utilization.begin(),N,tmp);

	for(int k=0;k!=N;k++)
	{
		vector<int> previous = SPT_all_previous[k];
		vector<int> MST_previous  = MST_all_previous[k];

		for(int t=0;t!=N;t++)
		{
			if(t==k)continue;

			if(!prev_path[k].at(t))			//Spanning Tree
			{
				int v = t, u = v;
				while(u!=k)
				{
					v = u;
					u = MST_previous[u];
					//ֻ���ǵ���
					Utilization[u].at(v) += TM[k].at(t);	
				}
			}
			else
			{
				int v = t, u = v;
				while(u!=k)
				{
					v = u;
					u = previous[u];
					Utilization[u].at(v) += TM[k].at(t);
				}
			}
		}
	}

	Find_Utilization(Utilization,Capacity,utility_info);
}

void Budget_Adjust(FLOAT_MATRIX &Utilization, FLOAT_MATRIX &TM, FLOAT_MATRIX &Capacity,
				   LINK_MAP &MST_map, LINK_MAP &SPT_map,INT_MATRIX &MST_pass,
				   BOOL_MATRIX &MST, BOOL_MATRIX &prev_path,
				   TH_MAP &dFlow_th)
{
	int N = TM.size();
	//��ÿһ������MST��·������node����Ӧ��ֵ��Ҫ����
	for(int u=0;u!=N;u++)
		for(int v=0;v!=N;v++)
		{
			if(u==v)continue;
			
			//�����Ƿ����budget
			if(MST[u].at(v))
			{
				pair<int,int> index(u,v);
				//ͳ�Ƴ���link�ϵ�SPT������Ȼ��������
				double SPT_traffic = 0, MST_traffic = 0, MST_traffic_max = 0;
				for(int i=0;i!=N;i++)
				{
					double node_traffic = 0, node_traffic2 = 0, node_traffic3 = 0;
					bool dirST = MST_map[i].find(index)!=MST_map[i].end();	//spanning tree link����һ��
					bool dirSP = SPT_map[i].find(index)!=SPT_map[i].end();	//shortest path link����һ��
					if(dirST)			
					{
						for(int j=0;j!=((MST_map[i])[index]).size();j++)
						{
							if(!prev_path[i].at(((MST_map[i])[index])[j])		//����MST
								&& (!dirSP || (dirSP && find(((SPT_map[i])[index]).begin(),((SPT_map[i])[index]).end(),((MST_map[i])[index])[j])==((SPT_map[i])[index]).end())))//SP��ST��������·�ϲ��غ�(�غ�ʱ��Ϊ��SP traffic��	
								node_traffic += TM[i].at(((MST_map[i])[index])[j]);

							if(!dirSP || (dirSP && find(((SPT_map[i])[index]).begin(),((SPT_map[i])[index]).end(),((MST_map[i])[index])[j])==((SPT_map[i])[index]).end()))
								node_traffic3 += TM[i].at(((MST_map[i])[index])[j]);
						}
					}

					if(dirSP)
					{
						for(int j=0;j!=((SPT_map[i])[index]).size();j++)
							if(prev_path[i].at(((SPT_map[i])[index])[j])			//����SPT
								|| (dirST && find(((MST_map[i])[index]).begin(),((MST_map[i])[index]).end(),((SPT_map[i])[index])[j])!=((MST_map[i])[index]).end()))//SP��ST��������·���غϣ���ʱ��Ϊ��SP 
								node_traffic2 += TM[i].at(((SPT_map[i])[index])[j]);
					}

					SPT_traffic += node_traffic2;
					MST_traffic += node_traffic;
					MST_traffic_max += node_traffic3;
				}
				/*if(SPT_traffic+MST_traffic!=Utilization[u][v]*Capacity[u][v]){
						cout<<u<<"-"<<v<<": "<<SPT_traffic<<"\t"<<MST_traffic<<"\t"<<SPT_traffic+MST_traffic<<" "<<Utilization[u][v]*Capacity[u][v]<<endl;
				}
				if(u==17 && v==25){
					cout<<"haha";
				}*/

				//ȷ����Ҫ���ٶ���MST����
				double rm_rate,delta_h;
				if(Utilization[u].at(v)>UTILITY_UB)
				{
					//������һ��������������SPT_traffic���Ѿ�������60%����Ӧ�ý�����budget����
					if(SPT_traffic<=Capacity[u].at(v)*UTILITY_UB){
						rm_rate = (Capacity[u].at(v)*UTILITY_UB-SPT_traffic)/MST_traffic;			//����budget
					}
					else		//�����н���budget����
					{
						for(int i=0;i!=N;i++)
						{
							(dFlow_th[i])[index] = 0;
						}
						continue;
					}

				}
				else if(Utilization[u].at(v)<UTILITY_LB){
					delta_h = (Capacity[u].at(v)*UTILITY_UB-SPT_traffic-MST_traffic)/MST_pass[u].at(v);		//����budget
				}
				else{;
				}
				
				//��ÿ����㣬���¸���·����ֵ
				for(int i=0;i!=N;i++)
				{
					bool dirST = MST_map[i].find(index)!=MST_map[i].end();
					bool dirSP = SPT_map[i].find(index)!=SPT_map[i].end();
					if(!dirST)continue;

					vector<int> node = (MST_map[i])[index];
					vector<double> tm_i;
					for(int j=0;j!=node.size();j++)
						tm_i.push_back(TM[i].at(node[j]));

					//ͳ�Ƴ���һ��·����MST�ϵ���Ӧ����
					double MST_i_traffic = 0;
					double MST_i_traffic_max = 0;
					for(int j=0;j!=node.size();j++)
					{
						if(!prev_path[i].at(node[j]))
							MST_i_traffic += tm_i[j];

						MST_i_traffic_max += tm_i[j];
					}

					double new_th;
					if(Utilization[u].at(v)>UTILITY_UB)
						new_th = MST_i_traffic*rm_rate;
					else if(Utilization[u].at(v)<UTILITY_LB)
						new_th = MST_i_traffic+delta_h;
					else
						new_th = MST_i_traffic;

					(dFlow_th[i])[index] = new_th;
				}
			}
			
		}
}

void double_sort(vector<int> &OD, vector<double> &TM, bool ascend)
{
	int n = TM.size();
	for(int i=0;i!=n;i++)
	{
		double min_TM = TM[i];
		int min_index = i;
		for(int j=i+1;j!=n;j++)
		{
			if(ascend)
			{
				if(min_TM>TM[j])
				{
					min_TM = TM[j];
					min_index = j;
				}

			}
			else
			{
				if(min_TM<TM[j])
				{
					min_TM = TM[j];
					min_index = j;
				}

			}
		}

		double tmp = TM[min_index];
		int tmp2 = OD[min_index];

		TM[min_index] = TM[i];
		OD[min_index] = OD[i];
		TM[i] = tmp;
		OD[i] = tmp2;
	}
}

void Scheduling(INT_MATRIX &SPT_all_previous, INT_MATRIX &MST_all_previous,
						 BOOL_MATRIX &MST, LINK_MAP &SPT_MAP, LINK_MAP &MST_MAP,
						 FLOAT_MATRIX &Capacity, FLOAT_MATRIX &TM,
						 BOOL_MATRIX &prev_path, BOOL_MATRIX &bOverlap, 
						 TH_MAP &dFlow_th)
{
	srand((int)time(NULL));
	int N = SPT_all_previous.size();

	//��ÿ��ingress������
	for(int s=0;s!=N;s++)
	{
		//������һ��round�Ľ�����з���
		vector<int> ST_OD;							//��һ��round��ST
		vector<int> SP_OD;							//��һ��round��SP
		vector<double> TM_ST_OD;					//ST����Ӧ����
		vector<double> TM_SP_OD;					//SP����Ӧ����

		LINK_BUDGET budget = dFlow_th[s];			//link budget for s
		LINK_BUDGET total_traffic;					//traffic already on link
		LINK_BUDGET link_remove;					//ȷ��ÿ��link��Ӧ�����(1)���Ǽ���(0)����
		for(LINK_BUDGET::iterator it=budget.begin();it!=budget.end();it++)
		{
			total_traffic[it->first] = 0 ;
		}

		vector<int> previous = SPT_all_previous[s];
		vector<int> MST_previous = MST_all_previous[s];

		//divide OD pairs according to last round's route
		for(int d=0;d!=N;d++)
		{
			if(s==d)continue;
			//���s->d��SP��ST��ȫ�غϣ���������Ͳ�����
			if(bOverlap[s].at(d))
			{
				prev_path[s].at(d)=true;
				continue;
			}

			if(prev_path[s].at(d))
			{
				SP_OD.push_back(d);
				TM_SP_OD.push_back(TM[s].at(d));
			}
			else
			{
				ST_OD.push_back(d);
				TM_ST_OD.push_back(TM[s].at(d));

			}
		}

		switch(SCHEDULING){
			case SILO:
			//����TM������ingress router��������
			double_sort(ST_OD,TM_ST_OD,false);			//����
			double_sort(SP_OD,TM_SP_OD,true);			//����
			break;
			case LISO:
			//����TM������ingress router��������
			double_sort(ST_OD,TM_ST_OD,true);			//����
			double_sort(SP_OD,TM_SP_OD,false);			//����
			break;
			case RS:
			//������ingress router����������
			//double_sort(ST_OD,TM_ST_OD,false);	
			//double_sort(SP_OD,TM_SP_OD,true);			
			break;
		}

		//����һ��round·����ÿ��link�ϵ�����д��total_traffic��ֻд��ST��������Ϊ���ж�budget
		for(int k=0;k!=ST_OD.size();k++)
		{
			int d = ST_OD[k];
			int v = d, u = v;
			while(u!=s)
			{
				v = u; u = MST_previous[u];
				pair<int,int> tmp(u,v);
				//���������link��,s->d��SP��ST�غϣ��򲻿��Ǹ�����
				//�������linkͬʱ��s->d��ST��SP������link���򲻼���ST����
				bool dir = SPT_MAP[s].find(tmp)!=SPT_MAP[s].end();
				if(dir)				//SP��ȷʵ�и���·
					if(find((SPT_MAP[s])[tmp].begin(),(SPT_MAP[s])[tmp].end(),d)!=(SPT_MAP[s])[tmp].end())		//�غ�
						continue;
				total_traffic[tmp] += TM[s].at(d);
			}
		}
		//ȷ��ÿ��link�����ӻ��Ǽ�������
		for(LINK_BUDGET::iterator it=total_traffic.begin();it!=total_traffic.end();it++)
		{
			if(it->second>budget[it->first])	//��������
				link_remove[it->first] = 0;
			else
				link_remove[it->first] = 1;
		}

		//migrate from ST to SP
		for(int k=0;k!=ST_OD.size();k++)
		{
			//��ÿ��OD pair��������·��
			int d = ST_OD[k];

			int v = d, u = v;
			double min_p = NO_PR;				//�������ĸ���
			bool brm = false;
			while(u!=s)
			{
				v = u; u = MST_previous[u];
				pair<int,int> link(u,v);

				//�������linkͬʱ��s->d��ST��SP������link���򲻽���������ж�
				bool dir = SPT_MAP[s].find(link)!=SPT_MAP[s].end();
				if(dir)				//SP��ȷʵ�и���·
					if(find((SPT_MAP[s])[link].begin(),(SPT_MAP[s])[link].end(),d)!=(SPT_MAP[s])[link].end())		//�غ�
						continue;

				//����spanning tree�����budget���������󣬾Ͳ���Ҫ�Ƴ���������Ҫ
				if(link_remove[link]==0)		//��Ҫ�������������Կ��Ǹ���·
				{
					if(total_traffic[link]>=budget[link])		//we should move traffic out of link
					{
						brm = true;
						if(total_traffic[link]-budget[link]>TM[s].at(d))	//(s,d) should be definitely removed
							min_p = 0;
						else		//probabilistic method
						{
							if(PROBABILISTIC == NO_PROB){
								min_p = 0;
								break;
							}
							else{
							if(min_p>1-(total_traffic[link]-budget[link])/TM[s].at(d))
								min_p = 1-(total_traffic[link]-budget[link])/TM[s].at(d);
							}
						}
					}
					else //else��Ӧ������ȻҪ���٣���֮ǰ�Ĳ����Ѿ��ɹ�����������
						min_p = min_p<1 ? min_p : 1;			//ֻҪ�����·�����Ϳ��Ա���������
				}
			}
			if(min_p!=NO_PR)			//should remove the TM with probability
			{
				double p = (double)rand()/(double)RAND_MAX;
				if((PROBABILISTIC == PROB && p>min_p)||(PROBABILISTIC == NO_PROB && min_p == 0))
				{
					prev_path[s].at(d) = true;
					int v=d, u = v;
					while(u!=s)
					{
						v=u;u=MST_previous[u];
						pair<int,int> index(u,v);

						//�������linkͬʱ��s->d��ST��SP������link���򲻼���ST����
						bool dir = SPT_MAP[s].find(index)!=SPT_MAP[s].end();
						if(dir)				//SP��ȷʵ�и���·
							if(find((SPT_MAP[s])[index].begin(),(SPT_MAP[s])[index].end(),d)!=(SPT_MAP[s])[index].end())		//�غ�
								continue;

						total_traffic[index] -= TM[s].at(d);
					}
				}
			}

			//pair<int,int> index(17,25);
			//if(s == 10)
			//	cout<<"\tST2SP: "<<total_traffic[index]/dFlow_th[10][index]<<endl;
		}

		//migrate from SP to ST
		for(int k=0;k!=SP_OD.size();k++)
		{
			//��ÿ��OD pair��������·��
			int d = SP_OD[k];

			int v = d, u = v;
			double min_p = NO_PR;		//��Ǩ�Ƶĸ���
#ifdef TM_SP_TRESHOLD
			if(TM[s].at(d)<TM_sp_threshold)
#endif
			while(u!=s)
			{
				v = u; u = MST_previous[u];
				pair<int,int> link(u,v);

				//�������linkͬʱ��s->d��ST��SP������link���򲻽���������ж�
				bool dir = SPT_MAP[s].find(link)!=SPT_MAP[s].end();
				if(dir)				//SP��ȷʵ�и���·
					if(find((SPT_MAP[s])[link].begin(),(SPT_MAP[s])[link].end(),d)!=(SPT_MAP[s])[link].end())		//�غ�
						continue;

				//����SP�����budget��������С�����޷����롣�����������ʼ�Ѿ���spanning tree���������˵���
				if(link_remove[link]==1)			//��Ҫ��������
				{
					if(total_traffic[link]<=budget[link])
					{
						if(budget[link]-total_traffic[link]>TM[s].at(d))
							min_p = min_p<1 ? min_p : 1; // equal 1 means OK for this link
						else
						{
							if(PROBABILISTIC == NO_PROB){
								min_p = 0;
								break;
							}
							else{
							if(min_p>(budget[link]-total_traffic[link])/TM[s].at(d))
								min_p = (budget[link]-total_traffic[link])/TM[s].at(d);
							}
						}

					}
					else				//else��Ӧ������Ȼ��Ҫ���ӣ���ǰ��Ĳ����Ѿ��ɹ�����������
						min_p = 0;					//��ʱ����������������
					
				}
				else{
					min_p = 0;//tend to be remove traffic
				}
				
			}

			if(min_p!=NO_PR)			//should add the TM with probability
			{
				double p = PROB_FACTOR * (double)rand()/(double)RAND_MAX;


				if((PROBABILISTIC == PROB && p<min_p)||(PROBABILISTIC == NO_PROB && min_p == 1))
				{
					prev_path[s].at(d) = false;
					int v=d, u = v;
					while(u!=s)
					{
						v=u;u=MST_previous[u];
						pair<int,int> index(u,v);

						//�������linkͬʱ��s->d��ST��SP������link���򲻼���ST����
						bool dir = SPT_MAP[s].find(index)!=SPT_MAP[s].end();
						if(dir)				//SP��ȷʵ�и���·
							if(find((SPT_MAP[s])[index].begin(),(SPT_MAP[s])[index].end(),d)!=(SPT_MAP[s])[index].end())		//�غ�
								continue;

						//if(index.first == 17 && index.second == 25){
						//	cout<<SP_OD[k]
						//}
						total_traffic[index] += TM[s].at(d);
					}
				}
			}
			//pair<int,int> index(17,25);
			//if(s == 10)
			//	cout<<"\tSP2ST: "<<total_traffic[index]/dFlow_th[10][index]<<endl;

		}
	}
}

extern vector<TEST_DESCRIPTOR> test_vector;
extern int current_i;
void log_report(string out){
	ofstream fout(LOG_FILE,ios::app);
	time_t itime;
	time(&itime);
	fout<<ctime(&itime)<<" "<<test_vector[current_i]._str<<" "<<out<<endl;
	fout.close();
}

void fail_report(){
	ofstream fout(FAIL_LIST,ios::app);
	fout<<test_vector[current_i]._str<<endl;
	fout.close();
}

extern int dumplinktime;

void DumpLinkUtilization(FLOAT_MATRIX& Topology, FLOAT_MATRIX& Util, BOOL_MATRIX& MST, bool MST_TRUE){

	string fileoutname;
	if(MST_TRUE){
		fileoutname = LOG_MST_UTIL_FILE;
		fileoutname += test_vector[current_i]._str;
	}
	else{
		fileoutname = LOG_ORDINARY_UTIL_FILE;
		fileoutname += test_vector[current_i]._str;
	}

	ofstream fout;
	if(dumplinktime < 2)
		fout.open(fileoutname.c_str());
	else
		fout.open(fileoutname.c_str(),ios::app);
	dumplinktime++;

	bool stop = false;
	for(int i=0; i<Topology.size();i++){
		for(int j=i+1; j<Topology.size();j++){
			if(Topology[i][j]==-1) continue;
			if(MST_TRUE && MST[i][j]) fout<<i<<" "<<j<<" "<<Util[i][j]<<endl;
			if(!MST_TRUE && !MST[i][j]) fout<<i<<" "<<j<<" "<<Util[i][j]<<endl;
			if(Util[i][j]>1){
				cout<<"    overload link:"<<i<<"-"<<j<<endl;
				stop = true;
			}
		}
	}
	fout.close();
	if(stop)
	system("pause");
}

void DumpPrev_path(FLOAT_MATRIX& Topology, BOOL_MATRIX& prev_path){

	string fileoutname;
		fileoutname = LOG_PREV_PATH_FILE;
		fileoutname += test_vector[current_i]._str;
	ofstream fout;
		fout.open(fileoutname.c_str(),ios::app);

	for(int i=0; i<Topology.size();i++){
		for(int j=0; j<Topology.size();j++){
			fout<<prev_path[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();
}

void DumpRouting(FLOAT_MATRIX Topology, INT_MATRIX& spt_routing,INT_MATRIX& mst_routing){
	string fileout = ROUTING_DIR;
	fileout += test_vector[current_i]._str;
	fileout += ".txt";
	ofstream fout(fileout.c_str());
	fout<<"{"<<endl;
	int N = spt_routing.size();
	for(int src = 0; src<N; src++){
		for(int dst = 0; dst<N; dst++){
			if(dst == src)continue;
			fout<<"\""<<src<<"_"<<dst<<"\": ["<<endl;

			//spf rouitng
			double weight_sum = 0;
			fout<<"{"<<endl<<"\"path\": ["<<endl;
			int v = src; int u = spt_routing[dst][v];
			weight_sum += Topology[v][u];
			fout<<v<<","<<endl;
			while(u!=dst){
				fout<<u<<","<<endl;
				v = u;
				u = spt_routing[dst][v];
				weight_sum += Topology[v][u];
			}
			fout<<u<<endl;
			fout<<"],"<<endl<<"\"wt\": "<<(int)weight_sum<<endl<<"},"<<endl;

			//mst routing
			weight_sum = 0;
			fout<<"{"<<endl<<"\"path\": ["<<endl;
			 v = src;  u = mst_routing[dst][v];
			weight_sum += Topology[v][u];
			fout<<v<<","<<endl;
			while(u!=dst){
				fout<<u<<","<<endl;
				v = u;
				u = mst_routing[dst][v];
				weight_sum += Topology[v][u];
			}
			fout<<u<<endl;
			fout<<"],"<<endl<<"\"wt\": "<<(int)weight_sum<<endl<<"}"<<endl<<"],"<<endl;
		}
	}
	fout<<"}";
	fout.close();
}

void DumpCoreLinks(BOOL_MATRIX MST){
	string fileout = CORELINKS_DIR;
	fileout += test_vector[current_i]._str;
	fileout += ".txt";
	ofstream fout(fileout.c_str());
	fout<<"["<<endl;
	int N = MST.size();
	bool comma = false;

	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			if(MST[i][j]==1){
				if(comma == false){comma = true;}
				else{fout<<"\n,";}
				fout<<"\""<<i<<"_"<<j<<"\"";
			}
		}
	}
	fout<<"\n]";
	fout.close();

	//for(int src = 0; src<N; src++){
	//	for(int dst = 0; dst<N; dst++){
	//		if(dst == src)continue;
	//		fout<<"\""<<src<<"_"<<dst<<"\": ["<<endl;

	//		//spf rouitng
	//		double weight_sum = 0;
	//		fout<<"{"<<endl<<"\"path\": ["<<endl;
	//		int v = src; int u = spt_routing[dst][v];
	//		weight_sum += Topology[v][u];
	//		fout<<v<<","<<endl;
	//		while(u!=dst){
	//			fout<<u<<","<<endl;
	//			v = u;
	//			u = spt_routing[dst][v];
	//			weight_sum += Topology[v][u];
	//		}
	//		fout<<u<<endl;
	//		fout<<"],"<<endl<<"\"wt\": "<<(int)weight_sum<<endl<<"},"<<endl;

	//		//mst routing
	//		weight_sum = 0;
	//		fout<<"{"<<endl<<"\"path\": ["<<endl;
	//		 v = src;  u = mst_routing[dst][v];
	//		weight_sum += Topology[v][u];
	//		fout<<v<<","<<endl;
	//		while(u!=dst){
	//			fout<<u<<","<<endl;
	//			v = u;
	//			u = mst_routing[dst][v];
	//			weight_sum += Topology[v][u];
	//		}
	//		fout<<u<<endl;
	//		fout<<"],"<<endl<<"\"wt\": "<<(int)weight_sum<<endl<<"}"<<endl<<"],"<<endl;
	//	}
	//}
	//fout<<"}";
	//fout.close();
}