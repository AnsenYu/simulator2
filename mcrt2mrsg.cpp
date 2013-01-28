#include "mcrt2mrsg.h"
#include <vector>
#include <stack>
#include <queue>

typedef struct{
    int _i;
	int _j;
	double _dis;
}triple;

void AddLink(int i, int j, int o_root, INT_MATRIX& parent, BOOL_MATRIX& MST);
void UpdateDis(FLOAT_MATRIX Topology,vector<triple>& ordinary_link, INT_MATRIX parent, int o_root);
int tree_path(int i, int j, int root, vector<int>& path, INT_MATRIX& parent);
int path_length(int i, int j,int root, INT_MATRIX parent);
void bfs_build_tree(int i, BOOL_MATRIX MST, INT_MATRIX& parent);


bool triple_cmp(const triple x, const triple y){
	if(x._dis > y._dis)
		return true;
	else
		return false;
}

void MRSG_from_MCRT(FLOAT_MATRIX Topology, BOOL_MATRIX& MST, int root,double mlu){
        int NO_LINK = -1;
	int N = Topology.size();
	INT_MATRIX topo_neighbor;
	topo_neighbor.clear();
        vector<int> tmpvector(2);
	INT_MATRIX parent(N,tmpvector);	
	bfs_build_tree(root, MST, parent);	
	vector<triple> ordinary_link;
	int total_link;
	int always_on_link;
	total_link = always_on_link = N-1;

	for(int i = 0; i<N;++i){
		for(int j = i+1; j<N;++j){
			if(Topology[i][j]!=NO_LINK && MST[i][j]==false){
			    triple tmp; tmp._i = i, tmp._j = j; tmp._dis = -1;
				ordinary_link.push_back(tmp);
				total_link++;
			}
		}
	}	

	while(true){
		UpdateDis(Topology,ordinary_link, parent, root);
		sort(ordinary_link.begin(), ordinary_link.end(), triple_cmp);// Descending
		if(ordinary_link.empty() || ordinary_link[0]._dis==0) break;//it is empty or all ordinary links are already in the same loop
		AddLink(ordinary_link[0]._i, ordinary_link[0]._j, root, parent, MST);
		always_on_link++;
	}
	if(ordinary_link.empty())return;

#ifdef RANDOM_ADD_LINK
#ifdef ALWAYS_ON_LINK_PERTENTAGE
	while(((double)always_on_link)/total_link < ALWAYS_ON_LINK_PERTENTAGE ){
		int rand_i = rand()%N;
		int rand_j = rand()%N;
		if(rand_i!=rand_j && Topology[rand_i][rand_j]!=-1 && MST[rand_i][rand_j]==false){
			MST[rand_i][rand_j] = true;
			always_on_link++;
		}
	}
#else
	int index = 0;
	while(mlu>0 && index<ordinary_link.size() && ((double)always_on_link)/total_link <  mlu/100){
		//if(MST[ordinary_link[index]._i][ordinary_link[index]._j] == false){
		//	MST[ordinary_link[index]._i][ordinary_link[index]._j] = true;
		//	always_on_link++;
		//}
		//index++;
		int rand_i = rand()%N;
		int rand_j = rand()%N;
		if(rand_i!=rand_j && Topology[rand_i][rand_j]!=-1 && MST[rand_i][rand_j]==false){
			MST[rand_i][rand_j] = true;
			always_on_link++;
		}
	}
#endif
#endif
        return;
}

int TopParent(INT_MATRIX parent,int i){
    while(parent[i][1]!=i) i=parent[i][1]; //to find the top parent of the same loop, initially itself
	return i;
}

void AddLink(int i, int j, int o_root, INT_MATRIX& parent, BOOL_MATRIX& MST){
    vector<int> path;
    int root = tree_path(i,j,o_root, path,parent);
    for(int k=0; k<path.size(); k++){
        parent[path[k]][1] = TopParent(parent, root);//some problem!!!
    }
    MST[i][j] = MST[j][i] = true;
}

void UpdateDis(FLOAT_MATRIX Topology, vector<triple>& ordinary_link, INT_MATRIX parent, int o_root){
    for(int k = 0; k<ordinary_link.size(); k++){
	    int li = ordinary_link[k]._i;
	    int lj = ordinary_link[k]._j;
	    if(TopParent(parent, li) == TopParent(parent, lj)){
			//they are already in the same circle
		    ordinary_link[k]._dis = 0;
		}
	    else{
#ifdef A0B1
			if(Topology[li][lj]!=0)
		        ordinary_link[k]._dis = 1/Topology[li][lj];
			else
				ordinary_link[k]._dis = 1/0.001;
#endif
#ifdef A1B1
			if(Topology[li][lj]!=0)
		        ordinary_link[k]._dis = path_length(li,lj,o_root, parent)/Topology[li][lj];
			else
				ordinary_link[k]._dis = path_length(li,lj,o_root, parent)/0.001;
#endif
#ifdef A1B0
		        ordinary_link[k]._dis = path_length(li,lj,o_root, parent);//o_root is the original root
#endif
	    }
    }
}

int tree_path(int i, int j, int root, vector<int>& path, INT_MATRIX& parent){
    //return the exact path length from i to the cycle of j
    stack<int> stacki,stackj;
	stacki.push(i); stackj.push(j);
	while(parent[i][0] != root){
	    stacki.push(parent[i][0]);
		i = parent[i][0];
	}
	while(parent[j][0] != root){
	    stackj.push(parent[j][0]);
		j = parent[j][0];
	}
	int top = root;
	while(stacki.top()==stackj.top()){
	    top = stacki.top();
		if(stacki.size()==1){//考虑i,j到root的路径重合
			path.clear();
			while(!stackj.empty()){
				path.push_back(stackj.top());stackj.pop();
			}
			return top;
		}
		if(stackj.size()==1){
			path.clear();
			while(!stacki.empty()){
				path.push_back(stacki.top());stacki.pop();
			}
			return top;
		}
	    stacki.pop(); stackj.pop();
	}
	stack<int> tmp;
	while(!stacki.empty()){
		tmp.push(stacki.top());
		stacki.pop();
	}
	path.clear();
	while(!tmp.empty()){
	    path.push_back(tmp.top()); tmp.pop();
	}
	path.push_back(top);
	while(!stackj.empty()){
	    path.push_back(stackj.top()); stackj.pop();
	}	
	return top;
}


int path_length(int i, int j,int root, INT_MATRIX parent){
    vector<int> path;
    tree_path(i,j,root,path,parent);
    int dis = 0;
    for(int k=1; k<path.size(); k++){
	if(TopParent(parent,path[k])!=TopParent(parent,path[k-1])) dis++; //if the current node is in the same loop of the previous node, no need to increase the dis
    }
    return dis;
}


void bfs_build_tree(int i, BOOL_MATRIX MST, INT_MATRIX& parent){
//for each in parent is a 2-tuple, the first is its actual direct parent, the next is the loops it belongs to.
	vector<bool> visited(MST.size(), false);
	queue<int> queue_i;
	queue_i.push(i);
	parent[i][0] = i; parent[i][1] = i;
	while(!queue_i.empty()){
		int top_node = queue_i.front();
		for(int k=0; k<MST.size(); ++k){
			if(k == parent[top_node][0]) continue;
			if(MST[top_node][k]!=false){
			    parent[k][0] = top_node; //set tree parent to top_node
			    parent[k][1] = k; //initialized belong to self-loops
			    queue_i.push(k);
			}
		}		
		queue_i.pop();
	}
    return;
}

