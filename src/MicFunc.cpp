#include "MicFunc.h"

/*
__attribute__((target(mic))) void funcheck(int i)
{
	#ifdef __MIC__
		printf("Index on MIC:%d\n",i);
	#else
		printf("Index on CPU:%d\n",i);
	#endif
}
*/

void cmsingleEdgeGraphMining(const Edge &e, vector<Edge> &single_edge_graph, int thread_num, int begin, int end, int mic_thread)
{
	
	int mic_num=2;
		
	// GS <- GS - es (make sure the results will not rebundant)
	for(int i=begin; i < end; i++)
	{
		Edge tmp_e = single_edge_graph[i];
		for (int j = 0; j < nr_graph; j++)  
			GS[j].removeEdge(tmp_e.x, tmp_e.a, tmp_e.y); 
	}
	
	GraphCode gc;
	gc.seq.push_back(&e);
	for (int j = 0; j < nr_graph; ++j)  
		if (GS[j].hasEdge(e.x, e.a, e.y))  
			gc.gs.push_back(j);
		
	vector<GraphCode> two_edges_child_gcs;
	vector<int> two_edges_nexts;
	vector<GraphCode> two_edges_child_gcs_cpu;
	vector<int> two_edges_nexts_cpu;
	vector<GraphCode> two_edges_child_gcs_mic;
	vector<int> two_edges_nexts_mic;
	
	//first expansion to get the children set(two edges frequent graph)
	one_edge_expansion(gc, 2, two_edges_child_gcs, two_edges_nexts);
	
	int n = two_edges_child_gcs.size();
	int cpu_len;
	int* index1 = new int[n];
	split_data_interval(n, mic_num+1, 0, index1, cpu_len);
	for(int i=0;i<cpu_len;i++)
	{
		two_edges_child_gcs_cpu.push_back(two_edges_child_gcs[index1[i]]);
		two_edges_nexts_cpu.push_back(two_edges_nexts[index1[i]]);
	}
	delete[] index1;
	
	//submining per level on CPU
	freqGraphMiningfrom2edgesOnCPU(two_edges_child_gcs_cpu, two_edges_nexts_cpu, thread_num);	
	
	int mic_len;
	int *index2 = new int[n];
	
	for(int i=1;i<=mic_num;i++)
	{
		two_edges_child_gcs_mic.clear();
		two_edges_nexts_mic.clear();
		
		split_data_interval(n, mic_num+1, i, index2, mic_len);
		for(int j=0;j<mic_len;j++)
		{
			two_edges_child_gcs_mic.push_back(two_edges_child_gcs[index2[j]]);
			two_edges_nexts_mic.push_back(two_edges_nexts[index2[j]]);
		}
		
		/* apply for space for preparing data to offload */
		int offload_len = two_edges_child_gcs_mic.size();
		int *ix = (int *)malloc((offload_len+1)*sizeof(int));
		int *iy = (int *)malloc((offload_len+1)*sizeof(int));
		int *x = (int *)malloc((offload_len+1)*sizeof(int));
		int *a = (int *)malloc((offload_len+1)*sizeof(int));
		int *y = (int *)malloc((offload_len+1)*sizeof(int));
		int gs_len=0;
		for(size_t j=0;j<offload_len;j++)
		{
			gs_len += two_edges_child_gcs_mic[j].gs.size();
		}
		gs_len += offload_len;
		int *gs = (int *)malloc((gs_len)*sizeof(int));
		int *nexts = (int *)malloc((offload_len)*sizeof(int));
		
		/* create bit-copyable data to offload */
		createDataForOffload(two_edges_child_gcs_mic, two_edges_nexts_mic, ix, iy, x, a, y, gs, nexts);
		
		/* offload to mic to execute */
		vector<GraphCode> two_edges_gcs;
		vector<int> two_edges_nexts;
		reconstructDataforMiningOnMIC(ix, iy, x, a, y, gs, nexts, offload_len, gs_len, two_edges_gcs, two_edges_nexts);
		freqGraphMiningfrom2edgesOnMIC(two_edges_gcs, two_edges_nexts, mic_thread);	
		
		
		delete[] ix;
		delete[] iy;
		delete[] x;
		delete[] a;
		delete[] y;
		delete[] gs;
		delete[] nexts;
	}
	delete[] index2;		
}

void createDataForOffload(vector<GraphCode> &two_edges_child_gcs_mic, vector<int> &two_edges_nexts_mic,   /* input paras */
				int *ix, int *iy, int *x, int *a, int *y,  /* split edge parameters(n+1)  */ 
				int *gs, /* combination of gss split by -1 */
				int *nexts /* nexts(n) */)
{
	
}

void reconstructDataforMiningOnMIC(int *ix, int *iy, int *x, int *a, int *y,  /* input split edge parameters(offload_len+1)  */ 
				int *gs, int *nexts, /* input combination of gss split by -1(gs_len) and nexts (offload_len) */ 
				int offload_len, int gs_len, /* length of input parameters  */
				vector<GraphCode> &two_edges_child_gcs_mic, vector<int> &two_edges_nexts_mic /* output  */)
{
					
}

void freqGraphMiningfrom2edgesOnMIC(vector<GraphCode> &two_edges_child_gcs, vector<int> &two_edges_nexts, int mic_thread)
{

	vector<GraphCode> current_global_child_gcs;
	vector<int> current_global_nexts;
	vector<GraphCode> tmp_global_child_gcs;
	vector<int> tmp_global_nexts;	
	
	int mic_len = two_edges_child_gcs.size();
	
	current_global_child_gcs.swap(two_edges_child_gcs);
	current_global_nexts.swap(two_edges_nexts);
	//submining per level in cpu
	while(mic_len != 0)
	{	
		#pragma omp parallel num_threads(mic_thread)
		{					
			vector<GraphCode> local_child_gcs;
			vector<int> local_nexts;
			vector<GraphCode> current_local_child_gcs;
			vector<int> current_local_nexts;
			
			#pragma omp for schedule(dynamic)
			for(size_t i=0; i<mic_len; i++)
			{
				//carry out current child's one edge expansion
				one_edge_expansion(current_global_child_gcs[i], current_global_nexts[i], local_child_gcs, local_nexts);
					
				//add the current child's expansion to local result
				current_local_child_gcs.insert(current_local_child_gcs.end(), local_child_gcs.begin(), local_child_gcs.end());
				current_local_nexts.insert(current_local_nexts.end(), local_nexts.begin(), local_nexts.end());
			
				//clear current child's expansion to carry out the next child's expansion
				local_child_gcs.clear();
				local_nexts.clear();
			}
					
			//merge local results into tmp global results  
			#pragma omp critical
			{
				tmp_global_child_gcs.insert(tmp_global_child_gcs.end(),current_local_child_gcs.begin(),current_local_child_gcs.end());
				tmp_global_nexts.insert(tmp_global_nexts.end(),current_local_nexts.begin(),current_local_nexts.end());
			}
		}
			
		//swap tmp global results and global results
		current_global_child_gcs.swap(tmp_global_child_gcs);
		current_global_nexts.swap(tmp_global_nexts);
		
		//clear tmp global results to carry out the next level mining
		tmp_global_child_gcs.clear();
		tmp_global_nexts.clear();
			
		//get the size of global results to judge whether continue next level mining
		mic_len = current_global_child_gcs.size();
	}		
}


void cmsingleEdgeGraphMining_simulationOnCPU(const Edge &e, vector<Edge> &single_edge_graph, int thread_num, int begin, int end, int mic_thread)
{
	
	int mic_num=2;
		
	// GS <- GS - es (make sure the results will not rebundant)
	for(int i=begin; i < end; i++)
	{
		Edge tmp_e = single_edge_graph[i];
		for (int j = 0; j < nr_graph; j++)  
			GS[j].removeEdge(tmp_e.x, tmp_e.a, tmp_e.y); 
	}
	
	GraphCode gc;
	gc.seq.push_back(&e);
	for (int j = 0; j < nr_graph; ++j)  
		if (GS[j].hasEdge(e.x, e.a, e.y))  
			gc.gs.push_back(j);
		
	vector<GraphCode> two_edges_child_gcs;
	vector<int> two_edges_nexts;
	vector<GraphCode> two_edges_child_gcs_cpu;
	vector<int> two_edges_nexts_cpu;
	vector<GraphCode> two_edges_child_gcs_mic;
	vector<int> two_edges_nexts_mic;
	
	//first expansion to get the children set(two edges frequent graph)
	one_edge_expansion(gc, 2, two_edges_child_gcs, two_edges_nexts);
		
	int n = two_edges_child_gcs.size();
	int cpu_len;
	int* index1 = new int[n];
	split_data_interval(n, mic_num+1, 0, index1, cpu_len);
	for(int i=0;i<cpu_len;i++)
	{
		two_edges_child_gcs_cpu.push_back(two_edges_child_gcs[index1[i]]);
		two_edges_nexts_cpu.push_back(two_edges_nexts[index1[i]]);
	}
	delete[] index1;
	
	//submining per level on CPU
	freqGraphMiningfrom2edgesOnCPU( two_edges_child_gcs_cpu, two_edges_nexts_cpu, thread_num);
	
	
	//submining per level on MIC
	vector<Graph *> MIC_S;
	
//	constructGraphSetOnMIC();
	
	int mic_len;
	int *index2 = new int[n];
	
	for(int i=1;i<=mic_num;i++)
	{
		MIC_S.clear();
		two_edges_child_gcs_mic.clear();
		two_edges_nexts_mic.clear();
		
		split_data_interval(n, mic_num+1, i, index2, mic_len);
		for(int j=0;j<mic_len;j++)
		{
			two_edges_child_gcs_mic.push_back(two_edges_child_gcs[index2[j]]);
			two_edges_nexts_mic.push_back(two_edges_nexts[index2[j]]);
		}
		
		freqGraphMiningfrom2edgesOnMIC_simulation(two_edges_child_gcs_mic, two_edges_nexts_mic, mic_thread, MIC_S);
		S.insert(S.end(), MIC_S.begin(), MIC_S.end());
	}
	delete[] index2;	
}

void freqGraphMiningfrom2edgesOnCPU(vector<GraphCode> &two_edges_child_gcs, vector<int> &two_edges_nexts, int thread_num)
{

	vector<GraphCode> current_global_child_gcs;
	vector<int> current_global_nexts;
	vector<GraphCode> tmp_global_child_gcs;
	vector<int> tmp_global_nexts;
	
	/*
	int n = two_edges_child_gcs.size();
	int cpu_len,iternum=0;
	
	if(n%2==0)
		cpu_len = n/2;
	else
		cpu_len = n/2+1;
	
	int* index = new int[cpu_len];
	split_data_interval(n, 0, index);
	*/
	
	int cpu_len = two_edges_child_gcs.size();
	
	current_global_child_gcs.swap(two_edges_child_gcs);
	current_global_nexts.swap(two_edges_nexts);
	//submining per level in cpu
	while(cpu_len != 0)
	{
		//iternum++;
		#pragma omp parallel num_threads(thread_num)
		{		
			
			vector<GraphCode> local_child_gcs;
			vector<int> local_nexts;
			vector<GraphCode> current_local_child_gcs;
			vector<int> current_local_nexts;
			
			#pragma omp for schedule(dynamic)
			for(size_t i=0; i<cpu_len; i++)
			{
				/*
				int pos;
				if(iternum == 1)
					pos = index[i];
				else
					pos = i;
				*/
				
				//carry out current child's one edge expansion
				one_edge_expansion(current_global_child_gcs[i], current_global_nexts[i], local_child_gcs, local_nexts);
					
				//add the current child's expansion to local result
				current_local_child_gcs.insert(current_local_child_gcs.end(), local_child_gcs.begin(), local_child_gcs.end());
				current_local_nexts.insert(current_local_nexts.end(), local_nexts.begin(), local_nexts.end());
			
				//clear current child's expansion to carry out the next child's expansion
				local_child_gcs.clear();
				local_nexts.clear();
			}
					
			//merge local results into tmp global results  
			#pragma omp critical
			{
				tmp_global_child_gcs.insert(tmp_global_child_gcs.end(),current_local_child_gcs.begin(),current_local_child_gcs.end());
				tmp_global_nexts.insert(tmp_global_nexts.end(),current_local_nexts.begin(),current_local_nexts.end());
			}
		}
		
		//swap tmp global results and global results
		current_global_child_gcs.swap(tmp_global_child_gcs);
		current_global_nexts.swap(tmp_global_nexts);
		
		//clear tmp global results to carry out the next level mining
		tmp_global_child_gcs.clear();
		tmp_global_nexts.clear();
			
		//get the size of global results to judge whether continue next level mining
		cpu_len = current_global_child_gcs.size();
	}	
	//delete index;	
}

void freqGraphMiningfrom2edgesOnMIC_simulation(vector<GraphCode> &two_edges_child_gcs, vector<int> &two_edges_nexts, int mic_thread, vector<Graph *> &MIC_S)
{

	vector<GraphCode> current_global_child_gcs;
	vector<int> current_global_nexts;
	vector<GraphCode> tmp_global_child_gcs;
	vector<int> tmp_global_nexts;
	
	/*
	int n = two_edges_child_gcs.size();
	int mic_len,iternum=0;
	
	mic_len = n/2;
	
	int* index = new int[mic_len];
	split_data_interval(n, 1, index);
	*/
	
	int mic_len = two_edges_child_gcs.size();
	
	current_global_child_gcs.swap(two_edges_child_gcs);
	current_global_nexts.swap(two_edges_nexts);
	//submining per level in cpu
	while(mic_len != 0)
	{
		//iternum++;
		#pragma omp parallel num_threads(mic_thread)
		{					
			vector<GraphCode> local_child_gcs;
			vector<int> local_nexts;
			vector<GraphCode> current_local_child_gcs;
			vector<int> current_local_nexts;
			
			#pragma omp for schedule(dynamic)
			for(size_t i=0; i<mic_len; i++)
			{
				/*
				int pos;
				if(iternum == 1)
					pos = index[i];
				else
					pos = i;
				*/
				
				//carry out current child's one edge expansion
				one_edge_expansion_mic_simulation(current_global_child_gcs[i], current_global_nexts[i], local_child_gcs, local_nexts, MIC_S);
					
				//add the current child's expansion to local result
				current_local_child_gcs.insert(current_local_child_gcs.end(), local_child_gcs.begin(), local_child_gcs.end());
				current_local_nexts.insert(current_local_nexts.end(), local_nexts.begin(), local_nexts.end());
			
				//clear current child's expansion to carry out the next child's expansion
				local_child_gcs.clear();
				local_nexts.clear();
			}
					
			//merge local results into tmp global results  
			#pragma omp critical
			{
				tmp_global_child_gcs.insert(tmp_global_child_gcs.end(),current_local_child_gcs.begin(),current_local_child_gcs.end());
				tmp_global_nexts.insert(tmp_global_nexts.end(),current_local_nexts.begin(),current_local_nexts.end());
			}
		}
			
		//swap tmp global results and global results
		current_global_child_gcs.swap(tmp_global_child_gcs);
		current_global_nexts.swap(tmp_global_nexts);
		
		//clear tmp global results to carry out the next level mining
		tmp_global_child_gcs.clear();
		tmp_global_nexts.clear();
			
		//get the size of global results to judge whether continue next level mining
		mic_len = current_global_child_gcs.size();
	}
	//delete index;		
}

void one_edge_expansion_mic_simulation(GraphCode &gc, int next, vector<GraphCode> &child_gcs, vector<int> &nexts, vector<Graph *> &MIC_S)
{
	/* construct graph from DFS code */  
    Graph *g = new Graph;  
    vector<const Edge *> &s = gc.seq;  
    g->node_label.resize(next);  
    g->edge_label = new vector<int>[next];  
    g->edge_next = new vector<int>[next];  
	
    for (size_t i = 0; i < s.size(); ++i)  
    {  
		int ix = s[i]->ix, iy = s[i]->iy, a = s[i]->a; 
        g->node_label[ix] = s[i]->x;  
        g->node_label[iy] = s[i]->y;	
        g->edge_label[ix].push_back(a);	
        g->edge_label[iy].push_back(a); 
		g->edge_next[ix].push_back(iy);  
		g->edge_next[iy].push_back(ix);	
    }  
  
    /* minimum DFS code pruning stage (4) */  
    Traveler t(s, *g);  
    t.travel();  
    if (!t.isMin())  
    {  
        delete[] g->edge_label;  
        delete[] g->edge_next;  
        delete g; 
        return;  
    }  
	
//	g->gs = gc.gs;	  //final line’s swap operation will finish this work，this will be work because g is a point.
	#pragma omp critical
    MIC_S.push_back(g);   //add a new result to mic (frequent subgraph) 
    
    /* enumerate potential children with 1-edge growth */  
    map<Edge, vector<int> > m;  
    for (size_t i = 0; i < gc.gs.size(); ++i)  
    {  
        set<Edge> c;  
        SubTraveler st(s, GS[gc.gs[i]], c, min_support, EF);  
        st.travel(next);  
        for (set<Edge>::const_iterator j = c.begin(); j != c.end(); ++j)  
            m[*j].push_back(gc.gs[i]);    //get the map(edge:gs) for every child's frequency
    }  
	
    /* mining subgraph children */  
    for (map<Edge, vector<int> >::iterator i = m.begin(); i != m.end(); ++i)  
    {  
		//must new the edge to make sure the spare will not free after function return
		//there is no need to delete e in mining stage,  
		//because ‘child_gc.seq = s‘ lead every child will use parent seqs.
		//if delete before program end, must error(delete one pointer more than once)
        Edge* e = new Edge(i->first.ix,i->first.iy,i->first.x,i->first.a,i->first.y);
		//const Edge *e = &i->first;  
        vector<int> &spp = i->second;   //gs for current extension
        if ((int)spp.size() < min_support)  
            continue;  
        GraphCode child_gc;  
        child_gc.gs.swap(spp);  
        child_gc.seq = s;  
        child_gc.seq.push_back(e);  //get the correct extend and carry out the next subminingW
        
		child_gcs.push_back(child_gc);
		if (e->iy == next)   //vextex add, backward extension, continue to next level
			nexts.push_back(next+1);
        else  //vextex invariant, forward extension, continue to next level
            nexts.push_back(next);   
    } 
    g->gs.swap(gc.gs); 
}

void split_data_interval(int n, int num, int order, int* index, int& local_n)
{
	
	local_n = 0;
	
	if(n==0)
		return;
	
	for(int i=0; i<n; i++)
	{
		if(i%num==order)
		{	 //cpu,mic0,mic1...,mic[num-1]
			index[local_n++] = i;
		}
	}
}

void constructGraphSetOnMIC(int mic_thread, const vector<GraphData *> &v_gd,  /* input paras */
				int *freq_node_label, int *freq_edge_label  /* input paras */ )
{
	int max_node_label=0, max_edge_label=0;
	/* sort labels of vertices and edges in GS by their frequency */  
    int rank_node_label[LABEL_MAX + 1], rank_edge_label[LABEL_MAX + 1];
    for (int i = 0; i <= LABEL_MAX; ++i)  
    {  
        rank_node_label[i] = i;  
        rank_edge_label[i] = i;  
    }  
    for (int i = LABEL_MAX; i > 0; --i)  
    {  
        for (int j = 0; j < i; ++j)  
        {  
            int tmp;  
            if (freq_node_label[rank_node_label[j]] < freq_node_label[rank_node_label[j + 1]])  
            {  
                tmp = rank_node_label[j];  
                rank_node_label[j] = rank_node_label[j + 1];  
                rank_node_label[j + 1] = tmp;  
            }  
            if (freq_edge_label[rank_edge_label[j]] < freq_edge_label[rank_edge_label[j + 1]])  
            {  
                tmp = rank_edge_label[j];  
                rank_edge_label[j] = rank_edge_label[j + 1];  
                rank_edge_label[j + 1] = tmp;  
            }  
        }  
    }  
  
    /* remove infrequent vertices and edges */  
    /* ralabel the remaining vertices and edges in descending frequency */    
    for (int i = 0; i <= LABEL_MAX; ++i)  
    {  
        if (freq_node_label[rank_node_label[i]] >= min_support)  
            max_node_label = i;  
        if (freq_edge_label[rank_edge_label[i]] >= min_support)  
            max_edge_label = i;  
    } 
  
    memcpy(rank_to_node_label, rank_node_label, sizeof(rank_node_label));  
    for (int i = 0; i <= LABEL_MAX; ++i)  
        rank_node_label[rank_to_node_label[i]] = i;  
    memcpy(rank_to_edge_label, rank_edge_label, sizeof(rank_edge_label));  
    for (int i = 0; i <= LABEL_MAX; ++i)  
        rank_edge_label[rank_to_edge_label[i]] = i;  
  
    for (size_t i = 0; i < v_gd.size(); ++i)  
    {  
        GraphData &gd = *v_gd[i];  
        for (size_t j = 0; j < gd.nodel.size(); ++j)  
        {  
            if (freq_node_label[gd.nodel[j]] < min_support)  
                gd.nodev[j] = false;  
            else  
                gd.nodel[j] = rank_node_label[gd.nodel[j]];  
        }  
        for (size_t j = 0; j < gd.edgel.size(); ++j)  
        {  
            if (!gd.nodev[gd.edgex[j]] || !gd.nodev[gd.edgey[j]])  
            {  
                gd.edgev[j] = false;  
                continue;  
            }  
            if (freq_edge_label[gd.edgel[j]] < min_support)  
                gd.edgev[j] = false;  
            else  
                gd.edgel[j] = rank_edge_label[gd.edgel[j]];  
        }  
  
        /* re-map vertex index */  
        map<int, int> m;  
        int cur = 0;  
        for (size_t j = 0; j < gd.nodel.size(); ++j)  
        {  
            if (!gd.nodev[j]) continue;  
            m[j] = cur++;  
        }  
        for (size_t j = 0; j < gd.edgel.size(); ++j)  
        {  
            if (!gd.edgev[j]) continue;  
            gd.edgex[j] = m[gd.edgex[j]];  
            gd.edgey[j] = m[gd.edgey[j]];  
        }  
    }  
  
    /* build graph set */  
    nr_graph = (int)v_gd.size();  
    GS = new Graph[nr_graph];  
    for (int i = 0; i < nr_graph; ++i)  
    {  
        Graph &g = GS[i];  
        GraphData &gd = *v_gd[i];  
        for (size_t j = 0; j < gd.nodel.size(); ++j)  
            if (gd.nodev[j])  
                g.node_label.push_back(gd.nodel[j]);  
        g.edge_next = new vector<int>[g.node_label.size()];  
        g.edge_label = new vector<int>[g.node_label.size()];  
        for (size_t j = 0; j < gd.edgel.size(); ++j)  
        {  
            if (!gd.edgev[j]) continue;  
            g.edge_label[gd.edgex[j]].push_back(gd.edgel[j]);  
            g.edge_label[gd.edgey[j]].push_back(gd.edgel[j]);  
            g.edge_next[gd.edgex[j]].push_back(gd.edgey[j]);  
            g.edge_next[gd.edgey[j]].push_back(gd.edgex[j]);  
        }  
    }  

    /* enumerate all frequent 1-edge graphs in GS */  
    EF.init(max_node_label, max_edge_label);
	
	#pragma omp parallel for num_threads(mic_thread) schedule(dynamic)	
    for (int x = 0; x <= max_node_label; ++x)  
    {  
        for (int a = 0; a <= max_edge_label; ++a)  
        {  
            for (int y = x; y <= max_node_label; ++y)  
            {  
                int count = 0;  
                for (int i = 0; i < nr_graph; ++i)  
                    if (GS[i].hasEdge(x, a, y))  
                        count++;  
                EF(x, a, y) = count;  
                EF(y, a, x) = count;  
            }  
        }  
    } 
}

void prepareDataOnMIC(char *input , int mic_thread , float min_Support_Rate)
{
	int freq_node_label[LABEL_MAX + 1], freq_edge_label[LABEL_MAX + 1];
	//we cannot memset the array in the function, 
	//because the array is transfered into the function as a copying pointer,
	//the sizeof will not work anymore
	memset(freq_node_label, 0, sizeof(freq_node_label));  
    memset(freq_edge_label, 0, sizeof(freq_edge_label));	
    vector<GraphData *> v_gd;
	/* load data from file */
	char Res[128];
	sprintf(Res,"/home/root/%s",input);  //add path prefix
	load_data(Res, v_gd, freq_node_label, freq_edge_label);
	
	min_support = (int) (v_gd.size() * min_Support_Rate);
	//min_support = v_gd.size() * 0.1;
    if (min_support < 1)  
        min_support = 1;
	
	constructGraphSetOnMIC(mic_thread, v_gd, freq_node_label, freq_edge_label);	
}

void write_resultsOnMIC(char *output)
{
	FILE *fp;
	char Res[128];
	/*
	char *file_path_getcwd;
    file_path_getcwd=(char *)malloc(80);
    getcwd(file_path_getcwd,80);
    printf("%s\n",file_path_getcwd);
	*/	
	//current path : "/var/volatile/tmp/coi_procs/1/xxx"	
	sprintf(Res,"../../../%s.txt",output);  //the data can only be written in tmp directory
    fp = fopen(Res, "w");  
    assert(fp);  
    for (int i = 0; i < (int)S.size(); ++i)  
    {  
        const Graph &g = *S[i];  
		fprintf(fp, "graph number: %d\n", i);  
        fprintf(fp, "frequent number: %d\n", (int)g.gs.size());  
		fprintf(fp, "oringal dataset graph list where it appeared:\n");  
        for (size_t j = 0; j < g.gs.size(); ++j)  
            fprintf(fp, " %d", g.gs[j]); 
        fprintf(fp, "\n"); 
		fprintf(fp, "node list:\n"); 		
        for (size_t j = 0; j < g.node_label.size(); ++j)  
            fprintf(fp, "v %d %d\n", (int)j, rank_to_node_label[g.node_label[j]]);  
        if (g.node_label.size() < 2)  
        {  
            fprintf(fp, "\n");  
            continue;  
        } 
		
		fprintf(fp, "edge list:\n");	
        for (int j = 0; j < (int)g.node_label.size(); ++j)  
        {  
            for (size_t k = 0; k < g.edge_label[j].size(); ++k)  
                if (j < g.edge_next[j][k])  
                    fprintf(fp, "e %d %d %d\n", j, g.edge_next[j][k], rank_to_edge_label[g.edge_label[j][k]]);  
        }  
        fprintf(fp, "\n");  
    }  
    fclose(fp);
}

