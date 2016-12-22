#include "MicFunc.h"

__attribute__((target(mic))) void funcheck(int i)
{
	#ifdef __MIC__
		printf("Index on MIC:%d\n",i);
	#else
		printf("Index on CPU:%d\n",i);
	#endif
}

void cmsingleEdgeGraphMining(int my_rank, const Edge &e, vector<Edge> &single_edge_graph, int thread_num, int begin, int end, int mic_thread)
{
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
	//first expansion to get the children set(two edges frequent graph)
	one_edge_expansion(gc, 2, two_edges_child_gcs, two_edges_nexts);
	
	if(my_rank == 0)
		printf("ALL LEN:%d\n", two_edges_child_gcs.size());
	
	//submining per level on MIC
	vector<Graph *> MIC_S;
	freqGraphMiningfrom2edgesOnMIC(my_rank, two_edges_child_gcs, two_edges_nexts, mic_thread, MIC_S);
	
	//submining per level on CPU
	freqGraphMiningfrom2edgesOnCPU(my_rank, two_edges_child_gcs, two_edges_nexts, thread_num);
	
	S.insert(S.end(), MIC_S.begin(), MIC_S.end());
}

void freqGraphMiningfrom2edgesOnCPU(int my_rank, vector<GraphCode> two_edges_child_gcs, vector<int> two_edges_nexts, int thread_num)
{

	vector<GraphCode> current_global_child_gcs;
	vector<int> current_global_nexts;
	vector<GraphCode> tmp_global_child_gcs;
	vector<int> tmp_global_nexts;
	
	int cpu_len =  two_edges_child_gcs.size()/2;
	
	current_global_child_gcs.swap(two_edges_child_gcs);
	current_global_nexts.swap(two_edges_nexts);
	//submining per level in cpu
	while(cpu_len != 0)
	{
		#pragma omp parallel num_threads(thread_num)
		{		
			
			vector<GraphCode> local_child_gcs;
			vector<int> local_nexts;
			vector<GraphCode> current_local_child_gcs;
			vector<int> current_local_nexts;
			
			#pragma omp for schedule(dynamic)
			for(size_t i=0; i<cpu_len; i++)
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
		cpu_len = current_global_child_gcs.size();
	}					
}

void freqGraphMiningfrom2edgesOnMIC(int my_rank, vector<GraphCode> two_edges_child_gcs, vector<int> two_edges_nexts, int mic_thread, vector<Graph *> &MIC_S)
{

	vector<GraphCode> current_global_child_gcs;
	vector<int> current_global_nexts;
	vector<GraphCode> tmp_global_child_gcs;
	vector<int> tmp_global_nexts;
	
	int mic_begin = two_edges_child_gcs.size()/2;
	int mic_end = two_edges_child_gcs.size();
	int mic_len = mic_end - mic_begin;
	
	if(my_rank == 0)
		printf("MIC LEN:%d\n", mic_len);
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
			for(size_t i=mic_begin; i<mic_end; i++)
			{
				//carry out current child's one edge expansion
				one_edge_expansion_mic(current_global_child_gcs[i], current_global_nexts[i], local_child_gcs, local_nexts, MIC_S);
					
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
		mic_begin = 0;
		mic_end = mic_len;
	}					
}

void one_edge_expansion_mic(GraphCode &gc, int next, vector<GraphCode> &child_gcs, vector<int> &nexts, vector<Graph *> &MIC_S)
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
