#include "SubTraveler.h"
    
void SubTraveler::DFS(int p)  
{  
	int x, y;  

	if (p >= (int)s.size())    //match to the size, success. to test the posible extends child edge.
	{  
		/* grow from rightmost path, pruning stage (3) */  
		int at = 0;  
		while (at >= 0)  
		{  
			x = s2g[at];  
			for (size_t i = 0; i < g.edge_label[x].size(); ++i)  
			{  
				y = g.edge_next[x][i];  
				if (f[x][y])  
					continue;  
				int ix = g2s[x], iy = g2s[y];  
				assert(ix == at);  
				Edge ne(0, 1, g.node_label[x], g.edge_label[x][i], g.node_label[y]);  
				/* pruning */  
				if (EF(ne.x, ne.a, ne.y) < min_support)  
					continue;  
				/* pruning stage (1) */  
				if (ne < *s[0])  
					continue;  
				/* pruning stage (2) */  
				if (iy >= 0)  //iy is already exist, judge if it is rightmost node to backward-extend
				{  
					if (ix <= iy)  
						continue;  
					bool flag = false;  
					for (size_t j = 0; j < g.edge_label[y].size(); ++j)  
					{  
						int z = g.edge_next[y][j], w = g.edge_label[y][j];  
						if (!f[y][z])  
							continue;  
						if (w > ne.a || (w == ne.a && g.node_label[z] > ne.x))  
						{  
							flag = true;  
							break;  
						}  
					}  
					if (flag)  
						continue;  
				}  
				else  //iy is not exist, backward-extend
					iy = next;  
				ne.ix = ix;  
				ne.iy = iy;  
				c.insert(ne);  

			}  
			at = rm[at];  
		}  
		return;  
	}  

	//not match success yet, find the current level, match it and carry out next level match work recursively。
	const Edge &e = *s[p];  
	x = s2g[e.ix];  
	assert(g.node_label[x] == e.x);  
	
	for (size_t i = 0; i < g.edge_label[x].size(); ++i)  //test every adjacent edges 
	{  
		if (g.edge_label[x][i] != e.a)  
			continue;  
		y = g.edge_next[x][i];  
		if (f[x][y])  
			continue;  
		if (g.node_label[y] != e.y)  
			continue;                   //label is matched, then to see the id
		 
		if (s2g[e.iy] < 0 && g2s[y] < 0)  //if id is not seted, set it and match successed for current level. continue next level
		{  
			s2g[e.iy] = y;  
			g2s[y] = e.iy;  
			f[x][y] = true;  
			f[y][x] = true;  
			DFS(p + 1);  
			f[x][y] = false;  
			f[y][x] = false;  
			g2s[y] = -1;  
			s2g[e.iy] = -1;  
		}  
		else  //if id is already seted, match it to decide whether successed and continue next level
		{  
			if (y != s2g[e.iy])  
				continue;  
			if (e.iy != g2s[y])  
				continue;  
			f[x][y] = true;  
			f[y][x] = true;  
			DFS(p + 1);  
			f[x][y] = false;  
			f[y][x] = false;  
		}  
	}  
}  
   
void SubTraveler::travel(int _next)  
{  
	next = _next;  
	s2g.resize(0);  
	s2g.resize(next, -1);  
	g2s.resize(0);  
	g2s.resize(g.node_label.size(), -1);  
	f.resize(g.node_label.size());  
	for (size_t i = 0; i < g.node_label.size(); ++i)  
	{  
		f[i].resize(0);  
		f[i].resize(g.node_label.size(), false);  
	}  

	/* find rightmost path */  
	rm.resize(0);  
	rm.resize(next, -1);  
	for (size_t i = 0; i < s.size(); ++i)  
		if (s[i]->iy > s[i]->ix && s[i]->iy > rm[s[i]->ix])  
			rm[s[i]->ix] = s[i]->iy;  

	//match every subgraph and explore the posible extends(travel every node for original graph g)
	for (size_t i = 0; i < g.node_label.size(); ++i)  
	{  
		if (g.node_label[i] != s[0]->x)  //match the first node
			continue;  
		s2g[s[0]->ix] = (int)i;  
		g2s[i] = s[0]->ix;  
		DFS(0);  
		g2s[i] = -1;  
		s2g[s[0]->ix] = -1;  
	}  
}