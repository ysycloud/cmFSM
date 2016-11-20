#include "Traveler.h"
 
void Traveler::DFS(vector<int> &v, int c, int next)  
{  
	if (c >= (int)s.size())  
		return;  
	vector<int> bak; 
	while (!v.empty())  
	{  
		bool flag = false;
		int x = v.back();  
		for (size_t i = 0; i < g.edge_next[x].size(); ++i)  
		{  
			int y = g.edge_next[x][i];  
			if (f[x][y]) continue;  
			flag = true;  
			if (g2s[y] < 0)  
			{  
				Edge e(g2s[x], next, g.node_label[x], g.edge_label[x][i], g.node_label[y]); 
				// if the edge of s is smaller, continue to use next node
				// if the edge of s is biger, return false directly
				// only if they are equal, continue to next level of DFS
				if (*s[c] < e) continue;  
				if (e < *s[c])  
				{  
					is_min = false;  
					return;  
				}  
				g2s[y] = next;  
				v.push_back(y);  
				f[x][y] = true;  
				f[y][x] = true;  
				DFS(v, c + 1, next + 1);  
				if (!is_min) return;  
				f[x][y] = false;  
				f[y][x] = false;  
				v.pop_back();  
				g2s[y] = -1;  
			}  
			else  
			{  
				Edge e(g2s[x], g2s[y], g.node_label[x], g.edge_label[x][i], g.node_label[y]);  
				if (*s[c] < e) continue;  
				if (e < *s[c])  
				{  
					is_min = false; 
					return;  
				}  
				f[x][y] = true;  
				f[y][x] = true;  
				DFS(v, c + 1, next);  
				if (!is_min) return;  
				f[x][y] = false;  
				f[y][x] = false;  
			}  
		}  
		if (flag) break;  
		bak.push_back(v.back());  
		v.pop_back();  
	}  
	while (!bak.empty())  
	{  
		v.push_back(bak.back());  
		bak.pop_back();  
	}  
}  
  
  
void Traveler::travel()  
{  
	g2s.resize(0);  
	g2s.resize(g.node_label.size(), -1);  
	f.resize(g.node_label.size());  
	for (size_t i = 0; i < g.node_label.size(); ++i)  
	{  
		f[i].resize(0);  
		f[i].resize(g.node_label.size(), false);  
	}  

	//from every node to DFS to judge whether s is minimum
	for (size_t i = 0; i < g.node_label.size(); ++i)  
	{  
		int x = g.node_label[i];  
		//the first x biger than first edge's x of s, current DFS from this node won't be minimum, continue using next node to DFS and judge
		if (x > s[0]->x) continue;  
//		assert(x == s[0]->x);  
		vector<int> v(1, i);  
		g2s[i] = 0;  
		DFS(v, 0, 1);  
		g2s[i] = -1;  
	}  
}  