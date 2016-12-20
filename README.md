# cmFSM
Frequent subgraph mining Tool in a parallel manner of cpu/mic cooperation 


### Download and setup

```shell
#git clone
git clone https://github.com/ysycloud/cmFSM.git
cd cmFSM
#make
make all
```

### Input Format:

(line of comment in the format: # (comment))

(transaction header in the format: t (comment))

(nodes and labels for each of the nodes in the format: v (node id) (label number))

(edges in the format: e (source node id) (target node id) (label number))

An example is:
```shell
# start
t # 1        (transactionid is 1)
v 0 1        (label of node 1 is 1)
v 1 2        (label of node 2 is 2)
e 0 1 3      (edge from node 1 to node 2, with label 3)
```

### Usage:
```shell
Usage:	paraGSpan [options]
	general options:
		-s --support: The minimal value of support rates(0-1).
		-d --division: the division strategy among processes[default:2].
			0: equality; 1: single; 2: increment; 3: circle; 4：dynamic
		-t --thread: the number of threads in per process_num[default:1].
		-m --micthread: the number of threads in per mic[default:0 not use mic]
	input/output options:
		-i --input: input file of graph set information.
		-o --output: the output file of frequent subgraph results.
```

As an example file Chemical_340 is provided. To run paraGSpan on this file,
type:
```shell
mpirun -n 1 ./paraGSpan -i data/Chemical_340 -o output -s 0.1 -d 2 -t 2
```

### Output Format:
the output at command line of the example above:  
```shell
loading file time: 0.000000 seconds
340 graphs with minSup = 34
single_vertex_graph_num: 16
prework for mining frequenct subgraph spend: 0.020000 seconds
single_edge_graph_num: 23
the division strategy among processes : increment
Mining frequenct subgraph spend: 0.410000 seconds
Found 860 frequent subgraphs
output mining results time: 0.010000 seconds
```

the results written in file of the example above:
```shell 
graph number: 0
frequent number: 321
oringal dataset graph list where it appeared:
 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 19 20 21 22 23 24 25 ......
node list:
v 0 1

graph number: 1
frequent number: 239
oringal dataset graph list where it appeared:
 3 4 5 6 8 9 10 11 12 15 16 19 20 22 23 25 26 27 28 29 30 32 33 ......
 node list:
v 0 9

......

graph number: 859
frequent number: 44
oringal dataset graph list where it appeared:
 2 3 7 8 20 25 33 35 46 49 62 63 68 93 102 103 115 119 120 123 136 ......
node list:
v 0 7
v 1 7
edge list:
e 0 1 1
```