#ifndef STRUCT_H
#define STRUCT_H

typedef struct {
	unsigned long s;
	unsigned long t;
} edge;

//edge list structure:
typedef struct {
  unsigned long n;//number of nodes
  unsigned long long e;//number of edges
  unsigned long long emax;//max number of edges
  edge *edges;//list of edges
  unsigned long long *cd;//cumulative degree cd[0]=0 length=n+1
  unsigned long *adj;//concatenated lists of neighbors of all nodes
  long double *weights;//concatenated lists of weights of neighbors of all nodes
  long double totalWeight;//total weight of the links
  unsigned long *map;//map[u]=original label of node u
} adjlist;

typedef struct {
  unsigned long n;//number of clusters
  adjlist** sg;//sg[i]=pointer to adjlist of cluster i
} clusters;


// louvain partition
typedef struct {
  // size of the partition
  unsigned long size;
  
  // community to which each node belongs
  unsigned long *node2Community;
  
  // in and tot values of each node to compute modularity 
  long double *in;
  long double *tot;

  // utility arrays to find communities adjacent to a node
  // communities are stored using three variables
  // - neighCommWeights: stores weights to communities
  // - neighCommPos: stores list of neighbor communities
  // - neighCommNb: stores the number of neighbor communities
  long double *neighCommWeights;
  unsigned long *neighCommPos;
  unsigned long neighCommNb;
} partition;

#endif
