#include "partition.h"

#define K 5
#define MIN_IMPROVEMENT 0.005

unsigned long randomPartition(adjlist *g,unsigned long *lab);
unsigned long louvainPartition(adjlist *g, unsigned long *lab);

partitionFunction choose_partition(char *c){
	printf("Chosen partition algorithm: ");
	if (strcmp(c,"0")==0){
		printf("Random partition\n");
		return randomPartition;
	}
	//if (strcmp(c,"1")==0){
	//	printf("Greedy Sparsest Cut\n");
	//	return greedySparsestcut;
	//}

	if (strcmp(c,"2")==0){
		printf("Louvain partition\n");
		return louvainPartition;
	}
	printf("unknown\n");
	exit(1);
}

// -----------------------------------------------------------
// START Random utility functions

//generating n random labels (boolean values)
unsigned long randomPartition(adjlist *g,unsigned long *lab){
	unsigned long i,n=g->n;
	unsigned long nlab=(K>n)?n:K;

	//random side for each node:
	for (i=0;i<n;i++){
		lab[i]=rand()%nlab;
	}

	return nlab;
}


// END Random utility functions
// -----------------------------------------------------------


// -----------------------------------------------------------
// START Louvain utility functions

/*
int myCompare (const void * a, const void * b, void * array2) {
  long diff = ((unsigned long *)array2)[*(unsigned long *)a] - ((unsigned long *)array2)[*(unsigned *)b];
  int res = (0 < diff) - (diff < 0);
  return  res;
}

unsigned long * mySort(unsigned long *part, unsigned long size) {
  unsigned long *nodes = (unsigned long *)malloc(size * sizeof(unsigned long));
  for (unsigned long i = 0; i < size; i++) {
    nodes[i]=i;
  }

  qsort_r(nodes, size, sizeof(unsigned long), myCompare, (void *)part);

  return nodes;
}
*/



void freeLouvainPartition(partition *p) {
  free(p->in);
  free(p->tot);
  free(p->neighCommWeights);
  free(p->neighCommPos);
  free(p->node2Community);
  free(p);
}


partition *createLouvainPartition(adjlist *g) {
  partition *p = (partition *)malloc(sizeof(partition));

  p->size = g->n;

  p->node2Community = (unsigned long *)malloc(p->size * sizeof(unsigned long));
  p->in = (long double *)malloc(p->size * sizeof(long double));
  p->tot = (long double *)malloc(p->size * sizeof(long double));

  p->neighCommWeights = (long double *)malloc(p->size * sizeof(long double));
  p->neighCommPos = (unsigned long *)malloc(p->size * sizeof(unsigned long));
  p->neighCommNb = 0;

  for (unsigned long node = 0; node < p->size; node++) {
    p->node2Community[node] = node;
    p->in[node]  = selfloopWeighted(g, node);
    p->tot[node] = degreeWeighted(g, node);
    p->neighCommWeights[node] = -1;
    p->neighCommPos[node] = 0;
  }

  return p;
}

long double modularity(partition *p, adjlist *g) {
  long double q  = 0.0L;
  long double m2 = g->totalWeight;

  for (unsigned long node = 0; node < p->size; node++) {
    if (p->tot[node] > 0.0L)
      q += p->in[node] - (p->tot[node] * p->tot[node]) / m2;
  }

  return q / m2;
}

void neighCommunitiesInit(partition *p) {
  for (unsigned long i = 0; i < p->neighCommNb; i++) {
    p->neighCommWeights[p->neighCommPos[i]] = -1;
  }
  p->neighCommNb = 0;
}

/*
Computes the set of neighbor communities of a given node (excluding self-loops)
*/
void neighCommunities(partition *p, adjlist *g, unsigned long node) {
  p->neighCommPos[0] = p->node2Community[node];
  p->neighCommWeights[p->neighCommPos[0]] = 0.;
  p->neighCommNb = 1;

  // for all neighbors of node, add weight to the corresponding community
  for (unsigned long long i = g->cd[node]; i < g->cd[node + 1]; i++) {
    unsigned long neigh  = g->adj[i];
    unsigned long neighComm = p->node2Community[neigh];
    long double neighW = (g->weights == NULL)?1.0:g->weights[i];

    // if not a self-loop
    if (neigh != node) {
      // if community is new (weight == -1)
      if (p->neighCommWeights[neighComm] == -1) {
	p->neighCommPos[p->neighCommNb] = neighComm;
	p->neighCommWeights[neighComm] = 0.;
	p->neighCommNb++;
      }
      p->neighCommWeights[neighComm] += neighW;
    }
  }
}

/*
Same behavior as neighCommunities except:
- self loop are counted
- data structure if not reinitialised
*/
void neighCommunitiesAll(partition *p, adjlist *g, unsigned long node) {
  for (unsigned long long i = g->cd[node]; i < g->cd[node + 1]; i++) {
    unsigned long neigh  = g->adj[i];
    unsigned long neighComm = p->node2Community[neigh];
    long double neighW = (g->weights == NULL)?1.0:g->weights[i];
    
    // if community is new
    if (p->neighCommWeights[neighComm] == -1) {
      p->neighCommNb++;
      p->neighCommPos[p->neighCommNb] = neighComm;
      p->neighCommWeights[neighComm] = 0.;      
    }
    p->neighCommWeights[neighComm] += neighW;
  }
}


unsigned long updatePartition(partition *p, unsigned long *part, unsigned long size) {
  // Renumber the communities in p
  unsigned long *renumber = (unsigned long *)calloc(p->size, sizeof(unsigned long));
  unsigned long last = 1;
  for (unsigned node = 0; node < p->size; node++) {
    if (renumber[p->node2Community[node]] == 0) {
      renumber[p->node2Community[node]] = last++;
    }
  }

  // Update part with the renumbered communities in p
  for (unsigned long i = 0; i < size; i++) {
    part[i] = renumber[p->node2Community[part[i]]] - 1;
  }

  free(renumber);
  return last-1;
}


// Compute one pass of Louvain and returns the improvement
long double louvainOneLevel(partition *p, adjlist *g) {
  unsigned long nbMoves;
  long double startModularity = modularity(p, g);
  long double newModularity = startModularity;
  long double curModularity;

  // generate a random order for nodes' movements
  /*
  unsigned long *randomOrder = (unsigned long *)malloc(p->size * sizeof(unsigned long));
  for (unsigned long i = 0; i < p->size; i++)
    randomOrder[i] = i;
  for (unsigned long i = 0; i < p->size - 1; i++) {
    unsigned long randPos = rand()%(p->size - i) + i;
    unsigned long  tmp = randomOrder[i];
    randomOrder[i] = randomOrder[randPos];
    randomOrder[randPos] = tmp;
  }
  */
  
  // repeat while 
  //   there are some nodes moving
  //   or there is an improvement of quality greater than a given epsilon 
  do {
    curModularity = newModularity;
    nbMoves = 0;

    // for each node:
    //   remove the node from its community
    //   compute the gain for its insertion in all neighboring communities
    //   insert it in the best community with the highest gain
    for (unsigned long i = 0; i < g->n; i++) {
      unsigned long node = i;//randomOrder[nodeTmp];
      unsigned long oldComm = p->node2Community[node];
      long double degreeW = degreeWeighted(g, node);


      // computation of all neighboring communities of current node
      neighCommunitiesInit(p);
      neighCommunities(p, g, node);

      // remove node from its current community
      removeNode(p, g, node, oldComm, p->neighCommWeights[oldComm]);

      // compute the gain for all neighboring communities
      // default choice is the former community
      unsigned long bestComm = oldComm;
      long double bestCommW  = 0.0L;
      long double bestGain = 0.0L;
      for (unsigned long i = 0; i < p->neighCommNb; i++) {
	unsigned long newComm = p->neighCommPos[i];
	long double newGain = gain(p, g, newComm, p->neighCommWeights[newComm], degreeW);

	if (newGain > bestGain) {
	  bestComm = newComm;
	  bestCommW = p->neighCommWeights[newComm];
	  bestGain = newGain;
	}
      }

      // insert node in the nearest community
      insertNode(p, g, node, bestComm, bestCommW);

      if (bestComm != oldComm) {
	nbMoves++;
      }
    }

    newModularity = modularity(p, g);
  } while (nbMoves>0 && 
    	   newModularity - curModularity > MIN_IMPROVEMENT);
  
  //  free(randomOrder);

  return newModularity - startModularity;
}

unsigned long louvainPartition(adjlist *g, unsigned long *lab) {
  // Initialize partition with trivial communities
  for (unsigned long i = 0; i < g->n; i++) {
    lab[i] = i;
  }
  
  // Execution of Louvain method
  //  while(1) {
    partition *gp = createLouvainPartition(g);
    long double improvement = louvainOneLevel(gp, g);
    //    if (improvement < MIN_IMPROVEMENT) {
    //      freeLouvainPartition(gp);
    //      break;
    //    }
    
    unsigned long n = updatePartition(gp, lab, g->n);
    

    //adjlist *g2 = partition2Graph(gp, g);
    //    if (g->n < size) {
    //      freeGraph(g);
    //    }
    freeLouvainPartition(gp);
    //    g = g2;
    //  }

  return n;
}



// END Louvain utility functions
// -----------------------------------------------------------
