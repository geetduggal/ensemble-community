// Simple program for annealing via the Massen-Doye scheme 
// Geet Duggal

#include "ht.h"
#include "avl.h"
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <signal.h>

char usage[] = "Usage: ./massendoye graph.edg\n\
  graph.edg: each edge is \"%%u\\t%%u\\n\"\n\
  outQ.txt will contain a history of modularity values to plot\n\
  out.par will contain the partition at maximal modularity before close\n";

typedef struct {
  avl info;     // AVL tree information
  uint32_t v;   // An ID for AVL tree comparison
  uint32_t k;   // Degree
  uint32_t c;   // Community ID
} node_t;

typedef struct {
  uint32_t v;
  uint32_t w; 
  uint32_t n;
} edge_t;

typedef struct {
  uint32_t a;             // Actually storing 2m*a and dividing later
  uint32_t num_nodes;     // The number of nodes for the community
  uint32_t g;             // Group ID
  avl_tree nodes;         // Each community has an AVL tree of associated nodes
} comm_t;

// Globals
FILE * outpar, * outQ;

uint32_t niter, liter, last_max,  m, n, nc, ranv, ranj;
double Q, maxQ, minQ, T;

ht_t * edges;
node_t * nodes;
comm_t * comms;
uint32_t * ks;

edge_t e;
uint8_t * found, found_modules;

const gsl_rng_type * Tr;
gsl_rng * r;

void catchsig(int signum) {
  signal(SIGINT, catchsig);
  found_modules = 1;
}

int cmpint(void * a, void * b) {
  return ((node_t *)a)->v - ((node_t *)b)->v;
}

// Calculate Ejv, traversing the AVL and hashing for edge existence 
void calcE(avl * a, int v, uint32_t * E) {
  if (a == 0)
    return;

  if (a->right) 
    calcE(a->right, v, E);
   
  // Create the current edge to hash against (obeying v < w in the table)
  if ( ((node_t *)a)->v < v) {
    e.v = ((node_t *)a)->v; 
    e.w = v;
  }
  else {
    e.v = v;
    e.w = ((node_t *)a)->v; 
  }

  // If this edge exists, increment Ejv
  found = ht_find(edges, (uint8_t *)&e);
  if (ht_occupied(edges, found))
    (*E)+= ((edge_t *)found)->n;
 
  if (a->left) 
    calcE(a->left, v, E);
}

double calcdQ(uint32_t v, uint32_t i) {
  uint32_t Eiv, Ejv;

  // Determine the number of new within community edges for both communities 
  Eiv = 0;
  calcE(comms[nodes[v].c].nodes.root, v, &Eiv);
  Ejv = 0;
  calcE(comms[i].nodes.root, v, &Ejv);

  // Efficient dQ calculation 
  return 1/(double)m*( ((double)Ejv - Eiv) - (double)nodes[v].k/(2*m)*
      (nodes[v].k + (double)comms[i].a - comms[nodes[v].c].a) );
}



void move_node(uint32_t v, uint32_t i) {
  avl_remove(&comms[nodes[v].c].nodes, (avl *)&nodes[v]);
  comms[nodes[v].c].num_nodes--; 

  avl_insert(&comms[i].nodes, (avl *)&nodes[v]);
  comms[i].num_nodes++;

  // If we create a new community upon acceptance, increment nc
  if (comms[i].num_nodes == 1)
    nc++;

  // Moving a node from size-1 community to any other one kills a community
  if (comms[nodes[v].c].num_nodes == 0)
    nc--;

  // Perform our state updates for the next Monte Carlo step
  comms[nodes[v].c].a -= nodes[v].k;
  comms[i].a += nodes[v].k;
  nodes[v].c = i;
}

uint32_t commidx (uint32_t j) {
  uint32_t i = 0;
  // If we're creating a new community find the first free slot (no free list)
  if (j == nc) {
    for (i=0; i<n; i++) {
      if (comms[i].num_nodes == 0)
        break;
    }
    j = i;
  }
  // Otherwise update the community index for the ranj'th non-empty community 
  else {
    j++;
    do {
      if (comms[i].num_nodes > 0)
        j--;

      i++;
    } while (j > 0 && i < n);
    j = i-1;
  }
  return j;
}

// A Metropolis Monte Carlo step
void monte_step() {
  double dQ, Qu, Qv, ranp;

  // Pick a node and a community at random
  ranv = gsl_rng_uniform_int(r, n);
  if (comms[nodes[ranv].c].num_nodes == 1)
    ranj = gsl_rng_uniform_int(r, nc);
  else
    ranj = gsl_rng_uniform_int(r, nc+1);

  ranj = commidx(ranj);

  // Transferring the node to its original community does nothing
  if (nodes[ranv].c == ranj)
    return;

  // Determine the number of new within community edges for both communities 
  dQ = calcdQ(ranv, ranj);

  // Populate Qu and Qv (normalized)
  Qu = (Q + 1)/2;
  Qv = (Q+dQ + 1)/2;

  ranp = gsl_rng_uniform(r);

  // If Metropolis-accepted, we change state
  if ( Qv > Qu || ranp < exp(dQ/T) ) {
    // Update Q
    Q += dQ;
    outQ = fopen("outQ.txt", "a");
    fprintf(outQ, "%g\n", Q);
    fclose(outQ);

    // Perform the move of ranv from the current community to ranj
    move_node(ranv, ranj);

    // Update if we have a new maximum Q and print it out to the log
    if (Q > maxQ) {
      maxQ = Q;
      last_max = niter;
    }
  }

}

int main (int argc, char ** argv) {
  FILE * infile; 
  uint32_t v, tmpv, sumsqk;
  ht_t * nodecnt;

  if (argc != 2 && argc != 3) {
    fprintf(stderr, usage);
    return 1;
  }

  infile = fopen(argv[1], "r");

  // Count the number of edges
  m = 0;
  while(fscanf(infile, "%u\t%u\n", &e.v, &e.w) == 2)
    m++; 

  // Store the edges in a newly allocated hash table
  rewind(infile);
  edges = ht_alloc(log(m)/log(2)+3, 2*sizeof(uint32_t),  sizeof(edge_t));
  nodecnt = ht_alloc(log(m)/log(2)+3, sizeof(uint32_t),  sizeof(uint32_t));
  while(fscanf(infile, "%u\t%u\n", &e.v, &e.w) == 2) { 
    e.n = 0;
    // Always insert an edge (v,w) such that v < w
    if (e.w < e.v) { 
      tmpv = e.v;
      e.v = e.w;
      e.w = tmpv;
    }

    // Insert in hash if not already found
    found = ht_find(edges, (uint8_t *)&e);  
    if (!ht_occupied(edges, found))
      found = ht_insert(edges, found, (uint8_t *)&e);

    ((edge_t *)found)->n++;

    // Start counting the nodes by hashing them
    found = ht_find(nodecnt, (uint8_t *)&e.v);
    if (!ht_occupied(nodecnt, found))
      ht_insert(nodecnt, found, (uint8_t *)&e.v);

    found = ht_find(nodecnt, (uint8_t *)&e.w);
    if (!ht_occupied(nodecnt, found))
      ht_insert(nodecnt, found, (uint8_t *)&e.w);
  }

  // Grab the node size and allocate node and community arrays
  nc = n = nodecnt->num_entries;
  ht_free(nodecnt);
  nodes = calloc(n, sizeof(node_t));
  ks = calloc(n, sizeof(uint32_t));
  comms = calloc(n, sizeof(comm_t));

  // Initialize node values
  for (v=0; v<n; v++) {
    nodes[v].v = v;
    nodes[v].k = 0;
    ks[v] = 0;
    nodes[v].c = v;
  }

  // Populate node degree
  rewind(infile);
  while(fscanf(infile, "%u\t%u\n", &e.v, &e.w) == 2) {
    nodes[e.v].k++;
    nodes[e.w].k++;
    ks[e.v]++;
    ks[e.w]++;
  }
  
  fclose(infile);

  // Initialize community values and Q including AVL tree inits
  sumsqk = Q = 0;
  for (v=0; v<n; v++) {
    comms[v].a = nodes[v].k;
    comms[v].nodes.compar = cmpint;
    comms[v].nodes.root = 0;

    // Each node should be placed in its own community
    avl_insert(&comms[v].nodes, (avl *)&nodes[v]);
    comms[v].num_nodes = 1;

    sumsqk += nodes[v].k*nodes[v].k;
  }
  minQ = maxQ = Q = -(double)sumsqk/(4*m*m);

  // Set up random number generator
  gsl_rng_env_setup();
  Tr = gsl_rng_default;
  r = gsl_rng_alloc(Tr);

  last_max = 0;
  unlink("outQ.txt");
  niter = 0;

  T = .3;
  found_modules = 0;
  signal(SIGINT, catchsig);
  if (argv[2]) {
    T = atof(argv[2]);
  }

  // MCMC for the modules
  while(1) {
    monte_step();

    if (((niter % (n)) == 0)) {
      uint32_t npars = 100;
      char dirname[256];
      char fname[256];
      double Qavg;

      //if (argv[2]) {
      while (npars) {
        if (npars == 100) Qavg = Q;
        if (niter % (5*n) == 0) { 
	  if(argv[2]) {
	    sprintf(dirname, "ensemble/T%f", T);
	    mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	    sprintf(fname, "ensemble/T%f/coarsened-T-%f-%03u.par", T, T, npars);
	    outpar = fopen(fname, "w");
	    { uint32_t i;
	      for (i=0; i<n; i++) 
		fprintf(outpar, "%u\t%u\n", nodes[i].v, nodes[i].c);
	    }
	    fclose(outpar);
	    sprintf(fname, "ensemble/coarsened-T-%f-VEC.mod", T);
	    outpar = fopen(fname, "a");
	    fprintf(outpar, "%g\t%g\n", T, Q);
	    fclose(outpar);
	  }
          else {
            Qavg += Q; 
          }
          npars--;
        }
        monte_step();
        niter++;
      }
      //}
      /*
      if (argv[2]) {
        printf("%f\n", Qavg/100); fflush(stdout);
      }
      */
      fprintf(stdout, "%8u\t%4u\t%5e\t%5e\t%5e\n", 
              niter, nc, T, Q, Qavg/100);
      fflush(stdout);

      //T -= .005;
      T *= .96;
    }

    if (found_modules || T<= 1e-22) {
	    outpar = fopen("out.par", "w");
	    { uint32_t i;
	      for (i=0; i<n; i++) 
		fprintf(outpar, "%u\t%u\n", nodes[i].v, nodes[i].c);
	    }
	    fclose(outpar);
      exit(0);
    }

    niter++;
  }

  ht_free(edges);
  free(nodes); 
  free(comms); 
  free(r);
  return 0;
}
