#ifndef __SIMILARITY__
#define __SIMILARITY__

typedef struct ProtNei //list of nearest neighbour of a protein "id"
{
  long id;
  VLong nei; //list of neighbour that is in core, store the core equivalence class numbers
  VLong level;//distance of nei to id 
} ProtNei;

void allocProtNei(ProtNei *&p)
{
  p = (ProtNei*) calloc(1, sizeof(ProtNei));
}

void allocProtNei(long id, ProtNei *&p)
{
  allocProtNei(p);
  p->id = id;
}

void deleteProtNei(ProtNei *&p)
{
  free(p);
}

void add_nei(ProtNei *&p, long p1, long lev)
{
  if (p == NULL) allocProtNei(p);
  if (isin(p1, p->nei) >=0) return;
  p->nei.push_back(p1);
  p->level.push_back(lev);
}


class CSimilarity
{
  ProtNei **prot_A, **prot_B;
  long n_protA, n_protB;
  VLong core_protA, core_protB, mask_sim;
  double** ABSim; //similarity matrix
  CGraph *GA, *GB; //list of graphs
  int mode;
 public:
  CSimilarity(double ** simmat, CGraph *G1, CGraph* G2)
    {
      mode =0;
      startup(simmat, G1, G2);
    };
  CSimilarity(double ** simmat, CGraph *G1, CGraph* G2, VLong mask)
    {
      mode =1;
      startup(simmat, G1, G2);
    };
  ~CSimilarity(){};

  void startup(double **simmat, CGraph *G1, CGraph *G2)
  {
      ABSim = simmat;

      GA = G1;
      GB = G2;
  }
  void init( VLong prot1s,VLong prot2s, VLong coreA, VLong coreB);
  void free_memory();
  void get_neighbours();
  void get_neighbours(ProtNei **&ps, long n_prot, CGraph *G, VLong core);
  void get_neighbours(ProtNei *&ps, CGraph *G, VLong core);
  double get_sim(long i, long j); //by local index
  double get_sim_prot(long protA, long protB); //by global index
  double **get_sim_mat();
};

void CSimilarity::init(VLong prot1s, VLong prot2s, VLong coreA, VLong coreB)
{
  n_protA = prot1s.size();
  n_protB = prot2s.size();

  prot_A = (ProtNei**) calloc(n_protA, sizeof(ProtNei*));
  prot_B = (ProtNei**) calloc(n_protB, sizeof(ProtNei*));

  for(long i =0; i!= n_protA; ++i) allocProtNei(prot1s[i], prot_A[i]);
  for(long i =0; i!= n_protB; ++i) allocProtNei(prot2s[i], prot_B[i]);

  core_protA =coreA;
  core_protB =coreB;
}

void CSimilarity::free_memory()
  {
    if (prot_A != NULL)
      {
	for(long i =0; i!= n_protA; ++i) deleteProtNei(prot_A[i]);
	free(prot_A);
      }
    if (prot_B != NULL)
      {
	for(long i =0; i!= n_protB; ++i) deleteProtNei(prot_B[i]);
	free(prot_B);
      }
    n_protA = n_protB = 0;
    
  }

void CSimilarity::get_neighbours()
  {
    //    printf("%f \n", ABSim[0][883]);
    get_neighbours(prot_A, n_protA, GA, core_protA);
    printf("Finished getting neighbours of candidate proteins in the first graph\n");
    get_neighbours(prot_B, n_protB, GB, core_protB);
    printf("Finished getting neighbours of candidate proteins in the second graph\n");
  }

void CSimilarity::get_neighbours(ProtNei **&ps, long n_prot, CGraph *G, VLong core)
{
  for(long i = 0; i!= n_prot; ++i)
    {
      get_neighbours(ps[i],  G, core);
      if (i%100==0) printf("Get neighbours of %d\n", i);
    }
}

void CSimilarity::get_neighbours(ProtNei *&p,  CGraph *G, VLong core)
{
  VLong fnei, snei;
  fnei  = G->GetNeigh(p->id); // get first neighbours
  snei  = G->Get2Neigh(p->id);//get second neighbours
  for (long i = 0; i!= fnei.size(); ++i)
    for(long j =0; j!= core.size(); ++j)
      if (fnei[i] == core[j])
	add_nei(p, j, 1);
  //  printf("n_nei = %d\n", p->nei.size());
  /*
  for (long i = 0; i!= snei.size(); ++i)
    for(long j =0; j!= core.size(); ++j)
      if (snei[i] == core[j])
	add_nei(p, j, 2);
  */
 //  printf("n_nei = %d\n", p->nei.size());
}

double **CSimilarity::get_sim_mat()
{
  printf("Updating similarity matrix\n");
  double **out = (double**) malloc(n_protA*sizeof(double*));
  for(long i =0; i!= n_protA; ++i)
    {
      out[i] = (double*) malloc(n_protB * sizeof(double));
      for(long j=0; j!= n_protB; ++j)
	{
	  out[i][j] = get_sim(i, j);
	  //  if (out[i][j]>1)
	  //	  if (i%100==0 && j%100==0)
	  // printf("updated sim %d %d %f\n", i, j, out[i][j]);
	}
    }
  return out;
}


double CSimilarity::get_sim(long i, long j)
{
  ProtNei *p1, *p2;
  p1 = prot_A[i];
  p2 = prot_B[j];
  double out = ABSim[p1->id][p2->id];
  double x = 0;
  long id1, id2, k=0;

  for(long ii=0; ii!= p1->nei.size(); ++ii)
    for(long jj=0; jj!= p2->nei.size(); ++jj)
      {
	if (p1->nei[ii] == p2->nei[jj]) 
	  {
	    id1 = core_protA[p1->nei[ii]];
	    id2 = core_protB[p2->nei[jj]];
	    x = ABSim[id1][id2];
	    out += 1.0/((p1->level[ii] +1) * (p2->level[jj] +1)) *x;
	    k++;
	  }
      }

  return out;
}


double CSimilarity::get_sim_prot(long protA, long protB)
{
  long i = -1;
  long j = -1;

  for(long ii =0; ii!= n_protA; ++i)
    if (prot_A[ii]->id == protA) {i = ii; break;}
  
  for(long ii =0; ii!= n_protB; ++i)
    if (prot_B[ii]->id == protA) { j = ii; break;}
  
  if (i == -1 || j == -1) return -1;
  return  get_sim(i, j);
}

#endif
