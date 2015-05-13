#ifndef __HUNGARIAN__
#define __HUNGARIAN__
#include "similarity.h"
//#define max(a, b) (((a) > (b)) ? (a) : (b)).
//#define min(a, b) (((a) < (b)) ? (a) : (b)).
#define INF 10000000
class CHungarian
{
  long n, m, maxN, max_match; //size of group 1, size of group 2
  long *xx, *yy;            //global index of proteins in group1 and proteins in group 2
  long mode;                // =0: match nodes, = 1 match clusters, will make difference when comes to CSim. manage this later
  CSimilarity  *ProtSim;
  long *xy, *yx;   //    
  double *lx, *ly; //label of x, label of y
  bool *S, *T;     //set S, T in the algorithm 
  double *slack;
  long *slackx; //slackx[y] such a vertex, that
                // l(slackx[y]) + l(y) - w(slackx[y],y) = slack[y]
  long *prev;   //array for memorizing alternating paths
  double score;

  double **cost;
  double ratio;

public:
  CHungarian(CSimilarity *sim)
    {
      ProtSim = sim;
      mode = 2;  //match using ProtSim
      xx= yy = NULL;
      prev =       xy = yx = slackx = NULL;
      S = T = NULL;
      lx = ly = slack = NULL;
      cost = NULL;
      ratio = 0.8;
    }

  CHungarian(double **costx)
    {
      cost = costx;
      mode = 1;
      xx= yy = NULL;
      prev =  xy = yx = slackx = NULL;
      S = T = NULL;
      lx = ly = slack = NULL;

    }
  ~CHungarian()
    {
      if (xx!= NULL) free(xx);
      if (yy!= NULL) free(yy);
      if (S!= NULL ) free(S);
      if (T!= NULL ) free(T);
      if (slack != NULL) free(slack);
      if (slackx!= NULL) free(slackx);
      if (prev  != NULL) free(prev);
      free(lx); free(ly);
      if (yx != NULL) free(yx);
      if (xy != NULL) free(xy);

      lx = ly = NULL;	 
      xx = yy = xy= prev = yx = NULL;
      S = NULL; T = NULL; slack = NULL; slackx = NULL;
    }
  
  void init(long *gr1, long *gr2, long size1, long size2);
  void init(VLong gr1, VLong gr2);
  double Cost(long i, long j); //get the similarity between 2 proteins or 2 members of 2 groups
  void hungarian(long *gr1, long *gr2, long size1, long size2);
  void hungarian(char **gr1, char **gr2, long size1, long size2, CNodes* allnodes);
  void hungarian(VLong gr1, VLong gr2);
  void init_labels();
  void update_labels();
  void add_to_tree(long x, long prevx);
  void augment();              

  void print_mapping(FILE *f, CNodes *allnodes);
  void print_mapping(FILE *f);

  void print_mapping(char *fn, CNodes *allnodes);
  void print_mapping(char *fn);

  double get_score_mapping()
  {
    return score;
  };
  long *get_map(){return xy;};
};

void CHungarian::init(VLong gr1, VLong gr2)
{
  n = gr1.size(); m = gr2.size();
  maxN = max(m, n);

  if (xx == NULL) xx = (long *) calloc(maxN, sizeof(long));
  if (yy == NULL) yy = (long *) calloc(maxN, sizeof(long));

  if (n == maxN) 
    {
      for(long i =0; i!= n; ++i) xx[i] = gr1[i];
      for(long i =0; i!= m; ++i) yy[i] = gr2[i];
      for(long i =m; i!= maxN; ++i) yy[i] = -1; //add dummy nodes
    }

  if (m == maxN) 
    {
      for(long i =0; i!= m; ++i) yy[i] = gr2[i];
      for(long i =0; i!= n; ++i) xx[i] = gr1[i];
      for(long i =n; i!= maxN; ++i) xx[i] = -1;
    }
  

  xy = (long *) calloc(maxN, sizeof(long));
  yx = (long *) calloc(maxN, sizeof(long));
  prev =  (long *) calloc(maxN, sizeof(long));
  slack=  (double *) calloc(maxN, sizeof(double));
  slackx =(long *) calloc(maxN, sizeof(long));
  lx  = (double*) calloc(maxN, sizeof(double));
  ly  = (double*) calloc(maxN, sizeof(double));
  S   = (bool * ) calloc(maxN, sizeof(bool));
  T   = (bool * ) calloc(maxN, sizeof(bool));
}

void CHungarian::init(long *gr1, long *gr2, long size1, long size2)
{

  n = size1; m = size2;
  maxN = max(m, n);

  if (xx == NULL) xx = (long *) calloc(maxN, sizeof(long));
  if (yy == NULL) yy = (long *) calloc(maxN, sizeof(long));

  if (n == maxN) 
    {
      for(long i =0; i!= n; ++i) xx[i] = gr1[i];
      for(long i =0; i!= m; ++i) yy[i] = gr2[i];
      for(long i =m; i!= maxN; ++i) yy[i] = -1; //add dummy nodes
    }

  if (m == maxN) 
    {
      for(long i =0; i!= m; ++i) yy[i] = gr2[i];
      for(long i =0; i!= n; ++i) xx[i] = gr1[i];
      for(long i =n; i!= maxN; ++i) xx[i] = -1;
    }
  

  xy = (long *) calloc(maxN, sizeof(long));
  yx = (long *) calloc(maxN, sizeof(long));
  prev =  (long *) calloc(maxN, sizeof(long));
  slack=  (double *) calloc(maxN, sizeof(double));
  slackx =(long *) calloc(maxN, sizeof(long));
  lx  = (double*) calloc(maxN, sizeof(double));
  ly  = (double*) calloc(maxN, sizeof(double));
  S   = (bool * ) calloc(maxN, sizeof(bool));
  T   = (bool * ) calloc(maxN, sizeof(bool));

  //  printf("finished hungarian initiation\n");
}

double CHungarian::Cost(long x, long y)
{
  if (xx[x] ==-1 || yy[y] == -1) return 0;
  if (mode ==1)
    {
      if (x >=n || y>=m) return 0;
      return cost[x][y]; //TODO: change later for the case of assigning clusters
    }
  
  if(mode == 2)
    {
      return ProtSim->get_sim(x, y);
    }
}

void CHungarian::hungarian(char **gr1, char **gr2, long size1, long size2, CNodes *allnodes)
{
  xx = (long*) calloc(size1, sizeof(long));
  yy = (long*) calloc(size2, sizeof(long));
  
  for(long i=0; i!= size1; ++i) xx[i] = allnodes->Search(gr1[i]);
  for(long i=0; i!= size2; ++i) yy[i] = allnodes->Search(gr2[i]);
  hungarian(xx, yy, size1, size2);
}

void CHungarian::hungarian(long *gr1, long *gr2, long size1, long size2)
{
  init(gr1, gr2, size1, size2);

  score = 0;                      //weight of the optimal matching
  max_match = 0;                    //number of vertices in current matching
  for(long i =0; i!= maxN; ++i) 
    xy[i] =  yx[i] = -1;

  init_labels();                    //step 0
  augment();                        //steps 1-3

  for (long x = 0; x < maxN; x++)       //forming answer there
    score += Cost(x, xy[x]);
}

void CHungarian::hungarian(VLong gr1, VLong gr2)
{
  init(gr1, gr2);
  score = 0;                      //weight of the optimal matching
  max_match = 0;                    //number of vertices in current matching
  for(long i =0; i!= maxN; ++i) 
    xy[i] =  yx[i] = -1;
  init_labels();                    //step 0
  augment();                        //steps 1-3
  for (long x = 0; x < maxN; x++)       //forming answer there
    score += Cost(x, xy[x]);
}

void CHungarian::print_mapping(char *fn, CNodes *allnodes)
{
  FILE *f= fopen(fn, "w");
  print_mapping(f, allnodes);
  fclose(f);
}

void CHungarian::print_mapping(FILE *f, CNodes *allnodes)
{
  fprintf(f, "Mapping score %f\n", score);
  for(long x = 0; x!= maxN; ++x)
    if (xx[x] != -1 && yy[xy[x]] != -1) 
      fprintf(f, "%s\t%s\t%f\n", allnodes->GetLabel(xx[x]), allnodes->GetLabel(yy[xy[x]]), Cost(x, xy[x]));
}

void CHungarian::print_mapping(char *fn)
{
  FILE *f= fopen(fn, "w");
  print_mapping(f);
  fclose(f);
}

void CHungarian::print_mapping(FILE *f)
{
  fprintf(f, "Mapping score %f\n", score);
  for(long x = 0; x!= maxN; ++x)
    if (xx[x] != -1 && yy[xy[x]] != -1) 
      fprintf(f, "%d\t%d\t%f\n", x, xy[x], Cost(x, xy[x]));
}


void CHungarian::init_labels()
{
  for(long i =0; i!= maxN; ++i) lx[i] =ly[i] = 0.0;

  for (long x = 0; x != maxN; x++)
    {    
      for (long y = 0; y != maxN; y++)
	{
	  lx[x] = max(lx[x], Cost(x, y));
	}
    }
}

void CHungarian::add_to_tree(long x, long prevx) 
//x - current vertex,prevx - vertex from X before x in the alternating path,
//so we add edges (prevx, xy[x]), (xy[x], x)
{
  S[x] = true;                    //add x to S
  prev[x] = prevx;                //we need this when augmenting
  for (long y = 0; y != maxN; y++)    //update slacks, because we add new vertex to S
    if (lx[x] + ly[y] - Cost(x, y) < slack[y])
      {
	slack[y] = lx[x] + ly[y] - Cost(x, y);
	slackx[y] = x;
      }
}


void CHungarian::update_labels()
{
  long x, y;
  double delta = INF;             //init delta as infinity
  for (y = 0; y != maxN; ++y)            //calculate delta using slack
    if (!T[y]) delta = min(delta, slack[y]);
  for (x = 0; x != maxN; ++x)            //update X labels
    if (S[x]) lx[x] -= delta;
  for (y = 0; y != maxN; ++y)            //update Y labels
    if (T[y]) ly[y] += delta;
  for (y = 0; y != maxN; y++)            //update slack array
    if (!T[y]) slack[y] -= delta;
}


void CHungarian::augment()                         //main function of the algorithm
{
  //printf("in augment %d\n", max_match);  
  if (max_match == maxN) return;     //check whether matching is already perfect
  long x, y, root;                   //just counters and root vertex
  long *q, wr = 0, rd = 0;      //q - queue for bfs, wr,rd - write and read pos in queue
  q = NULL;
  q = (long*) calloc(maxN, sizeof(long));
  if (q == NULL)
    { printf("not enough memory\n");
      return;
    }

  for(long i = 0; i!= maxN; ++i) 
    {
      S[i] = T[i] = false;
      prev[i] = -1;
    }


  for (x = 0; x != maxN; ++x)            //finding root of the tree
    {
      if (xy[x] == -1)
	{
	  q[wr++] = root = x;
	  prev[x] = -2;
	  S[x] = true;
	  break;
	  
	}
    }


  for (y = 0; y != maxN; ++y)            //initializing slack array
    {
      slack[y] = lx[root] + ly[y] - Cost(root, y);
      slackx[y] = root;
    }

//second part of augment() function
  while (true)                                                        //main cycle
    {

      while (rd < wr)                                                 //building tree with bfs cycle
        {

	  x = q[rd++];                                                //current vertex from X part
	  for (y = 0; y != maxN; ++y)                                     //iterate through all edges in equality graph
	    if (Cost(x, y) == lx[x] + ly[y] &&  !T[y])
	      {
		if (yx[y] == -1) break;                             //an exposed vertex in Y found, so
		//augmenting path exists!
		T[y] = true;                                        //else just add y to T,
		q[wr++] = yx[y];                                    //add vertex yx[y], which is matched
		//with y, to the queue
		add_to_tree(yx[y], x);                              //add edges (x,y) and (y,yx[y]) to the tree
	      }
	  if (y < maxN) break;                                           //augmenting path found!
        }
        if (y < maxN) break;                                               //augmenting path found!

        update_labels();                                                //augmenting path not found, so improve labeling
	
        wr = rd = 0;                
        for (y = 0; y != maxN; ++y)        
	  //in this cycle we add edges that were added to the equality graph as a
	  //result of improving the labeling, we add edge (slackx[y], y) to the tree if
	  //and only if !T[y] &&  slack[y] == 0, also with this edge we add another one
	  //(y, yx[y]) or augment the matching, if y was exposed
	  if (!T[y] &&  slack[y] == 0)
            {
	      if (yx[y] == -1)                                        //exposed vertex in Y found - augmenting path exists!
                {
		  x = slackx[y];
		  break;
                }
	      else
                {
		  T[y] = true;                                        //else just add y to T,
		  if (!S[yx[y]])    
                    {
		      q[wr++] = yx[y];                                //add vertex yx[y], which is matched with
		      //y, to the queue
		      add_to_tree(yx[y], slackx[y]);                  //and add edges (x,y) and (y,
		      //yx[y]) to the tree
                    }
                }
            }
        if (y < maxN) break;                                               //augmenting path found!
    }
  
    if (y < maxN)                                                          //we found augmenting path!
    {
      max_match++;                                                    //increment matching
      //in this cycle we inverse edges along augmenting path
      for (long cx = x, cy = y, ty; cx != -2; cx = prev[cx], cy = ty)
        {
	  ty = xy[cx];
	  yx[cy] = cx;
	  xy[cx] = cy;
        }
      augment();                                                      //recall function, go to step 1 of the algorithm
    }
    if (q != NULL) free(q);
    q = NULL;
}//end of augment() function


#endif




