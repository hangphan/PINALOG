#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<vector>
#include<algorithm>
#include <unistd.h>
using namespace std;

#include"TrieTree.cpp"
#include"TrieTree.h"
#include"nodes.h"
#include"nodes.cpp"
#include"graph_vector.h"
#include"graph_vector.cpp"
#include"go_sim.h"
#include"go_sim.cpp"
#include"similarity.h"
#include"Hungarian_vector.h"
#include"utility.h"


CNodes allnodes, *nodesA, *nodesB;
COntology_RW Ont;
CGraph A, B;
double **AA,**BB, **AB; //distance matrices of species A with itself, B with itself and A with B
long sizeA, sizeB, ppisizeA, ppisizeB;
double MaxScore;
double theta;

//clusters - for core mapping
VVLong grs1, grs2;
//core
VLong mark1, mark2,  gg1, gg2;;
VLong core_g1, core_g2; //proteins in core - being added after each extension
VLong gr1, gr2;
VLong mark, mask;  //mark cluster id for cross validation of function, use local id of graph 1
                   //mask global id of allnodes to get maskout function in function similarity calculation
//candidates for extension mapping
VLong aln1, aln2;  //list of candidates for mapping
long done1, done2; // count number of nodes finished mapping in both network
VDouble score, scores;

double ratio; //sim = ratio * seqsim + (1-ratio) * funcsim
double percentile = 0.85;
char communities_fn1[128], communities_fn2[128];

long fold = -1;


double initial_val =0;

FILE *fo;
CSimilarity Sim(AB, &A, &B);
void get_core();
void get_core_mapping();
int  map_clusters();
void make_communities_files(char *x1, char* x2);// Does not check if file exist, create anyway
void make_communities_files1(char *networkA, char *networkB); //Check if file exist then no need to create community files, otherwise run CFinder to find it
void core_mapping();
void extension(double x);

void mapping();

void filter_results();
void write_aligned_pairs(char *fn);
void write_conserved_edges(char *fon);

void filter_core()
{
  double mean = 0;
  for(long i=0; i!= score.size(); ++i)
    mean+= score[i];
  mean = mean / score.size();
  VDouble temp = score;
  sort(temp.begin(), temp.end());
  long xx =(long) (percentile*temp.size());
  double tempd = temp[xx];

  double thres1 = tempd;

  printf("%f %f %d %d\n", temp[xx], mean, xx, temp.size());

  if (tempd<0.0001) thres1 = mean;

  printf("score threshold  %f\n", thres1 );
  
  VLong temp1, temp2;
  temp = score;
  temp1 = core_g1;
  temp2 = core_g2;
  core_g1.clear();
  core_g2.clear();
  score.clear();
  for(long i =0; i!= temp1.size(); ++i)
    {
      if (temp[i]>= thres1)
	{
	  core_g1.push_back(temp1[i]);
	  core_g2.push_back(temp2[i]);
	  score.push_back(temp[i]);
	}
    }
  printf("Done filtering %d\n", core_g1.size());

}

void init_firstmark()
{

  for(long i =0; i!= sizeA; ++i) mark1.push_back(-1);
  for(long i =0; i!= sizeB; ++i) mark2.push_back(-1);
  done1 = done2 = 0;

  //if mark1[i] = -1: node not done anything,  =1: is core; =0 is candidate for mapping, remember i is local index in graph 
  for(long i = 0; i!= core_g1.size(); ++i)
    {
      if (mark1[core_g1[i]] == -1) {mark1[core_g1[i]] = 1; done1++;} //mark proteins that are in core
      if (mark2[core_g2[i]] == -1) {mark2[core_g2[i]] = 1; done2++;} //mark proteins that are in core
    }
  printf("done first mark\n");
}


void get_candidate_for_mapping()
{
  VLong nei1, nei2;
  gr1.clear();
  gr2.clear();
  //  FILE *f = fopen("temp.txt", "w");
  for(long i = 0; i!= core_g1.size(); ++i)
    {
      nei1 = A.GetNeigh(core_g1[i]);
      nei2 = B.GetNeigh(core_g2[i]);

      for(long j =0; j!= nei1.size(); ++j)
       	if(mark1[nei1[j]] == -1) 
	  {
	    mark1[nei1[j]] = 0;
	    gr1.push_back(nei1[j]);
	    //    fprintf(f, "%s %s\n", nodesA->GetLabel(core_g1[i]), nodesA->GetLabel(nei1[j]));
	  }

      for(long j =0; j!= nei2.size(); ++j)
	if(mark2[nei2[j]] == -1) 
	  {

	    mark2[nei2[j]] = 0;
	    gr2.push_back(nei2[j]);
	  }
      nei1.clear();
      nei2.clear();
    }
  //  fclose(f);
  printf("gr1, gr2 .size () %d %d \n",gr1.size(), gr2.size());
  //  exit(0);
  if (gr1.size()==0 && done1 <sizeA)
    for(long i =0; i!= sizeA; ++i)
      if (mark1[i] == -1) gr1.push_back(i);
  if (gr2.size() ==0 && done2 <sizeB)
    for(long i =0; i!= sizeB; ++i)
      if (mark2[i] == -1) gr2.push_back(i);
}


void read_clusters(char *fn, VVLong &cls, CNodes *nodes, int opt)
{
  FILE *fi = fopen(fn, "r");
  
  if (fi == NULL) 
    {
      printf("cannnot open file %s\n", fn);
      return;
      //      exit(0);
    }
  //read clusters from file : format : cluster_size protein1 protein2 ...
  if (opt ==0)
    while(1)
      {
	if (feof(fi)) break;
	long temp;
	fscanf(fi, "%d", &temp);
	VLong tempv;
	for(long i =0; i!= temp; ++i)
	  {
	    char temps[128];
	    fscanf(fi, "%s", temps);
	    long id = nodes->Search(temps);
	    if (id >=0) tempv.push_back(id);
	  }
	cls.push_back(tempv);
      }
  //read cluster from file: format: ignore lines start with #, other lines start with cluster name, then proteins in that cluster
  int buff_size = 1024 * 1000;  
  char buff[buff_size], st[128];
  
  if (opt ==1)
    while(fgets(buff, buff_size, fi) != NULL)
      {
	if (strlen(buff) <2) continue;
	if (buff[0] == '#') continue;

	VLong tempv;
	char *p = buff;
	
	sscanf(p, "%s", st);
	p += strlen(st) +1;
	while (strlen(p) > 2)
	  {
	    sscanf(p, "%s", st);
	    p += strlen(st) +1;
	    long id = nodes->Search(st);
	    if (id >=0) tempv.push_back(id);
	  }
	//	printf("%d\n", tempv.size());
	if (tempv.size() <=150 && tempv.size()>0)
	  cls.push_back(tempv);
	if (feof(fi))
	  break;
	
      }
  fclose(fi);
}

//get similarity matrix of a subset of nodes in graph A and graph B
//remember to free memory after using this function
double **get_sim_mat(VLong x, VLong y)
{
  double **out = (double**) malloc(x.size()*sizeof(double*));
  for(long i =0; i!= x.size(); ++i)
    {
      out[i] = (double*) malloc(y.size() * sizeof(double));
      for(long j=0; j!= y.size(); ++j)
	{
	  out[i][j] = AB[x[i]][y[j]];
	  //	  printf("%d %d %lf\n", i,j, out[i][j]);
	}
    }
  return out;
}



void mapping()
{
  core_mapping();
  extension(0.01);
  return ;
}

void make_communities_files(char *networkA, char *networkB)
{
  char cmdline[1024];
  sprintf(communities_fn1, "Communities_net1.txt");
  sprintf(communities_fn2, "Communities_net2.txt");
  printf("%s %s %s %s\n", networkA, networkB, communities_fn1, communities_fn2);

  sprintf(cmdline, "./CFinder_commandline -i %s -t 1 -o CFinder_net1  >/dev/null", networkA);  system(cmdline);
  printf("%s\n", cmdline);
  sprintf(cmdline, "cat ./CFinder_net1/k=*/communities >%s", communities_fn1);  system(cmdline);
  sprintf(cmdline, "rm -r ./CFinder_net1");  system(cmdline);
  
  sprintf(cmdline, "./CFinder_commandline -i %s -t 1 -o CFinder_net2 >/dev/null", networkB);  system(cmdline);
  sprintf(cmdline, "cat ./CFinder_net2/k=*/communities >%s", communities_fn2); 
  printf("%s\n", cmdline); system(cmdline);
  sprintf(cmdline, "rm -r ./CFinder_net2");  system(cmdline);
}

void make_communities_files1(char *networkA, char *networkB)
{
  char cmdline[1024];
  sprintf(communities_fn1, "Communities_net1.txt");
  sprintf(communities_fn2, "Communities_net2.txt");
  printf("%s %s %s %s\n", networkA, networkB, communities_fn1, communities_fn2);
  if( access( communities_fn1, F_OK ) == -1 )
    {
      sprintf(cmdline, "./CFinder_commandline -i %s -t 1 -o CFinder_net1", networkA);  system(cmdline);
      printf("%s\n", cmdline);
      sprintf(cmdline, "cat ./CFinder_net1/k=*/communities >%s", communities_fn1);  system(cmdline);
      sprintf(cmdline, "rm -r ./CFinder_net1");  system(cmdline);
    }
  if( access( communities_fn2, F_OK ) == -1 )
    {
      sprintf(cmdline, "./CFinder_commandline -i %s -t 1 -o CFinder_net2", networkB);  system(cmdline);
      sprintf(cmdline, "cat ./CFinder_net2/k=*/communities >%s", communities_fn2); 
      printf("%s\n", cmdline); system(cmdline);
      sprintf(cmdline, "rm -r ./CFinder_net2");  system(cmdline);
    }
}

void core_mapping()
{
  read_clusters(communities_fn1, grs1, nodesA, 1);
  read_clusters(communities_fn2, grs2, nodesB, 1);
  
  if (grs1.size() >0 && grs2.size() >0)
    map_clusters();
  else
    {
      core_g1.clear();
      core_g2.clear();
      
      VLong maxCol, maxRow;
      long i, j;
      double temp;
      
      
      for( i=0; i!= sizeA; ++i)
	{
	  temp =0;
	  maxRow.push_back(-1);
	  
	  for(j = 0; j!= sizeB; ++j)
	    {
	      if (AB[i][j] > temp) {temp = AB[i][j]; maxRow[i] = j; }
	    }
	  
	}
      for( j=0; j!= sizeB; ++j)
	  {
	    temp =0;
	    maxCol.push_back(-1);
	    for(i = 0; i!= sizeA; ++i)
	      if (AB[i][j] > temp) {temp = AB[i][j]; maxCol[j] = i;}
	  }
      long count = 0;
      for(long i=0; i!= sizeA; ++i)
	{
	    if (maxRow[i] >=0)
	      if (maxCol[maxRow[i]] ==i) 
		{
		  core_g1.push_back(i);
		  core_g2.push_back(maxRow[i]);
		}
	    
	}
    }
}


void ReadBlast(char *fn)
{
  FILE *f = fopen(fn, "r");
  char st1[128], st2[128];
  long id1, id2, Aid1, Aid2, Bid1, Bid2, td1, td2;
  double scorex=0;
  MaxScore= 0;
  
  while(!feof(f))
    {
      fscanf(f, "%s %s %lf", st1, st2, &scorex);
      Aid1 = nodesA->Search(st1);
      Aid2 = nodesA->Search(st2);
      if (Aid1<Aid2) {td1 = Aid1 ; Aid1=Aid2; Aid2 = td1;}

      Bid1 = nodesB->Search(st1);
      Bid2 = nodesB->Search(st2);

      if (Bid1<Bid2) {td1 = Bid1 ; Bid1=Bid2; Bid2 = td1;}
      //
      if (Aid1 >=0 && Aid2 >=0)
	{
	  //	  printf("%s %s %d %d %d %d\n",st1, st2, Aid1, Aid2, Bid1, Bid2);
	  if (scorex >MaxScore) MaxScore = scorex;	  
	  if (scorex > AA[Aid1][Aid2] ) 	 
	    AA[Aid1][Aid2] = scorex;
	  continue;
	}
      if (Aid1 >=0 && Bid1 >=0 )
	{
	  //	  printf("%s %s %d %d %d %d\n",st1, st2, Aid1, Aid2, Bid1, Bid2);
	  if (scorex >MaxScore) MaxScore = scorex;
	  if (scorex > AB[Aid1][Bid1])
	    AB[Aid1][Bid1] = scorex;
	  continue;
	}
      if (Bid1 >=0 && Bid2 >=0)
	{
	  //	  printf("%s %s %d %d %d %d\n",st1, st2, Aid1, Aid2, Bid1, Bid2);
	  if (scorex >MaxScore) MaxScore = scorex;
	  if (scorex > BB[Bid1][Bid2])
	    BB[Bid1][Bid2] = scorex;
	  continue;
	}
    }

  fclose(f);

  return ;
}

void ReadBlastEval(char *fn)
{
  FILE *f = fopen(fn, "r");
  char st1[128], st2[128];
  long id1, id2, Aid1, Aid2, Bid1, Bid2, td1, td2;
  double scorex=0;
  
  while(!feof(f))
    {
      fscanf(f, "%s %s %lf", st1, st2, &scorex);
      Aid1 = nodesA->Search(st1);
      Aid2 = nodesA->Search(st2);
      if (Aid1<Aid2) {td1 = Aid1 ; Aid1=Aid2; Aid2 = td1;}

      Bid1 = nodesB->Search(st1);
      Bid2 = nodesB->Search(st2);

      if (Bid1<Bid2) {td1 = Bid1 ; Bid1=Bid2; Bid2 = td1;}
      //      printf("%s %s %d %d %d %d\n",st1, st2, Aid1, Aid2, Bid1, Bid2);
      if (Aid1 >=0 && Aid2 >=0)
	{
	  if (scorex < AA[Aid1][Aid2] ) 	 
	    AA[Aid1][Aid2] = scorex;
	  continue;
	}
      if (Aid1 >=0 && Bid1 >=0)
	{
	  if (scorex < AB[Aid1][Bid1])
	    AB[Aid1][Bid1] = scorex;
	  continue;
	}
      if (Bid1 >=0 && Bid2 >=0)
	{
	  if (scorex < BB[Bid1][Bid2])
	    BB[Bid1][Bid2] = scorex;
	  continue;
	}
    }

  fclose(f);

  return ;
}

void AddSim(char *fn)
{
  FILE *f = fopen(fn, "r");
  char st1[128], st2[128];
  long id1, id2, Aid1, Aid2, Bid1, Bid2, td1, td2;
  double score;
  while(!feof(f))
    {
      fscanf(f, "%s %s %lf", st1, st2, &score);
      
      Aid1 = nodesA->Search(st1);
      Aid2 = nodesA->Search(st2);
      if (Aid1<Aid2) {td1 = Aid1 ; Aid1=Aid2; Aid2 = td1;}

      Bid1 = nodesB->Search(st1);
      Bid2 = nodesB->Search(st2);
      if (Bid1<Bid2) {td1 = Bid1 ; Bid1=Bid2; Bid2 = td1;}

      if (Aid1 >=0 && Aid2 >=0)
	{
	  AA[Aid1][Aid2] = score;
	  continue;
	}
      if (Aid1 >=0 && Bid1 >=0)
	{
	  AB[Aid1][Bid1] = score;
	  continue;
	}
      if (Bid1 >=0 && Bid2 >=0)
	{
	  BB[Bid1][Bid2] = score;
	  continue;
	}
    }

  fclose(f);
  return ;
}

void NormalizeDistMatrix()
{
  long i,j;
  printf("Start normalizing matrix %lf\n", MaxScore);
  
  for( i =1; i!= sizeA; ++i)
    {
      for (j = 0; j!= i; ++j)
	{
	  if (AA[i][i] ==0 || AA[j][j] ==0) AA[i][j] /= MaxScore;
	  else AA[i][j] /= sqrt(AA[i][i] * AA[j][j]);

	}
    }
  //  printf("done AA\n");
  for( i =1; i!= sizeB; ++i)
    {
      for (j = 0; j!= i; ++j)
	{
	  if (BB[i][i] ==0 || BB[j][j] ==0) BB[i][j] /= MaxScore;
	  else BB[i][j] /= sqrt(BB[i][i] * BB[j][j]);
	}
    }
  //  printf("done BB\n");
  for( i =0; i!= sizeA; ++i)
    {
      for(j =0; j!= sizeB; ++j)
	{
	  double temp = AB[i][j];
	  if (AA[i][i]==0 || BB[j][j] ==0) temp /=MaxScore;
	  else temp /= sqrt(AA[i][i] * BB[j][j]);
	  if (temp >1) AB[i][j] = AB[i][j] / max(AA[i][i], BB[j][j]);
	  else AB[i][j] = temp;
	  if (AB[i][j] >1) printf("%f\n", AB[i][j]);
	}
    }
  //  printf("xxx %f \nDone AB\n", AB[0][883]);
  printf("Done normalizing matrix\n");
}

double GetSeqFuncRatio()
{
  VLong maxCol, maxRow;
  long i, j;
  double temp;
 

  for( i=0; i!= sizeA; ++i)
    {
      temp =0;
      maxRow.push_back(-1);
      
      for(j = 0; j!= sizeB; ++j)
	{
	  if (AB[i][j] > temp) {temp = AB[i][j]; maxRow[i] = j; }
	}
      
    }
  for( j=0; j!= sizeB; ++j)
    {
      temp =0;
      maxCol.push_back(-1);
      for(i = 0; i!= sizeA; ++i)
	if (AB[i][j] > temp) {temp = AB[i][j]; maxCol[j] = i;}
    }
  long count = 0;
  for(long i=0; i!= sizeA; ++i)
    {
      if (maxRow[i] >=0)
	{
	  if (maxCol[maxRow[i]] ==i) count++;
	}
    }
  double tempf = (double) count *1.0  / (sizeA + sizeB);
  return 1-tempf;
}

void AddFuncSim()
{
  long i,j;
  long gid1, gid2;
  VLong gos1, gos2;
  printf("adding function similarity\n");

  double functhres = 0.01;
  
  for( i =1; i!= sizeA; ++i)
    {
      for (j = 0; j!= i; ++j)
	{
	  gos1 = allnodes.get_goids_from_node(gg1[i]);
	  gos2 = allnodes.get_goids_from_node(gg1[j]);
	  double x = Ont.function_similarity(gos1, gos2, SCHLICKER);
	  AA[i][j] = AA[i][j] * theta + (1-theta) * x;
	  //if (AA[i][j] >functhres) fprintf(fo, "%s\t%s\t%f\n", allnodes.GetLabel(gid1), allnodes.GetLabel(gid2), AA[i][j]);
	}
    }
  for( i =1; i!= sizeB; ++i)
    {
      for (j = 0; j!= i; ++j)
	{
	  gos1 = allnodes.get_goids_from_node(gg2[i]);
	  gos2 = allnodes.get_goids_from_node(gg2[j]);
	  double x = Ont.function_similarity(gos1, gos2, SCHLICKER);
	  BB[i][j] = BB[i][j] * theta + (1-theta) * x;
	  //	  if (BB[i][j] >functhres) fprintf(fo, "%s\t%s\t%f\n", allnodes.GetLabel(gid1), allnodes.GetLabel(gid2), BB[i][j]);
	}
    }

  for( i =1; i!= sizeA; ++i)
    {
      for(j =0; j!= sizeB; ++j)
	{
	  gos1 = allnodes.get_goids_from_node(gg1[i]);
	  gos2 = allnodes.get_goids_from_node(gg2[j]);
	  double x = Ont.function_similarity(gos1, gos2, SCHLICKER);
	  AB[i][j] = AB[i][j] * theta + (1-theta) *x;
	  //	  if (AB[i][j] >functhres) fprintf(fo, "%s\t%s\t%f\n", allnodes.GetLabel(gid1), allnodes.GetLabel(gid2), AB[i][j]);
	}
    }
}


double **Alloc1(long size)
{
  double **out;
  out = (double**) malloc(size * sizeof(double*)); 
  if (out == NULL) return NULL;
  long i, j;
  for(i =0; i!= size; ++i)
    {
      out[i] = (double *)malloc((i+1) * sizeof(double));
      for(long j=0; j!= i+1; ++j) out[i][j] = initial_val; //default =0
      if (out[i] ==NULL) break;
    }

  if (i < size)
    {
      j = i;
      for(i =1; i!=j; ++i) free(out[i]);
      return NULL; //not enough memory
    }
  return out;
}

double **Alloc2(long size1, long size2)
{
  double **out;
  out = (double**) malloc(size1 * sizeof(double*)); 
  if (out == NULL) return NULL;
  
  long mem = 0;
  long i, j;
  for(i =0; i!= size1; ++i)
    {
      out[i] = (double *)malloc(size2 * sizeof(double));
      for(long j=0; j!= size2; ++j) out[i][j] = initial_val; //default: = 0
      if (out[i] ==NULL){mem = 1;  break;}
    }

  if (mem == 1)
    {
      j = i;
      for(i =0; i!=j; ++i) free(out[i]);
      printf("not enough memory");
      return NULL; //not enough memory
    }
  return out;
}

void Free(long size, double **&x)
{
  for(long i=0; i!= size; ++i) 
    if (x[i] != NULL) {free(x[i]); x[i] = NULL;}
	
  free(x);
  x= NULL;
}

void extension(double thres)
{
  printf("Start extension\n");
  init_firstmark();

  long countx, toquit;
  toquit = 0;
  countx = -1;
  while(1)
    {
      printf("Proteins already aligned:  graph A %d , graph B %d\n", done1, done2);
      if (countx ==0) //no more node pairs can be added anymore, then map all remaining pairs of nodes
	{
	  gr1.clear(); gr2.clear();
	  for(long i =0; i!= sizeA; ++i)
	    if (mark1[i] == -1) gr1.push_back(i); //global id of proteins unprocessed in A
	  for(long i =0; i!= sizeB; ++i)
	    if (mark2[i] == -1) gr2.push_back(i); //global id of proteins unprocessed in B
	  toquit = 1;
	}
      else
	{
	  printf("in here - get candidate\n");

	  get_candidate_for_mapping();

	  if (gr1.size()==0 || gr2.size() ==0) break; //check quit condition
	}


      printf("Start getting neighbours of candidate  proteins (not in core) %d %d %d \n", core_g1.size(), gr1.size(), gr2.size());
      CSimilarity SimF(AB, &A, &B);
      SimF.init(gr1, gr2,  core_g1, core_g2);
      SimF.get_neighbours();
      
      double **costx = SimF.get_sim_mat();
      printf("Finish getting neighbour\n");
      
      CHungarian MappingP(costx);
      MappingP.hungarian(gr1, gr2); //global index
      printf("Finish Hungarian mapping\n");
      //add pairs with high scores to core
      long *mapped = MappingP.get_map();
      long maxs = max(gr1.size(), gr2.size());
      countx = 0;

      for(long i = 0; i!= maxs; ++i)
	{
	  if (i >=gr1.size()) {mark2[gr2[mapped[i]]] = -1; continue;}
	  if (mapped[i] >= gr2.size()) {mark1[gr1[i]] = -1; continue;}
	  double cost = costx[i][mapped[i]];//MappingP.Cost(i, mapped[i]);

	  if (cost > thres)
	    {
	      //add this pair to core list
	      countx ++; 	      
	      done1 ++;	      
	      done2 ++;
	      core_g1.push_back(gr1[i]);
	      core_g2.push_back(gr2[mapped[i]]);
	      score.push_back(cost);
	      mark1[gr1[i]]  = 1;
	      mark2[gr2[mapped[i]]]  = 1;
	      continue;
	    }

	  mark1[gr1[i]] = -1;
	  mark2[gr2[mapped[i]]] = -1;
	}

      printf("Size of mapped pairs %d (%d new)\n", core_g1.size(), countx);
      for(long i =0; i!= gr1.size(); ++i) free(costx[i]);
      free(costx);
      if(toquit ==1) break;
    }
  //print results
  printf("Finish extension mapping\n");
  return;
}



void get_core()
{

  long core_msize= min(sizeA, sizeB) *5/100;
  long mark;
  double temp_score,x, min_score_core;

  for(long i = 0; i!= sizeA; ++i)
    {
      x = 0;
      mark = 0;
      for(long j=0; j!= sizeB; ++j)
      {
	temp_score = AB[i][j];
	if (temp_score >x) 
	  {
	    mark = j;
	  }
      }

      if (AB[i][mark] ==0) continue;
	
      if (core_g1.size() < core_msize) 
	{
	  core_g1.push_back(i); 
	  core_g2.push_back(mark); 
	  score.push_back(AB[i][mark]);
	  continue;
	}
      

	//	continue;
      long min_core_mark = -1;
      double min_score_core = 1;
      
      for(long k = 0; k != score.size(); ++k)
	{
	  if (score[k] < min_score_core)
	    {
	      min_score_core = score[k];
	      min_core_mark = k;
	    }
	}
      
      if (AB[i][mark] >min_score_core && min_core_mark != -1)
	{
	  core_g1[min_core_mark] = i;
	  core_g2[min_core_mark] = mark;
	  score[min_core_mark] = AB[i][mark];
	}
    }
}



int map_clusters()
{
  printf("starting core mapping %d %d\n",  grs1.size(), grs2.size());  

  //Calculate similarity score of mapping clusters with each other  
  double **cl_sim;
  cl_sim = (double**) calloc(grs1.size(), sizeof(double*));
  for(long i=0; i!= grs1.size(); ++i)
    cl_sim[i] = (double*) calloc(grs2.size(), sizeof(double));

  for(long i = 0; i!= grs1.size(); ++i) 
    for(long j =0; j!= grs2.size(); ++j)
      {

	double **costx = get_sim_mat(grs1[i], grs2[j]);
	CHungarian Mapping(costx);

	Mapping.hungarian(grs1[i], grs2[j]);
	cl_sim[i][j] = Mapping.get_score_mapping();

	for(long ii =0; ii!= grs1[i].size(); ++ii) free(costx[ii]);
	free(costx);
      }
  
  CHungarian MappingC(cl_sim);
  VLong x, y;
  for(long i =0; i!= grs1.size(); ++i) x.push_back(i);
  for(long i =0; i!= grs2.size(); ++i) y.push_back(i);
  MappingC.hungarian(x, y);
  printf("finished mapping clusters\n start map details\n");


  long *mapc = MappingC.get_map();


  for(long i =0; i!= grs1.size(); ++i)
    {
        long cl1 = i; 
      long cl2 = mapc[i];
      if (cl2 >= grs2.size()) continue;
      printf("%d %d  %d %d %f\n", cl1, cl2, grs1[cl1].size(), grs2[cl2].size(), cl_sim[cl1][cl2]);
      double **costx = get_sim_mat(grs1[cl1], grs2[cl2]);
      CHungarian MappingP(costx);

      MappingP.hungarian(grs1[cl1], grs2[cl2]);
      long *mapp = MappingP.get_map();
      
      for(long ii=0; ii!= grs1[cl1].size(); ii++)
	{
	  long p1, p2; 
	  p1 = grs1[cl1][ii];
	  if(mapp[ii] >= grs2[cl2].size()) continue;
	 
	  p2 = grs2[cl2][mapp[ii]];
	  int isin = -1;
	  for(long k=0; k!= core_g1.size(); ++k)
	    if(core_g1[k] == p1 && core_g2[k] == p2) {isin =k; break;}
	  if(isin ==-1)
	    {
	      core_g1.push_back(p1);
	      core_g2.push_back(p2);
	      score.push_back(AB[p1][p2]);
	    }
	  
	}
      for(long ii =0; ii!= grs1[cl1].size(); ++ii) free(costx[ii]);
      free(costx);
    }

  for(long i =0; i!= grs1.size(); ++ i) free(cl_sim[i]);
  free(cl_sim);
  
  
  filter_core();
  return 0;
}


void filter_results()
{
  for(long i=0; i!= core_g1.size(); ++i)
    {
      long id1 = core_g1[i];
      long id2 = core_g2[i];
      double scorex = score[i];
      long isin = -1;
      long status = 0;
      for(long i=0; i!= aln1.size(); ++i)
	{
	  if (id1 == aln1[i] && id2 == aln2[i]) {status = 3; isin = i; break;}
	  if (id1 == aln1[i] && id2 != aln2[i]) {status = 1; isin = i; break;}
	  if (id1 != aln1[i] && id2 == aln2[i]) {status = 2; isin = i; break;}
	}

      if (status ==3 && scores[isin] < scorex) {scores[isin] = scorex; continue;}
      if (status ==1 && scores[isin] < scorex) {scores[isin] = scorex; aln2[isin] = id2;continue;}
      if (status ==2 && scores[isin] < scorex) {scores[isin] = scorex; aln1[isin] = id1;continue;}

      if (isin == -1)
	{
	  aln1.push_back(id1);
	  aln2.push_back(id2);
	  scores.push_back(scorex);
	  continue;
	}
    }
}


void write_aligned_pairs(char *fn)
{
  FILE *f = fopen(fn, "w");

  for(long i =0; i!= aln1.size(); ++i)
    fprintf(f, "%s\t%s\t%lf\n", nodesA->GetLabel(aln1[i]), nodesB->GetLabel(aln2[i]) , AB[aln1[i]][aln2[i]]);

  fclose(f);
} 



void write_conserved_edges(char *fon)
{
  FILE *f = fopen(fon, "w");

  for(long i =1; i!= aln1.size(); ++i)
    for(long j = 0; j!= i; ++j)
      if (A.IsNeighbour(aln1[i], aln1[j]) >= 0 && B.IsNeighbour(aln2[i], aln2[j]) >= 0)
	fprintf(f , "%s_%s\t%s_%s\n", nodesA->GetLabel(aln1[i]), nodesB->GetLabel(aln2[i]), nodesA->GetLabel(aln1[j]),  nodesB->GetLabel(aln2[j]));
  
  fclose(f);
}
