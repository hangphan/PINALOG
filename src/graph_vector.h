#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "nodes.h"
//#include "tree.h"
#ifndef __graph__
#define __graph__
#include "utility.h"
#define MaxNumEdge 200000
#define MaxNumNode 50000


class CGraph
{
  CNodes nodes;
  void GetProteinName(char *st, char *ptn);
  long num_edges;
  VDouble node_scores;
  //Dealing with edges
  //  VVLong Edges; //store local index of pairs of nodes
  VVLong adlist;  //adjacency list: NOTICE: id in adlist is local index in the "nodeids" array
                 //row i: list of nodes adjacent to i-th node in the graph (i-th is the local index in nodes) 
  void LineProcess(char *buffer, long option); //option = 0: read from original file (taxid taxid protein1 protein2 - protein names has database name attached)
                                              //option = 1: read from simpler file (protein1 protein2- simply take in protein names and do not care about databases)
 public:
  CGraph()
    {
      num_edges=0;
    }
  ~CGraph()
    {
    }

  //Nodes
  long get_num_nodes(){return nodes.GetNumNodes();}
  long get_num_edges();
  char *get_label(long gid) {return nodes.GetLabel(gid);}
  VVLong GetEdges() {return adlist;}

  //Edges and graphs
  //long get_num_edges(){return Edges[0].size();}
  void Read(char *fn, long option); //option = 0: read from original file (taxid taxid protein1 protein2 - protein names has database name attached)
                                   //option = 1: read from simpler file (protein1 protein2- simply take in protein names and do not care about databases)
  void Read(FILE *f, long option); 
  CNodes *get_cnodes(){return &nodes;}
  VLong GetNeigh(long lid) {return adlist[lid];}
  VLong Get2Neigh(long gid);
  long IsNeighbour(long id1, long id2);
  void assign_node_scores(VDouble x) {node_scores = x;}

  void Write2GraphML(char *fn);
  void Write2GraphML(FILE *f);

  void Write2HTML(char *fn);
  void Write2HTML(FILE *f);
};

#endif
