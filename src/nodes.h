#ifndef __NODES__
#define __NODES__
#define MaxLen 128
#define MaxLenBuffer 1024*10

#include"utility.h"
#include"go_sim.h"
#include "TrieTree.h"

class CNode{
  char label[100]; //TODO
  long  index;   
  VLong goids;
 public:
  CNode()
    {
    };
  ~CNode()
    {
    };

  void SetupNode(char *s, long id)   { strcpy(label, s); index = id;}
  char* GetLabel() {return label;}
  long  GetIndex(){return index;}
  VLong get_goids(){return goids;}
  void add_goid(long goid) ;
  void add_goid1(long goid, COntology_RW *x);
  void add_goid2(long goid, COntology_RW *x);
  void replace_goid(long goid, long i);
};


class CNodes
{
 protected:
  CNode **nodes;
  long n;
  CTrieTree Trie;
  long inc;
  long mode; // association of GO terms read in mode,  =0: add only bottom line terms, =1 add all ancestors
 public:
  CNodes(long Inc) ;
  CNodes();
  ~CNodes();
  
  long   Search(char *s) {return Trie.Search(s);}
  CNode* SearchNode(char* s) 
    {
      int id = Trie.Search(s); 
      if (id <0) return NULL;
      else return nodes[id];
    }
  long   Add(char *s);
  long   GetNumNodes(){return n;}
  CNode *GetNode(long id){return nodes[id];}
  void  Print(FILE *f);
  char *GetLabel(long i) {nodes[i]->GetLabel();}
  void  Read(FILE *f) ; //read of a list of nodes
  void  Read(char *fn) 
  {
    FILE *f = fopen(fn, "r");
    Read(f);
    fclose(f);
  } //read of a list of nodes

  void  ReadFromGraph(FILE *f) ; //read of a list of nodes
  void  ReadFromGraph(char *fn) 
  {
    FILE *f = fopen(fn, "r");
    ReadFromGraph(f);
    fclose(f);
  } //read of a list of nodes

  CTrieTree *get_trie(){return &Trie;};

  void add_gene_associations(FILE *f, COntology_RW *ontology);
  void add_gene_associations1(FILE *f, COntology_RW *ontology);
  void add_gene_associations_isoform(FILE *f, COntology_RW *ontology);

  void add_goid(long id, long goid) {nodes[id]->add_goid(goid);}
  void add_goid(long id, long goid, COntology_RW *ontology)
  {
    nodes[id]->add_goid2(goid, ontology);
  }
  
  VLong get_goids_from_node(long id);
  VLong get_associated_prots(long goid);
  void change_assoc_mode(long m) {mode =m;}
  void update_ontology_protein_count(COntology_RW *ontology);
};
#endif
