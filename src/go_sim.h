#ifndef __GO_SIM__
#define __GO_SIM__

#include<string.h>
#include "utility.h"

//#include "nodes.h"

//#include "nodes.cpp"
enum OntologyType {CC, MF, BP} ; //cellular component, molecular function, biological process

enum RelationType {ISA, PARTOF};

enum SemSimType {RESNIK, LIN, JIANG, RELEVANCE, WANG}; //semantic similarity type

enum FuncSimType {MAX, AVERAGE, SCHLICKER, WANG2}; //functional similarity type

class CGOTerm
{
protected:
  long id,  ontology_type;
  char goid[15];
  char name[1028]; //description of goterm (to make it easier to understand the term);
  VLong parents, ancestors, children, parent_relations, children_relations;
  long initx;  

public:
  CGOTerm();
  ~CGOTerm();

  void init_goterm(long iid, char* goiid, char *namex, long ont_type);
  long inited(){return initx;} // check if term already initiated 
  void add_parent(long iid, long r_type);  
  void add_ancestor(long iid);
  void add_child(long iid, long r_type);

  char *get_goid() {return goid;}
  long get_id() {return id;}
  long get_ont_type(){return ontology_type;}
  char *get_name(){return name;}
  
  VLong get_ancestors(){return ancestors;}
  long get_ancestor(long i){return ancestors[i];}
  long get_num_ancestor(){return ancestors.size();}

  VLong get_parents(){return parents;}
  long get_parent(long i){return parents[i];}
  long get_num_parents(){return parents.size();}
  VLong get_parent_relations(){return parent_relations;}  

  VLong get_children() {return children;}
  long get_child(long i){return children[i];}
  long get_num_children(){return children.size();}
  VLong get_children_relations() {return children_relations;}
  long get_child_relation(long i) {return children_relations[i];}

  void print_term(FILE *f);

};


/**
 * This class is used to compute the Resnik semantic similarity
 *
 *
 */

class CGOTerm_RW: public CGOTerm
{
  long n_prot; // number of protein being annotated to itself both directly and indirectly
  VLong prots;
  long prev_prot; //previous protein being added - used to compare with a candidate to be added


  //for Wang's method
  VDouble ancestor_scores;
  double semantic_value;
  long is_in_ancestors(long iid) ;
  void compute_ancestor_score(long anc, CGOTerm_RW **CG, double ww1, double ww2);

 public:
  CGOTerm_RW():CGOTerm()
    {
      n_prot =0;
      semantic_value = 0;
      prev_prot = -1;
  }
    
    ~CGOTerm_RW()
      {
      };
    
    void add_annotation(long prot_id); 
    void add_annotation1(long prot_id); 
    long get_num_prot(){return n_prot;}
    double compute_ancestor_scores(CGOTerm_RW **CG, double ww1, double ww2);
    void compute_semantic_value();
    double get_ancestor_score(long i){return ancestor_scores[i];}
    double get_semantic_value() {return semantic_value;}
    VLong get_prots(){return prots;}
    void update_protein_count(long count) {n_prot = count;}
};




class COntology_RW
{
  CGOTerm_RW **terms;
  long n_terms, n_type[3], n_protein; //n_protein: number of proteins being assigned, n_type[CC] number of proteins being assigned to CC ontology and so on
  long roots[3];
  VLong mark;

  CTrieTree Trie;  //to store and seach GOID quickly

  VDouble sem_sim; //Semantic similarity array, a n*(n+1)/2 array  with n = n_terms
  double w1, w2; //for Wang's method, w1 is for is_a relationship contribution, w2 is for part_of relationship contribution 

  double semantic_similarity_resnik(long i, long j);
  double semantic_similarity_jiang(long i, long j);
  double semantic_similarity_relevance(long i, long j);
  double semantic_similarity_lin(long i, long j);
  double semantic_similarity_wang(long i, long j);
  double term_information(long i);
  long   validate_goid(char *st);
  long   get_ont_type(char *st);
  void  induce_all_ancestors(long i);
  void  err()  {    printf("error reading file\n");    exit(0);   }


  double function_similarity_max(VLong, VLong);
  double function_similarity_average(VLong, VLong);
  double  GOScore(VLong, VLong);
  double function_similarity_schlicker(VLong, VLong);
  double function_similarity_wang(VLong, VLong);
  long   get_index(long i, long j) { if (i >j) return i * (i+1)/2 + j ;
    else return j * (j+1)/2 +i;}


 public:
  COntology_RW();
  ~COntology_RW();

  double get_sem_sim(long i, long j);
  long get_n_terms(){return n_terms;}
  void read_go_ids(FILE *f);
  void read_go_ids(char *fn);
  void read_relations(FILE *f);
  void read_relations(char *fn);
  void add_association(long i, long prot); 
  void add_protein_associations(FILE *f);
  void add_protein_associations(char *fn);
  void add_protein_associations2(char *fn);
  void add_protein_associations(FILE *f, CTrieTree *);
  void add_protein_associations(char *fn, CTrieTree *);
  void add_protein_associations1(char *fn, CTrieTree *);
  void add_protein_associations2(char *fn, CTrieTree *);
  void add_term_frequency(char *fn);
  void add_term_frequency(FILE *f);
  void print_frequency(char *fn);

  double semantic_similarity(long i, long j, long option);
  double semantic_similarity(char *goid1, char *goid2, long option);  //if 2 terms are from different ontologies, the distance is set to -1
  double function_similarity(char **, long size1, char **, long size2, long opt);
  double function_similarity(VLong , VLong, long opt);
  double *function_similarity_schlicker_vec(VLong, VLong);

  void calculate_semantic_similarities(long option);
  double get_semantic_similarity(char *, char*);

  //for other program to access 
  long Search(char *goterm) {return Trie.Search(goterm);};
  char *get_goterm(long goid) {return terms[goid]->get_goid();}
  char *get_goterm_name(long goid) {return terms[goid]->get_name();}
  long get_ont_type(long goid) {return terms[goid]->get_ont_type();}
  VLong get_ancestors(long goid) {return terms[goid]->get_ancestors();}
  long get_n_prot(long goid){return terms[goid]->get_num_prot();};
  long get_n_proteins(){return n_protein;};
  long get_n_proteins(long ont_type){return n_type[ont_type];};
  VLong get_assoc_prots(long goid){return terms[goid]->get_prots();}
  VLong get_parents(long goid) {return terms[goid]->get_parents();}
  VLong get_children(long goid) {return terms[goid]->get_children();}

  void update_protein_count(long goid, long count){terms[goid]->update_protein_count(count);}
  VLong get_to_second_level();

};





#endif
