#include "go_sim.h"

CGOTerm::CGOTerm()
{
}

CGOTerm::~CGOTerm()
{
  parents.clear();
  parent_relations.clear();
  children.clear();
  children_relations.clear();
  ancestors.clear();
}


void CGOTerm::init_goterm(long iid, char* goiid, char *namex, long ont_type)
{
  id = iid;
  strcpy(goid, goiid);
  strcpy(name, namex);
  ontology_type = ont_type;
  add_ancestor(iid);
  initx =1;
}


void CGOTerm::add_parent(long iid, long r_type)
{
  for(long i =0; i!= parents.size(); ++i)
    if (parents[i] == iid) return;
  parents.push_back(iid);
  parent_relations.push_back(r_type);
}

void CGOTerm::add_ancestor(long iid)
{
  for(long i = 0; i!= ancestors.size(); ++i)
    if (ancestors[i] == iid) return;
  ancestors.push_back(iid);
}

void CGOTerm::add_child(long iid, long r_type)
{
  for(long i=0; i!= children.size(); ++i)
    if (children[i] == iid) return;
  children.push_back(iid);
  children_relations.push_back(r_type);
}
  
void CGOTerm_RW::add_annotation(long prot_id)
{

  for(long i=n_prot-1; i!= -1; --i)
    if(prots[i] == prot_id) 
	return;
  prots.push_back(prot_id);
  n_prot++;
}

void CGOTerm_RW::add_annotation1(long prot_id)
{
  if (prev_prot <0) 
    {
      prev_prot = prot_id;
      n_prot ++;
      return;
    }
  if (prot_id == prev_prot) return;

  prev_prot = prot_id;
  n_prot++;
}


void CGOTerm::print_term(FILE *f)
{
  fprintf(f, "%s\t", goid);
}

long CGOTerm_RW::is_in_ancestors(long iid) 
  {
    for(long i = 0; i!= ancestors.size(); ++i)
      if (ancestors[i] == iid) return i;
    return -1;
  }



COntology_RW::COntology_RW()
{
  long MaxNumTerm = 50000;
  terms = (CGOTerm_RW **) calloc(MaxNumTerm, sizeof(CGOTerm_RW*));

  n_terms = 0;
  n_type[0] =n_type[1]=n_type[2] = 0;
  w1 = 0.8 ;//for is-a relationship
  w2 = 0.6 ;//for part-of relationship
}

COntology_RW::~COntology_RW()
{
  for(long i = 0; i!= n_terms; ++i)
    if (terms[i] != NULL) delete terms[i];
  free(terms);
}

long COntology_RW::validate_goid(char *st)
{
  if (st[0] == 'G' && st[1] == 'O' && st[2] == ':') return 1;
  else return -1;
}

long COntology_RW::get_ont_type(char *st)
{
  if (strcmp(st, "molecular_function") ==0) return MF;
  if (strcmp(st, "cellular_component") ==0) return CC;
  if (strcmp(st, "biological_process") ==0) return BP;
  return -1;
}


/**
 *  Read in all the GO ids in the input file to make a Trie Tree for further processing
 *  Might change this to read in different type of ontology type (Biological process, Molecular function or Cellular Component)
 *  Or keep it the same and extract types based on ont_type information later
 */
void COntology_RW::read_go_ids(char *fn)
{

  FILE *f = fopen(fn, "r");

  read_go_ids(f);
  fclose(f);
}
void COntology_RW::read_go_ids(FILE *f)
{
  //  printf("ok here\n");
  long buff_size = 10240;
  char buffer[buff_size];
  char goiid[15];
  long obsolete = 0;
  long ont_type;
  long count = 0;
  char namex[1024];
  
  while(!feof(f))
    {
      fgets(buffer, buff_size, f);
      char st1[128], st2[128], st3[128];
      sscanf(buffer, "%s", st1);
      if (strstr(st1, "is_obsolete") != NULL) 
	{
	  obsolete = 1;
	}
      if (strcmp(st1, "[Term]")==0)
	{
	  //add old term
	  if (count >0 && obsolete ==0)
	    {
	      Trie.Add(goiid); 
	      if (n_terms  == Trie.Size()) break;  // move on if goid already exists
	      terms[n_terms]= new CGOTerm_RW;
	      if (namex[strlen(namex)-1] == '\n') namex[strlen(namex) -1] = '\0';
	      terms[n_terms]->init_goterm(n_terms, goiid, namex, ont_type);  
	      n_terms++;
	    }
	  count++;
	  obsolete = 0;
	  //read new term: get GO id
	  fgets(buffer, buff_size, f);
	  sscanf(buffer, "%s %s", st2, st3);
	  if (strcmp(st2, "id:") != 0) break;
	  if (validate_goid(st3) == -1) break; //move on if not a goid
	  strcpy(goiid, st3);  

	  //get ontology type 
	  fgets(buffer, buff_size, f);
	  sscanf(buffer, "%s %s", st2, st3);
	  if (strcmp(st2, "name:") ==0)
	    {
	      char *p = buffer;
	      p += 6;
	      strcpy(namex, p);
	    }
	  if (strcmp(st2, "namespace:") ==0) ont_type = get_ont_type(st3);
	  else 
	    {
	      fgets(buffer, buff_size, f);
	      sscanf(buffer, "%s %s", st2, st3);
	      if (strcmp(st2, "namespace:") ==0) ont_type = get_ont_type(st3);
	    }
	}
    }
  roots[BP] = Trie.Search("GO:0008150"); //biological process root
  roots[MF] = Trie.Search("GO:0003674"); //moleculuar root
  roots[CC] = Trie.Search("GO:0005575"); //cellular componenet root
  //  printf("n term = %d\n", n_terms);
}

/**
 * Read the file containing list of GO terms and their relationship with each other, read into the DAG structure of the program
 */

void COntology_RW::read_relations(char *fn)
{
  FILE *f = fopen(fn, "r");
  read_relations(f);
  fclose(f);
}

void COntology_RW::read_relations(FILE *f) 
{
  long buff_size = 1024;
  char buffer[buff_size];
  //  printf("Trie size %d\n", Trie.Size());
  while(!feof(f))
    {
      fgets(buffer, buff_size, f);
      char st1[128], st2[128], st3[128], st4[128];
      sscanf(buffer, "%s %s", st1, st2);

      // [Term] is the sign of starting of a new term. Other than that, just read through and skip
      if (strcmp(st1, "[Term]")==0)
	{
	  fgets(buffer, buff_size, f); //Read next line to get GO ID 
	  sscanf(buffer, "%s %s", st2, st3);
	  if (strcmp(st2, "id:") != 0) break;	  
	  if (validate_goid(st3) == -1) break; //move on if not a goid
	  long id = Trie.Search(st3);
	  if (id <0) continue;  //move on if goid not exists in the list (i.e., list imported before reading this file)

	    
	  while (strlen(buffer) > 1)
	    {
	      long parent ;
	      fgets(buffer, buff_size, f);
	      sscanf(buffer, "%s %s %s", st2, st3, st4);

	      if (strcmp(st2, "is_a:") ==0) 
		{
		  parent = Trie.Search(st3);
		  if (parent >=0)
		    {
		      terms[id]->add_parent(parent, ISA);
		      terms[parent]->add_child(id, ISA);
		    }
		}
	      if (strcmp(st2, "relationship:") ==0)
		{
		  parent = Trie.Search(st4);
		  if (parent >=0)
		    {
		      terms[id]->add_parent(parent, PARTOF);
		      terms[parent]->add_child(id, PARTOF);
		    }
		}
	      
	    }
	}
      
    }
  for(long i = 0; i!= n_terms; ++i)
    mark.push_back(0);

  for(long i = 0; i!= n_terms; ++i)
    induce_all_ancestors(i);
}

/**
 *  Induce all the ancestors of one go-term including itselft
 */
void COntology_RW::induce_all_ancestors(long id)
{
  
  if (mark[id] ==1)
    return; // this term has been visited and updated

  terms[id]->add_ancestor(terms[id]->get_id());
  
  long n_parent = terms[id]->get_num_parents();

  if (n_parent ==0) 
    {
      mark[id] =1;
      return;
    }

  for(long i = 0; i!= n_parent; ++i)
    {
      long parent = terms[id]->get_parent(i);
      terms[id]->add_ancestor(parent);
      induce_all_ancestors(parent); // truy ho^`i 
    
      long npp = terms[parent]->get_num_ancestor();
      
      for(long j = 0; j!= npp; ++j)
	  terms[id]->add_ancestor(terms[parent]->get_ancestor(j));

    }
  mark[id] = 1;
}

/**
 *  Given a file of annotation of proteins in terms of GO terms, this procedure read in the association information
 *  and count the number of proteins begin associated to each term in the DAG graph
 */

void COntology_RW::add_protein_associations(char *fn)
{
  FILE *f = fopen(fn, "r");
  add_protein_associations(f);
  fclose(f);
}

void COntology_RW::add_protein_associations2(char *fn)
{
  FILE *f = fopen(fn, "r");
  long buff_size = 3000;
  char buffer[buff_size], previous_prot[128], new_protein[128], st1[128], st2[128], st3[128], st4[128], st5[128];
  char *p1, *p2;
  long mark[3]; 
  long iid, ot; //temporary variable

  n_protein = 0;
  strcpy(previous_prot, "");
  mark[CC] = 0; mark[BP]=0; mark[MF] = 0;
  
  while(1)
    {
      //printf("ok here\n");  
      if (feof(f)) break;
      
      fgets(buffer, buff_size, f);
      //      fputs(buffer, stdout);
      sscanf(buffer, "%s %s %s ", st1, st2, st3);
      iid= Trie.Search(st3);
      //      printf("%s\n", st3);
      if (iid <0) continue;
      //CHECK HERE
      strcpy(new_protein, st2);
      ot = terms[iid]->get_ont_type();

      if (strcmp(new_protein, previous_prot) !=0)
	{
	  strcpy(previous_prot, new_protein);
	  n_protein++;
	  mark[CC] = 0; mark[BP]=0; mark[MF] = 0;
	  if (n_protein % 5000 ==0) printf("%s %d\n", st2, n_protein);
	}

      if (mark[ot] ==0) 
	{
	  n_type[ot] ++; 
	  mark[ot] = 1;
	}
      
      add_association(iid, n_protein);
    }




  fclose(f);
}

/**
 *   Add gene association information of GO terms and  genes/proteins to the DAG graph (counting the number of proteins/ genes being annotated with terms)
 *   If a child term is assigned with a protein, then all its ancestors are automatically assigned with that protein as well 
 *   The gene-association file may contain continuous lines with same protein and associated term (for different evidence codes)
 */

void COntology_RW::add_protein_associations(FILE *f)
{
  long buff_size = 3000;
  char buffer[buff_size], previous_prot[128], new_protein[128], st1[128], st2[128], st3[128], st4[128], st5[128];
  char *p1, *p2;
  long mark[3]; 
  long iid, ot; //temporary variable

  n_protein = 0;
  strcpy(previous_prot, "");
  mark[CC] = 0; mark[BP]=0; mark[MF] = 0;
  fgets(buffer, buff_size, f);
  while(!feof(f))
    {
      fgets(buffer, buff_size, f);
      sscanf(buffer, "%s %s %s %s %s", st1, st2, st3, st4, st5);
      p1 = NULL;
      p2 = NULL;
      p1 = strstr(st4, "GO:");
      p2 = strstr(st5, "GO:");
      if (p1 != NULL)
	  iid= Trie.Search(st4);
      else if (p2 != NULL)
	iid= Trie.Search(st5);

      p1 = strstr(st4, "NOT");

      if (p1 != NULL)  continue;
      if (iid <0) continue;
      //CHECK HERE
      strcpy(new_protein, st2);
      ot = terms[iid]->get_ont_type();

      if (strcmp(new_protein, previous_prot) !=0)
	{
	  strcpy(previous_prot, new_protein);
	  n_protein++;
	  mark[CC] = 0; mark[BP]=0; mark[MF] = 0;
	  if (n_protein % 5000 ==0) printf("%s %d\n", st2, n_protein);
	}

      if (mark[ot] ==0) 
	{
	  n_type[ot] ++; 
	  mark[ot] = 1;
	}
      
      add_association(iid, n_protein);
    }
  /**
  for(long i = 0; i!= 3; ++i)
    {
      printf("Check point: n_protein %d,  n_prot at root %d \n", n_type[i], terms[roots[i]]->get_num_prot(), i );
    }

  */  
}

VLong COntology_RW::get_to_second_level()
  {
    VLong out;
    out.push_back(roots[BP]);
    out.push_back(roots[MF]);
    out.push_back(roots[CC]);
    CGOTerm_RW *p = terms[roots[BP]];
    long np = p->get_num_children();
    for(long i=0; i!= np; ++i) out.push_back(p->get_id());
    
    p = terms[roots[MF]];
    np = p->get_num_children();
    for(long i=0; i!= np; ++i) out.push_back(p->get_id());
    
    p = terms[roots[CC]];
    np = p->get_num_children();
    for(long i=0; i!= np; ++i) out.push_back(p->get_id());
    
    return out;
  }

void COntology_RW::add_protein_associations(char *fn, CTrieTree* proteins)
{
  FILE *f = fopen(fn, "r");
  add_protein_associations(f, proteins);
  fclose(f);
}

/**
 *   Add gene association information of GO terms and  genes/proteins to the DAG graph (counting the number of proteins/ genes being annotated with terms)
 *   If a child term is assigned with a protein, then all its ancestors are automatically assigned with that protein as well 
 *   The gene-association file may contain continuous lines with same protein and associated term (for different evidence codes)
 */
void COntology_RW::add_protein_associations(FILE *f, CTrieTree* proteins)
{
  long buff_size = 3000;
  char buffer[buff_size], previous_prot[128], new_protein[128], st1[128], st2[128], st3[128], st4[128], st5[128];
  long mark[3], iid, prot, ot; 
  char *p1, *p2;

  mark[CC] = 0; mark[BP]=0; mark[MF] = 0;
  n_protein = 0;
  strcpy(previous_prot, "");

  while(!feof(f))
    {
      fgets(buffer, buff_size, f);
      sscanf(buffer, "%s %s %s %s %s", st1, st2, st3, st4, st5);
      p1 = NULL;
      p2 = NULL;
      p1 = strstr(st4, "GO:");
      p2 = strstr(st5, "GO:");
      if (p1 != NULL)	iid= Trie.Search(st4);
      else if (p2 != NULL) iid= Trie.Search(st5);
      p1 =strstr(st4, "NOT"); 

      if (  p1!= NULL)  continue;
      if (iid <0) continue;
      prot = proteins->Search(st2);
      if (prot <0) continue;

      strcpy(new_protein, st2);
      ot = terms[iid]->get_ont_type();
      
      if (strcmp(new_protein, previous_prot) !=0)
	{
	  strcpy(previous_prot, new_protein);
	  mark[CC] = 0; mark[BP]=0; mark[MF] = 0;
	  n_protein ++;
	}

      if (mark[ot] ==0) {n_type[ot] ++; mark[ot] = 1;}

      add_association(iid, prot);
    }
  /**
  for(long i = 0; i!= 3; ++i)
    {
      printf("Check point: n_protein %d,  n_prot at root %d \n", n_type[i], terms[roots[i]]->get_num_prot(), i );
    }

  */  
}


//the list of ancestors of a node contains itself 
void COntology_RW::add_association(long id, long prot)
{
  long np = terms[id]->get_num_ancestor();
  for(long i = 0; i!= np; ++i)
    {
      long ancestor = terms[id]->get_ancestor(i);
      terms[ancestor]->add_annotation(prot);
    }
}

/**
 * This block is dedicated to calculating the semantic similarity 
 * between two terms in the DAG of GO terms
 */

double COntology_RW::term_information(long i)
{
  long ot = terms[i]->get_ont_type();
  double np = (double) terms[i]->get_num_prot();
  return - log (np/n_type[ot]);
}



double COntology_RW::semantic_similarity(char *goid1, char *goid2, long option)
{
  long id1, id2;
  id1 = Trie.Search(goid1);
  id2 = Trie.Search(goid2);

  if (id1 <0 || id2<0) return -1;
  return semantic_similarity(id1, id2, option);
}

double COntology_RW::semantic_similarity(long id1, long id2, long option)
{
  long ot1, ot2;
  ot1 = terms[id1]->get_ont_type();
  ot2 = terms[id2]->get_ont_type();

  if (ot1 != ot2) return -1;

  if (option == RESNIK) return semantic_similarity_resnik(id1, id2);
  if (option == LIN)    return semantic_similarity_lin(id1, id2);
  if (option == JIANG)  return semantic_similarity_jiang(id1, id2);
  if (option == RELEVANCE) return semantic_similarity_relevance(id1, id2);
  if (option == WANG) return semantic_similarity_wang(id1, id2);
}

double COntology_RW::semantic_similarity_resnik(long id1, long id2)
{
  long ot = terms[id1]->get_ont_type();
  long na1, na2;
  na1 = terms[id1]->get_num_ancestor();
  na2 = terms[id2]->get_num_ancestor();
  
  long min = 300000;
  for(long ii=0; ii!= na1; ++ii)
    for(long jj=0; jj!= na2; ++jj)
      {
	long ac1 = terms[id1]->get_ancestor(ii);
	long ac2 = terms[id2]->get_ancestor(jj) ;
	if ( ac1 == ac2 )
	  {
	    long np = terms[ac1]->get_num_prot() ;
	    if (np>= 0 && np <min)
	      min = np;
	    
	  }
    }
  if (min ==0) return 0;
  double minf = (double) min;
  return  -log(minf/n_type[ot]);
}

double COntology_RW::semantic_similarity_lin(long id1, long id2)
{
  double ica = semantic_similarity_resnik(id1, id2);
  long np1 = terms[id1]->get_num_prot();
  long np2 = terms[id2]->get_num_prot();
  if (ica ==0 || np1 ==0 ||np2 ==0 ) return 0;
  return (double) 2*ica / (term_information(id1) + term_information(id2));
}

double COntology_RW::semantic_similarity_jiang(long id1, long id2)
{
  double ica = semantic_similarity_resnik(id1, id2);
  long np1 = terms[id1]->get_num_prot();
  long np2 = terms[id2]->get_num_prot();
  if (ica ==0 || np1 ==0 ||np2 ==0 ) return 0;
  return (double) 1/ (- 2*ica  + term_information(id1) + term_information(id2)  +1);
}


double COntology_RW::semantic_similarity_relevance(long id1, long id2)
{
  long ot = terms[id1]->get_ont_type();
  
  long na1, na2;
  na1 = terms[id1]->get_num_ancestor();
  na2 = terms[id2]->get_num_ancestor();
  
  if (id1 == id2) 
    {
      double np =(double) terms[id1]->get_num_prot();
      return 1 - np/n_type[ot];
    }
  long  ca;
  long min = 300000;
  for(long ii=0; ii!= na1; ++ii)
    for(long jj=0; jj!= na2; ++jj)
      {
	if (terms[id1]->get_ancestor(ii) == terms[id2]->get_ancestor(jj))
	  {
	    long common_ancestor = terms[id1]->get_ancestor(ii);
	    long np = terms[common_ancestor]->get_num_prot();
	    if ( np >=0 && np< min)
	      {
		ca = common_ancestor; 
		min = np;
		
	    }
	  }
      }
  long np1 = terms[id1]->get_num_prot();
  long np2 = terms[id2]->get_num_prot();
  //  printf("common ancestor %s %d %d\n", terms[ca]->get_goid(), terms[ca]->get_num_prot(), n_type[ot]);
  if (min ==0 || np1 ==0 ||np2 ==0 ) return 0;
  double minf = (double) min;
  double ica= -log(minf/n_type[ot]);
  
  //  if (id1 ==id2){ printf("%f %d %f %f \n", minf, n_type[ot], 1- minf/n_type[ot]);}
  return (double)  2*ica * (1- minf/n_type[ot] ) /(term_information(id1) + term_information(id2) );
}


double COntology_RW::semantic_similarity_wang(long id1, long id2)
{
  if (id1 == id2) return 1;
  long ot = terms[id1]->get_ont_type();
  long na1, na2;
  na1 = terms[id1]->get_num_ancestor();
  na2 = terms[id2]->get_num_ancestor();
  double sum = 0;
  for(long ii=0; ii!= na1; ++ii)
    for(long jj=0; jj!= na2; ++jj)
      if (terms[id1]->get_ancestor(ii) == terms[id2]->get_ancestor(jj))
	sum += terms[id1]->get_ancestor_score(ii) + terms[id2]->get_ancestor_score(jj);
  
  return sum/(terms[id1]->get_semantic_value() + terms[id2]->get_semantic_value());
}

//get semantic similarity of 2 terms i and j from either an array of similarity score or from a file
double COntology_RW::get_sem_sim(long i, long j)
{
      return semantic_similarity(i, j, RELEVANCE);
}


double COntology_RW::get_semantic_similarity(char *st1, char* st2)
{
  long id1 = Trie.Search(st1);
  long id2 = Trie.Search(st2);
  if (id1 <0 || id2 <0) return -1;
  else return get_sem_sim(id1, id2);
}

/** 
 * The next procedures are only used for Wang's method of semantic similarity. 
 * Others don't use these functions.
 * This function compute the score of ancestors of a node, on the amount of information the ancestor contributes to the 
 * meaning of the node term. The score of an ancestor A (in relation with term t) 
 *  S_t(A) = 1 if  A = t
 *         = max (w_e(A, cA) S_t(cA)) {cA is a child of A, w_e = 0.8 if A-cA is a "is_a" relationship and = 0.6 if it is a "part-of" relationship}
 */
void CGOTerm_RW::compute_ancestor_score(long anc_index, CGOTerm_RW **CG, double ww1, double ww2)
{
  if (ancestor_scores[anc_index] >0) return ;
  long  iid  = ancestors[anc_index];

  if (iid == id) 
    {    
      ancestor_scores[anc_index] = 1;
      return ;
    }

  double max = 0;
  
  CGOTerm_RW *term = CG[iid];
  long n_child = term->get_num_children();
  for(long i = 0; i!= n_child; ++i)
    {
      long cid = is_in_ancestors(term->get_child(i));
      if (cid >= 0)
	{
	  long cid_relation = term->get_child_relation(i);
	 
	  if (ancestor_scores[cid] < 0) 
	    compute_ancestor_score(cid, CG, ww1, ww2);
	  
	  double x = ancestor_scores[cid];
	  if (cid_relation == ISA) x *= ww1;
	  else x *= ww2;
	  if ( x > max)
	    max = x;
	}
    }
  ancestor_scores[anc_index] = max;
}

double CGOTerm_RW::compute_ancestor_scores(CGOTerm_RW **CG, double ww1, double ww2)
{
  for(long i = 0; i!= ancestors.size(); ++i)
    ancestor_scores[i] = -1;
  for(long i = 0; i!= ancestors.size(); ++i)
    compute_ancestor_score(i, CG, ww1, ww2);
}

void CGOTerm_RW::compute_semantic_value()
{
  semantic_value = 0;
  for(long i = 0; i!= ancestors.size(); ++i)
    semantic_value += ancestor_scores[i];
  //  print_term(stdout);
}

/**
 * This block is dedicated to computing the functional similarity of 2 list of GO terms 
 * associated with 2 proteins (here represented only by 2 lists of GO terms)
 */

double COntology_RW::function_similarity(char **goids1, long size1, char **goids2, long size2, long option)
{
  if (size1 ==0 || size2 ==0) return 0;
  
  VLong ids1, ids2;
  for(long i = 0; i!= size1; ++i)
      ids1.push_back( Trie.Search(goids1[i]));

  for(long j = 0; j!= size2; ++j)   
    ids2.push_back( Trie.Search(goids2[j]));

  return function_similarity(ids1, ids2, option); 
}

double COntology_RW::function_similarity(VLong ids1,  VLong ids2, long option)
{
  if (ids1.size() ==0 || ids2.size() ==0) return 0;
  if (option == MAX)     return   function_similarity_max(ids1,  ids2);
  if (option == AVERAGE)  return   function_similarity_average(ids1, ids2);
  if (option == SCHLICKER)  return  function_similarity_schlicker(ids1, ids2);
  if (option == WANG2)  return  function_similarity_wang(ids1, ids2);

}

double COntology_RW::function_similarity_max(VLong ids1, VLong ids2)
{
  double max = -1;
  for(long i = 0; i!= ids1.size(); ++i)
    for(long j = 0; j!= ids2.size(); ++j)
      {
	double sim = get_sem_sim(ids1[i], ids2[j]);
	if (sim > max) max = sim;
      }
  return max;
}


double COntology_RW::function_similarity_average(VLong ids1, VLong ids2)
{
  long count = 0;
  double sum = 0;
  for(long i = 0; i!= ids1.size(); ++i)
    for(long j = 0; j!= ids2.size(); ++j)
      {
	double sim = get_sem_sim(ids1[i], ids2[j]);
	if (sim >=0) 
	  {
	    sum+= sim;
	    count++;
	  }
      }
  return sum/count;
}

double COntology_RW::GOScore(VLong ids1, VLong ids2)
{
  //max (rowScore, columnScore) ; rowScore = 1/n sum(max in row), columnScore = 1/m sum(max in column);
  if (ids1.size() ==0 || ids2.size() ==0) return 0;

  double sumMaxRow = 0;
  long count = 0;
  for(long i = 0; i!= ids1.size(); ++i)
    {
      double maxRow = -1;
      
      for(long j = 0; j!= ids2.size(); ++j)
	{
	  double sim = get_sem_sim(ids1[i], ids2[j]);
	  if (sim >1)
	    exit(0);
	  if (sim > maxRow) maxRow = sim;
	}
      if (maxRow >=0) 
	{
	  sumMaxRow += maxRow;
	  count ++;
	}
    }

  if (count >0) 
    {
      sumMaxRow = sumMaxRow/count; 
    }
  else sumMaxRow = -1;

  double sumMaxCol = 0;
  count = 0;
  for(long i = 0; i!= ids2.size(); ++i)
    {
      double maxCol = -1;
      for(long j = 0; j!= ids1.size(); ++j)
	{
	  double sim = get_sem_sim(ids2[i], ids1[j]);
	  if (sim > maxCol) maxCol = sim;
	}
      if (maxCol >=0) 
	{
	  sumMaxCol += maxCol;
	  count ++;
	}
    }

  if (count >0) 
    {
      sumMaxCol = sumMaxCol/count; 
      if (sumMaxCol >1) printf("maxCol %f\n", sumMaxCol);
    }
  else sumMaxCol = -1;

  if (sumMaxCol> sumMaxRow) return sumMaxCol;
  else return sumMaxRow;
}

double *COntology_RW::function_similarity_schlicker_vec(VLong ids1, VLong ids2)
{

  double *out = (double*)calloc(3, sizeof(double));
  VVLong gos1, gos2;
  long  size_gos1[3], size_gos2[3];

  for(long i=0; i!= 3; ++i)
    {
      VLong temp;
      gos1.push_back(temp);
      gos2.push_back(temp);
    }

  for(long i =0; i!= ids1.size(); ++i)
    {
      long ot = terms[ids1[i]] ->get_ont_type() ;
      gos1[ot].push_back(ids1[i]);
    }
  for(long i =0; i!= ids2.size(); ++i)
    {
      long ot = terms[ids2[i]] ->get_ont_type() ;
      gos2[ot].push_back(ids2[i]);
    }
  out[MF] = GOScore(gos1[MF], gos2[MF]);
  out[BP] = GOScore(gos1[BP], gos2[BP]);
  out[CC] = GOScore(gos1[CC], gos2[CC]);

  return out;
}

double COntology_RW::function_similarity_schlicker(VLong ids1, VLong ids2)
{
  double *x = function_similarity_schlicker_vec(ids1, ids2);
  double out = sqrt(x[MF] *x[MF] + x[BP]*x[BP])/sqrt(2);
  free(x);
  return out;
}



//list of go_terms are from the same ontology, usually molecular function, does not accept those ffrom different ontologies
double COntology_RW::function_similarity_wang(VLong ids1, VLong ids2)
{
  double sumMaxRow = 0;
  for(long i = 0; i!= ids1.size(); ++i)
    {
      double maxRow = 0;
      for(long j = 0; j!= ids2.size(); ++j)
	{
	  double sim = get_sem_sim(ids1[i], ids2[j]);
	  if (sim > maxRow) maxRow = sim;
	}
      sumMaxRow += maxRow;
    }
  
  double sumMaxCol = 0;
  for(long i = 0; i!= ids2.size(); ++i)
    {
      double maxCol = 0;
      for(long j = 0; j!= ids1.size(); ++j)
	{
	  double sim = get_sem_sim(ids2[i], ids1[j]);
	  if (sim > maxCol) maxCol = sim;
	}
      sumMaxCol += maxCol;
    }
  return (sumMaxRow + sumMaxCol)/(ids1.size() + ids2.size());
}

void COntology_RW::add_term_frequency(char *fn)
{
  FILE *f= fopen(fn, "r");
  add_term_frequency(f);
  fclose(f);

}

void COntology_RW::add_term_frequency(FILE *f)
{
  char tterm[128];
  long ffreq;
  while(!feof(f))
    {
      fscanf(f, "%s %ld", tterm, &ffreq);
      long tid = Trie.Search(tterm);
      if (tid <0) continue;
      terms[tid]->update_protein_count(ffreq);
    }
}

void  COntology_RW::print_frequency(char *fn)
{
  FILE *f = fopen(fn, "w");
  for(long i =0; i!= n_terms; ++i)
    {
      fprintf(f, "%s\t%d\n", terms[i]->get_goid(), terms[i]->get_num_prot());
    }
  fclose(f);
}

void COntology_RW::add_protein_associations1(char *fn, CTrieTree* proteins)
{
  FILE *f = fopen(fn, "r");
  long buff_size = 3000;
  char buffer[buff_size], previous_prot[128], new_protein[128], st1[128], st2[128], st3[128], st4[128], st5[128];
  n_protein = 0;
  strcpy(previous_prot, "");
 
  long mark[3]; 
  mark[CC] = 0; mark[BP]=0; mark[MF] = 0;
  while(!feof(f))
    {
      fgets(buffer, buff_size, f);
      sscanf(buffer, "%s %s %s %s %s", st1, st2, st3, st4, st5);
      //      fputs(buffer, stdout);
      //      if (strstr(st1, "UniProtKB") ==NULL) continue;
      long iid;
      if (strstr(st4, "GO:") != NULL)
	  iid= Trie.Search(st4);
      else if (strstr(st5, "GO:") != NULL)
	iid= Trie.Search(st5);
      if (strstr(st4, "NOT")  != NULL)  continue;
      if (iid <0) continue;

      long prot = proteins->Add(st2);
      //      if (prot <0) continue;

      strcpy(new_protein, st2);
      long ot = terms[iid]->get_ont_type();
      
      if (strcmp(new_protein, previous_prot) !=0)
	{
	  strcpy(previous_prot, new_protein);
	  mark[CC] = 0; mark[BP]=0; mark[MF] = 0;
	  n_protein ++;
	}

      if (mark[ot] ==0) {n_type[ot] ++; mark[ot] = 1;}

      add_association(iid, prot);
    }
  fclose(f);
}

void COntology_RW::add_protein_associations2(char *fn, CTrieTree* proteins)
{
  FILE *f = fopen(fn, "r");
  long buff_size = 3000;
  char buffer[buff_size], previous_prot[128], new_protein[128], st1[128], st2[128], st3[128], st4[128], st5[128];
  n_protein = 0;
  strcpy(previous_prot, "");
 
  long mark[3]; 
  mark[CC] = 0; mark[BP]=0; mark[MF] = 0;
  while(!feof(f))
    {
      fgets(buffer, buff_size, f);
      sscanf(buffer, "%s %s %s %s ", st1, st2, st3, st4);
      //      fputs(buffer, stdout);
      //      if (strstr(st1, "UniProtKB") ==NULL) continue;
      long iid;
      if (strstr(st3, "GO:") != NULL)
	  iid= Trie.Search(st3);
      else if (strstr(st4, "GO:") != NULL)
	iid= Trie.Search(st4);
      if (strstr(st4, "NOT")  != NULL)  continue;
      if (iid <0) continue;

      long prot = proteins->Add(st2);
      //      if (prot <0) continue;

      strcpy(new_protein, st2);
      long ot = terms[iid]->get_ont_type();
      
      if (strcmp(new_protein, previous_prot) !=0)
	{
	  strcpy(previous_prot, new_protein);
	  mark[CC] = 0; mark[BP]=0; mark[MF] = 0;
	  n_protein ++;
	}

      if (mark[ot] ==0) {n_type[ot] ++; mark[ot] = 1;}

      add_association(iid, prot);
    }
  fclose(f);
}
