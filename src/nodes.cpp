#include "nodes.h"

 //add only those at the bottom level 
void CNode::add_goid1(long goid, COntology_RW *ontology)
{

 VLong gos = ontology->get_ancestors(goid);

 for(long k=0; k!= gos.size(); ++k)
   {
     for(long l = 0; l!= goids.size(); ++l)
       if (gos[k] == goids[l])
	 {
	   replace_goid(goid, l);
	   return;
	 }
   }
 
 add_goid(goid);
 return;
}

  // add all ancestors
void CNode::add_goid2(long goid, COntology_RW *ontology)
{

 VLong gos = ontology->get_ancestors(goid);

 for(long k=0; k!= gos.size(); ++k)
   add_goid(gos[k]);
 
 add_goid(goid);
 return;
}

void CNode::replace_goid(long goid, long id)
{
  goids[id] = goid;
}

//add goid as it is
void CNode::add_goid(long goid)
{
  for(long i =0; i!= goids.size(); ++i)
    if (goids[i] == goid) return;
  goids.push_back(goid);
}


CNodes::CNodes(long M)
{
  mode = M; //mode when get association GO terms : =0 : add only those at the bottom, =1 add all ancestors of the term
  n = 0;
}

CNodes::CNodes()
{

  mode = 0;
  nodes = NULL;
  n = 0;
}

CNodes::~CNodes()
{
  for(long i = 0; i != n; ++i) 
   if (nodes[i] != NULL) delete nodes[i];
  if (nodes!= NULL)  {free(nodes); nodes = NULL;}

}

long CNodes::Add(char *lab)
{
  long id = Search(lab);
  if (id > -1) return id ;
 
  Trie.Add(lab); //Add new label to TrieTree
 
 //Add new node to nodes

  n++;
  nodes = (CNode**) realloc(nodes, n* sizeof(CNode*));

  nodes[n-1] = new CNode;	     
  nodes[n-1]->SetupNode(lab, n-1);
  return (n-1);
}


void CNodes::Read(FILE *f)
{
  char temps[1024];
  while(!feof(f))
    {
      fscanf(f, "%s", temps);
      Add(temps);
    }
}

void CNodes::ReadFromGraph(FILE *f)
{
  char st1[1024], st2[1024];
  while(!feof(f))
    {
      fscanf(f, "%s %s", st1, st2);
      Add(st1);
      Add(st2);
      fgets(st1, 1024, f);
    }
}

void CNodes::add_gene_associations(FILE *f, COntology_RW *ontology)
{

  long size_buff = 50000;
  char buffer[size_buff];
  char st1[128], st2[128], st3[128], st4[128], st5[128];
  while (!feof(f))
    {
      fgets(buffer, size_buff, f);
      sscanf(buffer, "%s %s %s %s %s", st1, st2, st3, st4, st5);
      long node_id = Trie.Search(st2);
      if (node_id <0) continue;

      //st4 can either be the GO term itself or the qualifier of the GO term, which can be 
      //"NOT", "colocolizes_with", "contributes_to"
      long goid;
      if (strstr(st4, "GO:") != NULL)
	  goid= ontology->Search(st4);
      else if (strstr(st5, "GO:") != NULL)
	goid= ontology->Search(st5);
      
      if (strcmp(st4, "NOT") ==0) continue;
      if (goid <0) continue;

      if (mode ==0)
	nodes[node_id]->add_goid1(goid, ontology);
      if (mode ==1) 
	nodes[node_id]->add_goid2(goid, ontology);
      continue;
    }
}
void CNodes::add_gene_associations1(FILE *f, COntology_RW *ontology)
{

  long size_buff = 50000;
  char buffer[size_buff];
  char st1[128], st2[128], st3[128], st4[128], st5[128];
  while (!feof(f))
    {
      fgets(buffer, size_buff, f);
      sscanf(buffer, "%s %s %s %s", st1, st2, st3, st4);
      long node_id = Trie.Search(st2);
      if (node_id <0) continue;

      //st4 can either be the GO term itself or the qualifier of the GO term, which can be 
      //"NOT", "colocolizes_with", "contributes_to"
      long goid;
      if (strstr(st3, "GO:") != NULL)
	  goid= ontology->Search(st3);
      else if (strstr(st4, "GO:") != NULL)
	goid= ontology->Search(st3);
      
      if (goid <0) continue;

      if (mode ==0)
	nodes[node_id]->add_goid1(goid, ontology);
      if (mode ==1) 
	nodes[node_id]->add_goid2(goid, ontology);
      continue;
    }
}

void CNodes::add_gene_associations_isoform(FILE *f, COntology_RW *ontology)
{
  long size_buff = 50000;
  char buffer[size_buff];
  char st1[128], st2[128], st3[128], st4[128], st5[128];
  char *p1, *p2;
  long goid, node_id;

  while (!feof(f))
    {
      fgets(buffer, size_buff, f);
      sscanf(buffer, "%s %s %s %s", st1, st2, st3, st4);
      p1 = p2 =  NULL;
      p1 = strstr(st4, "GO:"); 
      p2 = strstr(st5, "GO:");

      if ( p1 != NULL)	  goid= ontology->Search(st4);
      else 
	if (p2!= NULL)	goid= ontology->Search(st5);

      p1 = strstr(st4, "NOT"); 
      if (p1  != NULL)  continue;
      if (goid <0) continue;

      node_id = Trie.Search(st2);
      if (node_id >= 0) 
	{
	  if (mode ==0)   nodes[node_id]->add_goid1(goid, ontology);
	  if (mode ==1)   nodes[node_id]->add_goid2(goid, ontology);
	}
      for(long i =0; i!= 20; ++i)
	{
	  char temps[128];
	  sprintf(temps, "%s-%d", st2, i);
	  node_id = Trie.Search(temps);
	  //	  printf("%s\t%d\n", temps, node_id);
	  if (node_id >=0) 
	    {
	      if (mode ==0)
		nodes[node_id]->add_goid1(goid, ontology);
	      if (mode ==1) 
		nodes[node_id]->add_goid2(goid, ontology);
	    }
	  
	}


    }
}

void CNodes::update_ontology_protein_count(COntology_RW *ontology)
{
  long nterms = ontology->get_n_terms();
  for(long k=0; k!= nterms; ++k)
    {
      long count =0;
      for(long i =0; i!= n; ++i)
	{

	  VLong goids = get_goids_from_node(i);
	  for(long j=0; j!= goids.size(); ++j)
	    if (goids[j] == k)
	      {
		count++;
		break;
	      }
	}
      ontology->update_protein_count(k, count);
    }

}



VLong CNodes::get_goids_from_node(long id)
{
  return nodes[id]->get_goids();
}

VLong CNodes::get_associated_prots(long goid)
{
  VLong temp;
  for(long i =0; i!= n; ++i)
    {

      VLong goids = get_goids_from_node(i);
      for(long j=0; j!= goids.size(); ++j)
	if (goids[j] == goid)
	  {
	    temp.push_back(i);
	    break;
	  }
    }
  return temp;
}
