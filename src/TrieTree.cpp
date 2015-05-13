#include "TrieTree.h"

CTrieNode::CTrieNode()
{
  nc = 0;
  label = '.';
  index = -1;
  id = -1;
  ns = 0;
  ExtendSize();
}

CTrieNode::~CTrieNode()
{
  //  printf("Out of TrieNode\n");
  
  if (child!= NULL) 
    {
      for (long i = 0; i != nc; ++i) if (child[i] != NULL) delete child[i];
      free(child);
      child = NULL;
    }
}

CTrieNode *CTrieNode::Search(long xid) 
{
  if (id == xid)   {return this;}
  for(long i =0; i!= nc; ++i)  
    {
      CTrieNode *p = child[i]->Search(xid);
      if ( p != NULL) return p;
    }
  return NULL;
}
 

void CTrieNode :: ExtendSize(void)
{
  /*
  child =(CTrieNode**) realloc((CTrieNode**) child, (ns+10) * sizeof(CTrieNode*));
  for (long i = ns; i != ns+10; ++i) child[i] = NULL;
  ns+= 10;
  */
  CTrieNode **p = (CTrieNode**)calloc(ns+10, sizeof(CTrieNode *));
  for (long i = 0; i != ns+10; ++i) p[i] = NULL;
  
  for (long i = 0; i != nc; ++i)
    p[i] = child[i];
  
  if (ns >0)
    {
      free(child);
    }
  child = p;
  ns+= 10;
  
}


void CTrieNode::Add(long nid, char lab, long numchild, long *children)
{
  Add(nid, lab, numchild);
  child = (CTrieNode **) calloc(nc, sizeof(CTrieNode *));
  //while(nc >ns+1) ExtendSize();

  for(long i = 0; i!= nc; ++i)
    child[i] = (CTrieNode *) calloc(1, sizeof(CTrieNode));

  for(long i = 0; i!= nc; ++i)
    {
      child[i]->SetUpID(children[i]);
      //      printf("%d\t", child[i]->ID());
    }

  //  printf("\n");
}

void  CTrieNode::Add(char *s, long nid, long *xid) //nid =   index  of new node ; xid = node id 
{
  if (strlen(s) == 0) {return;}


  if (s[0] != label)
    {
      printf("error: %c %c\n", s[0], label);
      exit(1);
    }

  
  if (id == -1) 
    {
      id = *xid;
      *xid  = (*xid)+1;
    }
  if (strlen(s) == 1)
    {
      if (index != -1) //
	{
	  printf("error in searching or something wrong\n");
	  exit(1);
	}

      index = nid;
      return;
    }

  for (long i = 0; i != nc; ++i)
    if (child[i]->Label() == s[1]) 
      {
	child[i]->Add(&s[1], nid, xid);
	return;
      }

  if (nc +1 >= ns)  ExtendSize();
  child[nc] = new CTrieNode;

  child[nc]->SetUp(s[1]);
  
  if (strlen(s) == 1) child[nc]->SetUpIndex(nid);
  else
    {
      child[nc]->Add(&s[1], nid, xid);
    }
  nc++;
}


long CTrieNode::Exist(char* s)
{
  
  if (strlen(s) == 0) return SearchError;
  if (s[0] != label) return SearchNonExist;

  if (strlen(s) == 1)
    return index;

  for (long i = 0; i != nc; ++i)
    {
      long x = child[i]->Exist(&s[1]);
      if (x >= 0) return x;
    }

  return SearchNonExist;
}


void CTrieNode::printNode(FILE *f)
{
  fprintf(f, "%d\t%d\t%c\t%d\t", id, index, label, nc);
  for(long i = 0; i!= nc; ++i) fprintf(f, "%d\t", child[i]->ID());
  fprintf(f, "\n");
  if (nc ==0) return;
  for(long i = 0; i!= nc; ++i) child[i]->printNode(f);
}

void CTrieTree::Read(FILE *f)
{
  //input file contain protein names only
  char buffer[MaxLen];
  while (!feof(f))
    {
      fscanf(f, "%s", buffer); 
      Add(buffer);
    }
}

void CTrieTree::LoadTrie(FILE *f)
{
  long start = 0;
  R->SetUpID(0);

  while (!feof(f))
    {
      long xid, nid;
      char lab;
      long nc;
      fscanf(f, "%d %d %c %d", &xid, &nid, &lab, &nc );
      
       //      printf("solving node %d\n", xid);

      CTrieNode *temp =  R->Search(xid);

      if (nc >0)
	{
	  long child[nc];
	  for(long i = 0; i!= nc; ++i) fscanf(f, "%d", &child[i]);
	  //add info to available node
	  temp->Add(nid, lab, nc, child);
	}
      if (nc == 0)
	temp->Add(nid, lab, nc); //add info to available node
      start++;
    }
}

long CTrieTree::Search(char *s)
{
  adddot(s);
  return R->Exist(buffer);
}

long CTrieTree :: Add(char *s)
{
  adddot(s);
  long x = Search(s);
  
  if (x < 0) //not exist
    R->Add(buffer, n++, &nn);
  return n;
}

void CTrieTree::adddot(char *s)
{
  if (s[0] != '.') sprintf(buffer, ".%s", s);
  else strcpy(buffer, s);
}

void CTrieTree::PrintTrie(FILE *f)
{
  R->printNode(f);
}


CTrieTree::CTrieTree()
{
  R = new CTrieNode;
  n = 0;
  nn = 0;
}

CTrieTree::~CTrieTree()
{
  //  printf("out of TrieTree\n");      
  if (R!= NULL) delete R;
  R = NULL;
}

