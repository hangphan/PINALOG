#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "utility.h"

#ifndef __TrieTree__
#define __TrieTree__


#define SearchError -2
#define SearchNonExist -5

class CTrieNode{

  char label; //label of the node
  long  index; //index of the node; -1: no index -> index of the protein in the database
  long  id;    // id of the node in the prefix tree
  long  nc; //number of children
  long  ns; //size of child
  CTrieNode **child;

public:
  CTrieNode();
  ~CTrieNode();
  long Exist(char *s);
  char Label() { return label;}
  long  ID() {return id;}
  void SetUp(char c) { label = c;}
  void SetUpIndex(long nid) { index = nid;}
  void SetUpID(long xid) {id = xid;}
  void ExtendSize(void); //Extend size of children set by 10 
  void Add(char *s, long nid, long *xid);
  void Add(long nid, char lab, long numchild, long  *children);
  void Add(long nid, char lab, long numchild) {index = nid; label = lab; nc = numchild; }
  void printNode(FILE *f);
  CTrieNode *Search(long xid) ;
} ;

class CTrieTree{

protected:
  #define MaxLen 128
  char buffer[MaxLen];
  long n; //number of node having index
  long nn; // number of nodes
  CTrieNode *R;

  void adddot(char *s);

public:
  CTrieTree();
  ~CTrieTree();
  void Read(FILE *f);     //*Read a list of protein name and construct a prefix tree from this
  void LoadTrie(FILE *f); //*Load a TrieTree that has been written down by PrintTrie();
  long Search(char *s);
  long Add(char *s);
  void PrintTrie(FILE *f);
  long  Size() {return n;} //get number of protein encoded by the tree
};


#endif
