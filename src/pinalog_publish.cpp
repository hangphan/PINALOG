#include"pinalog_publish.h"
char ontology_file[] = "gene_ontology.1_2.obo";
char gene_association_file[] = "intact.gene_assocs.txt";
//Option =1: mapping using only protein sequence similarity
//Option =2: mapping using both protein function and protein sequence similarity

int main(int ArgI, char** ArgC)
{
  long option =1;

  if (ArgI > 4) option =2;
  if (option  ==1) printf("PINALOG alignment using protein sequence similarity \n");
  if (option  ==2) printf("PINALOG alignment using protein sequence and function similarity \n");

  A.Read(ArgC[1], 1);
  B.Read(ArgC[2], 1);

  nodesA = A.get_cnodes();
  nodesB = B.get_cnodes();
  sizeA = A.get_num_nodes();
  sizeB = B.get_num_nodes();
  printf("option = %d\n",option);
  printf("Finish reading graph %s   %d proteins  %d interactions \n",  ArgC[1], sizeA, A.get_num_edges());
  printf("Finish reading graph %s   %d proteins  %d interactions \n",  ArgC[2], sizeB, B.get_num_edges());
  
  make_communities_files1(ArgC[1], ArgC[2]);

  long templ;

  for(long i =0; i!= sizeA; ++i) {  templ = allnodes.Add(nodesA->GetLabel(i)); gg1.push_back(templ);}
  for(long i =0; i!= sizeB; ++i) {  templ = allnodes.Add(nodesB->GetLabel(i)); gg2.push_back(templ);}
  
  AA = Alloc1(sizeA);
  if (AA ==NULL) return 1;
  BB = Alloc1(sizeB);
  if (BB ==NULL) return 1;
  AB = Alloc2(sizeA, sizeB);
  if (AB ==NULL) return 1;
  printf("Finish allocating space for  similarity matrix \n");

  ReadBlast(ArgC[3]);
  printf("Finish reading blast score. Calculating sequence/function ratio ...\n");

  theta = GetSeqFuncRatio();
  printf("Ratio = %f\n Normalizing similarity matrix \n", theta);
  NormalizeDistMatrix();

  if (option == 2)
    {
      Ont.read_go_ids(ontology_file);
      Ont.read_relations(ontology_file);
      printf("Finish reading ontology\n");
      
      Ont.add_protein_associations(gene_association_file, allnodes.get_trie());
      FILE *fi = fopen(ArgC[4], "r");
      allnodes.add_gene_associations_isoform(fi, &Ont);
      fclose(fi);
      printf("Finish adding protein associations\n");
      //    fo = fopen("net1_net2.dist", "w");
      AddFuncSim();
      //      fclose(fo);

    }
  
  char outfile1[128], outfile2[128];
  sprintf(outfile1, "net1_net2.pinalog.nodes_algn.txt");
  sprintf(outfile2, "net1_net2.pinalog.edges_algn.txt");

  if (option ==1 ) mapping();
  if (option ==2 ) mapping();

  filter_results();
  write_aligned_pairs(outfile1);
  write_conserved_edges(outfile2);
  

  Free(sizeA, AA);
  Free(sizeB, BB);
  Free(sizeA, AB);
  return 0;
}
