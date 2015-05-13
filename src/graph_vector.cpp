#include "graph_vector.h"

void CGraph::Read(char *fn, long option) 
{
  FILE * fi = fopen(fn, "r");
  Read(fi, option);
  fclose(fi);
}

 void CGraph::Read(FILE *f, long option)
 {
   long buff_size = 15000;
   char buffer[buff_size];

   while(!feof(f))
     {
       strcpy(buffer, "");
       fgets(buffer, buff_size, f);
       if (strlen(buffer) <1) break;
       LineProcess(buffer, option);
     }
 }

 void CGraph::LineProcess(char *buffer, long option)
 {
   char s1[MaxLen], s2[MaxLen], s3[MaxLen], s4[MaxLen], ns1[MaxLen], ns2[MaxLen];


   if (option ==0)
     {
       sscanf(buffer, "%s %s %s %s", s1, s2, s3, s4);//node info is in s3 and s4
       GetProteinName(s3, ns1);
       GetProteinName(s4, ns2);
     }
   else 
     {
       sscanf(buffer, "%s %s ", s1, s2);//node info is in s1 and s2
       strcpy(ns1, s1);
       strcpy(ns2, s2);
     }

   long id1 =  nodes.Add(ns1);
   long id2 =  nodes.Add(ns2);
   long size1 = max(id1, id2) - adlist.size() +1;
   if (size1 >0) 
     {
       VLong tempV;
       for(long i =0; i!= size1; ++i)	 adlist.push_back(tempV);
     }

   if (isIn(id2, adlist[id1]) <0) {adlist[id1].push_back(id2); sort(adlist[id1].begin(), adlist[id1].end());}
   if (isIn(id1, adlist[id2]) <0) {adlist[id2].push_back(id1); sort(adlist[id2].begin(), adlist[id2].end());}
 }

long CGraph::get_num_edges()
{
  if (num_edges>0) return num_edges;

  for(long i=0; i!= adlist.size(); ++i)
    {
      for(long j=0; j!= adlist[i].size(); ++j)
	if (adlist[i][j] >= i) num_edges++;
      
    }
  printf("in here\n");
  return num_edges;
}

void CGraph::GetProteinName(char *st, char *ptn)
{
  char *p1 = st;
  char temp1[MaxLen], temp2[MaxLen];
  sscanf(p1, temp1, '|');        //get the first accession number in the list of them
  
  char *p2 = strstr(temp1, ":"); // find the ":" character, if there is, get the accession number after it, else it is the accession number with no database info
  
  long length=0;                  //length of protein name
  
  if (p2 == NULL) p2 = temp1;
  else p2++;

  for(char *p = p2; p!= temp1 + strlen(temp1); ++p) ptn[length++] = p[0];
  ptn[length] = '\0';
  if (length > 30 || length == 0) 
    printf("%s length = %d\n", ptn, length);
}

long CGraph::IsNeighbour(long idx1, long idx2)//local ids
{
  if (adlist[idx1].size() < adlist[idx2].size())
    return isIn(idx2, adlist[idx1]);
  else
    return isIn(idx1, adlist[idx2]);
  
}



VLong CGraph::Get2Neigh(long gid)
{
  VLong out;
  VLong firstNei = GetNeigh(gid);

  if (firstNei.empty()) return out;
  
  for(long i =0; i!= firstNei.size(); ++i)
    {
      VLong temp = GetNeigh(firstNei[i]);
      for(long j = 0; j != temp.size(); ++j) 
	if (temp[j] != gid && IsNeighbour(gid, temp[j]) <0)
	  out.push_back(temp[j]);
    }
  return out;
}




void CGraph::Write2GraphML(char *fn)
{
  FILE *f = fopen(fn, "w");
  Write2GraphML(f);
  fclose(f);

}
void CGraph::Write2GraphML(FILE *f)
{
  fprintf(f, "<graphml>\n");
  fprintf(f, "  <graph>\n");
  
  long size = nodes.GetNumNodes();

  for(long i=0; i!= size; ++i)
    fprintf(f, "    <node id = \"%s\"/>\n", nodes.GetLabel(i));
  
  for(long i=0; i!= adlist.size(); ++i)
    for(long j =0; j!= adlist[i].size(); ++j)
      if (adlist[i][j] >= i)
	fprintf(f, "    <edge source=\"%s\" target=\"%s\" />\n", nodes.GetLabel(i), nodes.GetLabel(adlist[i][j]));
    fprintf(f, "  </graph>\n");
  fprintf(f, "</graphml>\n");
}


void CGraph::Write2HTML(char *fn)
{
  FILE *f = fopen(fn, "w");
  Write2HTML(f);
  fclose(f);



}


void CGraph::Write2HTML(FILE *f)
{
  fprintf(f, "<html>\n");
  fprintf(f, "  <head>\n");
  fprintf(f, "      <script type=\"text/javascript\" src=\"AC_OETags.min.js\"></script>\n");
  fprintf(f, "      <script type=\"text/javascript\" src=\"json2.min.js\"></script>\n");
  fprintf(f, "      <script type=\"text/javascript\" src=\"cytoscapeweb.min.js\"></script>\n");
  fprintf(f, "  </head>\n\n");
        
  fprintf(f, "  <body>\n");
  fprintf(f, "     <div id=\"Conservednet\" style=\"width:800px;height:600px;\"></div>\n");
  fprintf(f, "  </body>\n\n");
   
  fprintf(f, "  <script type=\"text/javascript\">\n");
  fprintf(f, " // network data could alternatively be grabbed via ajax\n");
  fprintf(f, "       var xml = '\\\n");
  fprintf(f, "         <graphml>\\\n");
  fprintf(f, "            <graph>\\\n");
  
  long size = nodes.GetNumNodes();
  if (node_scores.size() != 0)
    for(long i=0; i!= size; ++i)
      fprintf(f, "                <node id = \"%s\" score %lf/>\\\n", nodes.GetLabel(i), node_scores[i]);
  else
    for(long i=0; i!= size; ++i)
      fprintf(f, "                <node id = \"%s\"/>\\\n", nodes.GetLabel(i));

  for(long i=0; i!= adlist.size(); ++i)
    for(long j =0; j!= adlist[i].size(); ++j)
      if (adlist[i][j] >= i)
	fprintf(f, "            <edge source=\"%s\" target=\"%s\" />\\\n", nodes.GetLabel(i), nodes.GetLabel(adlist[i][j]));
  fprintf(f, "             </graph>\\\n");
  fprintf(f, "          </graphml>\\\n");
  fprintf(f, "         ';\n");

  fprintf(f, "      // init and draw\n");
  fprintf(f, "      var vis = new org.cytoscapeweb.Visualization(\"Conservednet\");\n");;
  fprintf(f, "      vis.draw({ network: xml });\n");
  fprintf(f, "  </script>\n");
  fprintf(f, "</html>\n");

}


