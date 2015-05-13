#ifndef __UTILITY__
#define __UTILITY__

#define NUM_DB 11
#include"TrieTree.h"
#include <math.h>

#define PI 3.14159265358979323846264338327950288419716939937510

#ifndef min
#define min(x, y)       ((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define max(x, y)       ((x) > (y) ? (x) : (y))
#endif

typedef std::vector<long> VLong;
typedef std::vector<VLong> VVLong;
typedef std::vector<double> VDouble;
typedef std::vector<VDouble> VVDouble;

long isin(long val,  VLong x)
{
  for(long i=0; i!= x.size(); ++i)
    if (val == x[i]) return  i;
  return -1;
}

//binary search
long isIn(long id, VLong x)
{
  if (x.size() ==0) return -1;
  if (x.size() ==1 && id != x[0]) return -1;
  if (x.size() ==1 && id == x[0]) return 0;

  long low = 0;
  long high = x.size();
  long mid;
  while (low < high) 
    {
      mid = (low + high) / 2;
      if (x[mid] < id) low = mid + 1;
      else high = mid;
    }
  //  printf("is In\n");
  if (x[low] == id) 
    return low;
  else   return -1;
}


void init_vector(VLong &x, long size, long init_val)
{
  if (!x.empty()) x.clear();
  for (long i=0; i!= size; ++i) x.push_back(init_val);
}

void init_vector(VVLong &x, long size)
{
  if (x.size() >0) x.clear();
  for(long i =0; i!= size; ++i) 
    {
      VLong templ;
      x.push_back(templ);
    }
}


int cmp(const void *x, const void *y)
{
  const long *px  = (const long *)x;
  const long *py  = (const long *)y;
  return (*px) > (*py);
}

#endif
