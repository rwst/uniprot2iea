#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <search.h>

#define BUFLEN 32768

enum Fields {AC = 0, DEF, DEC, DES, N_FIELDS};

void dummy() {}

typedef struct backend
{
  (*void)(record_start)();
  (*void)(record_end)();
  (*void)(field_start)();
  (*void)(field_end)();
} BACKEND;

BACKEND flat={dummy, dummy, dummy, dummy};

typedef struct node
{
  char *str;
  struct node *ptr;
} NODE;

typedef struct my_entry
{
  char key[16];
  NODE *list;
} MENTRY;

FILE *fp;
char *bufadr, *ptr;
MENTRY *entries;
int ecount, eco=0;

int main (int n, char **argv[])
{
  int count = 0, nr = 0, np = 0, cseq = 0, i,k,j;
  char *buf, *str[N_FIELDS], *a, *b;

  count = 0;
  if (n!=2 || argv[1][0]=='-')
  {
    printf ("uniprot2iea (c) R Stephan <ralf@ark.in-berlin.de>\n");
    printf ("Released under Artistic/GPL, see http://dev.perl.org/licenses\n");
    printf ("Home: http://code.google.com/p/uniprot2iea/\n");
    printf ("Usage: uniprot2iea filename\n\n");
  }
  printf ("Trying to read %s\n", argv[1]);
  fp = fopen (argv[1], "r");
  while ((buf = fgets (bufadr, BUFLEN-1, fp)) >0)
  {
    if (*buf == '#')
      continue;
    if (buf[0]=='A' && buf[1]=='C')
    {
      ++count;
      for (i=0; i<N_FIELDS; ++i) str[i]=NULL;
      a = buf+5;
      ptr = a;
      while (*ptr && *ptr!=';') ++ptr;
      str[AC] = stralloc (a, ptr);
    }
    if (buf[0]=='D' && buf[1]=='E')
    {
      if ((a=strstr(buf,"RecName Full="))!=NULL)
      {
         a += 13;
         ptr = a;
         while (*ptr && *ptr!=';') ++ptr;
         str[DEF] = stralloc (a, ptr);
      }
      if ((a=strstr(buf,"EC="))!=NULL)
      {
         a += 3;
         ptr = a;
         while (*ptr && *ptr!=';') ++ptr;
         str[DEC] = stralloc (a, ptr);
      }
      if ((a=strstr(buf,"Short="))!=NULL)
      {
         a += 6;
         ptr = a;
         while (*ptr && *ptr!=';') ++ptr;
         str[DES] = stralloc (a, ptr);
      }
    }
    if (buf[0]=='G' && buf[1]=='N')
    {
      if ((a=strstr(buf,"Name="))!=NULL)
      {
         a += 5;
         ptr = a;
         while (*ptr && *ptr!=';') ++ptr;
         str[GNN] = stralloc (a, ptr);
      }
      if ((a=strstr(buf,"OrderedLocusName="))!=NULL)
      {
         a += 17;
         ptr = a;
         while (*ptr && *ptr!=';') ++ptr;
         str[DEC] = stralloc (a, ptr);
      }
    }
    if (buf[0]=='D' && buf[1]=='R')
    {
      if ((a=strstr(buf,"GO:"))!=NULL)
      {
         a += 3;
         ptr = a;
         while (*ptr && *ptr!=';') ++ptr;
         ++ptr;                            /* add description */
         while (*ptr && *ptr!=';') ++ptr;
         add2list (golist, stralloc (a, ptr));
      }
      if ((a=strstr(buf,"InterPro;"))!=NULL)
      {
         a += 10;
         ptr = a;
         while (*ptr && *ptr!=';') ++ptr;
         add2list (iplist, stralloc (a, ptr));
      }
      if ((a=strstr(buf,"HAMAP;"))!=NULL)
      {
         a += 7;
         ptr = a;
         while (*ptr && *ptr!=';') ++ptr;
         add2list (halist, stralloc (a, ptr);
      }
    }
  }
  printf ("Number of groups: %d\n", count);
  printf ("Number of groups with MYCT only: %d\n", np);
  fclose (fp);
  }
