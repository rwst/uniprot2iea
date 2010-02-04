#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <search.h>

#define BUFLEN 32768
#define WM_GOONLY 1
#define WM_WITHGOA 2

char bufadr[BUFLEN];
enum Fields {AC = 0, DEF, DEC, DES, GNN, GNL, OS, N_FIELDS};

void dummy() {}
void fpf (const char* str) { fprintf (stdout, "%s\n", str); } 

typedef struct backend
{
  void (*record_start)();
  void (*record_end)();
  void (*field_start)();
  void (*field_end)();
  void (*lname)(const char*);
  void (*fname)(const char*);
  void (*sname)(const char*);
  void (*gname)(const char*);
  void (*oname)(const char*);
  void (*go)(const char*);     /* writes GO no+desc, type (C,P,F) */
  void (*goo)(const char*);    /* writes ann. type, reason */
  void (*upid)(const char*);
  void (*ec)(const char*);
} BACKEND;

BACKEND flat={dummy, dummy, dummy, dummy, fpf, fpf, fpf, fpf, fpf, fpf, fpf, fpf, fpf};

typedef struct node
{
  char *str;
  struct node *ptr;
} NODE;

NODE *golist = NULL, *goolist = NULL, *iplist = NULL, *halist = NULL;

char *stralloc (const char* p1, const char* p2)
{
  char *p = malloc (p2-p1+1);
  if (p==NULL) { fprintf (stderr, "Memory full.\n"); exit(1); }
  strncpy (p, p1, p2-p1);
  p[p2-p1] = '\0';
  return p;
}

void add2list (NODE* np, char *str)
{
  if (np == NULL)
  {
     np = malloc (sizeof(NODE));
     if (np==NULL) { fprintf (stderr, "Memory full.\n"); exit(1); }
     np->ptr = NULL;
     np->str = str;
  }
  else
  {
     while (np->ptr != NULL)
       np = np->ptr;
     np->ptr = malloc (sizeof(NODE));
     if (np->ptr==NULL) { fprintf (stderr, "Memory full.\n"); exit(1); }
     np = np->ptr;
     np->ptr = NULL;
     np->str = str;
  }
}

void freelist (NODE *np)
{
  NODE *p;
  if (np == NULL) return;
  do {
    if (np->str != NULL) free (np->str);
    p = np;
    np = np->ptr;
    free (p);
  }
  while (np != NULL);
}

int main (int n, char **argv)
{
  int count = 0, ac = 0, i,k,j, mode=WM_GOONLY;
  char *buf, *str[N_FIELDS], *a, *b, *ptr;
  FILE *fp;
  BACKEND be = flat;

  count = 0;
  if (n!=2 || argv[1][0]=='-')
  {
    printf ("uniprot2iea (c) R Stephan <ralf@ark.in-berlin.de>\n");
    printf ("Released under Artistic/GPL, see http://dev.perl.org/licenses\n");
    printf ("Home: http://code.google.com/p/uniprot2iea/\n");
    printf ("Usage: uniprot2iea filename\n\n");
    exit(1);
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
      for (i=0; i<N_FIELDS; ++i)
      {
        if (str[i]!=NULL)
          free (str[i]);
        str[i]=NULL;
      }
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
         str[GNL] = stralloc (a, ptr);
      }
    }
    if (buf[0]=='O' && buf[1]=='S')
    {
       a = buf+3;
       ptr = a;
       while (*ptr && *ptr!=';') ++ptr;
       str[OS] = stralloc (a, ptr);
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
         a = ++ptr;
         while (*ptr && *ptr!=';') ++ptr;
         add2list (goolist, stralloc (a, ptr));
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
         add2list (halist, stralloc (a, ptr));
      }
    }
    /* all fields for the entry are collected now, write them out */
    be.record_start();
    if (mode == WM_GOONLY)
    {
      NODE *np, *gnp;
      for (np = golist, gnp = goolist;
                (np!=NULL && (np->str)!=NULL);
                np = np->ptr, gnp = gnp->ptr)
      {
        be.lname (str[GNL]);
        be.upid (str[AC]);
        be.go (np->str);   /* writes GO no+desc, type (C,P,F) */
        be.goo (gnp->str); /* writes ann. type, reason */
        be.ec (str[DEC]);
        be.fname (str[DEF]);
        be.sname (str[DES]);
        be.gname (str[GNN]);
        be.oname (str[OS]);
        ++ac;
      }
    }
    be.record_end();

    freelist (golist);
    freelist (goolist);
    freelist (iplist);
    freelist (halist);
  }
  printf ("Number of proteins read: %d\n", count);
  printf ("Number of annotation records written: %d\n", ac);
  fclose (fp);
}
