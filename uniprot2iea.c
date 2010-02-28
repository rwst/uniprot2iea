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
void gaf1start() { fprintf (stdout, "!gaf-version: 1.0\n"); }
void gaf1rstart() { fprintf (stdout, "UniProtKB/Swiss-Prot\t"); }
void gaf1rend() { fprintf (stdout, "\n"); }
void gaf1fend() { fprintf (stdout, "\t"); }
void fpf (const char* str) { fprintf (stdout, "%s", str); gaf1fend(); }

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
  void (*goref)(const char*);
  void (*ref)(const char*);
  void (*gotype)(const char*);
  void (*goann)(const char*);    /* writes ann. type, reason */
  void (*upid)(const char*);
  void (*ec)(const char*);
} BACKEND;

BACKEND gaf1={gaf1rstart, gaf1rend, dummy, gaf1fend, fpf, fpf, fpf, fpf, fpf, fpf, fpf, fpf, fpf, fpf, fpf, fpf};

typedef struct node
{
  char *str;
  struct node *ptr;
} NODE;

NODE *iplist = NULL, *halist = NULL;

char *stralloc (const char* p1, const char* p2)
{
  char *p = malloc (p2-p1+1);
  if (p==NULL) { fprintf (stderr, "Memory full.\n"); exit(1); }
  strncpy (p, p1, p2-p1);
  p[p2-p1] = '\0';
  return p;
}

void add2list (NODE** np, char *str)
{
/*fprintf(stderr, "add2list %s\n", str);*/
  if (*np == NULL)
  {
     *np = malloc (sizeof(NODE));
     if (*np==NULL) { fprintf (stderr, "Memory full.\n"); exit(1); }
     (*np)->ptr = NULL;
     (*np)->str = str;
  }
  else
  {
     while ((*np)->ptr != NULL)
       (*np) = (*np)->ptr;
     (*np)->ptr = malloc (sizeof(NODE));
     if ((*np)->ptr==NULL) { fprintf (stderr, "Memory full.\n"); exit(1); }
     *np = (*np)->ptr;
     (*np)->ptr = NULL;
     (*np)->str = str;
  }
}

void freelist (NODE **np)
{
  NODE *p;
  if (*np == NULL) return;
  do {
    if ((*np)->str != NULL) free ((*np)->str);
    p = *np;
    *np = (*np)->ptr;
    free (p);
  }
  while (*np != NULL);
  *np = NULL;
}

const char *conv_goref (const char *str)
{
  static char *tab[] = {"EC", "HAMAP", "InterPro", "UniProtKB-KW", "UniProtKB-SubCell"};
  static char *ref[] = {"GO_REF:0000003", "GO_REF:0000020", "GO_REF:0000002", "GO_REF:0000004", "GO_REF:0000023" };
  int i;
  for (i=0; i<sizeof(tab)/sizeof(char *); ++i)
    if (!strcmp (str, tab[i]))
      return ref[i];
  fprintf (stderr, "error: unknown GO_REF: %s\n", str);
  return NULL;
}

int main (int n, char **argv)
{
  int count = 0, ac = 0, i,k,j, mode=WM_GOONLY;
  char *buf, *str[N_FIELDS], *a, *b, *ptr;
  FILE *fp;
  BACKEND be = gaf1;

  count = 0;
  for (i=0; i<N_FIELDS; ++i)
    str[i]=NULL;
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
      if ((a=strstr(buf,"OrderedLocusNames="))!=NULL)
      {
         a += 18;
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
    if (!strncmp (buf, "DR   GO;", 8))
    {
      char *go, *got, *gor, *goa;
      char *a = buf+9, *ptr = buf+19;
      go = stralloc (a, ptr);
      ++ptr; ++ptr;
      got = stralloc (ptr, ptr+1);
      ++ptr; ++ptr;
      while (*ptr && *ptr!=';') ++ptr;
      ++ptr;
      a = ++ptr;
      while (*ptr && *ptr!=':') ++ptr;
      goa = stralloc (a, ptr);
      a = ++ptr;
      while (*ptr && *ptr!='.') ++ptr;
      gor = stralloc (a, ptr);
/*      if ((a=strstr(buf,"InterPro;"))!=NULL)
      {
         a += 10;
         ptr = a;
         while (*ptr && *ptr!=';') ++ptr;
         add2list (&iplist, stralloc (a, ptr));
      }
      if ((a=strstr(buf,"HAMAP;"))!=NULL)
      {
         a += 7;
         ptr = a;
         while (*ptr && *ptr!=';') ++ptr;
         add2list (&halist, stralloc (a, ptr));
      }*/
    /* all fields for the entry are collected now, write them out */
      if (mode == WM_GOONLY)
      {
        be.record_start();
        be.upid (str[AC]);
        be.lname (str[GNN]!=NULL?str[GNN]:str[GNL]);
        be.go (go);
        be.goref (conv_goref (gor));
        be.goann (goa);
        //be.ref ();
        be.gotype (got);
        be.fname (str[DEF]);
        be.sname (str[DES]);
        be.oname (str[OS]);
        be.record_end();
        ++ac;
      }
      free (go);
      free (gor);
      free (got);
      free (goa);
    }

  }
  printf ("Number of proteins read: %d\n", count);
  printf ("Number of annotation records written: %d\n", ac);
  fclose (fp);
}
