/* dimacs2g.c  version 1.1; B D McKay, October 2022. */

#define USAGE "dimacs2g [-n#:#] [infile...]"

#define HELPTEXT \
" Read files of graphs in Dimacs format and write them to stdout.\n\
\n\
  -d     Use dreadnaut format (default is sparse6)\n\
  -n#:#  Specify a range of n values for output\n\
  -a\"string\"  A string to write before each graph.\n\
  -b\"string\"  A string to write after eacg graph.\n\
        -a and -b only operate for dreadnaut output;\n\
        and should be given in separate arguments.\n\
  -c     Don't copy \"c\" comments from the input.\n\
\n\
  Input files with name *.gz are ungzipped.\n"

/*************************************************************************/

#include "gtools.h" 

#if HAVE_GUNZIP && HAVE_POPEN
#if !POPEN_DEC
extern FILE *popen(const char *command, const char *type);
extern int pclose(FILE *stream);
#endif
#endif

#define GUNZIP "gunzip -c"

#define MAXCOMMENT 200
static char comment[MAXCOMMENT+3];
static int commentlen;

typedef struct 
{
   int v,w;
} vpair;

static int
nextchar(FILE *f)
{
    char s[2];

    if (fscanf(f,"%1s",s) != 1) return EOF;
    else                        return s[0];
}

static boolean
readdimacsgraph(FILE *f, sparsegraph *g)
/* Reads a graph from Bliss format into a sparse graph */
{
    int n,c;
    unsigned long ne,j;
    int haven;
    int i,v,w;
    int haveptn;
    DYNALLSTAT(vpair,elist,elist_sz);

    commentlen = 0;
    haven = 0;
    j = 0;
    while ((c = nextchar(f)) >= 0)
    {
        switch (c)
        {
        case 'c':
            commentlen = 0;
            while ((c = getc(f)) != '\n' && c != EOF)
            {
                if (commentlen < MAXCOMMENT)
                    comment[commentlen++] = c;
            }
            comment[commentlen] = '\0';
            break;

        case 'p':
            if (haven)
            {
                fprintf(stderr,"Duplicate p line\n");
                exit(1);
            }
            if (fscanf(f," edge %d %lu",&n,&ne) != 2)
            {
                fprintf(stderr,"Bad p line\n");
                return FALSE;
            }
            haven = 1;
            DYNALLOC1(vpair,elist,elist_sz,ne,"Alloc vpair");
            break;

        case 'n':
            if (!haven)
            {
                fprintf(stderr,"Missing p line\n");
                return FALSE;
            }  
            if (fscanf(f,"%d%d",&w,&v) != 2 || w < 1 || w > n)
            {
                fprintf(stderr,"Bad n line\n");
                return FALSE;
            }
            break;

        case 'e':
            if (!haven || j == ne)
            {
                fprintf(stderr,"Missing p line or too many e lines\n");
                return FALSE;
            }
            if (fscanf(f,"%d%d",&v,&w) != 2 || v < 1 || w < 1 || v > n || w > n)
            {
                fprintf(stderr,"Bad e line\n");
                return FALSE;
            }
            elist[j].v = v-1; elist[j].w = w-1;
            ++j;
            break;

        default:
            fprintf(stderr,"Unknown line %c\n",c);
            return FALSE;
        }
    }

    if (j != ne)
    {
        fprintf(stderr,"Wrong number of e lines\n");
        exit(1);
    }

    SG_ALLOC(*g,n,2*ne,"SG_ALLOC");
    g->nv = n;
    g->nde = 2*ne;

    for (i = 0; i < n; ++i) g->d[i] = 0;
    for (j = 0; j < ne; ++j) 
    {
        ++(g->d[elist[j].v]);
        ++(g->d[elist[j].w]);
    }
    g->v[0] = 0;
    for (i = 1; i < n; ++i) g->v[i] = g->v[i-1] + g->d[i-1];
    for (i = 0; i < n; ++i) g->d[i] = 0;

    for (j = 0; j < ne; ++j) 
    {
        v = elist[j].v;
        w = elist[j].w;
        g->e[g->v[v]+(g->d[v])++] = w;
        g->e[g->v[w]+(g->d[w])++] = v;
    }

    return TRUE;
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    FILE *infile;
    int j,firstarg;
    SG_DECL(g);
    size_t flen;
    boolean nocomment,nofiles,nswitch,iszip,dreadnaut,badargs;
    long nmin,nmax;
    char zcmd[550];
    char *arg,sw;
    char *prestring,*poststring;

    HELP; PUTVERSION;

    nswitch = dreadnaut = nocomment = FALSE;
    prestring = poststring = NULL;

    badargs = FALSE;
    firstarg = argc;

    for (j = 1; !badargs && j < argc; ++j)
    {
        arg = argv[j];
        if (arg[0] == '-' && arg[1] != '\0')
        {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                     SWBOOLEAN('d',dreadnaut)
                else SWBOOLEAN('c',nocomment)
                else SWRANGE('n',":-",nswitch,nmin,nmax,">E dimacs2g -n")
                else if (sw == 'a') { prestring = arg; break; }
                else if (sw == 'b') { poststring = arg; break; }
                else badargs = TRUE;
            }
        }
        else
        {
            firstarg = j;
            break;
        }
    }

    if (badargs) { fprintf(stderr,"Usage: %s\n",USAGE); exit(0); }

    if (!dreadnaut && (prestring != NULL || poststring != NULL))
        gt_abort(">E dimacs2g: -a and -b require -d\n");

    if (!nswitch) { nmin = 0; nmax = NAUTY_INFINITY-2; }

    nofiles = (argc == firstarg);
    for (j = firstarg; j < argc || nofiles; ++j)
    {
        if (nofiles)
            infile = stdin;
        else
        {
            flen = strlen(argv[j]);
            if (flen >= 3 && strcmp(argv[j]+flen-3,".gz") == 0)
            {
#if HAVE_GUNZIP && HAVE_POPEN
                if (strlen(argv[j]) > 500) gt_abort(">E file name too long\n");
                sprintf(zcmd,"%s \"%s\"",GUNZIP,argv[j]);
                if ((infile = popen(zcmd,"r")) == NULL)
                {
                    fprintf(stderr,
                       ">E dimacs2g: cannot open gunzip pipe for \"%s\"\n",
                       argv[j]);
                    gt_abort(NULL);
                }
                iszip = TRUE;
#else
                gt_abort(">E dimacs2g is not compiled with gunzip support\n");
#endif
            }
            else
            {
                if ((infile = fopen(argv[j],"r")) == NULL)
                {
                    fprintf(stderr,">E Can't open file %s\n",argv[j]);
                    gt_abort(NULL);
                }
                iszip = FALSE;
            }
        }

        if (!readdimacsgraph(infile,&g))
        {
            fprintf(stderr,">E Dimacs error in file %s\n",argv[j]);
            gt_abort(NULL);
        }
        else if (g.nv >= nmin && g.nv <= nmax)
        {
            if (dreadnaut)
            {
                if (commentlen > 0 && !nocomment)
                    printf("\"%s\"\n",comment);
                if (prestring) printf("%s\n",prestring);
                printf("n=%d $=0 g\n",g.nv);
                putgraph_sg(stdout,&g,80);
                printf("$$\n");
                if (poststring) printf("%s\n",poststring);
            }
            else
            {
                sortlists_sg(&g);
                writes6_sg(stdout,&g);
            }
        }

        if (nofiles)
            nofiles = FALSE;
        else if (iszip)
        {
#if HAVE_GUNZIP && HAVE_POPEN
            pclose(infile);
#endif
        }
        else
            fclose(infile);
    }

    exit(0);
}    
