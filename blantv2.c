//This file takes snippets of nauty to brute force
//a Traces canonical operation on two customly created graphs

#include "gtools.h"    /* which includes nauty.h, which includes stdio.h */
#include "nautinv.h"
#include "nauty.h"
#include "schreier.h"
#include "traces.h"
#include "binaryutil.h"

#define USAGE "dreadnaut [-o options]"

#define HELPTEXT \
" Enter nauty+traces test program.\n\
\n\
  -o options  - set initial options.  The parameter value is a string of\n\
                dreadnaut commands from the following set:\n\
                a,c,d,m,p,l,G,P,w,y,$,A,V,M\n\
                The effect is the same as if these commands are entered\n\
                at the beginning of the standard input.\n\
  For help within dreadnaut, use the h command.\n"

#define PM(x) ((x) ? '+' : '-')
#define SS(n,sing,plur)  (n),((n)==1?(sing):(plur))
#define WORKSIZE 60
#define FLUSHANDPROMPT do { flushline(INFILE); if (prompt) fprintf(PROMPTFILE,"> "); } while (0)


#include "sorttemplates.c"   /* define sort2ints(a,b,n) */

#define INFILE fileptr[curfile]
#define SCHREIER_DEFAULT 10

static long seed;

#if !MAXN
DYNALLSTAT(graph,g,g_sz);
DYNALLSTAT(graph,canong,canong_sz);
DYNALLSTAT(graph,savedg,savedg_sz);
DYNALLSTAT(setword,workspace,workspace_sz);
DYNALLSTAT(int,lab,lab_sz);
DYNALLSTAT(int,ptn,ptn_sz);
DYNALLSTAT(int,orbits,orbits_sz);
DYNALLSTAT(int,templab,templab_sz);
DYNALLSTAT(int,tempptn,tempptn_sz);
DYNALLSTAT(int,perm,perm_sz);
DYNALLSTAT(int,savedlab,savedlab_sz);
DYNALLSTAT(int,savedptn,savedptn_sz);
DYNALLSTAT(set,tempactive,tempactive_sz);
DYNALLSTAT(set,active,active_sz);
#else
static graph g[MAXM*1L*MAXN];
static graph canong[MAXM*1L*MAXN];
static graph savedg[MAXM*1L*MAXN];
static setword workspace[MAXM*2L*WORKSIZE];
static int lab[MAXN];
static int ptn[MAXN];
static int orbits[MAXN];
static int savedlab[MAXN],savedptn[MAXN];
static int perm[MAXN];
static int templab[MAXN];
static int tempptn[MAXN];
static int tempactive[MAXM];
static set active[MAXM];
#endif

static sparsegraph g_sg;
static sparsegraph canong_sg;
static sparsegraph savedg_sg;

static DEFAULTOPTIONS_GRAPH(options);
static DEFAULTOPTIONS_SPARSEGRAPH(options_sg);
static statsblk stats;
static int curfile;
static FILE *fileptr[MAXIFILES];
static FILE *outfile;
static char def_ext[] = DEFEXT;
static boolean firstpath;       /* used in usernode() */

DEFAULTOPTIONS_TRACES(traces_opts);
static TracesStats traces_stats;

#define TMP

#define DENSE_MODE  0
#define SPARSE_MODE 1
#define TRACES_MODE 2
#define SPARSEREP(mode) ((mode)==1||(mode)==2)
#define NOSPARSEYET(c) else if (SPARSEREP(mode)) { fprintf(ERRFILE,\
              "command %s is not implemented in the sparse case\n",c); }
#define NODENSEYET else if (!SPARSEREP(mode)) { fprintf(ERRFILE,\
              "command %c is not implemented in the dense case\n",c); }
#define NOTRACESYET if (mode==TRACES_MODE) {  fprintf(ERRFILE,\
              "command %c is not implemented for Traces\n",c); }

static int mode;

#define U_NODE  1               /* masks for u values */
#define U_AUTOM 2
#define U_LEVEL 4
#define U_TCELL 8     /* At version 2.4, usertcellproc() is gone */
#define U_REF  16
#define U_CANON 32

#ifndef  NODEPROC
#define NODEPROC usernode
#else
extern void NODEPROC(graph*,int*,int*,int,int,int,int,int,int);
#endif

#ifndef  AUTOMPROC
#define AUTOMPROC userautom
#else
extern void AUTOMPROC(int,int*,int*,int,int,int);
#endif

#ifndef  LEVELPROC
#define LEVELPROC userlevel
#else
extern void LEVELPROC(int*,int*,int,int*,statsblk*,int,int,int,int,int,int);
#endif

#ifndef  REFPROC
#define REFPROC NULL
#else
extern void REFPROC(graph*,int*,int*,int,int*,int*,set*,int*,int,int);
#endif

#ifndef  CANONPROC
#define CANONPROC usercanon
#else
extern int CANONPROC(graph*,int*,graph*,int,int,int,int);
#endif

#ifndef  INVARPROC
#define INVARPROC NULL
#define INVARPROCNAME "none"
#else
extern void INVARPROC(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
#define INVARPROCNAME "user-defined"
#endif

#ifndef  INVARPROC_SG
#define INVARPROC_SG NULL
#define INVARPROCNAME_SG "none"
#else
extern void INVARPROC_SG(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
#define INVARPROCNAME_SG "user-defined"
#endif

static struct invarrec
{
    void (*entrypoint)(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
    char *name;
    void (*entrypoint_sg)(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
    char *name_sg;
} invarproc[]
    = {{INVARPROC, INVARPROCNAME, INVARPROC_SG, INVARPROCNAME_SG},
       {NULL,        "none",        NULL,           "none"},
       {twopaths,    "twopaths",    NULL,           "unavailable"},
       {adjtriang,   "adjtriang",   NULL,           "unavailable"},
       {triples,     "triples",     NULL,           "unavailable"},
       {quadruples,  "quadruples",  NULL,           "unavailable"},
       {celltrips,   "celltrips",   NULL,           "unavailable"},
       {cellquads,   "cellquads",   NULL,           "unavailable"},
       {cellquins,   "cellquins",   NULL,           "unavailable"},
       {distances,   "distances",   distances_sg,   "distances_sg"},
       {indsets,     "indsets",     NULL,           "unavailable"},
       {cliques,     "cliques",     NULL,           "unavailable"},
       {cellcliq,    "cellcliq",    NULL,           "unavailable"},
       {cellind,     "cellind",     NULL,           "unavailable"},
       {adjacencies, "adjacencies", adjacencies_sg, "adjacencies_sg"},
       {cellfano,    "cellfano",    NULL,           "unavailable"},
       {cellfano2,   "cellfano2",   NULL,           "unavailable"},
       {refinvar,    "refinvar",    NULL,           "unavailable"}
      };
#define NUMINVARS ((int)(sizeof(invarproc)/sizeof(struct invarrec)))

static void help(FILE*, int);
static void userautom(int,int*,int*,int,int,int);
static void usernode(graph*,int*,int*,int,int,int,int,int,int);
static void userlevel(int*,int*,int,int*,statsblk*,int,int,int,int,int,int);
static int usercanon(graph*,int*,graph*,int,int,int,int);

static boolean options_writeautoms,options_writemarkers,
            options_digraph,options_getcanon,options_linelength;
static int options_invarproc,options_mininvarlevel,options_maxinvarlevel,
	    options_invararg,options_tc_level,options_cartesian;
static int options_schreier,options_keepgroup,options_verbosity,
	   options_strategy;

#if USE_ANSICONTROLS && !DREADTEST
#define PUTORBITS putorbitsplus
#else
#define PUTORBITS putorbits
#endif

#ifdef  EXTRADECLS
EXTRADECLS
#endif

#if !HAVE_SIGACTION
#undef ALLOW_INTERRUPT
#define ALLOW_INTERRUPT 0
#endif

#if ALLOW_INTERRUPT
/*****************************************************************************
*                                                                            *
*  Routines for catching SIGINT                                              *
*                                                                            *
*****************************************************************************/

void
sigintcatcher(int sig)
/* This is the routine called on SIGINT receipt. */
{
    struct sigaction ss;

    nauty_kill_request = 1;
    ss.sa_handler = SIG_DFL;
    sigemptyset(&ss.sa_mask);
    ss.sa_flags = 0;
    sigaction(SIGINT,&ss,0);
}

#endif
int main(){
int m,n,newm,newn;
boolean gvalid,ovalid,cvalid,pvalid,minus,prompt,doquot;
boolean gvalid_sg,cvalid_sg;
int i,j,k,worksize,numcells,savednc,refcode,umask,qinvar;
int oldorg,oldmode;
int maxsize,cell1,cell2;
boolean ranreg,same;
char *s1,*s2;
int c,d;
unsigned long uli;
size_t sli;
set *gp;
double timebefore,timeafter,mintime;
char filename[515];
int sgn,sgorg,nperm;
int multiplicity,actmult;
long zseed;
permnode *generators;
char *ap,*parameters;
boolean flushing;

    options_writeautoms = options_writemarkers = TRUE;
    options_digraph = FALSE;
    options_getcanon = options.getcanon;
    options_mininvarlevel = options.mininvarlevel;
    options_maxinvarlevel = options.maxinvarlevel;
    options_invararg = options.invararg;
    options_invarproc = 1; /* index into invarproc[] */
    options_tc_level = options.tc_level;
    options_cartesian = options.cartesian;
    options_linelength = options.linelength;
    options_schreier = SCHREIER_DEFAULT;
    options_keepgroup = FALSE;
    generators = NULL;
    options_verbosity = 1;
    options_strategy = 0;

    n = m = 1;
    worksize = WORKSIZE;

#if !MAXN
    n = WORDSIZE;
    DYNALLOC2(graph,g,g_sz,n,m,"dreadnaut");
    DYNALLOC1(int,lab,lab_sz,n,"dreadnaut");
    DYNALLOC1(int,ptn,ptn_sz,n,"dreadnaut");
    DYNALLOC1(int,orbits,orbits_sz,n,"dreadnaut");
    DYNALLOC1(int,perm,perm_sz,n,"dreadnaut");
    DYNALLOC1(set,active,active_sz,m,"dreadnaut");
    n = 1;
#endif

#ifdef DREADTEST
    seed = 1;
    ran_init(seed);
#else
#ifdef  INITSEED
    INITSEED;
    ran_init(seed);
#endif
#endif

    umask = 0;
    pvalid = FALSE;
    ovalid = FALSE;
    gvalid = gvalid_sg = FALSE;  /* at most one valid */
    cvalid = cvalid_sg = FALSE;  /* at most one valid */
    sgorg = labelorg = oldorg = 0;
    sgn = 0;
    multiplicity = 1;
    mintime = 0.0;
    flushing = FALSE;

#ifdef  INITIALIZE
    INITIALIZE;
#endif

if (prompt)
{
    fprintf(PROMPTFILE,"Dreadnaut version %s.\n",NAUTYVERSION);
    fprintf(PROMPTFILE,"> ");
}

nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);
nautinv_check(WORDSIZE,1,1,NAUTYVERSIONID);
nautil_check(WORDSIZE,1,1,NAUTYVERSIONID);
naututil_check(WORDSIZE,1,1,NAUTYVERSIONID);
nausparse_check(WORDSIZE,1,1,NAUTYVERSIONID);

SG_INIT(g_sg);
SG_INIT(canong_sg);
SG_INIT(savedg_sg);

minus = FALSE;

#if !MAXN
DYNALLOC2(graph,g,g_sz,n,m,"dreadnaut");

#endif
readgraph(INFILE,g,options_digraph,prompt,FALSE,
          options_linelength,m,n);

gvalid = TRUE;
cvalid = FALSE;

ovalid = FALSE;


int number = 4;
int numNodes = 3;
/***********************************************************
case 'n'
***********************************************************/
minus = FALSE;
i = numNodes;/**Some int we read either from user or from argv**/

gvalid = FALSE;
cvalid = FALSE;
gvalid_sg = FALSE;
cvalid_sg = FALSE;
pvalid = FALSE;
ovalid = FALSE;
n = i;
m = SETWORDSNEEDED(n);
freeschreier(NULL,&generators);
#if !MAXN
DYNALLOC1(int,lab,lab_sz,n,"dreadnaut");
DYNALLOC1(int,ptn,ptn_sz,n,"dreadnaut");
DYNALLOC1(int,orbits,orbits_sz,n,"dreadnaut");
DYNALLOC1(int,perm,perm_sz,n,"dreadnaut");
DYNALLOC1(set,active,active_sz,m,"dreadnaut");
#endif

/***********************************************************
case 'g'
***********************************************************/
minus = FALSE;
if (SPARSEREP(mode))
{
    readgraph_sg(INFILE,&g_sg,options_digraph,prompt,
                 options_linelength,n);
    gvalid_sg = TRUE;
    cvalid_sg = FALSE;
}
else
{
#if !MAXN
DYNALLOC2(graph,g,g_sz,n,m,"dreadnaut");
#endif
readgraph(INFILE,g,options_digraph,prompt,FALSE,
          options_linelength,m,n);
gvalid = TRUE;
cvalid = FALSE;
}
ovalid = FALSE;

int *edge = g_sg.e;

int *mat = numToLowMat(number, numNodes);
for(int row = 0; row < numNodes; row++){
    for(int col = 0; col < numNodes; col++){
        if(row == col) continue;
        if(mat[numNodes * row + col] == 1){
            *(edge++) = row;
            *(edge++) = col;
        }
    }
}

/***********************************************************
case 'x'
***********************************************************/

minus = FALSE;
            if (mode == TRACES_MODE)
	    {
                ovalid = FALSE;
                cvalid_sg = FALSE;
                if (!gvalid_sg)
                {
                    fprintf(ERRFILE,"g is not defined\n");
		    FLUSHANDPROMPT;
                    //break;
                }

	        traces_opts.getcanon = options_getcanon;
	        traces_opts.writeautoms = options_writeautoms;
	        traces_opts.cartesian = options_cartesian;
	        traces_opts.linelength = options_linelength;
	        traces_opts.digraph = options_digraph;
	        traces_opts.outfile = outfile;
		traces_opts.verbosity = options_verbosity;
		traces_opts.strategy = options_strategy;
		if (options_keepgroup)
		    traces_opts.generators = &generators;
		else
		    traces_opts.generators = NULL;

#if !MAXN
		DYNALLOC1(int,tempptn,tempptn_sz,n,"dreadnaut");
#endif
		if (!pvalid) unitptn(lab,ptn,&numcells,n);
		memcpy(tempptn,ptn,n*sizeof(int));
		savednc = numcells;
		if (options_invarproc != 1 && options_maxinvarlevel > 0)
                {
                    if (options_maxinvarlevel > 1) fprintf(ERRFILE,
		    "Warning: Traces only uses invariants at the top level\n");
                    if (invarproc[options_invarproc].entrypoint_sg)
		    {
#ifdef  CPUTIME
                        timebefore = CPUTIME;
#endif
			cellstarts(tempptn,0,active,m,n);
			doref((graph*)&g_sg,lab,tempptn,0,&savednc,&qinvar,perm,
			    active,&refcode,
                    	    options_sg.userrefproc ? options_sg.userrefproc :
                    	    refine_sg,
                    	    invarproc[options_invarproc].entrypoint_sg,0,0,
                    	    options_invararg,options_digraph,m,n);
                        fprintf(outfile,"Invariant %s %s; %d cell%s",
			    invarproc[options_invarproc].name_sg,
                            (qinvar == 2 ? "worked" : "failed"),
			    SS(savednc,"","s"));
#ifdef  CPUTIME
                        timeafter = CPUTIME;
			fprintf(outfile,"; cpu time = %.2f seconds\n",
                    	    timeafter-timebefore);
#else
		 	fprintf(outfile,"\n");
#endif
		    }
                }

#ifdef  CPUTIME
                timebefore = CPUTIME;
#endif
		actmult = 0;
                //setsigcatcher();
                for (;;)
                {
                    traces_opts.defaultptn = !pvalid;
                    Traces(&g_sg,lab,tempptn,orbits,&traces_opts,&traces_stats,
		       &canong_sg);
		    if (traces_stats.errstatus) break;
                    traces_opts.writeautoms = FALSE;
		    traces_opts.verbosity = 0;
		    ++actmult;
		    if (multiplicity > 0 && actmult >= multiplicity) break;
#ifdef  CPUTIME
		    if (mintime > 0.0 && (actmult < 20 || !(actmult&7))
		       	   && CPUTIME >= timebefore+mintime)
			break;
#endif
                }
#ifdef  CPUTIME
                timeafter = CPUTIME;
#endif
            if (traces_stats.errstatus)
            {
                if (traces_stats.errstatus == NAUABORTED)
                    fprintf(ERRFILE,"Traces aborted\n");
                else if (traces_stats.errstatus == NAUKILLED)
                    fprintf(ERRFILE,"Traces interrupted\n");
                else
                    fprintf(ERRFILE,
                      "Traces returned error status %d [this can't happen]\n",
                      traces_stats.errstatus);
                cvalid = cvalid_sg = ovalid = FALSE;
            }
            else
            {
                fprintf(outfile,"%d orbit%s",
				SS(traces_stats.numorbits,"","s"));
		fprintf(outfile,"; grpsize=");
		writegroupsize(outfile,
			       traces_stats.grpsize1,traces_stats.grpsize2);
                fprintf(outfile,"; %d gen%s",
                     SS(traces_stats.numgenerators,"","s"));
                fprintf(outfile,
		     "; %lu node%s ", SS(traces_stats.numnodes,"","s"));
		if (traces_stats.interrupted)
		    fprintf(outfile,
			    "(%lu interrupted, ",traces_stats.interrupted);
		else
		    fprintf(outfile,"(");
                fprintf(outfile,"%lu peak); maxlev=%d\n",
		    traces_stats.peaknodes,traces_stats.treedepth);
                if (options_getcanon)
                    fprintf(outfile,
		            "canupdates=%d; ",traces_stats.canupdates);
#ifdef  CPUTIME
                fprintf(outfile,actmult == 1 ?
                              "cpu time = %.2f seconds\n" :
                              "cpu time = %.7f seconds\n",
			      (timeafter-timebefore)/actmult);
#else
		fprintf(outfile,"\n");
#endif
		if (options_getcanon) cvalid_sg = TRUE;
		ovalid = TRUE;
            }
	    }
            else
            {
                ovalid = FALSE;
                cvalid = cvalid_sg = FALSE;
                if (!gvalid && !gvalid_sg)
                {
                    fprintf(ERRFILE,"g is not defined\n");
		    FLUSHANDPROMPT;
                    //break;
                }
		if (mode == DENSE_MODE)
		{
                    if (pvalid)
                    {
                        fprintf(outfile,"[fixing partition]\n");
                        options.defaultptn = FALSE;
                    }
                    else
                        options.defaultptn = TRUE;

                    options.outfile = outfile;
                    options.digraph = options_digraph;
                    options.cartesian = options_cartesian;
                    options.schreier = (options_schreier > 0);
                    options.getcanon = options_getcanon;
                    options.tc_level = options_tc_level;
                    options.linelength = options_linelength;
		    options.invarproc
				= invarproc[options_invarproc].entrypoint;
		    options.mininvarlevel = options_mininvarlevel;
		    if (options.invarproc)
		       options.maxinvarlevel = options_maxinvarlevel;
		    else
		       options.maxinvarlevel = 0;
		    options.invararg = options_invararg;
		    if (options_schreier > 0)
			schreier_fails(options_schreier);
/*
                    if (umask & U_NODE)  options.usernodeproc = NODEPROC;
                    else                 options.usernodeproc = NULL;
                    if (umask & U_AUTOM) options.userautomproc = AUTOMPROC;
                    else                 options.userautomproc = NULL;
                    if (umask & U_LEVEL) options.userlevelproc = LEVELPROC;
                    else                 options.userlevelproc = NULL;
                    if (umask & U_REF)   options.userrefproc = REFPROC;
                    else                 options.userrefproc = NULL;
                    if (umask & U_CANON) options.usercanonproc = CANONPROC;
                    else                 options.usercanonproc = NULL;
*/
#if !MAXN
                    if (options_getcanon)
                        DYNALLOC2(graph,canong,canong_sz,n,m,"dreadnaut");
                    DYNALLOC1(setword,workspace,workspace_sz,2*m*worksize,
								"dreadnaut");
#endif
                    firstpath = TRUE;
                    options.writeautoms = options_writeautoms;
                    options.writemarkers = options_writemarkers;
#ifdef  CPUTIME
                    timebefore = CPUTIME;
#endif
		    actmult = 0;
		    //setsigcatcher();
                    for (;;)
                    {
                        nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,
                             2*m*worksize,m,n,canong);
			if (stats.errstatus) break;
                        options.writeautoms = FALSE;
                        options.writemarkers = FALSE;
			++actmult;
			if (multiplicity > 0 && actmult >= multiplicity)
			    break;
#ifdef  CPUTIME
			if (mintime > 0.0 && (actmult < 20 || !(actmult&7))
	  		        && CPUTIME >= timebefore+mintime)
			    break;
#endif
                    }

#ifdef  CPUTIME
                    timeafter = CPUTIME;
#endif
		}
		if (mode == SPARSE_MODE)
		{
                    if (pvalid)
                    {
                        fprintf(outfile,"[fixing partition]\n");
                        options_sg.defaultptn = FALSE;
                    }
                    else
                        options_sg.defaultptn = TRUE;

                    options_sg.outfile = outfile;
                    options_sg.digraph = options_digraph;
                    options_sg.cartesian = options_cartesian;
                    options_sg.schreier = (options_schreier > 0);
                    options_sg.getcanon = options_getcanon;
                    options_sg.linelength = options_linelength;
		    options_sg.invarproc
                         = invarproc[options_invarproc].entrypoint_sg;
		    options_sg.mininvarlevel = options_mininvarlevel;
		    if (options_sg.invarproc)
		       options_sg.maxinvarlevel = options_maxinvarlevel;
		    else
		       options_sg.maxinvarlevel = 0;
		    options_sg.invararg = options_invararg;
                    options_sg.tc_level = options_tc_level;
		    if (options_schreier > 0)
			schreier_fails(options_schreier);

/*
                    if (umask & U_NODE)  options_sg.usernodeproc = NODEPROC;
                    else                 options_sg.usernodeproc = NULL;
                    if (umask & U_AUTOM) options_sg.userautomproc = AUTOMPROC;
                    else                 options_sg.userautomproc = NULL;
                    if (umask & U_LEVEL) options_sg.userlevelproc = LEVELPROC;
                    else                 options_sg.userlevelproc = NULL;
                    if (umask & U_REF)   options_sg.userrefproc = REFPROC;
                    else                 options_sg.userrefproc = NULL;
                    if (umask & U_CANON) options_sg.usercanonproc = CANONPROC;
                    else                 options_sg.usercanonproc = NULL;
*/
#if !MAXN
                    DYNALLOC1(setword,workspace,workspace_sz,2*m*worksize,
								"dreadnaut");
#endif

                    firstpath = TRUE;
                    options_sg.writeautoms = options_writeautoms;
                    options_sg.writemarkers = options_writemarkers;
#ifdef  CPUTIME
                    timebefore = CPUTIME;
#endif
		    actmult = 0;
		    //setsigcatcher();
                    for (;;)
                    {
                        nauty((graph*)&g_sg,lab,ptn,NULL,orbits,&options_sg,
                         &stats,workspace,2*m*worksize,m,n,(graph*)&canong_sg);
			if (stats.errstatus) break;
                        options_sg.writeautoms = FALSE;
                        options_sg.writemarkers = FALSE;
			++actmult;
			if (multiplicity > 0 && actmult >= multiplicity)
			    break;
#ifdef  CPUTIME
			if (mintime > 0.0 && (actmult < 20 || !(actmult&7))
	  		        && CPUTIME >= timebefore+mintime)
			    break;
#endif
                    }

#ifdef  CPUTIME
                    timeafter = CPUTIME;
#endif
		}

		if (stats.errstatus)
		{
		    if (stats.errstatus == NAUABORTED)
		        fprintf(ERRFILE,"nauty aborted\n");
		    else if (stats.errstatus == NAUKILLED)
		        fprintf(ERRFILE,"nauty interrupted\n");
                    else
                        fprintf(ERRFILE,
                        "nauty returned error status %d [this can't happen]\n",
                           stats.errstatus);
		    cvalid = cvalid_sg = ovalid = FALSE;
		}
                else
                {
                    if (options_getcanon)
		    {
			if (mode == DENSE_MODE) cvalid = TRUE;
			else                    cvalid_sg = TRUE;
		    }
                    ovalid = TRUE;
                    fprintf(outfile,"%d orbit%s",SS(stats.numorbits,"","s"));
		    fprintf(outfile,"; grpsize=");
		    writegroupsize(outfile,stats.grpsize1,stats.grpsize2);
                    fprintf(outfile,"; %d gen%s",
                            SS(stats.numgenerators,"","s"));
                    fprintf(outfile,"; %lu node%s",SS(stats.numnodes,"","s"));
                    if (stats.numbadleaves)
                        fprintf(outfile," (%lu bad lea%s)",
                            SS(stats.numbadleaves,"f","ves"));
                    fprintf(outfile,"; maxlev=%d\n", stats.maxlevel);
                    /* fprintf(outfile,"tctotal=%lu",stats.tctotal); */
                    if (options_getcanon)
                        fprintf(outfile,"canupdates=%lu; ",stats.canupdates);
#ifdef  CPUTIME
                    fprintf(outfile,actmult == 1 ?
                              "cpu time = %.2f seconds\n" :
                              "cpu time = %.7f seconds\n",
                            (timeafter-timebefore)/actmult);
#else
                    fprintf(outfile,"\n");
#endif
		    if (mode == DENSE_MODE && options_maxinvarlevel != 0
		       && invarproc[options_invarproc].entrypoint)
                    {
                        fprintf(outfile,"invarproc \"%s\" succeeded %lu/%lu",
                            invarproc[options_invarproc].name,
			    stats.invsuccesses,stats.invapplics);
                        if (stats.invarsuclevel > 0)
                            fprintf(outfile," beginning at level %d.\n",
                                    stats.invarsuclevel);
                        else
                            fprintf(outfile,".\n");
                    }
		    if (mode == SPARSE_MODE && options_maxinvarlevel != 0
		       && invarproc[options_invarproc].entrypoint_sg)
                    {
                        fprintf(outfile,"invarproc \"%s\" succeeded %lu/%lu",
                            invarproc[options_invarproc].name_sg,
			    stats.invsuccesses,stats.invapplics);
                        if (stats.invarsuclevel > 0)
                            fprintf(outfile," beginning at level %d.\n",
                                    stats.invarsuclevel);
                        else
                            fprintf(outfile,".\n");
                    }
                }
	    }

/***********************************************************
case '@'
***********************************************************/

            minus = FALSE;
            if (cvalid)
            {
#if !MAXN
                DYNALLOC2(graph,savedg,savedg_sz,n,m,"dreadnaut");
                DYNALLOC1(int,savedlab,savedlab_sz,n,"dreadnaut");
                DYNALLOC1(int,savedptn,savedptn_sz,n,"dreadnaut");
#endif
                sgn = n;
		memcpy(savedg,canong,m*(size_t)n*sizeof(setword));
                for (i = n; --i >= 0;)
		{
		    savedlab[i] = lab[i];
		    savedptn[i] = ptn[i];
		}
                sgorg = labelorg;
            }
            else if (cvalid_sg)
	    {
#if !MAXN
                DYNALLOC1(int,savedlab,savedlab_sz,n,"dreadnaut");
                DYNALLOC1(int,savedptn,savedptn_sz,n,"dreadnaut");
#endif
		sgn = n;
		copy_sg(&canong_sg,&savedg_sg);
                for (i = n; --i >= 0;)
		{
		    savedlab[i] = lab[i];
		    savedptn[i] = ptn[i];
		}
                sgorg = labelorg;
	    }
            else
	    {
                fprintf(ERRFILE,"h is not defined\n");
		FLUSHANDPROMPT;
	    }


/***********************************************************
case 'n'
***********************************************************/
number = 1;
numNodes = 3;

minus = FALSE;
i = numNodes;/**Some int we read either from user or from argv**/

gvalid = FALSE;
cvalid = FALSE;
gvalid_sg = FALSE;
cvalid_sg = FALSE;
pvalid = FALSE;
ovalid = FALSE;
n = i;
m = SETWORDSNEEDED(n);
freeschreier(NULL,&generators);
#if !MAXN
DYNALLOC1(int,lab,lab_sz,n,"dreadnaut");
DYNALLOC1(int,ptn,ptn_sz,n,"dreadnaut");
DYNALLOC1(int,orbits,orbits_sz,n,"dreadnaut");
DYNALLOC1(int,perm,perm_sz,n,"dreadnaut");
DYNALLOC1(set,active,active_sz,m,"dreadnaut");
#endif

/***********************************************************
case 'g'
***********************************************************/
minus = FALSE;
if (SPARSEREP(mode))
{
    readgraph_sg(INFILE,&g_sg,options_digraph,prompt,
                 options_linelength,n);
    gvalid_sg = TRUE;
    cvalid_sg = FALSE;
}
else
{
#if !MAXN
DYNALLOC2(graph,g,g_sz,n,m,"dreadnaut");
#endif
readgraph(INFILE,g,options_digraph,prompt,FALSE,
          options_linelength,m,n);
gvalid = TRUE;
cvalid = FALSE;
}
ovalid = FALSE;

edge = g_sg.e;

mat = numToLowMat(number, numNodes);
for(int row = 0; row < numNodes; row++){
    for(int col = 0; col < numNodes; col++){
        if(row == col) continue;
        if(mat[numNodes * row + col] == 1){
            *(edge++) = row;
            *(edge++) = col;
        }
    }
}

/***********************************************************
case 'x'
***********************************************************/

minus = FALSE;
            if (mode == TRACES_MODE)
	    {
                ovalid = FALSE;
                cvalid_sg = FALSE;
                if (!gvalid_sg)
                {
                    fprintf(ERRFILE,"g is not defined\n");
		    FLUSHANDPROMPT;
                    //break;
                }

	        traces_opts.getcanon = options_getcanon;
	        traces_opts.writeautoms = options_writeautoms;
	        traces_opts.cartesian = options_cartesian;
	        traces_opts.linelength = options_linelength;
	        traces_opts.digraph = options_digraph;
	        traces_opts.outfile = outfile;
		traces_opts.verbosity = options_verbosity;
		traces_opts.strategy = options_strategy;
		if (options_keepgroup)
		    traces_opts.generators = &generators;
		else
		    traces_opts.generators = NULL;

#if !MAXN
		DYNALLOC1(int,tempptn,tempptn_sz,n,"dreadnaut");
#endif
		if (!pvalid) unitptn(lab,ptn,&numcells,n);
		memcpy(tempptn,ptn,n*sizeof(int));
		savednc = numcells;
		if (options_invarproc != 1 && options_maxinvarlevel > 0)
                {
                    if (options_maxinvarlevel > 1) fprintf(ERRFILE,
		    "Warning: Traces only uses invariants at the top level\n");
                    if (invarproc[options_invarproc].entrypoint_sg)
		    {
#ifdef  CPUTIME
                        timebefore = CPUTIME;
#endif
			cellstarts(tempptn,0,active,m,n);
			doref((graph*)&g_sg,lab,tempptn,0,&savednc,&qinvar,perm,
			    active,&refcode,
                    	    options_sg.userrefproc ? options_sg.userrefproc :
                    	    refine_sg,
                    	    invarproc[options_invarproc].entrypoint_sg,0,0,
                    	    options_invararg,options_digraph,m,n);
                        fprintf(outfile,"Invariant %s %s; %d cell%s",
			    invarproc[options_invarproc].name_sg,
                            (qinvar == 2 ? "worked" : "failed"),
			    SS(savednc,"","s"));
#ifdef  CPUTIME
                        timeafter = CPUTIME;
			fprintf(outfile,"; cpu time = %.2f seconds\n",
                    	    timeafter-timebefore);
#else
		 	fprintf(outfile,"\n");
#endif
		    }
                }

#ifdef  CPUTIME
                timebefore = CPUTIME;
#endif
		actmult = 0;
                //setsigcatcher();
                for (;;)
                {
                    traces_opts.defaultptn = !pvalid;
                    Traces(&g_sg,lab,tempptn,orbits,&traces_opts,&traces_stats,
		       &canong_sg);
		    if (traces_stats.errstatus) break;
                    traces_opts.writeautoms = FALSE;
		    traces_opts.verbosity = 0;
		    ++actmult;
		    if (multiplicity > 0 && actmult >= multiplicity) break;
#ifdef  CPUTIME
		    if (mintime > 0.0 && (actmult < 20 || !(actmult&7))
		       	   && CPUTIME >= timebefore+mintime)
			break;
#endif
                }

#ifdef  CPUTIME
                timeafter = CPUTIME;
#endif
            if (traces_stats.errstatus)
            {
                if (traces_stats.errstatus == NAUABORTED)
                    fprintf(ERRFILE,"Traces aborted\n");
                else if (traces_stats.errstatus == NAUKILLED)
                    fprintf(ERRFILE,"Traces interrupted\n");
                else
                    fprintf(ERRFILE,
                      "Traces returned error status %d [this can't happen]\n",
                      traces_stats.errstatus);
                cvalid = cvalid_sg = ovalid = FALSE;
            }
            else
            {
                fprintf(outfile,"%d orbit%s",
				SS(traces_stats.numorbits,"","s"));
		fprintf(outfile,"; grpsize=");
		writegroupsize(outfile,
			       traces_stats.grpsize1,traces_stats.grpsize2);
                fprintf(outfile,"; %d gen%s",
                     SS(traces_stats.numgenerators,"","s"));
                fprintf(outfile,
		     "; %lu node%s ", SS(traces_stats.numnodes,"","s"));
		if (traces_stats.interrupted)
		    fprintf(outfile,
			    "(%lu interrupted, ",traces_stats.interrupted);
		else
		    fprintf(outfile,"(");
                fprintf(outfile,"%lu peak); maxlev=%d\n",
		    traces_stats.peaknodes,traces_stats.treedepth);
                if (options_getcanon)
                    fprintf(outfile,
		            "canupdates=%d; ",traces_stats.canupdates);
#ifdef  CPUTIME
                fprintf(outfile,actmult == 1 ?
                              "cpu time = %.2f seconds\n" :
                              "cpu time = %.7f seconds\n",
			      (timeafter-timebefore)/actmult);
#else
		fprintf(outfile,"\n");
#endif
		if (options_getcanon) cvalid_sg = TRUE;
		ovalid = TRUE;
            }
	    }
            else
            {
                ovalid = FALSE;
                cvalid = cvalid_sg = FALSE;
                if (!gvalid && !gvalid_sg)
                {
                    fprintf(ERRFILE,"g is not defined\n");
		    FLUSHANDPROMPT;
                    //break;
                }
		if (mode == DENSE_MODE)
		{
                    if (pvalid)
                    {
                        fprintf(outfile,"[fixing partition]\n");
                        options.defaultptn = FALSE;
                    }
                    else
                        options.defaultptn = TRUE;

                    options.outfile = outfile;
                    options.digraph = options_digraph;
                    options.cartesian = options_cartesian;
                    options.schreier = (options_schreier > 0);
                    options.getcanon = options_getcanon;
                    options.tc_level = options_tc_level;
                    options.linelength = options_linelength;
		    options.invarproc
				= invarproc[options_invarproc].entrypoint;
		    options.mininvarlevel = options_mininvarlevel;
		    if (options.invarproc)
		       options.maxinvarlevel = options_maxinvarlevel;
		    else
		       options.maxinvarlevel = 0;
		    options.invararg = options_invararg;
		    if (options_schreier > 0)
			schreier_fails(options_schreier);
/*
                    if (umask & U_NODE)  options.usernodeproc = NODEPROC;
                    else                 options.usernodeproc = NULL;
                    if (umask & U_AUTOM) options.userautomproc = AUTOMPROC;
                    else                 options.userautomproc = NULL;
                    if (umask & U_LEVEL) options.userlevelproc = LEVELPROC;
                    else                 options.userlevelproc = NULL;
                    if (umask & U_REF)   options.userrefproc = REFPROC;
                    else                 options.userrefproc = NULL;
                    if (umask & U_CANON) options.usercanonproc = CANONPROC;
                    else                 options.usercanonproc = NULL;
*/
#if !MAXN
                    if (options_getcanon)
                        DYNALLOC2(graph,canong,canong_sz,n,m,"dreadnaut");
                    DYNALLOC1(setword,workspace,workspace_sz,2*m*worksize,
								"dreadnaut");
#endif
                    firstpath = TRUE;
                    options.writeautoms = options_writeautoms;
                    options.writemarkers = options_writemarkers;
#ifdef  CPUTIME
                    timebefore = CPUTIME;
#endif
		    actmult = 0;
		    //setsigcatcher();
                    for (;;)
                    {
                        nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,
                             2*m*worksize,m,n,canong);
			if (stats.errstatus) break;
                        options.writeautoms = FALSE;
                        options.writemarkers = FALSE;
			++actmult;
			if (multiplicity > 0 && actmult >= multiplicity)
			    break;
#ifdef  CPUTIME
			if (mintime > 0.0 && (actmult < 20 || !(actmult&7))
	  		        && CPUTIME >= timebefore+mintime)
			    break;
#endif
                    }

#ifdef  CPUTIME
                    timeafter = CPUTIME;
#endif
		}
		if (mode == SPARSE_MODE)
		{
                    if (pvalid)
                    {
                        fprintf(outfile,"[fixing partition]\n");
                        options_sg.defaultptn = FALSE;
                    }
                    else
                        options_sg.defaultptn = TRUE;

                    options_sg.outfile = outfile;
                    options_sg.digraph = options_digraph;
                    options_sg.cartesian = options_cartesian;
                    options_sg.schreier = (options_schreier > 0);
                    options_sg.getcanon = options_getcanon;
                    options_sg.linelength = options_linelength;
		    options_sg.invarproc
                         = invarproc[options_invarproc].entrypoint_sg;
		    options_sg.mininvarlevel = options_mininvarlevel;
		    if (options_sg.invarproc)
		       options_sg.maxinvarlevel = options_maxinvarlevel;
		    else
		       options_sg.maxinvarlevel = 0;
		    options_sg.invararg = options_invararg;
                    options_sg.tc_level = options_tc_level;
		    if (options_schreier > 0)
			schreier_fails(options_schreier);

/*
                    if (umask & U_NODE)  options_sg.usernodeproc = NODEPROC;
                    else                 options_sg.usernodeproc = NULL;
                    if (umask & U_AUTOM) options_sg.userautomproc = AUTOMPROC;
                    else                 options_sg.userautomproc = NULL;
                    if (umask & U_LEVEL) options_sg.userlevelproc = LEVELPROC;
                    else                 options_sg.userlevelproc = NULL;
                    if (umask & U_REF)   options_sg.userrefproc = REFPROC;
                    else                 options_sg.userrefproc = NULL;
                    if (umask & U_CANON) options_sg.usercanonproc = CANONPROC;
                    else                 options_sg.usercanonproc = NULL;
8/
#if !MAXN
                    DYNALLOC1(setword,workspace,workspace_sz,2*m*worksize,
								"dreadnaut");
#endif

                    firstpath = TRUE;
                    options_sg.writeautoms = options_writeautoms;
                    options_sg.writemarkers = options_writemarkers;
#ifdef  CPUTIME
                    timebefore = CPUTIME;
#endif
		    actmult = 0;
		    setsigcatcher();
                    for (;;)
                    {
                        nauty((graph*)&g_sg,lab,ptn,NULL,orbits,&options_sg,
                         &stats,workspace,2*m*worksize,m,n,(graph*)&canong_sg);
			if (stats.errstatus) break;
                        options_sg.writeautoms = FALSE;
                        options_sg.writemarkers = FALSE;
			++actmult;
			if (multiplicity > 0 && actmult >= multiplicity)
			    break;
#ifdef  CPUTIME
			if (mintime > 0.0 && (actmult < 20 || !(actmult&7))
	  		        && CPUTIME >= timebefore+mintime)
			    break;
#endif
                    }
		    unsetsigcatcher();
#ifdef  CPUTIME
                    timeafter = CPUTIME;
#endif
		}

		if (stats.errstatus)
		{
		    if (stats.errstatus == NAUABORTED)
		        fprintf(ERRFILE,"nauty aborted\n");
		    else if (stats.errstatus == NAUKILLED)
		        fprintf(ERRFILE,"nauty interrupted\n");
                    else
                        fprintf(ERRFILE,
                        "nauty returned error status %d [this can't happen]\n",
                           stats.errstatus);
		    cvalid = cvalid_sg = ovalid = FALSE;
		}
                else
                {
                    if (options_getcanon)
		    {
			if (mode == DENSE_MODE) cvalid = TRUE;
			else                    cvalid_sg = TRUE;
		    }
                    ovalid = TRUE;
                    fprintf(outfile,"%d orbit%s",SS(stats.numorbits,"","s"));
		    fprintf(outfile,"; grpsize=");
		    writegroupsize(outfile,stats.grpsize1,stats.grpsize2);
                    fprintf(outfile,"; %d gen%s",
                            SS(stats.numgenerators,"","s"));
                    fprintf(outfile,"; %lu node%s",SS(stats.numnodes,"","s"));
                    if (stats.numbadleaves)
                        fprintf(outfile," (%lu bad lea%s)",
                            SS(stats.numbadleaves,"f","ves"));
                    fprintf(outfile,"; maxlev=%d\n", stats.maxlevel);
                    /* fprintf(outfile,"tctotal=%lu",stats.tctotal); */
                    if (options_getcanon)
                        fprintf(outfile,"canupdates=%lu; ",stats.canupdates);
#ifdef  CPUTIME
                    fprintf(outfile,actmult == 1 ?
                              "cpu time = %.2f seconds\n" :
                              "cpu time = %.7f seconds\n",
                            (timeafter-timebefore)/actmult);
#else
                    fprintf(outfile,"\n");
#endif
		    if (mode == DENSE_MODE && options_maxinvarlevel != 0
		       && invarproc[options_invarproc].entrypoint)
                    {
                        fprintf(outfile,"invarproc \"%s\" succeeded %lu/%lu",
                            invarproc[options_invarproc].name,
			    stats.invsuccesses,stats.invapplics);
                        if (stats.invarsuclevel > 0)
                            fprintf(outfile," beginning at level %d.\n",
                                    stats.invarsuclevel);
                        else
                            fprintf(outfile,".\n");
                    }
		    if (mode == SPARSE_MODE && options_maxinvarlevel != 0
		       && invarproc[options_invarproc].entrypoint_sg)
                    {
                        fprintf(outfile,"invarproc \"%s\" succeeded %lu/%lu",
                            invarproc[options_invarproc].name_sg,
			    stats.invsuccesses,stats.invapplics);
                        if (stats.invarsuclevel > 0)
                            fprintf(outfile," beginning at level %d.\n",
                                    stats.invarsuclevel);
                        else
                            fprintf(outfile,".\n");
                    }
                }
	    }

/***********************************************************
case '##'
***********************************************************/

if (cvalid || cvalid_sg)
            {
                if (sgn > 0)
                {
                    if (sgn != n)
                        fprintf(outfile,
                              "h and h' have different sizes.\n");
                    else
                    {
			if (cvalid)
		 	{
                            for (sli = 0; sli < m*(size_t)n; ++sli)
                                if (savedg[sli] != canong[sli]) break;
			    same = (sli == m*(size_t)n);
			}
			else
			    same = aresame_sg(&canong_sg,&savedg_sg);

                        if (!same)
                            fprintf(outfile,"h and h' are different.\n");
                        else
                        {
			    for (i = 0; i < n; ++i)
				if ((ptn[i] == 0) != (savedptn[i] == 0))
				    break;
			    if (i < n)
                                fprintf(outfile,
                 "h and h' are identical but have incompatible colourings.\n");
			    else
                                fprintf(outfile,
                                    "h and h' are identical.\n");

                                putmapping(outfile,savedlab,sgorg,
                                       lab,labelorg,options_linelength,n);
                        }
                    }
                }
                else
		{
                    fprintf(ERRFILE,"h' is not defined\n");
		    FLUSHANDPROMPT;
		}
            }
            else
	    {
                fprintf(ERRFILE,"h is not defined\n");
		FLUSHANDPROMPT;
	    }
    return 0;
}
