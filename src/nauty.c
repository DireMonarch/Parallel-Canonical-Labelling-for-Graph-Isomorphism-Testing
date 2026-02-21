/*****************************************************************************
 *                                                                            *
 *  Main source file for version 2.7 of nauty.                                *
 *                                                                            *
 *   Copyright (1984-2018) Brendan McKay.  All rights reserved.  Permission   *
 *   Subject to the waivers and disclaimers in nauty.h.                       *
 *                                                                            *
 *   CHANGE HISTORY                                                           *
 *       10-Nov-87 : final changes for version 1.2                            *
 *        5-Dec-87 : renamed to version 1.3 (no changes to this file)         *
 *       28-Sep-88 : renamed to version 1.4 (no changes to this file)         *
 *       23-Mar-89 : changes for version 1.5 :                                *
 *                   - add use of refine1 instead of refine for m==1          *
 *                   - changes for new optionblk syntax                       *
 *                   - disable tc_level use for digraphs                      *
 *                   - interposed doref() interface to refine() so that       *
 *                        options.invarproc can be supported                  *
 *                   - declared local routines static                         *
 *       28-Mar-89 : - implemented mininvarlevel/maxinvarlevel < 0 options    *
 *        2-Apr-89 : - added invarproc fields in stats                        *
 *        5-Apr-89 : - modified error returns from nauty()                    *
 *                   - added error message to ERRFILE                         *
 *                   - changed MAKEEMPTY uses to EMPTYSET                     *
 *       18-Apr-89 : - added MTOOBIG and CANONGNIL                            *
 *        8-May-89 : - changed firstcode[] and canoncode[] to short           *
 *       10-Nov-90 : changes for version 1.6 :                                *
 *                   - added dummy routine nauty_null (see dreadnaut.c)       *
 *        2-Sep-91 : changes for version 1.7 :                                *
 *                   - moved MULTIPLY into nauty.h                            *
 *       27-Mar-92 : - changed 'n' into 'm' in error message in nauty()       *
 *        5-Jun-93 : renamed to version 1.7+ (no changes to this file)        *
 *       18-Aug-93 : renamed to version 1.8 (no changes to this file)         *
 *       17-Sep-93 : renamed to version 1.9 (no changes to this file)         *
 *       13-Jul-96 : changes for version 2.0 :                                *
 *                   - added dynamic allocation                               *
 *       21-Oct-98 : - made short into shortish for BIGNAUTY as needed        *
 *        7-Jan-00 : - allowed n=0                                            *
 *                   - added nauty_check() and a call to it                   *
 *       12-Feb-00 : - used better method for target cell memory allocation   *
 *                   - did a little formating of the code                     *
 *       27-May-00 : - fixed error introduced on Feb 12.                      *
 *                   - dynamic allocations in nauty() are now deallocated     *
 *                     before return if n >= 320.                             *
 *       16-Nov-00 : - use function prototypes, change UPROC to void.         *
 *                   - added argument to tcellproc(), removed nvector         *
 *                   - now use options.dispatch, options.groupopts is gone.   *
 *       22-Apr-01 : - Added code for compilation into Magma                  *
 *                   - Removed nauty_null() and EXTDEFS                       *
 *        2-Oct-01 : - Improved error message for bad dispatch vector         *
 *       21-Nov-01 : - use NAUTYREQUIRED in nauty_check()                     *
 *       20-Dec-02 : changes for version 2.2:                                 *
 *                   - made tcnode0 global                                    *
 *                   - added nauty_freedyn()                                  *
 *       17-Nov-03 : changed INFINITY to NAUTY_INFINITY                       *
 *       14-Sep-04 : extended prototypes even to recursive functions          *
 *       16-Oct-04 : disallow NULL dispatch vector                            *
 *       11-Nov-05 : changes for version 2.3:                                 *
 *                   - init() and cleanup() optional calls                    *
 *       23-Nov-06 : changes for version 2.4:                                 *
 *                   - use maketargetcell() instead of tcellproc()            *
 *       29-Nov-06 : add extra_autom, extra_level, extra_options              *
 *       10-Dec-06 : remove BIGNAUTY                                          *
 *       10-Nov-09 : remove shortish and permutation types                    *
 *       16-Nov-11 : added Shreier option                                     *
 *       15-Jan-12 : added TLS_ATTR to static declarations                    *
 *       18-Jan-13 : added signal aborting                                    *
 *       19-Jan-13 : added usercanonproc()                                    *
 *       14-Oct-17 : corrected code for n=0                                   *
 *                                                                            *
 *****************************************************************************/
// #define MPI

/* Print debug messages, if TRUE */
#define __DEBUG__ FALSE
#define __DEBUG_VERBOSE__ FALSE
#define __DEBUG_REF__ FALSE
#define __DEBUG_ACTIVE_PRUNE__ FALSE
#define __DEBUG_CANON__ FALSE

#define ONE_WORD_SETS
#include "nauty.h"
#include "schreier.h"
#include "node.h"
#include "nsb_stack.h"
#include "path.h"
#include "code_path.h"
#include "util.h"

#ifdef MPI
#include "mpi.h"
#include "mpi_routines.h"
#endif /* if MPI */

#ifdef NAUTY_IN_MAGMA
#include "cleanup.e"
#endif

#define NAUTY_ABORTED (-11)
#define NAUTY_KILLED (-12)

typedef struct tcnode_struct
{
    struct tcnode_struct *next;
    set *tcellptr;
} tcnode;

/* aproto: header new_nauty_protos.h */

#ifndef NAUTY_IN_MAGMA
#if !MAXN
static int firstpathnode0(Node *, Stack *, tcnode *);
static int othernode0(Node *, Stack *, tcnode *);
static int firstnode0(Node *, Stack *, tcnode *);
#else
static int firstnode(Node *, Stack *);
static int firstpathnode(Node *, Stack *);
static int othernode(Node *, Stack *);
#endif
static void firstterminal(int *, int, Node *);
static int processnode(int *, int *, int, int, Node *);
static void recover(int *, int, Node *);
static void writemarker(int, int, int, int, int, int);
static void new_cannon_found(Node *node, boolean first);
static void active_prune_all_children(Path *, int, Stack *);

#endif

#if MAXM == 1
#define M 1
#else
#define M m
#endif

#define OPTCALL(proc) \
    if (proc != NULL) \
    (*proc)

/* copies of some of the options: */
static TLS_ATTR
    boolean getcanon,
    digraph, writeautoms, domarkers, cartesian, doschreier;
static TLS_ATTR int linelength, tc_level, mininvarlevel, maxinvarlevel, invararg;
static TLS_ATTR void (*usernodeproc)(graph *, int *, int *, int, int, int, int, int, int);
static TLS_ATTR void (*userautomproc)(int, int *, int *, int, int, int);
static TLS_ATTR void (*userlevelproc)(int *, int *, int, int *, statsblk *, int, int, int, int, int, int);
static TLS_ATTR int (*usercanonproc)(graph *, int *, graph *, unsigned long, int, int, int);
static TLS_ATTR void (*invarproc)(graph *, int *, int *, int, int, int, int *, int, boolean, int, int);
static TLS_ATTR FILE *outfile;
static TLS_ATTR dispatchvec dispatch;

/* local versions of some of the arguments: */
static TLS_ATTR int m, n;
static TLS_ATTR graph *g, *canong;
static TLS_ATTR int *orbits;
static TLS_ATTR statsblk *stats;
/* temporary versions of some stats: */
static TLS_ATTR unsigned long invapplics, invsuccesses;
static TLS_ATTR int invarsuclevel;

/* working variables: <the "bsf leaf" is the leaf which is best guess so
                           far at the canonical leaf>  */

static TLS_ATTR int

    /**
     * Seem to be truly global, not per node
     */
    stglb_allsamelevel, /* level of least ancestor of first leaf for which all descendant leaves are known to be equivalent */
    stglb_samerows,     /* number of rows of canong which are correct for the bsf leaf  BDM:correct description? */
    stglb_canonlevel,   /* level of bsf leaf */
    stglb_stabvertex;
/* point fixed in ancestor of first leaf at level gca_canon */ /* Seems to only be used in userautomproc */

static TLS_ATTR boolean needshortprune; /* used to flag calls to shortprune */

#if !MAXN
DYNALLSTAT(set, defltwork, defltwork_sz);                           /* workspace in case none provided */
DYNALLSTAT(int, workperm, workperm_sz);                             /* various scratch uses */
DYNALLSTAT(int, firstlab, firstlab_sz); /* label from first leaf */ /* once set int he first leaf, shouldn't change */
DYNALLSTAT(int, canonlab, canonlab_sz); /* label from bsf leaf */   /* set at first leaf, and at any time a better Canonica Label is found */
// DYNALLSTAT(short,firstcode,firstcode_sz);   /* codes for first leaf */
// DYNALLSTAT(short,canoncode,canoncode_sz);   /* codes for bsf leaf */
DYNALLSTAT(int, firsttc, firsttc_sz); /* index of target cell for left path */

/* In the dynamically allocated case (MAXN=0), each level of recursion
   needs one set (tcell) to represent the target cell.  This is
   implemented by using a linked list of tcnode anchored at the root
   of the search tree.  Each node points to its child (if any) and to
   the dynamically allocated tcell.  Apart from the first node of
   the list, each node always has a tcell good for m up to alloc_m.
   tcnodes and tcells are kept between calls to nauty, except that
   they are freed and reallocated if m gets bigger than alloc_m.  */

static TLS_ATTR tcnode tcnode0 = {NULL, NULL};
static TLS_ATTR int alloc_m = 0;

#else
static TLS_ATTR set defltwork[2 * MAXM]; /* workspace in case none provided */
static TLS_ATTR int workperm[MAXN];      /* various scratch uses */
static TLS_ATTR int firstlab[MAXN],      /* label from first leaf */
    canonlab[MAXN];                      /* label from bsf leaf */

static TLS_ATTR int firsttc[MAXN + 2]; /* index of target cell for left path */
#endif

static TLS_ATTR set *workspace, *worktop; /* first and just-after-last
                     addresses of work area to hold automorphism data */
static TLS_ATTR set *fmptr;               /* pointer into workspace */

static TLS_ATTR schreier *gp; /* These two for Schreier computations */
static TLS_ATTR permnode *gens;

/**
 * Stuff I added that shouldn't be here, but I haven't moved out of global yet
 */
static TLS_ATTR Path *firstpath; /* path of first leaf */ /* I added this, maybe shouldn't be static global */
static TLS_ATTR Path *canonpath; /* path of bsf leaf */   /* I added this, maybe shouldn't be static global */
static TLS_ATTR CodePath *firstcodes;                     /* codes for the first leaf */
static TLS_ATTR CodePath *canoncodes;                     /* codes for the canon (bsf) leaf */

#ifdef MPI
static TLS_ATTR MPIState *mpi_state;
static TLS_ATTR CanonState canonstate;

void mpi_received_active_prune(MPIState *mpi_state, Path *path, int to_level, Stack *stack); /* forward declaring function */
#endif                                                                                       /* ifdef MPI */

// TODO: Delete This
// FILE *tracefile;
// TO HERE

/**
 * End of stuff I added that shouldn't be here
 */

/*****************************************************************************
 *                                                                            *
 *  This procedure finds generators for the automorphism group of a           *
 *  vertex-coloured graph and optionally finds a canonically labelled         *
 *  isomorph.  A description of the data structures can be found in           *
 *  nauty.h and in the "nauty User's Guide".  The Guide also gives            *
 *  many more details about its use, and implementation notes.                *
 *                                                                            *
 *  Parameters - <r> means read-only, <w> means write-only, <wr> means both:  *
 *           g <r>  - the graph                                               *
 *     lab,ptn <rw> - used for the partition nest which defines the colouring *
 *                  of g.  The initial colouring will be set by the program,  *
 *                  using the same colour for every vertex, if                *
 *                  options->defaultptn!=FALSE.  Otherwise, you must set it   *
 *                  yourself (see the Guide). If options->getcanon!=FALSE,    *
 *                  the contents of lab on return give the labelling of g     *
 *                  corresponding to canong.  This does not change the        *
 *                  initial colouring of g as defined by (lab,ptn), since     *
 *                  the labelling is consistent with the colouring.           *
 *     active  <r>  - If this is not NULL and options->defaultptn==FALSE,     *
 *                  it is a set indicating the initial set of active colours. *
 *                  See the Guide for details.                                *
 *     orbits  <w>  - On return, orbits[i] contains the number of the         *
 *                  least-numbered vertex in the same orbit as i, for         *
 *                  i=0,1,...,n-1.                                            *
 *    options  <r>  - A list of options.  See nauty.h and/or the Guide        *
 *                  for details.                                              *
 *      stats  <w>  - A list of statistics produced by the procedure.  See    *
 *                  nauty.h and/or the Guide for details.                     *
 *  workspace  <w>  - A chunk of memory for working storage.                  *
 *  worksize   <r>  - The number of setwords in workspace.  See the Guide     *
 *                  for guidance.                                             *
 *          m  <r>  - The number of setwords in sets.  This must be at        *
 *                  least ceil(n / WORDSIZE) and at most MAXM.                *
 *          n  <r>  - The number of vertices.  This must be at least 1 and    *
 *                  at most MAXN.                                             *
 *     canong  <w>  - The canononically labelled isomorph of g.  This is      *
 *                  only produced if options->getcanon!=FALSE, and can be     *
 *                  given as NULL otherwise.                                  *
 *                                                                            *
 *  FUNCTIONS CALLED: firstpathnode(),updatecan()                             *
 *                                                                            *
 *****************************************************************************/

#ifdef MPI
void nauty(graph *g_arg, int *lab, int *ptn, set *active_arg,
           int *orbits_arg, optionblk *options, statsblk *stats_arg,
           set *ws_arg, int worksize, int m_arg, int n_arg, graph *canong_arg, int max_work_size_to_send, int work_split, int nodes_per_poll,
           int argc, char **argv, int *mpi_details)
#else  /* if MPI */
void nauty(graph *g_arg, int *lab, int *ptn, set *active_arg,
           int *orbits_arg, optionblk *options, statsblk *stats_arg,
           set *ws_arg, int worksize, int m_arg, int n_arg, graph *canong_arg)
#endif /* if MPI */
{
#ifdef MPI

    /** Set up MPI */
    mpi_state = (MPIState *)malloc(sizeof(MPIState));
    MPI_Init(&argc, &argv);                                     /* Initialize MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_state->my_rank);         /* Fetch rank */
    MPI_Comm_size(MPI_COMM_WORLD, &(mpi_state->num_processes)); /* Fetch number of processes */

    if (mpi_state->my_rank == 0)
        printf("MPI Active with %d processes\n\n", mpi_state->num_processes);
    printf("MPI: Process %d: Initialized\n", mpi_state->my_rank);

    /** */
    // double mpi_start_time = MPI_Wtime();  /* Start time for MPI operations */
    // int *processes_running = NULL;  /* only used by process 0 for tracking which processes have finished */
    boolean run_first_path_node = FALSE;
    canonstate.queue_sizes = NULL;
    canonstate.process_waiting_for_work = NULL;
    if (mpi_state->my_rank == 0)
    {
        run_first_path_node = TRUE;
        canonstate.queue_sizes = (int *)malloc(mpi_state->num_processes * sizeof(int));
        canonstate.process_waiting_for_work = (int *)malloc(mpi_state->num_processes * sizeof(int));
        for (int i = 0; i < mpi_state->num_processes; ++i) {
            canonstate.queue_sizes[i] = MPI_QUEUE_SIZE_UNDEFINED;
            canonstate.process_waiting_for_work[i] = FALSE;
        }
    }

    if (max_work_size_to_send < 0)
    {
        max_work_size_to_send = MPI_CONST_MAX_WORK_SIZE_TO_SEND; /* default value, can be overridden by command line argument */
    }
    if (work_split < 1)
    {
        work_split = MPI_CONST_SPLIT_RATIO; /* default value, can be overridden by command line argument */
    }

    if (nodes_per_poll < 1)
    {
        nodes_per_poll = MPI_CONST_NODES_BETWEEN_COMM_POLLS;
    }

    /* update mpi_details for return to main */
    mpi_details[0] = mpi_state->my_rank;
    mpi_details[1] = mpi_state->num_processes;
    mpi_details[2] = max_work_size_to_send;
    mpi_details[3] = work_split;
    mpi_details[4] = nodes_per_poll;

#else /* if MPI */

    printf("MPI Not active, running standalone\n\n");
    boolean run_first_path_node = TRUE;

#endif /* if MPI */

    int i;
    int numcells;
    int retval = 0;
    int initstatus;

#if !MAXN
    tcnode *tcp, *tcq;
#endif /* !MAXN */

    /* determine dispatch vector */

    if (options->dispatch == NULL)
    {
        fprintf(ERRFILE, ">E nauty: null dispatch vector\n");
        fprintf(ERRFILE, "Maybe you need to recompile\n");
        exit(1);
    }
    else
        dispatch = *(options->dispatch);

    if (options->userrefproc)
        dispatch.refine = options->userrefproc;
    else if (dispatch.refine1 && m_arg == 1)
        dispatch.refine = dispatch.refine1;

    if (dispatch.refine == NULL || dispatch.updatecan == NULL || dispatch.targetcell == NULL || dispatch.cheapautom == NULL)
    {
        fprintf(ERRFILE, ">E bad dispatch vector\n");
        exit(1);
    }
    /* take copies of some args, and options: */
    m = m_arg;
    n = n_arg;

    /**
     * Create the stack
     */
    Stack *stack;
    DYNALLOCSTACK(stack, "Create Stack");

    Node *root; // needs to be declared for all ranks, or it won't compile

    if (run_first_path_node)
    {
        /**
         * Allocate root node
         *
         * This should only be done standalone, or on MPI rank 0
         */
        DYNALLOCNODE(root, "Allocate Root");
        root->tcell = NULL;
        root->level = 1;
        root->target_cell = -1;
        root->target_vertex = -1;
        root->numcells = -1;
        root->gca_first = -1;
        root->gca_canon = -1;
        root->eqlev_first = -1;
        root->eqlev_canon = -1;
        root->comp_canon = -1;
        root->noncheaplevel = -1;

        root->lab = lab;
        root->ptn = ptn;

        root->active = (set *)ALLOCS(m, sizeof(set));
        root->fixedpts = (set *)ALLOCS(m, sizeof(set));

        DYNALLOCPATH(root->path, 0, "Root Path");
        DYNALLOCCODEPATH(root->codes, 1, "Root Codes");
        root->codes->data[0] = 0;
        root->cosetindex = 0; // <---------------------------------------------  WTH does this even do?

        /* check for excessive sizes: */

#if !MAXN
        if (m_arg > NAUTY_INFINITY / WORDSIZE + 1)
        {
            stats_arg->errstatus = MTOOBIG;
            fprintf(ERRFILE, "nauty: need m <= %d, but m=%d\n\n",
                    NAUTY_INFINITY / WORDSIZE + 1, m_arg);
            return;
        }
        if (n_arg > NAUTY_INFINITY - 2 || n_arg > WORDSIZE * m_arg)
        {
            stats_arg->errstatus = NTOOBIG;
            fprintf(ERRFILE, "nauty: need n <= min(%d,%d*m), but n=%d\n\n",
                    NAUTY_INFINITY - 2, WORDSIZE, n_arg);
            return;
        }
#else                   /* !MAXN */
        if (m_arg > MAXM)
        {
            stats_arg->errstatus = MTOOBIG;
            fprintf(ERRFILE, "nauty: need m <= %d\n\n", MAXM);
            return;
        }
        if (n_arg > MAXN || n_arg > WORDSIZE * m_arg)
        {
            stats_arg->errstatus = NTOOBIG;
            fprintf(ERRFILE,
                    "nauty: need n <= min(%d,%d*m)\n\n", MAXM, WORDSIZE);
            return;
        }
#endif                  /* !MAXN */
        if (n_arg == 0) /* Special code for zero-sized graph */
        {
            stats_arg->grpsize1 = 1.0;
            stats_arg->grpsize2 = 0;
            stats_arg->numorbits = 0;
            stats_arg->numgenerators = 0;
            stats_arg->errstatus = 0;
            stats_arg->numnodes = 1;
            stats_arg->numbadleaves = 0;
            stats_arg->maxlevel = 1;
            stats_arg->tctotal = 0;
            stats_arg->canupdates = (options->getcanon != 0);
            stats_arg->invapplics = 0;
            stats_arg->invsuccesses = 0;
            stats_arg->invarsuclevel = 0;

            g = canong = NULL;
            initstatus = 0;
            OPTCALL(dispatch.init)(g_arg, &g, canong_arg, &canong,
                                root->lab, root->ptn, root->active, options, &initstatus, m, n);
            if (initstatus)
                stats->errstatus = initstatus;

            if (g == NULL)
                g = g_arg;
            if (canong == NULL)
                canong = canong_arg;
            OPTCALL(dispatch.cleanup)(g_arg, &g, canong_arg, &canong,
                                    root->lab, root->ptn, options, stats_arg, m, n);
            return;
        }
    } // if (run_first_path_node)

    nautil_check(WORDSIZE, m, n, NAUTYVERSIONID);
    OPTCALL(dispatch.check)(WORDSIZE, m, n, NAUTYVERSIONID);

#if !MAXN
    DYNALLOC1(set, defltwork, defltwork_sz, 2 * m, "nauty");
    DYNALLOC1(int, workperm, workperm_sz, n, "nauty");
    DYNALLOC1(int, firstlab, firstlab_sz, n, "nauty");
    DYNALLOC1(int, canonlab, canonlab_sz, n, "nauty");
    // DYNALLOC1(short,firstcode,firstcode_sz,n+2,"nauty");
    // DYNALLOC1(short,canoncode,canoncode_sz,n+2,"nauty");
    DYNALLOC1(int, firsttc, firsttc_sz, n + 2, "nauty");
    firsttc[0] = 0;
    if (m > alloc_m)
    {
        tcp = tcnode0.next;
        while (tcp != NULL)
        {
            tcq = tcp->next;
            FREES(tcp->tcellptr);
            FREES(tcp);
            tcp = tcq;
        }
        alloc_m = m;
        tcnode0.next = NULL;
    }
#endif /* !MAXN */

    orbits = orbits_arg;
    stats = stats_arg;

    getcanon = options->getcanon;
    digraph = options->digraph;
    writeautoms = options->writeautoms;
    domarkers = options->writemarkers;
    cartesian = options->cartesian;
    doschreier = options->schreier;
    if (doschreier)
        schreier_check(WORDSIZE, m, n, NAUTYVERSIONID);
    linelength = options->linelength;
    if (digraph)
        tc_level = 0;
    else
        tc_level = options->tc_level;
    outfile = (options->outfile == NULL ? stdout : options->outfile);
    usernodeproc = options->usernodeproc;
    userautomproc = options->userautomproc;
    userlevelproc = options->userlevelproc;
    usercanonproc = options->usercanonproc;

    invarproc = options->invarproc;
    if (options->mininvarlevel < 0 && options->getcanon)
        mininvarlevel = -options->mininvarlevel;
    else
        mininvarlevel = options->mininvarlevel;
    if (options->maxinvarlevel < 0 && options->getcanon)
        maxinvarlevel = -options->maxinvarlevel;
    else
        maxinvarlevel = options->maxinvarlevel;
    invararg = options->invararg;

    initstatus = 0;

    if (run_first_path_node)
    {

        if (getcanon)
            if (canong_arg == NULL)
            {
                stats_arg->errstatus = CANONGNIL;
                fprintf(ERRFILE,
                        "nauty: canong=NULL but options.getcanon=TRUE\n\n");
                return;
            }

        /* initialize everything: */

        if (options->defaultptn)
        {
            for (i = 0; i < n; ++i) /* give all verts same colour */
            {
                root->lab[i] = i;
                root->ptn[i] = NAUTY_INFINITY;
            }
            root->ptn[n - 1] = 0;
            EMPTYSET(root->active, m);
            ADDELEMENT(root->active, 0);
            numcells = 1;
        }
        else
        {
            root->ptn[n - 1] = 0;
            numcells = 0;
            for (i = 0; i < n; ++i)
                if (root->ptn[i] != 0)
                    root->ptn[i] = NAUTY_INFINITY;
                else
                    ++numcells;
            if (active_arg == NULL)
            {
                EMPTYSET(root->active, m);
                for (i = 0; i < n; ++i)
                {
                    ADDELEMENT(root->active, i);
                    while (root->ptn[i])
                        ++i;
                }
            }
            else
                for (i = 0; i < M; ++i)
                    root->active[i] = active_arg[i];
        }

        g = canong = NULL;
        OPTCALL(dispatch.init)(g_arg, &g, canong_arg, &canong, root->lab, root->ptn, root->active, options, &initstatus, m, n);

        EMPTYSET(root->fixedpts, m);

        root->noncheaplevel = 1;
        root->eqlev_canon = -1; /* needed even if !getcanon */
        root->numcells = numcells;
        root->level = 1;
    } // if (run_first_path_node) {

    if (initstatus)
    {
        stats->errstatus = initstatus;
        return;
    }

    if (g == NULL)
        g = g_arg;
    if (canong == NULL)
        canong = canong_arg;

    if (doschreier)
        newgroup(&gp, &gens, n);

    for (i = 0; i < n; ++i)
        orbits[i] = i;
    stats->grpsize1 = 1.0;
    stats->grpsize2 = 0;
    stats->numgenerators = 0;
    stats->numnodes = 0;
    stats->numbadleaves = 0;
    stats->tctotal = 0;
    stats->canupdates = 0;
    stats->numorbits = n;

    if (worksize >= 2 * m)
        workspace = ws_arg;
    else
    {
        workspace = defltwork;
        worksize = 2 * m;
    }
    worktop = workspace + (worksize - worksize % (2 * m));
    fmptr = workspace;

    printf("workspace size: %lu   %lu\n", worktop - workspace, (worktop - workspace) / (2 * m));

    /* here goes: */
    stats->errstatus = 0;
    needshortprune = FALSE;
    invarsuclevel = NAUTY_INFINITY;
    invapplics = invsuccesses = 0;

    // TODO: Delete this!!
    // char filename[128];
    // sprintf(filename, "trace_%02d.txt", mpi_state->my_rank);
    // tracefile = fopen(filename, "w");
    // TO HERE

    Node *curr = NULL; /* pointer to node currently being acted upon */

    if (run_first_path_node)
    {
#if !MAXN
        retval = firstnode0(root, stack, &tcnode0);
#else  /* if !MAXN */
        retval = firstnode(root, stack);
#endif /* if !MAXN */

        if (__DEBUG__)
        {
            printf("numnodes: %lu\n", stats->numnodes);
            if (__DEBUG_VERBOSE__)
            {
                path_visualize(root->path);
                printf("  fpchild: %d", root->fpchild);
                printf("\n");
                node_visualize(root, n);
                printf("noncheaplevel: %d\n", root->noncheaplevel);
                printf("\n");
            }
        }
        int found_first_leaf = 0;

        /**
         * This loop runs the down the left hand path of the tree and finds the first terminal node (leaf)
         */
        while (!found_first_leaf)
        {
            FREENODE(curr);
            curr = stack_pop(stack);

            if (orbits[curr->target_vertex] == curr->target_vertex) /* ie, not equiv to previous child */
            {
#if !MAXN
                found_first_leaf = firstpathnode0(curr, stack, &tcnode0);
#else  /* if !MAXN */
                found_first_leaf = firstpathnode(curr, stack);
#endif /* if !MAXN */
            }
        }
    } // if (run_first_path_node) {

    /**
     * This is where we would put logic to determine if the problem is small enough to skip running in parallel
     */

#ifdef MPI

    /**
     * Send start or abort message as needed.  If determeined just above here that we want to skip running parallel, send abort.
     */
    // int mpi_return_status = 1;

    int firstleaflevel = -1;
    if (mpi_state->my_rank == 0)
        firstleaflevel = curr->level;
    mpi_start_messages(mpi_state, &firstpath, firstlab, &firstcodes, firsttc, firstleaflevel, &stglb_canonlevel, &stglb_allsamelevel, &canonpath, canonlab, &canoncodes, n);
    // if (mpi_state->my_rank == 0) {
    //     if (stack_size(stack) > 0) {
    //         mpi_send_start_message(mpi_state, firstpath, firstlab, firstcodes, firsttc, curr->level, stglb_canonlevel, stglb_allsamelevel, n);
    //     } else {
    //         mpi_send_start_abort_message(mpi_state);
    //     }
    // } else {
    //     mpi_return_status = mpi_wait_for_start(mpi_state, &firstpath, firstlab, &firstcodes, &canonpath, canonlab, &canoncodes, firsttc, &stglb_canonlevel, &stglb_allsamelevel, &stglb_samerows, n);
    // }
    if (__DEBUG_MPI__)
        printf("MPI: Process %d: Begining main loop\n", mpi_state->my_rank);
#endif /* ifdef MPI */

#ifdef MPI
    /**
     * if running MPI, we need an outer loop to handle asking for more work.
     */
    int last_comm_check = 0;
    mpi_state->state = MPI_STATE_WORKING;

    /** set up state struct for passing to MPI functions */
    canonstate.stack = stack;
    canonstate.canonpath = &canonpath;
    canonstate.canoncodes = &canoncodes;
    canonstate.canonlab = canonlab;
    canonstate.canonlevel = &stglb_canonlevel;
    canonstate.samerows = &stglb_samerows;
    canonstate.worktop = worktop;
    canonstate.workspace = workspace;
    canonstate.fmptr = &fmptr;
    canonstate.mpi_received_active_prune = mpi_received_active_prune;
    canonstate.max_work_size_to_send = max_work_size_to_send;
    canonstate.work_split = work_split;
    canonstate.numnodes = stats->numnodes;
    canonstate.n = n;
    canonstate.m = m;
    // TODO: DELETE THIS
    // canonstate.tracefile = tracefile;
    // TO HERE


    /** This is a special loop for MPI, as if the queue goes empty, we aren't done, we need to ask for more work */
    mpi_state->mpi_last_peroidic_update = 0; /* initialize periodic update timer to 0 */
    
    while (mpi_state->state == MPI_STATE_WORKING || mpi_state->state == MPI_STATE_ASKING_FOR_WORK || mpi_state->state == MPI_STATE_WAIT_FOR_QUEUE_SIZE_RESPONSES || mpi_state->state == MPI_STATE_WAIT_FOR_WORK_RETURN)
    {
#endif /* ifdef MPI */

        /**
         * This is the main loop that processes the stack.  When the stack is empty, the work is done.
         */
        while (stack_size(stack) > 0)
        {
            FREENODE(curr);
            curr = stack_pop(stack);

            curr->comp_canon = codepath_find_eqlev_level(curr->codes, canoncodes, &curr->eqlev_canon);
            codepath_find_eqlev_level(curr->codes, firstcodes, &curr->eqlev_first);

            if ((curr->fpchild == FALSE || orbits[curr->target_vertex] == curr->target_vertex) && curr->comp_canon > -1) /* ie, not equiv to previous child */
            {
#if !MAXN
                int rtnlevel = othernode0(curr, stack, &tcnode0);
#else
            int rtnlevel = othernode(curr, stack);
#endif

                if (rtnlevel < curr->level - 1)
                {
#ifdef MPI
                    canonstate.canonpath = &canonpath;
                    canonstate.canoncodes = &canoncodes;
                    canonstate.canonlab = canonlab;
                    canonstate.canonlevel = &stglb_canonlevel;
                    canonstate.samerows = &stglb_samerows;
                    canonstate.worktop = worktop;
                    canonstate.fmptr = &fmptr;
                    canonstate.numnodes = stats->numnodes;
                    mpi_send_active_prune_message(mpi_state, curr->path, rtnlevel + 1, &canonstate);

#endif /* ifdef MPI */
                    active_prune_all_children(curr->path, rtnlevel + 1, stack);
                }
            }
#ifdef MPI
            if (stats->numnodes > last_comm_check + nodes_per_poll || mpi_state->my_rank == 0)
            {
                last_comm_check = stats->numnodes;
                canonstate.canonpath = &canonpath;
                canonstate.canoncodes = &canoncodes;
                canonstate.canonlab = canonlab;
                canonstate.canonlevel = &stglb_canonlevel;
                canonstate.samerows = &stglb_samerows;
                canonstate.worktop = worktop;
                canonstate.fmptr = &fmptr;
                canonstate.numnodes = stats->numnodes;
                mpi_poll_for_messages(mpi_state, &canonstate);
                mpi_report_stack_size(mpi_state, &canonstate);
            }
#endif /* if MPI */
        } // while (stack_size(stack) > 0)   <--- Main Loop

#ifdef MPI

        /* We are here, and out of work, need to ask for more, if we are not process 0 */
        if (mpi_state->my_rank != 0)
        {
            canonstate.canonpath = &canonpath;
            canonstate.canoncodes = &canoncodes;
            canonstate.canonlab = canonlab;
            canonstate.canonlevel = &stglb_canonlevel;
            canonstate.samerows = &stglb_samerows;
            canonstate.worktop = worktop;
            canonstate.fmptr = &fmptr;
            canonstate.numnodes = stats->numnodes;
            mpi_ask_for_work(mpi_state, &canonstate);
        }
        else
        {
            /**
             * IF we are process 0, our queue is empty, first we will check if others are still workig, if so we are going to ask for more work, and re-enter the running state.
             */
        
            int number_running = 0;
            for(int i = 1; i < mpi_state->num_processes; ++i) if (canonstate.process_waiting_for_work[i] != TRUE) ++number_running;
            if (number_running > 0) {
                /* HERE because there is at least one process still running */
                last_comm_check = stats->numnodes;
                canonstate.canonpath = &canonpath;
                canonstate.canoncodes = &canoncodes;
                canonstate.canonlab = canonlab;
                canonstate.canonlevel = &stglb_canonlevel;
                canonstate.samerows = &stglb_samerows;
                canonstate.worktop = worktop;
                canonstate.fmptr = &fmptr;
                canonstate.numnodes = stats->numnodes;
                mpi_poll_for_messages(mpi_state, &canonstate);            
                
                mpi_ask_for_work_return(mpi_state, &canonstate);
            } else {
                /* HERE because all other processes have stopped */
                if (__DEBUG_MPI_PEROIDIC_STATUS__) {
                    printf("%011f MPI: Process %d: WORK END STATUS: [0] %d, -1", MPI_Wtime(), mpi_state->my_rank, stack_size(canonstate.stack));
                    for(int i = 1; i < mpi_state->num_processes; ++i) printf("   [%d] %d, %d", i, canonstate.queue_sizes[i], canonstate.process_waiting_for_work[i]);
                    printf("\n");
                }
                mpi_broadcast_work_request_reject(mpi_state, &canonstate);  /* send reject to all processes, which will put them in the MPI_STATE_WAIT_FOR_WORK_END state */
                mpi_state->state = MPI_STATE_WAIT_FOR_WORK_END; /* put ourselves in the MPI_STATE_WAIT_FOR_WORK_END state, so we can exit the main loop */
            }
             /* OLD LOGIC HERE */
            // /* if we are process 0, we enter the work end state */
            // mpi_state->state = MPI_STATE_WAIT_FOR_WORK_END;
        }

    } /* while (mpi_state->state == MPI_STATE_WORKING || mpi_state->state == MPI_STATE_ASKING_FOR_WORK || mpi_state->state == MPI_STATE_WAIT_FOR_QUEUE_SIZE_RESPONSES) */

    /* safety check that we are actually in work end state */
    if (mpi_state->state != MPI_STATE_WAIT_FOR_WORK_END)
    {
        /* if we are in the wrong state, print error message and exit with a non zero state, so we don't miss the error */
        printf("MPI: Process %d: Fell out of the main work loop while not in MPI_STATE_WAIT_FOR_WORK_END state, in state is %d\n", mpi_state->my_rank, mpi_state->state);
        exit(1);
    }

    /**
     * Process zero still has to direct the rest of the operation.  Needs to wait for all processes to finish, and continue forwarding broadcast messages.
     */
    if (mpi_state->my_rank == 0) {
        /* this shouldn't be needed, as process 0 no longer exist the main work loop, until all other processes have completed */
        // recalculate_processes_running(canonstate.queue_sizes, mpi_state->num_processes);
        // while (canonstate.queue_sizes[0]) {
        //     mpi_poll_for_messages(mpi_state, &canonstate);            
        //     recalculate_processes_running(canonstate.queue_sizes, mpi_state->num_processes);
        // } /* while (canonstate.queue_sizes[0]) */
        
        FREES(canonstate.queue_sizes);
    }

    printf("MPI: Process %d: finished main work loop, refined: %zu\n", mpi_state->my_rank, stats->numnodes);
// TODO: DELETE THIS
// fprintf(tracefile, "MPI: Process %d: %.10f Finished\n", mpi_state->my_rank, MPI_Wtime()); 
// TO HERE 

// TODO: Delete This
// fclose(tracefile);
// TO HERE

    /** visualize automorphisms */
    // dump_automorphisms(mpi_state->my_rank, workspace, fmptr, m);

    /** collect statistics */
    unsigned long total_numnodes;
    MPI_Reduce(&stats->numnodes, &total_numnodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_state->my_rank == 0)
    {
        stats->numnodes = total_numnodes;
    }
    /** */
#endif /* ifdef MPI */

    if (retval == NAUTY_ABORTED)
        stats->errstatus = NAUABORTED;
    else if (retval == NAUTY_KILLED)
        stats->errstatus = NAUKILLED;
    else
    {
        if (getcanon && run_first_path_node)
        {
            // (*dispatch.updatecan)(g,canong,canonlab,stglb_samerows,M,n);
            (*dispatch.updatecan)(g, canong, canonlab, 0, M, n);
            for (i = 0; i < n; ++i)
                root->lab[i] = canonlab[i];
        }
        stats->invarsuclevel =
            (invarsuclevel == NAUTY_INFINITY ? 0 : invarsuclevel);
        stats->invapplics = invapplics;
        stats->invsuccesses = invsuccesses;
    }

#if !MAXN
#ifndef NAUTY_IN_MAGMA
    if (n >= 320)
#endif
    {
        nautil_freedyn();
        OPTCALL(dispatch.freedyn)();
        nauty_freedyn();
    }
#endif
    if (run_first_path_node)
    {
        OPTCALL(dispatch.cleanup)(g_arg, &g, canong_arg, &canong, root->lab, root->ptn, options, stats, m, n);
    } // if (run_first_path_node)

    if (doschreier)
    {
        freeschreier(&gp, &gens);
        if (n >= 320)
            schreier_freedyn();
    }

    // printf("MPI: Process %d: finished canonpath: ", mpi_state->my_rank); path_visualize(canonpath); printf("\n");
    // printf("MPI: Process %d: finished canoncodes: ", mpi_state->my_rank); codepath_visualize(canoncodes); printf("\n");
    // printf("MPI: Process %d: finished canonlab: ", mpi_state->my_rank); for (i = 0; i < n; ++i) printf("%d ", canonlab[i]); printf("\n");

    /**
     * Free space used by root node
     *
     * Not using FREENODE here, because FREENODE would FREE the lab and ptn as well,
     * but at this level, those are owned by the calling function.
     */
    if (run_first_path_node)
    {
        // FREENODE(root);  /* oops */

        if (root)
        { // FREES(root->lab);
            // FREES(root->ptn);
            FREECODEPATH(root->codes);
            FREEPATH(root->path);
            FREES(root->fixedpts);
            FREES(root->active);
            // FREES(root->tcell);
            FREES(root);
            root = NULL;
        }
    }
    FREEPATH(firstpath);
    FREEPATH(canonpath);
    while (stack->top != NULL)
    {
        StackNode *sn = stack->top;
        stack->top = sn->next;
        FREENODE(sn->node);
        FREES(sn);
    }

#ifdef MPI
    // double mpi_end_time = MPI_Wtime();  /* End time for MPI operations */
    // printf("MPI: Process %d: finished nauty in %.3f seconds\n", mpi_state->my_rank, mpi_end_time - mpi_start_time);

    if (mpi_state->my_rank != 0)
        stats->errstatus = -1; /* set errstatus to -1 to indicate to the calling code, that non zero rank processes don't report out results */
    /** Shut down MPI and exit */
    printf("MPI: Process %d: shutting down normally\n", mpi_state->my_rank);
    MPI_Finalize();
    FREES(mpi_state);
    /** */
#endif /* if MPI */
}

/**
 * First Node Here
 */

static int
#if !MAXN
firstnode0(Node *node, Stack *stack, tcnode *tcnode_parent)
#else
firstnode(Node *node, Stack *stack)
#endif
{
    int target_vertex;
    int first_target_vertex, index, rtnlevel, tcellsize, target_cell, childcount, qinvar, refcode;

#if !MAXN
    set *tcell;
    tcnode *tcnode_this;

    tcnode_this = tcnode_parent->next;
    if (tcnode_this == NULL)
    {
        if ((tcnode_this = (tcnode *)ALLOCS(1, sizeof(tcnode))) == NULL ||
            (tcnode_this->tcellptr = (set *)ALLOCS(alloc_m, sizeof(set))) == NULL)
            alloc_error("tcell");
        tcnode_parent->next = tcnode_this;
        tcnode_this->next = NULL;
    }
    tcell = tcnode_this->tcellptr;
#else
    set tcell[MAXM];
#endif

    if (__DEBUG_REF__)
    {
#ifdef MPI
        printf("MPI: Process %d: ", mpi_state->my_rank);
#endif
        printf("REFINE ");
        path_visualize(node->path);
        printf("\n");
    }
    ++stats->numnodes;
    /* refine partition : */
    doref(g, node->lab, node->ptn, node->level, &(node->numcells), &qinvar, workperm,
          node->active, &refcode, dispatch.refine, invarproc,
          mininvarlevel, maxinvarlevel, invararg, digraph, M, n);

    FREECODEPATH(firstcodes);
    DYNALLOCCODEPATH(firstcodes, node->level + 1, "firstcodes in firstnode");
    firstcodes->data[0] = 0;
    firstcodes->data[node->level] = (short)refcode;

    if (qinvar > 0)
    {
        ++invapplics;
        if (qinvar == 2)
        {
            ++invsuccesses;
            if (mininvarlevel < 0)
                mininvarlevel = node->level;
            if (maxinvarlevel < 0)
                maxinvarlevel = node->level;
            if (node->level < invarsuclevel)
                invarsuclevel = node->level;
        }
    }

    target_cell = -1;
    if (node->numcells != n)
    {
        /* locate new target cell, setting tc to its position in lab, tcell
                         to its contents, and tcellsize to its size: */
        maketargetcell(g, node->lab, node->ptn, node->level, tcell, &tcellsize,
                       &target_cell, tc_level, digraph, -1, dispatch.targetcell, M, n);
        stats->tctotal += tcellsize;
    }
    firsttc[node->level] = target_cell;

    /* optionally call user-defined node examination procedure: */
    OPTCALL(usernodeproc)
    (g, node->lab, node->ptn, node->level, node->numcells, target_cell, (int)firstcodes->data[node->level], M, n);

    if (node->numcells == n) /* found first leaf? */
    {
        firstterminal(node->lab, node->level, node);
        OPTCALL(userlevelproc)(node->lab, node->ptn, node->level, orbits, stats, 0, 1, 1, n, 0, n);
        if (getcanon && usercanonproc != NULL)
        {
            (*dispatch.updatecan)(g, canong, canonlab, stglb_samerows, M, n);
            stglb_samerows = n;
            if ((*usercanonproc)(g, canonlab, canong, stats->canupdates,
                                 (int)canoncodes->data[node->level], M, n))
                return NAUTY_ABORTED;
        }
        return node->level - 1; /* Not used for anything, right?? */
    }

#ifdef NAUTY_IN_MAGMA
    if (main_seen_interrupt)
        return NAUTY_KILLED;
#else
    if (nauty_kill_request)
        return NAUTY_KILLED;
#endif

    /* use the elements of the target cell to produce the children: */
    index = 0;
    int temp_queue[tcellsize];
    int i = tcellsize - 1;
    /**
     * Putting the child elements in an array in reverse order, as we are going to add them to
     * the stack, and stack is LIFO, but we want to process the first one first!
     */
    for (first_target_vertex = target_vertex = nextelement(tcell, M, -1); target_vertex >= 0;
         target_vertex = nextelement(tcell, M, target_vertex))
    {
        temp_queue[i] = target_vertex;
        --i;
    }
    Node *child_node;
    for (i = 0; i < tcellsize; ++i)
    {
        child_node = node_make_child(node, tcell, n, m);
        child_node->fpchild = TRUE;
        child_node->target_cell = target_cell;
        child_node->target_vertex = temp_queue[i];
        child_node->path->data[child_node->path->sz - 1] = child_node->target_vertex;
        child_node->codes->data[child_node->codes->sz - 1] = (short)refcode;
        stack_push(stack, child_node);
    }

    stglb_stabvertex = first_target_vertex;

    return node->level - 1; /* Not used for anything, right?? */
}

/**
 *
 *  firstpathnode produces a node on the leftmost
 *  path down the tree.  The parameters describe the level and the current
 *  colour partition.  The set of active cells is taken from the global set
 *  'active'.  If the refined partition is not discrete, the leftmost child
 *  is produced by calling firstpathnode, and the other children by calling
 *  othernode.
 *  For MAXN=0 there is an extra parameter: the address of the parent tcell
 *  structure.
 *  The value returned is the level to return to.
 *
 *  FUNCTIONS CALLED: (*usernodeproc)(),doref(),cheapautom(),
 *                    firstterminal(),nextelement(),breakout(),
 *                    firstpathnode(),othernode(),recover(),writestats(),
 *                    (*userlevelproc)(),(*tcellproc)(),shortprune()
 *
 */

static int
#if !MAXN
firstpathnode0(Node *node, Stack *stack, tcnode *tcnode_parent)
#else
firstpathnode(Node *node, Stack *stack)
#endif
{
    int target_cell = node->target_cell;
    int target_vertex = node->target_vertex;
    int index, rtnlevel, tcellsize, childcount, qinvar, refcode;

    /**
     * Moved this here so the noncheaplevel is computed for the current node
     * before the breakout happens.  This is essentially the order the process
     * happens in when this is recursive, only in the code, it is recalculated
     * before the breakout at the end.
     */
    if (node->noncheaplevel >= node->level - 1 && !(*dispatch.cheapautom)(node->ptn, node->level, digraph, n))
        node->noncheaplevel = node->level;

    /**
     * breakout individualizes the target_vertex, and splits it from the target_cell.
     * This is done so that we can further refine the partition.  node->active is updated
     * as appropriate.
     */
    breakout(node->lab, node->ptn, node->level, target_cell, target_vertex, node->active, M);
    ++node->numcells;
    ADDELEMENT(node->fixedpts, target_vertex);
    node->cosetindex = target_vertex;

#if !MAXN
    set *tcell;
    tcnode *tcnode_this;

    tcnode_this = tcnode_parent->next;
    if (tcnode_this == NULL)
    {
        if ((tcnode_this = (tcnode *)ALLOCS(1, sizeof(tcnode))) == NULL ||
            (tcnode_this->tcellptr = (set *)ALLOCS(alloc_m, sizeof(set))) == NULL)
            alloc_error("tcell");
        tcnode_parent->next = tcnode_this;
        tcnode_this->next = NULL;
    }
    tcell = tcnode_this->tcellptr;
#else
    set tcell[MAXM];
#endif

    if (__DEBUG_REF__)
    {
#ifdef MPI
        printf("MPI: Process %d: ", mpi_state->my_rank);
#endif
        printf("REFINE ");
        path_visualize(node->path);
        printf("\n");
    }
    ++stats->numnodes;
    /* refine partition : */
    doref(g, node->lab, node->ptn, node->level, &node->numcells, &qinvar, workperm,
          node->active, &refcode, dispatch.refine, invarproc,
          mininvarlevel, maxinvarlevel, invararg, digraph, M, n);

    /**
     * Adding this in here to resize firstcodes as needed.  This was causing an issue with memory overflow, as the original size
     * was way too small.   The addition of 10 is completely arbitrary, and could probably be optimized.
     */
    if (firstcodes->sz <= node->level)
    {
        firstcodes = codepath_resize_to_new(firstcodes, firstcodes->sz + 10);
    }
    firstcodes->data[node->level] = (short)refcode;

    if (__DEBUG__)
        printf("firstpathnode called at level %d  refcode %d   numcells %d\n", node->level, (short)refcode, node->numcells);
    if (__DEBUG__)
    {
        printf("numnodes: %lu\n", stats->numnodes);
        if (__DEBUG_VERBOSE__)
        {
            path_visualize(node->path);
            printf("  fpchild: %d", node->fpchild);
            printf("\n");
            node_visualize(node, n);
            printf("noncheaplevel: %d\n", node->noncheaplevel);
            printf("\n");
        }
    }

    if (qinvar > 0)
    {
        ++invapplics;
        if (qinvar == 2)
        {
            ++invsuccesses;
            if (mininvarlevel < 0)
                mininvarlevel = node->level;
            if (maxinvarlevel < 0)
                maxinvarlevel = node->level;
            if (node->level < invarsuclevel)
                invarsuclevel = node->level;
        }
    }

    target_cell = -1;
    if (node->numcells != n)
    {
        /* locate new target cell, setting target_cell to its position in lab, tcell
                        to its contents, and tcellsize to its size: */
        maketargetcell(g, node->lab, node->ptn, node->level, tcell, &tcellsize,
                       &target_cell, tc_level, digraph, -1, dispatch.targetcell, M, n);
        stats->tctotal += tcellsize;
    }
    firsttc[node->level] = target_cell;

    /* optionally call user-defined node examination procedure: */
    OPTCALL(usernodeproc)
    (g, node->lab, node->ptn, node->level, node->numcells, target_cell, (int)firstcodes->data[node->level], M, n);

    if (node->numcells == n) /* found first leaf? */
    {
        firstterminal(node->lab, node->level, node);
        OPTCALL(userlevelproc)(node->lab, node->ptn, node->level, orbits, stats, 0, 1, 1, n, 0, n);
        if (getcanon && usercanonproc != NULL)
        {
            (*dispatch.updatecan)(g, canong, canonlab, stglb_samerows, M, n);
            stglb_samerows = n;
            if ((*usercanonproc)(g, canonlab, canong, stats->canupdates,
                                 (int)canoncodes->data[node->level], M, n))
                return NAUTY_ABORTED;
        }
        return 1; /* return 0 indicating we found the first leaf node */
    }

#ifdef NAUTY_IN_MAGMA
    if (main_seen_interrupt)
        return NAUTY_KILLED;
#else
    if (nauty_kill_request)
        return NAUTY_KILLED;
#endif

    /* use the elements of the target cell to produce the children: */
    index = 0;

    int temp_queue[tcellsize];
    int i = tcellsize - 1;

    /**
     * Putting the child elements in an array in reverse order, as we are going to add them to
     * the stack, and stack is LIFO, but we want to process the first one first!
     */
    for (target_vertex = nextelement(tcell, M, -1); target_vertex >= 0;
         target_vertex = nextelement(tcell, M, target_vertex))
    {
        temp_queue[i] = target_vertex;
        --i;
    }
    Node *child_node;
    for (i = 0; i < tcellsize; ++i)
    {
        child_node = node_make_child(node, tcell, n, m);
        child_node->fpchild = TRUE;
        child_node->target_cell = target_cell;
        child_node->target_vertex = temp_queue[i];
        child_node->path->data[child_node->path->sz - 1] = child_node->target_vertex;
        child_node->codes->data[child_node->codes->sz - 1] = (short)refcode;
        stack_push(stack, child_node);
    }

    MULTIPLY(stats->grpsize1, stats->grpsize2, index);

    return 0; /* return 0 to indicate we haven't found a leaf node, keep looking to the left */
}

/**
 *
 *  othernode produces a node other than an ancestor
 *  of the first leaf.  The parameters describe the level and the colour
 *  partition.  The list of active cells is found in the global set 'active'.
 *  The value returned is the level to return to.
 *
 *  FUNCTIONS CALLED: (*usernodeproc)(),doref(),refine(),recover(),
 *                    processnode(),cheapautom(),(*tcellproc)(),shortprune(),
 *                    nextelement(),breakout(),othernode(),longprune()
 *
 */
static int
#if !MAXN
othernode0(Node *node, Stack *stack, tcnode *tcnode_parent)
#else
othernode(Node *node, Stack *stack)
#endif
{
    int target_cell = node->target_cell;
    int target_vertex = node->target_vertex;
    int first_target_vertex, refcode, rtnlevel, tcellsize, qinvar;
    short code;
#if !MAXN
    set *tcell;
    tcnode *tcnode_this;

    tcnode_this = tcnode_parent->next;
    if (tcnode_this == NULL)
    {
        if ((tcnode_this = (tcnode *)ALLOCS(1, sizeof(tcnode))) == NULL ||
            (tcnode_this->tcellptr = (set *)ALLOCS(alloc_m, sizeof(set))) == NULL)
            alloc_error("tcell");
        tcnode_parent->next = tcnode_this;
        tcnode_this->next = NULL;
    }
    tcell = tcnode_this->tcellptr;
#else
    set tcell[MAXM];
#endif

#ifdef NAUTY_IN_MAGMA
    if (main_seen_interrupt)
        return NAUTY_KILLED;
#else
    if (nauty_kill_request)
        return NAUTY_KILLED;
#endif

    /**
     * Moved this here so the noncheaplevel is computed for the current node
     * before the breakout happens.  This is essentially the order the process
     * happens in when this is recursive, only in the code, it is recalculated
     * before the breakout at the end.
     */
    if (!(*dispatch.cheapautom)(node->ptn, node->level, digraph, n))
        node->noncheaplevel = node->level;

    /**
     * Longprune
     *
     * This is the traditional pruning of branches based on automorphisms.
     */
    if (!node->fpchild && node->level > 2)
    {
        /**
         * Need to generate tfix, which is the list of points in the trivial cells of
         * the partition.  It is used by longprune below.
         */
        set *tfix, *tmcr;
        tfix = (set *)malloc(M * sizeof(set));
        tmcr = (set *)malloc(M * sizeof(set));

        // TODO: Change this up, this is just for testing right now

        set *tmp_tcell = (set *)malloc(M * sizeof(set));

        // for (int i = node->path->sz-1; i >= 0; --i) 
        for (int i = 0; i < node->path->sz; ++i)
        {
            FILLSET(tmp_tcell, m, n);
            fmptn(node->lab, node->ptn, i + 1, tfix, tmcr, m, n); /* level-1 is i+1 because the zero index of the path is tree level 2 */
            longprune(tmp_tcell, tfix, workspace, fmptr, m);
            if (!ISELEMENT(tmp_tcell, node->path->data[i]))
            {

                // fprintf(tracefile, "MPI: Process %d: %.10f newlongprune ",mpi_state->my_rank, MPI_Wtime()); fpath_visualize(tracefile, node->path); fprintf(tracefile, "  tv: %d node-lvl: %d lvl: %d  vert: %d\n", node->target_vertex, node->level, i+2, node->path->data[i]);

                FREES(tmp_tcell);
                FREES(tfix);
                FREES(tmcr);
                return i + 1;
            }
        }


        FREES(tmp_tcell);


        // To HERE

        /**
         * trying to comment this out, as it's a repeat of the "new" method.
         */
        // fmptn(node->lab, node->ptn, node->level - 1, tfix, tmcr, m, n);

        // /**
        //  * longprune removes, if appropriate, vertices from target cell (node->tcell)
        //  * based on tfix and the automorphism details stored in the chunk of memory between
        //  * &workspace and &fmptr.
        //  */
        // longprune(node->tcell, tfix, workspace, fmptr, m);

        FREES(tfix);
        FREES(tmcr);

        /**
         * If the target vertex is no longer part of the target cell, after pruning,
         * then we can stop processing, and prune this node.
         */
        if (!ISELEMENT(node->tcell, node->target_vertex))
        {
            if (__DEBUG__)
            {
                printf("longprune ");
                path_visualize(node->path);
                printf("\n");
            }
            // TODO: Remove
            // fprintf(tracefile, "MPI: Process %d: %.10f longprune ",mpi_state->my_rank, MPI_Wtime()); fpath_visualize(tracefile, node->path); fprintf(tracefile, "  tv: %d  level: %d\n", node->target_vertex, node->level);
            // To Here
            // Not sure we need to run active_prune here, as no children of this node should have been generated, and there's no need to communicate that, as ther are no children!  <---------------------------------------------------
            // active_prune_all_children(node->path, node->level, stack);

            return node->level - 1;
        }
    }

    breakout(node->lab, node->ptn, node->level, target_cell, target_vertex, node->active, M);
    ++node->numcells;
    ADDELEMENT(node->fixedpts, target_vertex);

    ++stats->numnodes;
    if (__DEBUG_REF__)
    {
#ifdef MPI
        printf("MPI: Process %d: ", mpi_state->my_rank);
#endif
        printf("REFINE ");
        path_visualize(node->path);
        printf("\n");
    }
    /* refine partition : */
    doref(g, node->lab, node->ptn, node->level, &(node->numcells), &qinvar, workperm, node->active,
          &refcode, dispatch.refine, invarproc, mininvarlevel, maxinvarlevel,
          invararg, digraph, M, n);
    code = (short)refcode;

    // TODO: Remove
    // fprintf(tracefile, "MPI: Process %d: %.10f refine ",mpi_state->my_rank, MPI_Wtime()); fpath_visualize(tracefile, node->path); fprintf(tracefile, "\n");
    // To Here

    if (__DEBUG__)
    {
        printf("numnodes: %lu\n", stats->numnodes);
        if (__DEBUG_VERBOSE__)
        {
            path_visualize(node->path);
            printf("  cpath: ");
            path_visualize(canonpath);
            printf("  fpchild: %d  OTHER", node->fpchild);
            printf("\n");
            node_visualize(node, n);
            printf("noncheaplevel: %d\n", node->noncheaplevel);
            printf("\n");
        }
    }

    if (qinvar > 0)
    {
        ++invapplics;
        if (qinvar == 2)
        {
            ++invsuccesses;
            if (node->level < invarsuclevel)
                invarsuclevel = node->level;
        }
    }

    target_cell = -1;
    /* If children will be required, find new target cell and set tc to its
    position in lab, tcell to its contents, and tcellsize to its size: */

    if (node->numcells < n && (node->eqlev_first == node->level ||
                               (getcanon && node->comp_canon >= 0)))
    {
        if (!getcanon || node->comp_canon < 0)
        {
            maketargetcell(g, node->lab, node->ptn, node->level, tcell, &tcellsize, &target_cell,
                           tc_level, digraph, firsttc[node->level], dispatch.targetcell, M, n);
            if (target_cell != firsttc[node->level])
                node->eqlev_first = node->level - 1;
        }
        else
            maketargetcell(g, node->lab, node->ptn, node->level, tcell, &tcellsize, &target_cell,
                           tc_level, digraph, -1, dispatch.targetcell, M, n);
        stats->tctotal += tcellsize;
    }

    /* optionally call user-defined node examination procedure: */
    OPTCALL(usernodeproc)(g, node->lab, node->ptn, node->level, node->numcells, target_cell, (int)code, M, n);

    /* call processnode to classify the type of this node: */
    rtnlevel = processnode(node->lab, node->ptn, node->level, node->numcells, node);

    if (rtnlevel < node->level)
    {

        return rtnlevel;
    }

    if (needshortprune)
    {
        needshortprune = FALSE;
        shortprune(tcell, fmptr - M, M);
    }

    /* use the elements of the target cell to produce the children: */

    int temp_queue[tcellsize];
    int i = tcellsize - 1;

    /**
     * Putting the child elements in an array in reverse order, as we are going to add them to
     * the stack, and stack is LIFO, but we want to process the first one first!
     */
    for (first_target_vertex = target_vertex = nextelement(tcell, M, -1); target_vertex >= 0;
         target_vertex = nextelement(tcell, M, target_vertex))
    {

        temp_queue[i] = target_vertex;
        --i;
    }


    Node *child_node;
    for (++i; i < tcellsize; ++i)
    {
        child_node = node_make_child(node, tcell, n, m);
        child_node->fpchild = FALSE;
        child_node->target_cell = target_cell;
        child_node->target_vertex = temp_queue[i];
        child_node->path->data[child_node->path->sz - 1] = child_node->target_vertex;
        child_node->codes->data[child_node->codes->sz - 1] = code;
        stack_push(stack, child_node);
    }

    return node->level - 1;
}

/*****************************************************************************
 *                                                                            *
 *  Process the first leaf of the tree.                                       *
 *                                                                            *
 *  FUNCTIONS CALLED: NONE                                                    *
 *                                                                            *
 *****************************************************************************/

static void
firstterminal(int *lab, int level, Node *node)
{
    int i;

    stats->maxlevel = level;
    node->gca_first = stglb_allsamelevel = node->eqlev_first = level;

    firsttc[level + 1] = -1;

    for (i = 0; i < n; ++i)
        firstlab[i] = lab[i];

    /* my addition */
    FREEPATH(firstpath);
    FREEPATH(canonpath);
    firstpath = path_deep_copy(node->path);
    canonpath = path_deep_copy(node->path);
    /* end of my addition */

    if (getcanon)
    {
        stglb_canonlevel = node->eqlev_canon = node->gca_canon = level;
        node->comp_canon = 0;
        stglb_samerows = 0;
        for (i = 0; i < n; ++i)
            canonlab[i] = lab[i];
        new_cannon_found(node, TRUE);

        stats->canupdates = 1;
    }
}

/*****************************************************************************
 *                                                                            *
 *  Process a node other than the first leaf or its ancestors.  It is first   *
 *  classified into one of five types and then action is taken appropriate    *
 *  to that type.  The types are                                              *
 *                                                                            *
 *  0:   Nothing unusual.  This is just a node internal to the tree whose     *
 *         children need to be generated sometime.                            *
 *  1:   This is a leaf equivalent to the first leaf.  The mapping from       *
 *         firstlab to lab is thus an automorphism.  After processing the     *
 *         automorphism, we can return all the way to the closest invocation  *
 *         of firstpathnode.                                                  *
 *  2:   This is a leaf equivalent to the bsf leaf.  Again, we have found an  *
 *         automorphism, but it may or may not be as useful as one from a     *
 *         type-1 node.  Return as far up the tree as possible.               *
 *  3:   This is a new bsf node, provably better than the previous bsf node.  *
 *         After updating canonlab etc., treat it the same as type 4.         *
 *  4:   This is a leaf for which we can prove that no descendant is          *
 *         equivalent to the first or bsf leaf or better than the bsf leaf.   *
 *         Return up the tree as far as possible, but this may only be by     *
 *         one level.                                                         *
 *                                                                            *
 *  Types 2 and 3 can't occur if getcanon==FALSE.                             *
 *  The value returned is the level in the tree to return to, which can be    *
 *  anywhere up to the closest invocation of firstpathnode.                   *
 *                                                                            *
 *  FUNCTIONS CALLED:    isautom(),updatecan(),testcanlab(),fmperm(),         *
 *                       writeperm(),(*userautomproc)(),orbjoin(),            *
 *                       shortprune(),fmptn()                                 *
 *                                                                            *
 *****************************************************************************/

static int processnode(int *lab, int *ptn, int level, int numcells, Node *node)
{
    int i, code, save, newlevel;
    boolean ispruneok;
    int sr;

    node->gca_canon = path_greatest_common_ancestor(node->path, canonpath);
    node->gca_first = path_greatest_common_ancestor(node->path, firstpath);

    code = 0;
    if (node->eqlev_first != level && (!getcanon || node->comp_canon < 0))
        code = 4;
    else if (numcells == n)
    {
        if (node->eqlev_first == level)
        {
            for (i = 0; i < n; ++i)
                workperm[firstlab[i]] = lab[i];

            if (node->gca_first >= node->noncheaplevel ||
                (*dispatch.isautom)(g, workperm, digraph, M, n))
                code = 1;
        }
        if (code == 0)
        {
            if (getcanon)
            {
                sr = 0;
                if (node->comp_canon == 0)
                {
                    if (level < stglb_canonlevel)
                        node->comp_canon = 1;
                    else
                    {
                        (*dispatch.updatecan)(g, canong, canonlab, stglb_samerows, M, n);
                        stglb_samerows = n;
                        node->comp_canon = (*dispatch.testcanlab)(g, canong, lab, &sr, M, n);
                    }
                }
                if (node->comp_canon == 0)
                {
                    for (i = 0; i < n; ++i)
                        workperm[canonlab[i]] = lab[i];
                    code = 2;
                }
                else if (node->comp_canon > 0)
                    code = 3;
                else
                    code = 4;
            }
            else
                code = 4;
        }
    }

    if (code != 0 && level > stats->maxlevel)
        stats->maxlevel = level;

    switch (code)
    {
    case 0: /* nothing unusual noticed */
        return level;

    case 1: /* lab is equivalent to firstlab */
        if (fmptr == worktop)
            fmptr -= 2 * M;
        fmperm(workperm, fmptr, fmptr + M, M, n);

        /**
         * Moved this before the send automorphism logic, as timing was goobering things up.
         *
         * The send automorphism script can also receive a new automorphism while trying to send this one.
         * It seems this was happening, and the receipt of a new automorphism, while sending this one would
         * cause this one to be overwritten (because fmptr is pointing to this one still).  Then this will
         * increment fmptr again, which gives us the same number of autos, but one of them is uninitialized garbage.
         */
        fmptr += 2 * M;
#ifdef MPI
        canonstate.canonpath = &canonpath;
        canonstate.canoncodes = &canoncodes;
        canonstate.canonlab = canonlab;
        canonstate.canonlevel = &stglb_canonlevel;
        canonstate.samerows = &stglb_samerows;
        canonstate.worktop = worktop;
        canonstate.fmptr = &fmptr;
        canonstate.numnodes = stats->numnodes;

        mpi_send_new_automorphism(mpi_state, fmptr - (2 * M), fmptr - M, M, &canonstate); /* had to subtract 2m and m here to account for moving the fmptr earlier */
#endif                                                                                    /* ifdef MPI */
        if (writeautoms)
            writeperm(outfile, workperm, cartesian, linelength, n);
        stats->numorbits = orbjoin(orbits, workperm, n);
        ++stats->numgenerators;
        OPTCALL(userautomproc)(stats->numgenerators, workperm, orbits,
                               stats->numorbits, stglb_stabvertex, n);
        if (doschreier)
            addgenerator(&gp, &gens, workperm, n);
        return node->gca_first;

    case 2: /* lab is equivalent to canonlab */
        if (fmptr == worktop)
            fmptr -= 2 * M;
        fmperm(workperm, fmptr, fmptr + M, M, n);
        /**
         * Moved this before the send automorphism logic, as timing was goobering things up.
         *
         * The send automorphism script can also receive a new automorphism while trying to send this one.
         * It seems this was happening, and the receipt of a new automorphism, while sending this one would
         * cause this one to be overwritten (because fmptr is pointing to this one still).  Then this will
         * increment fmptr again, which gives us the same number of autos, but one of them is uninitialized garbage.
         */
        fmptr += 2 * M;
#ifdef MPI
        canonstate.canonpath = &canonpath;
        canonstate.canoncodes = &canoncodes;
        canonstate.canonlab = canonlab;
        canonstate.canonlevel = &stglb_canonlevel;
        canonstate.samerows = &stglb_samerows;
        canonstate.worktop = worktop;
        canonstate.fmptr = &fmptr;
        canonstate.numnodes = stats->numnodes;

        mpi_send_new_automorphism(mpi_state, fmptr - (2 * M), fmptr - M, M, &canonstate); /* had to subtract 2m and m here to account for moving the fmptr earlier */
#endif                                                                                    /* ifdef MPI */
        save = stats->numorbits;
        stats->numorbits = orbjoin(orbits, workperm, n);
        if (stats->numorbits == save)
        {
            if (node->gca_canon != node->gca_first)
                needshortprune = TRUE;
            return node->gca_canon;
        }
        if (writeautoms)
            writeperm(outfile, workperm, cartesian, linelength, n);
        ++stats->numgenerators;
        OPTCALL(userautomproc)(stats->numgenerators, workperm, orbits,
                               stats->numorbits, stglb_stabvertex, n);
        if (doschreier)
            addgenerator(&gp, &gens, workperm, n);
        if (orbits[node->cosetindex] < node->cosetindex)
            return node->gca_first;
        if (node->gca_canon != node->gca_first)
            needshortprune = TRUE;
        return node->gca_canon;

    case 3: /* lab is better than canonlab */
        /* my addition */
        FREEPATH(canonpath);
        canonpath = path_deep_copy(node->path);
        /* end of my addition */

        ++stats->canupdates;
        for (i = 0; i < n; ++i)
            canonlab[i] = lab[i];
        stglb_canonlevel = node->eqlev_canon = node->gca_canon = level;
        new_cannon_found(node, FALSE);
        node->comp_canon = 0;

        stglb_samerows = sr;
        if (getcanon && usercanonproc != NULL)
        {
            (*dispatch.updatecan)(g, canong, canonlab, stglb_samerows, M, n);
            stglb_samerows = n;

            if ((*usercanonproc)(g, canonlab, canong, stats->canupdates, (int)canoncodes->data[level], M, n))
                return NAUTY_ABORTED;
        }
        break;

    case 4: /* non-automorphism terminal node */
        ++stats->numbadleaves;
        break;
    } /* end of switch statement */

    /* only cases 3 and 4 get this far: */
    if (level != node->noncheaplevel)
    {
        ispruneok = TRUE;
        if (fmptr == worktop)
            fmptr -= 2 * M;
        fmptn(lab, ptn, node->noncheaplevel, fmptr, fmptr + M, M, n);
        fmptr += 2 * M;
    }
    else
        ispruneok = FALSE;

    save = (stglb_allsamelevel > node->eqlev_canon ? stglb_allsamelevel - 1 : node->eqlev_canon);
    newlevel = (node->noncheaplevel <= save ? node->noncheaplevel - 1 : save);

    if (ispruneok && newlevel != node->gca_first)
        needshortprune = TRUE;
    return newlevel;
}

/*****************************************************************************
 *                                                                            *
 *  Write statistics concerning an ancestor of the first leaf.                *
 *                                                                            *
 *  level = its level                                                         *
 *  tv = the vertex fixed to get the first child = the smallest-numbered      *
 *               vertex in the target cell                                    *
 *  cellsize = the size of the target cell                                    *
 *  index = the number of vertices in the target cell which were equivalent   *
 *               to tv = the index of the stabiliser of tv in the group       *
 *               fixing the colour partition at this level                    *
 *                                                                            *
 *  numorbits = the number of orbits of the group generated by all the        *
 *               automorphisms so far discovered                              *
 *                                                                            *
 *  numcells = the total number of cells in the equitable partition at this   *
 *               level                                                        *
 *                                                                            *
 *  FUNCTIONS CALLED: itos(),putstring()                                      *
 *                                                                            *
 *****************************************************************************/

static void
writemarker(int level, int tv, int index, int tcellsize,
            int numorbits, int numcells)
{
    char s[30];

#define PUTINT(i) \
    itos(i, s);   \
    putstring(outfile, s)
#define PUTSTR(x) putstring(outfile, x)

    PUTSTR("level ");
    PUTINT(level);
    PUTSTR(":  ");
    if (numcells != numorbits)
    {
        PUTINT(numcells);
        PUTSTR(" cell");
        if (numcells == 1)
            PUTSTR("; ");
        else
            PUTSTR("s; ");
    }
    PUTINT(numorbits);
    PUTSTR(" orbit");
    if (numorbits == 1)
        PUTSTR("; ");
    else
        PUTSTR("s; ");
    PUTINT(tv + labelorg);
    PUTSTR(" fixed; index ");
    PUTINT(index);
    if (tcellsize != index)
    {
        PUTSTR("/");
        PUTINT(tcellsize);
    }
    PUTSTR("\n");
}

/*****************************************************************************
 *                                                                            *
 *  nauty_check() checks that this file is compiled compatibly with the       *
 *  given parameters.   If not, call exit(1).                                 *
 *                                                                            *
 *****************************************************************************/

void nauty_check(int wordsize, int m, int n, int version)
{
    if (wordsize != WORDSIZE)
    {
        fprintf(ERRFILE, "Error: WORDSIZE mismatch in nauty.c\n");
        exit(1);
    }

#if MAXN
    if (m > MAXM)
    {
        fprintf(ERRFILE, "Error: MAXM inadequate in nauty.c\n");
        exit(1);
    }

    if (n > MAXN)
    {
        fprintf(ERRFILE, "Error: MAXN inadequate in nauty.c\n");
        exit(1);
    }
#endif

    if (version < NAUTYREQUIRED)
    {
        fprintf(ERRFILE, "Error: nauty.c version mismatch\n");
        exit(1);
    }

#if !HAVE_TLS
    if ((version & 1))
    {
        fprintf(ERRFILE,
                "*** Warning: program with TLS calling nauty without TLS ***\n");
    }
#endif
}

/*****************************************************************************
 *                                                                            *
 *  extra_autom(p,n)  - add an extra automorphism, hard to do correctly       *
 *                                                                            *
 *****************************************************************************/

void extra_autom(int *p, int n)
{
    if (writeautoms)
        writeperm(outfile, p, cartesian, linelength, n);
    stats->numorbits = orbjoin(orbits, p, n);
    ++stats->numgenerators;
    OPTCALL(userautomproc)(stats->numgenerators, p, orbits,
                           stats->numorbits, stglb_stabvertex, n);
}

/*****************************************************************************
 *                                                                            *
 *  extra_level(level,lab,ptn,numcells,tv1,index,tcellsize,childcount)        *
 *     creates an artificial level in the search.  This is dangerous.         *
 *                                                                            *
 *****************************************************************************/

void extra_level(int level, int *lab, int *ptn, int numcells, int tv1, int index,
                 int tcellsize, int childcount, int n)
{
    MULTIPLY(stats->grpsize1, stats->grpsize2, index);
    if (domarkers)
        writemarker(level, tv1, index, tcellsize, stats->numorbits, numcells);
    OPTCALL(userlevelproc)(lab, ptn, level, orbits, stats, tv1, index, tcellsize,
                           numcells, childcount, n);
}

/*****************************************************************************
 *                                                                            *
 *  nauty_freedyn() frees all the dynamic memory used in this module.         *
 *                                                                            *
 *****************************************************************************/

void nauty_freedyn(void)
{
#if !MAXN
    tcnode *tcp, *tcq;

    tcp = tcnode0.next;
    while (tcp != NULL)
    {
        tcq = tcp->next;
        FREES(tcp->tcellptr);
        FREES(tcp);
        tcp = tcq;
    }
    alloc_m = 0;
    tcnode0.next = NULL;
    DYNFREE(firsttc, firsttc_sz);
    // DYNFREE(canoncode,canoncode_sz);
    // DYNFREE(firstcode,firstcode_sz);
    DYNFREE(workperm, workperm_sz);
    DYNFREE(canonlab, canonlab_sz);
    DYNFREE(firstlab, firstlab_sz);
    DYNFREE(defltwork, defltwork_sz);
#endif
}

/**
 * Function called when a new best so far node has been found.
 */
static void new_cannon_found(Node *node, boolean first)
{
    if (__DEBUG_CANON__)
    {
#ifdef MPI
        printf("MPI: Process %d: ", mpi_state->my_rank);
#endif
        printf("New Cannon Found: ");
        path_visualize(node->path);
        printf("  =?=  ");
        path_visualize(canonpath);
        printf("\n");
    }
    FREECODEPATH(canoncodes);
    canoncodes = codepath_deep_copy(node->codes);
#ifdef MPI
    canonstate.canonpath = &canonpath;
    canonstate.canoncodes = &canoncodes;
    canonstate.canonlab = canonlab;
    canonstate.canonlevel = &stglb_canonlevel;
    canonstate.samerows = &stglb_samerows;
    canonstate.worktop = worktop;
    canonstate.fmptr = &fmptr;
    canonstate.numnodes = stats->numnodes;
    if (!first)
    {
        /* we don't want to send broadcast on first terminal, as the processes aren't ready to receive.  we'll send the CL as part of the initial broadcast */
        mpi_send_new_best_cl(mpi_state, canonpath, canoncodes, canonlab, stglb_canonlevel, stglb_samerows, n, &canonstate);
        // TODO: DELETE THIS
        // fprintf(tracefile, "MPI: Process %d: %.10f Send newCL ",mpi_state->my_rank, MPI_Wtime()); fpath_visualize(tracefile, canonpath); fprintf(tracefile, "\n");
        // TO HERE
    }
#endif /* ifdef MPI*/
}

/**
 * Function called when we discover a branch that needs to be pruned.
 *
 * This should remove all children from the stack, and notify other processes
 * to prune these children as well.
 */
static void active_prune_all_children(Path *path, int to_level, Stack *stack)
{

    StackNode *parent = NULL;
    StackNode *curr = stack->top;
    while (curr != NULL)
    {
        if (path_compare_to_level(curr->node->path, path, to_level))
        {
            if (__DEBUG_ACTIVE_PRUNE__)
            {
                printf("Pruning to level %d  ", to_level);
                path_visualize(curr->node->path);
                printf("  fpchild: %d", curr->node->fpchild);
                printf("  old_stack_size: %d", stack->size);
            }
            if (parent != NULL)
            {
                parent->next = curr->next;
                FREENODE(curr->node);
                FREESTACKNODE(curr);
                curr = parent->next;
                stack->size--;
            }
            else
            {
                stack->top = curr->next;
                FREENODE(curr->node);
                FREESTACKNODE(curr);
                curr = stack->top;
                stack->size--;
            }
            if (__DEBUG_ACTIVE_PRUNE__)
            {
                printf("  new_stack_size: %d", stack->size);
                printf("\n");
            }
        }
        else
        {
            parent = curr;
            curr = curr->next;
        }
    }
}

#ifdef MPI
void mpi_received_active_prune(MPIState *mpi_state, Path *path, int to_level, Stack *stack)
{
    // active_prune_all_children(path, to_level, stack);

    /**  Duplicated active_prune_all_children so we can troubleshoot during MPI without modifying the original function calls */
    StackNode *parent = NULL;
    StackNode *curr = stack->top;

    while (curr != NULL)
    {
        if (path_compare_to_level(curr->node->path, path, to_level))
        {
            // TODO: DELETE THIS
            // fprintf(tracefile, "MPI: Process %d: %.10f Recv Active Prune: level: %d path: ",mpi_state->my_rank, MPI_Wtime(), to_level); fpath_visualize(tracefile, path); fprintf(tracefile, "\n");
            // TO HERE

            if (parent != NULL)
            {
                parent->next = curr->next;
                FREENODE(curr->node);
                FREESTACKNODE(curr);
                curr = parent->next;
                stack->size--;
            }
            else
            {
                stack->top = curr->next;
                FREENODE(curr->node);
                FREESTACKNODE(curr);
                curr = stack->top;
                stack->size--;
            }
        }
        else
        {
            parent = curr;
            curr = curr->next;
        }
    }
}
#endif /* ifdef MPI */
