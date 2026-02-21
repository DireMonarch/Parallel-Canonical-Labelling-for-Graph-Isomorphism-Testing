/*****************************************************************************
*                                                                            *
*  Graph-specific auxiliary source file for version 2.8 of nauty.            *
*                                                                            *
*   Copyright (1984-2019) Brendan McKay.  All rights reserved.               *
*   Subject to waivers and disclaimers in nauty.h.                           *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       16-Nov-00 : initial creation out of nautil.c                         *
*       22-Apr-01 : added aproto line for Magma                              *
*                   EXTDEFS is no longer required                            *
*                   removed dynamic allocation from refine1()                *
*       21-Nov-01 : use NAUTYREQUIRED in naugraph_check()                    *
*       23-Nov-06 : add targetcell(); make bestcell() local                  *
*       10-Dec-06 : remove BIGNAUTY                                          *
*       10-Nov-09 : remove shortish and permutation types                    *
*       23-May-10 : add densenauty()                                         *
*       15-Jan-12 : add TLS_ATTR attributes                                  *
*       23-Jan-13 : add some parens to make icc happy                        *
*       15-Oct-19 : fix default size of dnwork[] to match densenauty()       *
*        6-Apr-21 : increase work space in densenauty()                      *
*                                                                            *
*****************************************************************************/

#define ONE_WORD_SETS
#include "nauty.h"

    /* macros for hash-codes: */
#define MASH(l,i) ((((l) ^ 065435) + (i)) & 077777)
    /* : expression whose long value depends only on long l and int/long i.
     Anything goes, preferably non-commutative. */

#define CLEANUP(l) ((int)((l) % 077777))
    /* : expression whose value depends on long l and is less than 077777
     when converted to int then short.  Anything goes. */

#if  MAXM==1
#define M 1
#else
#define M m
#endif

/* aproto: header new_nauty_protos.h */

dispatchvec dispatch_graph =
  {isautom,testcanlab,updatecan,refine,refine1,cheapautom,targetcell,
   naugraph_freedyn,naugraph_check,NULL,NULL};

#if !MAXN
DYNALLSTAT(set,workset,workset_sz);
DYNALLSTAT(int,workperm,workperm_sz);
DYNALLSTAT(int,bucket,bucket_sz);
DYNALLSTAT(set,dnwork,dnwork_sz);
#else
static TLS_ATTR set workset[MAXM];   /* used for scratch work */
static TLS_ATTR int workperm[MAXN];
static TLS_ATTR int bucket[MAXN+2];
static TLS_ATTR set dnwork[2*500*MAXM];
#endif

/*****************************************************************************
*                                                                            *
*  isautom(g,perm,digraph,m,n) = TRUE iff perm is an automorphism of g       *
*  (i.e., g^perm = g).  Symmetry is assumed unless digraph = TRUE.           *
*                                                                            *
*****************************************************************************/

boolean
isautom(graph *g, int *perm, boolean digraph, int m, int n)
{
    set *pg;
    int pos;
    set *pgp;
    int posp,i;

    for (pg = g, i = 0; i < n; pg += M, ++i)
    {
        pgp = GRAPHROW(g,perm[i],M);
        pos = (digraph ? -1 : i);

        while ((pos = nextelement(pg,M,pos)) >= 0)
        {
            posp = perm[pos];
            if (!ISELEMENT(pgp,posp)) return FALSE;
        }
    }
    return TRUE;
}

/*****************************************************************************
*                                                                            *
*  testcanlab(g,canong,lab,samerows,m,n) compares g^lab to canong,           *
*  using an ordering which is immaterial since it's only used here.  The     *
*  value returned is -1,0,1 if g^lab <,=,> canong.  *samerows is set to      *
*  the number of rows (0..n) of canong which are the same as those of g^lab. *
*                                                                            *
*  GLOBALS ACCESSED: workset<rw>,permset(),workperm<rw>                      *
*                                                                            *
*****************************************************************************/

int
testcanlab(graph *g, graph *canong, int *lab, int *samerows, int m, int n)
{
    int i,j;
    set *ph;

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n,"testcanlab");
    DYNALLOC1(set,workset,workset_sz,m,"testcanlab");
#endif

    for (i = 0; i < n; ++i) workperm[lab[i]] = i;

    for (i = 0, ph = canong; i < n; ++i, ph += M)
    {
        permset(GRAPHROW(g,lab[i],M),workset,M,workperm);
        for (j = 0; j < M; ++j)
            if (workset[j] < ph[j])
            {
                *samerows = i;
                return -1;
            }
            else if (workset[j] > ph[j])
            {
                *samerows = i;
                return 1;
            }
    }

    *samerows = n;
    return 0;
}

/*****************************************************************************
*                                                                            *
*  updatecan(g,canong,lab,samerows,m,n) sets canong = g^lab, assuming        *
*  the first samerows of canong are ok already.                              *
*                                                                            *
*  GLOBALS ACCESSED: permset(),workperm<rw>                                  *
*                                                                            *
*****************************************************************************/

void
updatecan(graph *g, graph *canong, int *lab, int samerows, int m, int n)
{
    int i;
    set *ph;

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n,"updatecan");
#endif

    for (i = 0; i < n; ++i) workperm[lab[i]] = i;

    for (i = samerows, ph = GRAPHROW(canong,samerows,M);
                                               i < n; ++i, ph += M)
        permset(GRAPHROW(g,lab[i],M),ph,M,workperm);
}

/*****************************************************************************
*                                                                            *
*  refine(g,lab,ptn,level,numcells,count,active,code,m,n) performs a         *
*  refinement operation on the partition at the specified level of the       *
*  partition nest (lab,ptn).  *numcells is assumed to contain the number of  *
*  cells on input, and is updated.  The initial set of active cells (alpha   *
*  in the paper) is specified in the set active.  Precisely, x is in active  *
*  iff the cell starting at index x in lab is active.                        *
*  The resulting partition is equitable if active is correct (see the paper  *
*  and the Guide).                                                           *
*  *code is set to a value which depends on the fine detail of the           *
*  algorithm, but which is independent of the labelling of the graph.        *
*  count is used for work space.                                             *
*                                                                            *
*  GLOBALS ACCESSED:  workset<w>,bit<r>,nextelement(),bucket<w>,workperm<w>  *
*                                                                            *
*****************************************************************************/

void
refine(graph *g, int *lab, int *ptn, int level, int *numcells,
       int *count, set *active, int *code, int m, int n)
{

#if MAXM==1
    refine1(g,lab,ptn,level,numcells,count,active,code,m,n);
}
#else

    int i,c1,c2,labc1;
    setword x;
    set *set1,*set2;
    int split1,split2,cell1,cell2;
    int cnt,bmin,bmax;
    long longcode;
    set *gptr;
    int maxcell,maxpos,hint;

#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n,"refine");
    DYNALLOC1(set,workset,workset_sz,m,"refine");
    DYNALLOC1(int,bucket,bucket_sz,n+2,"refine");
#endif

    longcode = *numcells;
    split1 = -1;
    hint = 0;
    while (*numcells < n && ((split1 = hint, ISELEMENT(active,split1))
                         || (split1 = nextelement(active,M,split1)) >= 0
                         || (split1 = nextelement(active,M,-1)) >= 0))
    {
        DELELEMENT(active,split1);
        for (split2 = split1; ptn[split2] > level; ++split2) {}
        longcode = MASH(longcode,split1+split2);
        if (split1 == split2)       /* trivial splitting cell */
        {
            gptr = GRAPHROW(g,lab[split1],M);
            for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
            {
                for (cell2 = cell1; ptn[cell2] > level; ++cell2) {}
                if (cell1 == cell2) continue;
                c1 = cell1;
                c2 = cell2;
                while (c1 <= c2)
                {
                    labc1 = lab[c1];
                    if (ISELEMENT(gptr,labc1))
                        ++c1;
                    else
                    {
                        lab[c1] = lab[c2];
                        lab[c2] = labc1;
                        --c2;
                    }
                }
                if (c2 >= cell1 && c1 <= cell2)
                {
                    ptn[c2] = level;
                    longcode = MASH(longcode,c2);
                    ++*numcells;
                    if (ISELEMENT(active,cell1) || c2-cell1 >= cell2-c1)
                    {
                        ADDELEMENT(active,c1);
                        if (c1 == cell2) hint = c1;
                    }
                    else
                    {
                        ADDELEMENT(active,cell1);
                        if (c2 == cell1) hint = cell1;
                    }
                }
            }
        }

        else        /* nontrivial splitting cell */
        {
            EMPTYSET(workset,m);
            for (i = split1; i <= split2; ++i)
                ADDELEMENT(workset,lab[i]);
            longcode = MASH(longcode,split2-split1+1);

            for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
            {
                for (cell2 = cell1; ptn[cell2] > level; ++cell2) {}
                if (cell1 == cell2) continue;
                i = cell1;
                set1 = workset;
                set2 = GRAPHROW(g,lab[i],m);
                cnt = 0;
                for (c1 = m; --c1 >= 0;)
                    if ((x = ((*set1++) & (*set2++))) != 0)
                        cnt += POPCOUNT(x);

                count[i] = bmin = bmax = cnt;
                bucket[cnt] = 1;
                while (++i <= cell2)
                {
                    set1 = workset;
                    set2 = GRAPHROW(g,lab[i],m);
                    cnt = 0;
                    for (c1 = m; --c1 >= 0;)
                        if ((x = ((*set1++) & (*set2++))) != 0)
                            cnt += POPCOUNT(x);

                    while (bmin > cnt) bucket[--bmin] = 0;
                    while (bmax < cnt) bucket[++bmax] = 0;
                    ++bucket[cnt];
                    count[i] = cnt;
                }
                if (bmin == bmax)
                {
                    longcode = MASH(longcode,bmin+cell1);
                    continue;
                }
                c1 = cell1;
                maxcell = -1;
                for (i = bmin; i <= bmax; ++i)
                    if (bucket[i])
                    {
                        c2 = c1 + bucket[i];
                        bucket[i] = c1;
                        longcode = MASH(longcode,i+c1);
                        if (c2-c1 > maxcell)
                        {
                            maxcell = c2-c1;
                            maxpos = c1;
                        }
                        if (c1 != cell1)
                        {
                            ADDELEMENT(active,c1);
                            if (c2-c1 == 1) hint = c1;
                            ++*numcells;
                        }
                        if (c2 <= cell2) ptn[c2-1] = level;
                        c1 = c2;
                    }
                for (i = cell1; i <= cell2; ++i)
                    workperm[bucket[count[i]]++] = lab[i];
                for (i = cell1; i <= cell2; ++i) lab[i] = workperm[i];
                if (!ISELEMENT(active,cell1))
                {
                    ADDELEMENT(active,cell1);
                    DELELEMENT(active,maxpos);
                }
            }
        }
    }

    longcode = MASH(longcode,*numcells);
    *code = CLEANUP(longcode);
}
#endif /* else case of MAXM==1 */

/*****************************************************************************
*                                                                            *
*  refine1(g,lab,ptn,level,numcells,count,active,code,m,n) is the same as    *
*  refine(g,lab,ptn,level,numcells,count,active,code,m,n), except that       *
*  m==1 is assumed for greater efficiency.  The results are identical in all *
*  respects.  See refine (above) for the specs.                              *
*                                                                            *
*****************************************************************************/

void
refine1(graph *g, int *lab, int *ptn, int level, int *numcells,
       int *count, set *active, int *code, int m, int n)
{
    int i,c1,c2,labc1;
    setword x;
    int split1,split2,cell1,cell2;
    int cnt,bmin,bmax;
    long longcode;
    set *gptr,workset0;
    int maxcell,maxpos,hint;

#if !MAXN 
    DYNALLOC1(int,workperm,workperm_sz,n,"refine1"); 
    DYNALLOC1(int,bucket,bucket_sz,n+2,"refine1"); 
#endif

    longcode = *numcells;
    split1 = -1;

    hint = 0;
    while (*numcells < n && ((split1 = hint, ISELEMENT1(active,split1))
                         || (split1 = nextelement(active,1,split1)) >= 0
                         || (split1 = nextelement(active,1,-1)) >= 0))
    {
        /** 
         * split1 and split2 define the current cell of the partition that
         * the algorithm is refining against.  When calculating the degree
         * for each vertex, split1 and split2 define the bounds of the cell
         * that degree is calculated from.  This is the scope of the 
         * scoped degree.
         */
        DELELEMENT1(active,split1);  /* remove the current cell from the active group, as it's being refined against now */
        for (split2 = split1; ptn[split2] > level; ++split2) {} /* Find the end of the current cell that starts at split1 */
        longcode = MASH(longcode,split1+split2);

        /**
         * The "trivial splitting cell" here is a discrete cell (one vertex) that
         * is being refined against (it's the scope of the scoped degree).
         */
        if (split1 == split2)       /* trivial splitting cell */
        {
            /**
             * gptr is the row of the adjacency matrix associated with the vertex
             * that is the cell we are refining against.  With this, we can just
             * look at a given element for each vertex in each cell we refine.
             * The degree can only be 1 or 0.
             */
            gptr = GRAPHROW(g,lab[split1],1); 

            /**
             * Loop through all cells in the partition. cell1 is the first position
             * in the cell, cell2 becaues the last position.
             */
            for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
            {
                for (cell2 = cell1; ptn[cell2] > level; ++cell2) {} /* finds the end of the cell*/
                if (cell1 == cell2) continue; /* only look at non-discrete cells */
                
                /**
                 * this section, with the while loop is sorting the lab portion
                 * of the partition with the vertices that have a degree 1 to the
                 * scope (current cell being refined against) at the lower index
                 * end of the array.
                 */
                c1 = cell1;
                c2 = cell2;
                while (c1 <= c2)
                {
                    labc1 = lab[c1];
                    if (ISELEMENT1(gptr,labc1))
                        ++c1;
                    else
                    {
                        lab[c1] = lab[c2];
                        lab[c2] = labc1;
                        --c2;
                    }
                } /* while (c1 <= c2) */

                /**
                 * Pretty sure this is checking that the cell being refined is 
                 * not degree of all zeros or all ones.  which would mean the partition
                 * is ordered, the degree one vertices at the start, the last 
                 * vertex with degree one is at c2, and that's where the split
                 * needs to be made.
                 */
                if (c2 >= cell1 && c1 <= cell2)
                {
                    ptn[c2] = level; /* split at c2 */
                    longcode = MASH(longcode,c2);
                    
                    ++*numcells; /* increment numcells (number of cells) (also dereferenced)*/

                    /**
                     * checking to see if the current cell (cell1 to cell2) is
                     * in the active set or lower part of the newly split cell is larger
                     * then the upper part of the newly split cell.
                     */
                    if (ISELEMENT1(active,cell1) || c2-cell1 >= cell2-c1)
                    {
                        ADDELEMENT1(active,c1); /* make the c1 part (the higher part) active */
                        /**
                         * I **think** this is saying if the split all ones except
                         * for one single zero degree, make the zero degree vertex
                         * the next hint.
                         */
                        if (c1 == cell2) hint = c1;
                    }
                    else /* if (ISELEMENT1(active,cell1) || c2-cell1 >= cell2-c1) */
                    {
                        ADDELEMENT1(active,cell1); /* make the lower part active */
                        /**
                         * I **think** this is saying if the split all zeros except
                         * for one single one degree, make the one degree vertex
                         * the next hint.
                         */                        
                        if (c2 == cell1) hint = cell1;
                    }
                }
            } /* for (cell1 = 0; cell1 < n; cell1 = cell2 + 1) */
        }

        else        /* nontrivial splitting cell */
        {
            /**
             * workset0 is used as a mask.  The bit positions in workset0 are set
             * to 1 for each vertex that is in the the current cell being refined
             * against, as bounded by split1 and split2.  workset0 is bitwise ANDed
             * with the line of the adjacency marix for any give cell to make 
             * sure only the vertices who are members of the cell being refined 
             * against are considered for the degree calculation.
             */
            workset0 = 0; /* clear workset0 */
            for (i = split1; i <= split2; ++i)
                ADDELEMENT1(&workset0,lab[i]);
            longcode = MASH(longcode,split2-split1+1);

            /**
             * Loop through all cells in the partition. cell1 is the first position
             * in the cell, cell2 becaues the last position.
             */
            for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
            {
                for (cell2 = cell1; ptn[cell2] > level; ++cell2) {} /* finds the end of the cell*/

                if (cell1 == cell2) continue;  /* only look at non-discrete cells */

                i = cell1;  /* i is used as an indexer to walk through each 
                                position in the cell*/

                /**
                 * x becomes the row out of the adjacency matrix for the current
                 * vertex in the current cell (lab[i]) bitwise anded with workset0
                 * which is a mask to only look at the current cell being refined 
                 * against.  POPCOUNT simply counts the population of x, which is
                 * simply returning the number of 1 bits in the number, which is 
                 * really representative of the degree of the given vertex (lab[i])
                 * with respect tot he current cell being refined against.
                 */
                if ((x = workset0 & g[lab[i]]) != 0)
                    cnt = POPCOUNT(x);
                else
                    cnt = 0;

                /**
                 * count[] is an array that tracks the current degree (cnt) of 
                 * the vertex at index i in the partition, with resepct to the 
                 * cell that is being refined against.  
                 * 
                 * Note that count is not indexed on the vertex, but on the index
                 * of the vertex in the partition!
                 */
                 /**
                  * bmin and bmax are min and max index bounds for bucket.  They
                  * are tracked along with the algorithm so that later in the 
                  * function, they can be used as bounds for iterating over the
                  * bucket[] array.
                  */
                 /**
                  * This is the first vertex inspected in this cell, so bmin,
                  * bmax, and count[i] should all be set to cnt
                  */
                count[i] = bmin = bmax = cnt;

                /**
                 * bucket[] is an array that tracks the distribution of each 
                 * degree (cnt).  Each time a given degree value is encountered,
                 * bucket[cnt] is incremented.
                 */
                /**
                 * because this is the first vertex inspected in this cell, 
                 * bucket[cnt] can simply be set to 1, as any previous value is
                 * useless, and should be overwritten.
                 */
                bucket[cnt] = 1;

                /* Iterate over all remaining vertices in the current cell */
                while (++i <= cell2)
                {
                    /**
                     * same as above ^^^^^
                     */
                    if ((x = workset0 & g[lab[i]]) != 0)
                        cnt = POPCOUNT(x);
                    else
                        cnt = 0;
                    
                    /**
                     * The bucket[] array is never initialized to zero.  Therefore
                     * as cnt values outside the bmin to bmax range are found, bmin
                     * or bmax need to be adjusted.  While doing so, the bucket[]
                     * array is also walked and elements initialized to zero along
                     * the way.  This reduces work needed to initialize the array
                     * to all zeros for each cell.  Only the bucket[] elements
                     * betwen bmin and bmax are used, so only they need to be 
                     * initalized.
                     */
                    while (bmin > cnt) bucket[--bmin] = 0;
                    while (bmax < cnt) bucket[++bmax] = 0;

                    /* finally increment bucket[cnt] and set the value of count[i]*/
                    ++bucket[cnt];
                    count[i] = cnt;
                } /* while (++i <= cell2) */

                /**
                 * If this is true then all elements of the current cell have the
                 * same degree.  There is no need to do any further processing
                 * of this cell at this time.
                 */
                if (bmin == bmax)
                {
                    longcode = MASH(longcode,bmin+cell1);
                    continue;
                }

                /**
                 * The next three for loops split the current cell (cell1 to cell2) 
                 * of the partion (lab/ptn) at this level, based on the 
                 * distributions in bucket[].
                 * 
                 * New cells are formed such that the vertices in each cell are 
                 * the same scoped degree to the cell being refined against. The
                 * new cells are in order of increasing degree.                  
                 */

                /**
                 * c1 and c2 are used to iterate over the current cell (cell1 to 
                 * cell2).
                 * 
                 * maxcell is used to track the length of the largest new cell as
                 * the current cell is split.  The value of maxcell is only used
                 * for determining which cell is larges, not for anything else.
                 * 
                 * maxpos (below) is used to track the starting position (in the 
                 * partition) of the first element of the largest cell in the new
                 * split of the current cell (cell1 to cell2).
                 */
                c1 = cell1;
                maxcell = -1;

                /** 
                 * Iterate over bucket[bmin] to bucket[bmax]. 
                 * 
                 * This for loop does the work to split ptn up for the new cell
                 * splits of the current cell (cell1 to cell2).  It does not
                 * adjust the positions of the verices in lab, that is done later.
                 */
                for (i = bmin; i <= bmax; ++i)
                    if (bucket[i]) /* only do anything if bucket[i] is not zero */
                    {
                        /** 
                         * c2 is set to one past the end of the newly formed cell, which 
                         * starts at c1, and is bucket[i] long.
                         */
                        c2 = c1 + bucket[i]; 
                        bucket[i] = c1; /* this is used later to rearange the vertices */
                        longcode = MASH(longcode,i+c1);
                        /**
                         * check to see if this new cell is longer than the 
                         * current maxcell.  If it is set maxcell and maxpos
                         * as appropriate.
                         */
                        if (c2-c1 > maxcell)
                        {
                            maxcell = c2-c1;
                            maxpos = c1;
                        }
                        /**
                         * if the new cell is not the first cell in the split of
                         * the current cell (cell1 to cell2), then add it to the
                         * active set.  
                         */
                        if (c1 != cell1)
                        {
                            ADDELEMENT1(active,c1);
                            /**
                             * if this new cell (c1 to c2) is discrete, set hint
                             * to be the start of the new cell.  This will help
                             * determine which cell to refine against on the next
                             * iteration.  See the while loop at the top of the
                             * function.
                             */
                            if (c2-c1 == 1) hint = c1;
                            ++*numcells; /* increment numcells (number of cells) (also dereferenced)*/
                        }
                        /**
                         * if c2, which is one past the last vertex in the new cell
                         * being created, is less than or equal to cell2, the end
                         * of the current cell (cell1 to cell2), that means this 
                         * new cell is not the last one being created, so we set
                         * the ptn[c2-1] (where c2-1 is teh actual last vertex
                         * in the new cell) equal to the level.
                         * 
                         * The reason we don't do this for the last new cell, is
                         * that the ptn[cell2] is already equal to or less than 
                         * level, as it was partitioned before, so we don't want
                         * to change that value.
                         */
                        if (c2 <= cell2) ptn[c2-1] = level;
                        /**
                         * since c2 is really one past the last vertex in the new
                         * cell being created, it is also the first vertex of the
                         * next cell to be created, so set c1 = c2.
                         */
                        c1 = c2;
                    } /* if (bucket[i]) */ /* but also */ /* for (i = bmin; i <= bmax; ++i) */
                
                /**
                 * There is a lot going on here, but this basically reoders the 
                 * vertices of the current cell (cell1 to cell2) into workperm[]
                 * so that they all end up in the right cell at the end.
                 * 
                 * This uses some trickery with the values left in bucket[] after
                 * the last for loop split the ptn[] up.  It's a bit challenging
                 * to follow, but it does work.
                 */
                for (i = cell1; i <= cell2; ++i)
                    workperm[bucket[count[i]]++] = lab[i];

                /**
                 * This for loop simply copies the values within the bounds of
                 * the current partition (cell1 to cell2) from workperm into lab
                 */
                for (i = cell1; i <= cell2; ++i) lab[i] = workperm[i];

                /**
                 * This is algoritm specific, but it check if the current cell 
                 * (cell1 to cell2) is not part of the active set, it makes it
                 * part of the active set and makes sure the first largest cell
                 * of the new split (starting at maxpos) is not in the active set.
                 * 
                 * I'm not entirely sure why that is dependent on the current cell
                 * being not in active, but perhaps it's just avoiding the largest
                 * part of the cell being partitioned against again?
                 */
                if (!ISELEMENT1(active,cell1))
                {
                    ADDELEMENT1(active,cell1);
                    DELELEMENT1(active,maxpos);
                }
            } /* for (cell1 = 0; cell1 < n; cell1 = cell2 + 1) */
        } /* else */ /* if (split1 == split2) */
    } /* ugly while */

    longcode = MASH(longcode,*numcells);
    *code = CLEANUP(longcode);
}

/*****************************************************************************
*                                                                            *
*  cheapautom(ptn,level,digraph,n) returns TRUE if the partition at the      *
*  specified level in the partition nest (lab,ptn) {lab is not needed here}  *
*  satisfies a simple sufficient condition for its cells to be the orbits of *
*  some subgroup of the automorphism group.  Otherwise it returns FALSE.     *
*  It always returns FALSE if digraph!=FALSE.                                *
*                                                                            *
*  nauty assumes that this function will always return TRUE for any          *
*  partition finer than one for which it returns TRUE.                       *
*                                                                            *
*****************************************************************************/

boolean
cheapautom(int *ptn, int level, boolean digraph, int n)
{
    int i,k,nnt;

    if (digraph) return FALSE;

    k = n;
    nnt = 0;
    for (i = 0; i < n; ++i)
    {
        --k;
        if (ptn[i] > level)
        {
            ++nnt;
            while (ptn[++i] > level) {}
        }
    }

    return (k <= nnt + 1 || k <= 4);
}

/*****************************************************************************
*                                                                            *
*  bestcell(g,lab,ptn,level,tc_level,m,n) returns the index in lab of the    *
*  start of the "best non-singleton cell" for fixing.  If there is no        *
*  non-singleton cell it returns n.                                          *
*  This implementation finds the first cell which is non-trivially joined    *
*  to the greatest number of other cells.                                    *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>,workperm<rw>,workset<rw>,bucket<rw>              *
*                                                                            *
*****************************************************************************/

static int
bestcell(const graph *g, const int *lab, const int *ptn, int level,
         int tc_level, int m, int n)
{
    int i;
    set *gp;
    setword setword1,setword2;
    int v1,v2,nnt;

#if !MAXN 
    DYNALLOC1(int,workperm,workperm_sz,n,"bestcell"); 
    DYNALLOC1(set,workset,workset_sz,m,"bestcell"); 
    DYNALLOC1(int,bucket,bucket_sz,n+2,"bestcell"); 
#endif

   /* find non-singleton cells: put starts in workperm[0..nnt-1] */

    i = nnt = 0;

    while (i < n)
    {
        if (ptn[i] > level)
        {
            workperm[nnt++] = i;
            while (ptn[i] > level) ++i;
        }
        ++i;
    }

    if (nnt == 0) return n;

    /* set bucket[i] to # non-trivial neighbours of n.s. cell i */

    for (i = nnt; --i >= 0;) bucket[i] = 0;

    for (v2 = 1; v2 < nnt; ++v2)
    {
        EMPTYSET(workset,m);
        i = workperm[v2] - 1;
        do
        {
            ++i;
            ADDELEMENT(workset,lab[i]);
        }
        while (ptn[i] > level);
        for (v1 = 0; v1 < v2; ++v1)
        {
            gp = GRAPHROW(g,lab[workperm[v1]],m);
#if  MAXM==1
            setword1 = *workset & *gp;
            setword2 = *workset & ~*gp;
#else
            setword1 = setword2 = 0;
            for (i = m; --i >= 0;)
            {
                setword1 |= workset[i] & gp[i];
                setword2 |= workset[i] & ~gp[i];
            }
#endif
            if (setword1 != 0 && setword2 != 0)
            {
                ++bucket[v1];
                ++bucket[v2];
            }
        }
    }

    /* find first greatest bucket value */

    v1 = 0;
    v2 = bucket[0];
    for (i = 1; i < nnt; ++i)
        if (bucket[i] > v2)
        {
            v1 = i;
            v2 = bucket[i];
        }

    return (int)workperm[v1];
}

/*****************************************************************************
*                                                                            *
*  targetcell(g,lab,ptn,level,tc_level,digraph,hint,m,n) returns the index   *
*  in lab of the next cell to split.                                         *
*  hint is a suggestion for the answer, which is obeyed if it is valid.      *
*  Otherwise we use bestcell() up to tc_level and the first non-trivial      *
*  cell after that.                                                          *
*                                                                            *
*****************************************************************************/

int
targetcell(graph *g, int *lab, int *ptn, int level,
       int tc_level, boolean digraph, int hint, int m, int n)
{
    int i;

    if (hint >= 0 && ptn[hint] > level &&
                     (hint == 0 || ptn[hint-1] <= level))
        return hint;
    else if (level <= tc_level)
        return bestcell(g,lab,ptn,level,tc_level,m,n);
    else
    {
        for (i = 0; i < n && ptn[i] <= level; ++i) {}
        return (i == n ? 0 : i);
    }
}

/*****************************************************************************
*                                                                            *
*  densenauty(g,lab,ptn,orbits,&options,&stats,m,n,h)                        *
*  is a slightly simplified interface to nauty().  It allocates enough       *
*  workspace for 500 automorphisms and checks that the densegraph dispatch   *
*  vector is in use.                                                         *
*                                                                            *
*****************************************************************************/
#ifdef MPI
void
densenauty(graph *g, int *lab, int *ptn, int *orbits,
           optionblk *options, statsblk *stats, int m, int n, graph *h, int max_work_size_to_send, int work_split, int nodes_per_poll,
           int argc, char** argv, int *mpi_details)
#else /* ifdef MPI */
void
densenauty(graph *g, int *lab, int *ptn, int *orbits,
    optionblk *options, statsblk *stats, int m, int n, graph *h)
#endif /* ifdef MPI */
{
    if (options->dispatch != &dispatch_graph)
    {
        fprintf(ERRFILE,"Error: densenauty() needs standard options block\n");
        exit(1);
    }

#if !MAXN
    /* Don't increase 2*500*m in the next line unless you also increase
       the default declaration of dnwork[] earlier. */
    DYNALLOC1(set,dnwork,dnwork_sz,2*500*m,"densenauty malloc");
#endif

#ifdef MPI
    nauty(g,lab,ptn,NULL,orbits,options,stats,dnwork,2*500*m,m,n,h, max_work_size_to_send, work_split, nodes_per_poll, argc, argv, mpi_details);
#else /* ifdef MPI */
    nauty(g,lab,ptn,NULL,orbits,options,stats,dnwork,2*500*m,m,n,h);
#endif /* ifdef MPI */
}

/*****************************************************************************
*                                                                            *
*  naugraph_check() checks that this file is compiled compatibly with the    *
*  given parameters.   If not, call exit(1).                                 *
*                                                                            *
*****************************************************************************/

void
naugraph_check(int wordsize, int m, int n, int version)
{
    if (wordsize != WORDSIZE)
    {
        fprintf(ERRFILE,"Error: WORDSIZE mismatch in naugraph.c\n");
        exit(1);
    }

#if MAXN
    if (m > MAXM)
    {
        fprintf(ERRFILE,"Error: MAXM inadequate in naugraph.c\n");
        exit(1);
    }

    if (n > MAXN)
    {
        fprintf(ERRFILE,"Error: MAXN inadequate in naugraph.c\n");
        exit(1);
    }
#endif

    if (version < NAUTYREQUIRED)
    {
        fprintf(ERRFILE,"Error: naugraph.c version mismatch\n");
        exit(1);
    }
}

/*****************************************************************************
*                                                                            *
*  naugraph_freedyn() - free the dynamic memory in this module               *
*                                                                            *
*****************************************************************************/

void
naugraph_freedyn(void)
{
#if !MAXN
    DYNFREE(workset,workset_sz);
    DYNFREE(workperm,workperm_sz);
    DYNFREE(bucket,bucket_sz);
    DYNFREE(dnwork,dnwork_sz);
#endif
}
