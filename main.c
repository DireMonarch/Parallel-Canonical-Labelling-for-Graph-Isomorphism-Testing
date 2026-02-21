#include "src/nauty2_8_9/nauty.h"
#include "src/nauty2_8_9/gtools.h"
#include "util.h"
#include <time.h>
#include <libgen.h>

#ifdef MPI
#include "mpi_routines.h"
#endif /* ifdef MPI */


static void
putam(FILE *f, graph *g, int linelength, boolean space, boolean triang,
      int m, int n)   /* write adjacency matrix */
{
    set *gi;
    int i,j;
    boolean first;

    for (i = 0, gi = (set*)g; i < n - (triang!=0); ++i, gi += m)
    {
        first = TRUE;
        for (j = triang ? i+1 : 0; j < n; ++j)
        {
            if (!first && space) putc(' ',f);
            else                 first = FALSE;
            if (ISELEMENT(gi,j)) putc('1',f);
            else                 putc('0',f);
        }
        putc('\n',f);
    }
}


void usage(int argc, char *argv[])
{
    printf("Usage:\n");
    printf("%s <input graph filename> [OPTIONS]\n", argv[0]);
    printf("\nOPTIONS:\n");
    printf("\t-p\t\tWrite canonical label as adjacency matrix to file\n");
#ifdef MPI    
    printf("\t-w <int>\tSpecify work size to send, default is %d\n", MPI_CONST_MAX_WORK_SIZE_TO_SEND);
    printf("\t-s <int>\tSpecify the work split value, default is %d\n", MPI_CONST_SPLIT_RATIO);
    printf("\t-n <int>\tSpecify the number of nodes to process between communication polls, default is %d\n", MPI_CONST_NODES_BETWEEN_COMM_POLLS);
#endif /* ifdef MPI */
    printf("\n\n");

}


int main(int argc, char *argv[]) {

    FILE *infile;
    int codetype;
    char * infilename;

    /**
     * Required CLI argument, the input graph filename
     */
    if (argc < 2){
        usage(argc, argv);
        exit(1);
    }
    infilename = argv[1];

    /**
     * Parse optional CLI arguments
     */
    boolean print_cl = FALSE;
    int max_work_size_to_send = -1;
    int work_split = -1;
    int nodes_per_poll = -1;

    for (int i = 2; i < argc; ++i){
        if (argv[i][0] != '-' || argv[i][2] != 0) {
            printf("Invalid option: %s\n\n", argv[i]);
            usage(argc, argv);
            exit(1);
        }
        switch (argv[i][1])
        {
        case 'p':
            print_cl = TRUE;
            break;
        case 'w':
            if (i+1 >= argc){
                printf("Missing value for option: %s\n\n", argv[i]);
                usage(argc, argv);
                exit(1);                
            }
            max_work_size_to_send = atoi(argv[++i]);
            if (max_work_size_to_send < 0) {
                printf("Invalid work size to send %d, must be >= 0\n", max_work_size_to_send);
                exit(1);
            }
            break;          
        case 's':
            if (i+1 >= argc){
                printf("Missing value for option: %s\n\n", argv[i]);
                usage(argc, argv);
                exit(1);                
            }
            work_split = atoi(argv[++i]);
            if (work_split < 0) {
                printf("Invalid work split value %d, must be >= 0\n", work_split);
                exit(1);
            }
            break;         
        case 'n':
            if (i+1 >= argc){
                printf("Missing value for option: %s\n\n", argv[i]);
                usage(argc, argv);
                exit(1);                
            }
            nodes_per_poll = atoi(argv[++i]);
            if (nodes_per_poll < 0) {
                printf("Invalid number of nodes to process between communication polls value %d, must be >= 0\n", nodes_per_poll);
                exit(1);
            }
            break;
        default:
            printf("Unrecognized option: %s\n\n", argv[i]);
            usage(argc, argv);
            exit(1);        
            break;
        }
    }


    infile = opengraphfile(infilename,&codetype,FALSE,1);
    if (codetype != GRAPH6 && codetype != (GRAPH6+HAS_HEADER) && codetype != SPARSE6 && codetype != (SPARSE6+HAS_HEADER)){
        printf("Unsupported graph type %d encoutered.\n", codetype);
        exit(-1);
    }
    int m, n;
    graph *g = readg(infile, NULL, 0, &m, &n);
    fclose(infile);


    /* initialize and set options to default options */
    static DEFAULTOPTIONS_GRAPH(options);

    options.getcanon = TRUE;

    /* The following optional call verifies that we are linking
    to compatible versions of the nauty routines. */
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    int *lab, *ptn, *orbits;
    int lab_sz=0, ptn_sz=0, orbits_sz=0;

    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    statsblk stats;
    stats.maxlevel = 0;

    graph *h;
    int h_sz;
    h = malloc(1);
    DYNALLOC2(graph, h, h_sz, m, n, "malloc");

    double starttime = wtime();
    int mpi_details[5] = {0, 0, -2, -3, -4};
#ifdef MPI
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, h, max_work_size_to_send, work_split, nodes_per_poll, argc, argv, mpi_details);
#else /* ifdef MPI */
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, h);
#endif /* ifdef MPI */
    double endtime = wtime();

    if (mpi_details[0] == 0) {
        if (print_cl) {
            printf("::%d\n", max_work_size_to_send);
            char filename[256];
            // if (max_work_size_to_send > 0) {
            //     sprintf(filename, "%s-%d.txt", basename(infilename), max_work_size_to_send);
            // } else {
                sprintf(filename, "%s.txt", basename(infilename));
            // }
            FILE *f = fopen(filename, "w");
            putam(f, h, 0, TRUE, FALSE, m, n);
            // putam(stdout, h, 0, TRUE, FALSE, m, n);
            fclose(f);
        }

        FILE *f = fopen("statistics.csv", "a");

        fprintf(f, "%s,%d,%d,%d,%d,%d,%lu,%f\n", infilename, n, mpi_details[1], mpi_details[2], mpi_details[3], mpi_details[4], stats.numnodes, endtime - starttime);
        fclose(f);

        printf("\n");
        printf("Errstatus: %d\n", stats.errstatus);
        printf("Verticies: %d\n", n);
        printf("Number of nodes: %lu\n", stats.numnodes);
        printf("Time of Execution: %f s\n", endtime - starttime);
        printf("\n");
    }

    FREES(lab);
    FREES(ptn);
    FREES(orbits);
    FREES(g); 

    return 0;
}