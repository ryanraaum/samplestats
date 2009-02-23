#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "simple_getopt.h"
#include "fs.h"
#include "r2.h"
#include "tajd.h"

#define PACKAGE "sample_stats3"
#define VERSION "0.0.1"

/* String containing name the program is called with. */
const char *program_name;

/*  Calculates the frequency of an allele at a particular site
 *    across all chromosomes
 *
 *      allele      - the allele to count, "A", "G", "C"
 *      site        - the position to count, in the range 0, number of sites - 1 
 *      nsam        - the number of samples in the dataset (number of rows in list)
 *      list        - a matrix containing the data, samples in rows, positions in 
 *                    columns
 *
 *  Returns an integer
 */
int frequency(char allele, int site, int nsam, char **list) 
{
    int     i,                  /* iterator */
            count;              /* running count of allele */

    count = 0;

    for (i=0; i<nsam; i++) {
        count += (list[i][site] == allele ? 1 : 0);
    }

    return count;
} 

/*  Allocates space for a sample by positions list as
 *  a two dimensional matrix of chars
 *
 *      nsam        - the number of samples (rows)
 *      len         - the number of positions on the chromosome (columns)
 *
 *  Returns a two dimensional matrix of chars  
 */
char **create_list(int nsam, int len) 
{
    int     i;                  /* iterator */
    char    **list;             /* pointer to the matrix we are creating here */

    /* Try to allocate enough memory for an array of char pointers
     * of length "nsam". Throw an error if it fails.  */
    if (!(list = (char **)malloc((unsigned)(nsam*sizeof(char*))))) {
        perror("alloc error in create_list");
    }

    /* Now walk through that intial array and try to allocate enough
     * memory for a char array of length "len" for each row (sample).
     * Throw an error if it fails. */
    for (i=0; i<nsam; i++) {
        if (!(list[i] = (char *)malloc((unsigned)(len*sizeof(char))))) {
            perror("alloc error in cmatric. 2");
        }
    }
    
    return list;
}       

/*  Allocate more space for positions in each entry of the data list
 *
 *      nsam        - total number of samples in data list
 *      nmax        - new maximum number of positions
 *      list        - the data ( samples by positions matrix of chars )
 *
 *  Returns nothing
 */
void biggerlist(int nsam, unsigned nmax, char** list ) 
{
    int     i;                  /* iterator */

    /* run through the list and attempt to expand the size of each 
     * positions array */
    for (i=0; i<nsam; i++) {
        list[i] = (char *)realloc(list[i], nmax*sizeof(char)) ;
        if (list[i] == NULL) {
            perror("realloc error. bigger");
        }
    }
}      

/* Calculate the frequencies of 'A', 'G', 'C', and 'T' at each site
 * in the data list. For each site, the site_frequencies array has
 * an array with four positions. The convention in this program is 
 * that 0 -> 'A', 1 -> 'G', 2 -> 'C', 3 -> 'T'.
 *
 *      nsam        - total number of samples in data list
 *      nsites      - total number of positions
 *      list        - the data ( samples by positions matrix of chars )
 *      site_freqs  - 2D matrix with enough rows for each site and 4 int columns
 * 
 * Returns nothing (updates the passed-by-references site_frequencies 
 */
void calculate_site_frequencies(int   nsam, 
                                int   nsites, 
                                char  **list, 
                                int   **site_freqs)
{
    int     i;                  /* iterator */

    for (i=0; i<nsites; i++) {
        site_freqs[i][0] = frequency( 'A', i, nsam, list );
        site_freqs[i][1] = frequency( 'G', i, nsam, list );
        site_freqs[i][2] = frequency( 'C', i, nsam, list );
        site_freqs[i][3] = frequency( 'T', i, nsam, list );
    }
}

/*  Calculate the number of segregating sites
 *
 *      nsam            - total number of samples
 *      nsites          - total number of sites 
 *      site_freqs      - array with nucleotide counts per site 
 *
 *  Returns an integer
 */
int num_segregating_sites(int nsam, int nsites, int **site_freqs) 
{
    int     i,j,                /* iterators */
            addone,             /* flag to say if we should increment the 
                                 * count */
            count;              /* running count of segregating sites */

    count = 0;

    /* run through site_freqs, upping the segregating site count every
     * time we find a site where one of the individual nucleotide counts
     * is neither 0 nor the max (which means there must be variation!) */
    for (i=0; i<nsites; i++) {
        addone = 0;
        for (j=0; j<4; j++) {
            if (site_freqs[i][j] != 0 && site_freqs[i][j] != nsam) {
                addone = 1;
            }
        }
        if (addone) {
            count++;
        }
    }
    
    return count;
}

/*  Calculate pi (nucleotide diversity)
 *
 *      nsam            - total number of samples
 *      nsites          - total number of sites 
 *      site_freqs      - array with counts of '1' per site 
 *
 *  Returns a double 
 */
double theta_pi(int nsam, int nsites, int **site_freqs) 
{
    int     i,j;                /* iterators */
    double  pi,                 /* what we are calculating */
            ssh,                /* sum of site homozygosity */
            nd,                 /* nsam cast to double value */
            denom;              /* nsam * (nsam - 1) */

    pi = 0.0 ;

    /* create a double with the value of nsam */
    nd = nsam;
    denom = nd * (nd - 1.0);

    for (i=0; i<nsites; i++) {
        ssh = 0.0;
        for (j=0; j<4; j++) {
            ssh += (site_freqs[i][j] * (site_freqs[i][j] - 1.0))/denom;
        }
        pi += 1.0 - ssh;
    }
    
    return pi;
}

/*  Calculates Watterson's theta
 *
 *      nsam            - total number of samples in data list
 *      segsites        - total number of segregating sites 
 *
 *  Returns a double with the calculated value
 */
double theta_w(int nsam, int segsites) 
{
    int     i;                  /* iterator */
    double  dsegsites,          /* segregating sites as double */
            denom;              /* denominator of Watterson's theta */

    dsegsites = segsites;

    denom = 0.0;
    for (i=1; i<nsam; i++) {
        denom += 1.0/i;
    }
    
    return segsites/denom;
}

/*  Count up the haplotype frequencies in the data
 *
 *      nsam            - total number of samples in data list
 *      list            - the data ( samples by positions matrix of chars )
 *      hap_freqs       - the (initialized) array of integers to fill (length nsam)
 *
 *  Returns nothing (fills in the array given)
 */
void count_haplotype_frequencies(int nsam, char **list, int *hap_freqs) 
{
    int     i, j;               /* iterators */

    /* first initialize all haplotype counts to 0 */
    for (i=0; i<nsam; i++) {
        hap_freqs[i] = 0;
    }

    /* step through the data list, comparing each haplotype
     * to all those following it */
    for (i=0; i<nsam; i++) {
        /* if a haplotype has not been seen, it's count in the 
         * hap_freqs array will be 0 */
        if (hap_freqs[i] == 0) {
            /* we automatically have one of this haplotype */
            hap_freqs[i] = 1;
            /* compare the current haplotype to all following in the list */
            for (j=i+1; j<nsam; j++) {
                if (! strcmp(list[i], list[j])) {
                    /* if we have a match, up the current haplotype's count
                     * and set the count of the matching position to -9 so
                     * we know we've already counted it */
                    hap_freqs[i] += 1;
                    hap_freqs[j] = -9;
                }
            }
        }

    }
}

/*  Count the total number of haplotypes
 *
 *      nsam            - total number of samples in data list
 *      hap_freqs       - an array with counts per haplotype
 *                         (each entry corresponds to a row in the
 *                          original data; all entries should have 
 *                          an integer value; a value of -9 indicates 
 *                          that the entry at that row has already been
 *                          counted)
 *
 *  Returns an integer
 */
int num_haplotypes(int nsam, int* hap_freqs) 
{
    int     i,                  /* iterator */
            count;              /* running total */

    count = 0;

    for (i=0; i<nsam; i++) {
        if (hap_freqs[i] > 0) {
            count++;
        }
    }
    
    return count;
}

/*  Count the total number of singletons (haplotypes with only
 *    one occurence.
 *
 *      nsam            - total number of samples in data list
 *      hap_freqs       - an array with counts per haplotype
 *                         (each entry corresponds to a row in the
 *                          original data; all entries should have 
 *                          an integer value; a value of -9 indicates 
 *                          that the entry at that row has already been
 *                          counted)
 *
 *  Returns an integer
 */
int num_singletons(int nsam, int* hap_freqs) 
{
    int     i,                  /* iterator */
            count;              /* running total */

    count = 0;

    for (i=0; i<nsam; i++) {
        if (hap_freqs[i] == 1) {
            count++;
        }
    }
    
    return count;
}

/*  Count the total number of singleton sites in the data.
 *    This is the number of variant sites appearing once
 *    and only once in the entire dataset.
 *
 *      nsites          - total number of sites 
 *      site_freqs      - array with counts of '1' per site 
 *
 *  Returns an integer
 */
int num_singleton_sites(int nsites, int **site_freqs) 
{
    int     i,j,                /* iterators */
            num_ones,           /* how many sites have a frequency of 1 */
            num_gt_zero,        /* how many sites have a frequency > 0 */
            count;              /* running total */

    count = 0;

    for (i=0; i<nsites; i++) {
        num_ones = 0;
        num_gt_zero = 0;
        for (j=0; j<4; j++) {
            if (site_freqs[i][j] == 1)
                num_ones++;
            if (site_freqs[i][j] > 0)
                num_gt_zero++;
        }
        if (num_ones == 1 && num_gt_zero == 2)
          count++;
    }
    
    return count;
}

/*  Calculate homozygosity.  This is done by summing up the squares
 *    of haplotype proportions.
 *
 *      nsam            - total number of samples in data list
 *      hap_freqs       - an array with counts per haplotype
 *                         (each entry corresponds to a row in the
 *                          original data; all entries should have 
 *                          an integer value; a value of -9 indicates 
 *                          that the entry at that row has already been
 *                          counted)
 *
 *  Returns a double
 */
double homozygosity(int nsam, int* hap_freqs) 
{
    int     i;              /* iterator */
    double  total,          /* double value for the total number of samples */
            count,          /* temporary holder for each haplotype count */
            proportion,     /* temporary holder for each haplotype proportion */
            ho;             /* calculated as we go 
                             * (sum of per site homozygosities) */

    ho = 0.0;
    total = nsam;

    for (i=0; i<nsam; i++) {
        if (hap_freqs[i] > 0) {
            count = hap_freqs[i];
            proportion = count/total;
            ho += proportion*proportion;
        }
    }
    
    return ho;
}

/* Print help info. */
static void print_help (void) 
{
  printf ("Usage: %s [OPTIONS]\n", program_name);

  puts ("");
  fputs ("\
    -h        display this help and exit\n\
    -v        display version information and exit\n", stdout);

  puts ("");
  fputs ("\
  Use these options to calculate:\n\
    -p        pi:     nucleotide diversity\n\
    -S        ss:     number of segregating sites\n\
    -W        thetaW: Watterson's theta\n\
    -D        D:      Tajima's D\n\
    -H        ho:     Homozygosity\n\
    -n        nh:     number of haplotypes\n\
    -s        ns:     number of singletons\n\
    -N        nss:    number of singleton sites\n\
    -R        r2:     Ramos-Onsins and Rozas' R2\n\
    -U        Fs:     Fu's Fs\n", stdout);

  puts ("");
  fputs ("\
Calculate summary statistics from 'seq-gen' output. With no options specified,\n\
the core set calculated by Hudson's original 'sample_stats' program will be\n\
written to the screen, minus Fay's H and the difference between nucleotide\n\
diversity (pi) and Fay's H.  Calculating Fay's H requires us to know the\n\
ancestral and derived states at each site. While this is easy for 'ms'\n\
generated data, where the derived state is known (1), for sequence data, we\n\
don't know the derived state. Thus, I've dropped that calculation. Otherwise,\n\
the core set calculated with no options given is: pi, ss, D\n\
\n\
Examples:\n\
      ms 10 1 -T > treefile\n\
      seq-gen -mHKY -l 40 < treefile > seqdata\n\
\n\
      sample_stats3 < seqdata\n\
      sample_stats3 -SW < seqdata\n\
      sample_stats3 -pSWDHns < seqdata", stdout);
  printf ("\n");
}

/* Print version and copyright information. */
static void print_version (void) 
{
  printf ("(%s) version %s\n", PACKAGE, VERSION);
}

int main(int argc, char *argv[]) {
    int     i,                  /* iterator */
            maxline,            /* size of the line buffer */
            nsam,               /* number of samples in the dataset */
            nsites,             /* number of sites in the dataset */
            segsites,           /* number of segregating sites in the dataset */
            nextsam,            /* number of samples in the next replicate */
            nextsites,          /* number of sites in the next replicate */
            nh,                 /* number of haplotypes */
            intbuf,             /* integer buffer used for input processing */
            throwaway;          /* integer buffer used for input processing */
    
    char    **list,             /* a matrix containing the data, 
                                 *   samples in rows, positions in columns*/
            smallbuf[100],      /* small character buffer used to read in the 
                                 * first line of the phylip formatted data */
            *line;              /* temporary string to hold each line as it gets
                                 *   pulled from stdin */
    
    double  pi;                 /* nucleotide diversity */
    
    int     **site_frequencies, /* Array holding the nucleotide counts per site. 
                                 *   This is a 2D array with rows corresponding
                                 *   to sites and 4 columns, indexed in the normal
                                 *   C-style ( 0, 1, 2, 3 ) corresponding to
                                 *   'A', 'G', 'C', 'T' respectively. */
            *hap_frequencies,   /* array holding unique haplotypes counts */
            *unic_frequencies;  /* array holding count of unique sites per 
                                 *   sequence */

    int     ss_flag,            /* 0 or 1; output the number of segregating 
                                 *         sites */
            pi_flag,            /* 0 or 1; output nucleotide diversity */
            tw_flag,            /* 0 or 1; output Watterson's theta (thetaW) */
            nh_flag,            /* 0 or 1; output the number of haplotypes */
            ns_flag,            /* 0 or 1; output the number of singleton 
                                 *         haplotypes */
            ho_flag,            /* 0 or 1; output homozygosity */
            nss_flag,           /* 0 or 1; output the number of singleton 
                                 *         sites */
            td_flag,            /* 0 or 1; output Tajima's D */
            r2_flag,            /* 0 or 1; output Romas-Onsins & Rozas' R2 */
            fs_flag;            /* 0 or 1; output Fu's Fs */
    
    char    ch;                 /* current character iterator for getopt 
                                 *   option parsing */

    program_name = argv[0];

    /* by default, we print a similar set to the original sample_stats 
     * (Fay's H requires us to know the ancestral and derived states, and
     *  while it is easy to know this in ms output - 1 is the derived state - 
     *  we would have to estimate it in some way for seq-gen generated 
     *  sequence data. So we skip Fay's H and, by extention, H, the difference
     *  between nucleotide diversity (pi) and Fay's H.) */
    ss_flag = 1;
    pi_flag = 1;
    tw_flag = 0;
    nh_flag = 0;
    ns_flag = 0;
    ho_flag = 0;
    nss_flag= 0;
    td_flag = 1;
    r2_flag = 0;
    fs_flag = 0;

    /* However, if there are command line options, zero out all flags and
     * print only those specified in the options */
    if (argc > 1) {
      ss_flag = 0;
      pi_flag = 0;
      td_flag = 0;
    }

    /* Use getopt to parse the following flags:
     *      S - number of segregating sites
     *      p - pi
     *      W - Watterson's theta
     *      D - Tajima's D
     *      H - Homozygosity
     *      n - number of haplotypes
     *      s - number of singletons
     *      N - number of singleton sites
     *      R - Ramos-Onsins & Rozas' R2
     *      U - Fu's Fs
     *      */
    while ((ch = getopt(argc, argv, "SpWDHnshvNRU")) != -1) {
        switch (ch) {
	        case 'S':
		        ss_flag = 1;
		        break;
            case 'p':
                pi_flag = 1;
                break;
            case 'W':
                tw_flag = 1;
                break;
            case 'H':
                ho_flag = 1;
                break;
            case 'n':
                nh_flag = 1;
                break;
            case 's':
                ns_flag = 1;
                break;
            case 'D':
                td_flag = 1;
                break;
            case 'N':
                nss_flag = 1;
                break;
            case 'R':
                r2_flag = 1;
                break;
            case 'U':
                fs_flag = 1;
                break;
            case 'h':
                print_help();
                exit (EXIT_SUCCESS);
                break;
            case 'v':
                print_version();
                exit (EXIT_SUCCESS);
                break;
            case '?':
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                exit (EXIT_FAILURE);
		        break;
        }
    }

    /* read in the first line, bail out if no data */
    if (fgets(smallbuf, 1000, stdin) == NULL)
        exit(EXIT_FAILURE);

    /* first line should be phylip-style " NSAM NSITES", so pull that data out */
    sscanf(smallbuf, "%d %d", &nsam, &nsites);
    /* make sure our nsam and segsites variables are at least greater than zero */
    if (nsam <= 0 || nsites <= 0)
        exit(EXIT_FAILURE);

    /* initialize the two dimensional char matrix <list> that will hold our data */
    list = create_list(nsam,nsites+1);

    maxline = nsites + 100;

    /* initialize the line buffer */
    line = (char *)malloc((maxline+1)*sizeof(char));

    /* allocate enough space for the number of sites in the data */
    site_frequencies = (int **)malloc(nsites*sizeof(int*));

    /* initialize each site array */
    for (i=0; i<nsites; i++)
        site_frequencies[i] = (int *)malloc((4*sizeof(int)));

    /* allocate enough space for the number of samples in the data */
    hap_frequencies = (int *)malloc(nsam*sizeof(int));

    /* allocate enough space for the number of samples in the data */
    unic_frequencies = (int *)malloc(nsam*sizeof(int));

    /* repeat while we still find data (see the end of the loop) */
    while (nsam > 0 && nsites > 0) {
        
        /* for the number of samples, read in each line, first
         * pulling off the sample identifier (which is important
         * because we want to be able to match up with groups from
         * ms), then pulling the sequence and putting it in the
         * data list indexed by identifier - 1 (to work with C-style
         * zero-based indexing). */
        for (i=0; i<nsam; i++) {
            if (fgets(line, maxline, stdin) == NULL)
                exit(EXIT_FAILURE);
            sscanf(line, "%d", &intbuf);
            sscanf(line, "%d %s", &throwaway, list[intbuf-1]);
        }

        /* only perform calculations we need to */
        
        if (pi_flag || td_flag || tw_flag || ss_flag || nss_flag || r2_flag || fs_flag) 
            calculate_site_frequencies(nsam, nsites, list, site_frequencies);

        if (nh_flag || ns_flag || ho_flag || fs_flag) 
            count_haplotype_frequencies(nsam, list, hap_frequencies);
        
        /* fill in the unic_frequencies array if necessary */
        if (r2_flag)
            count_agct_unic_frequencies(nsam, nsites, list, 
                                        site_frequencies, unic_frequencies);

        if (pi_flag || td_flag || r2_flag || fs_flag) 
            pi = theta_pi(nsam, nsites, site_frequencies);

        if (ss_flag || tw_flag || td_flag) 
            segsites = num_segregating_sites(nsam, nsites, site_frequencies);

        if (nh_flag || fs_flag)
            nh = num_haplotypes(nsam, hap_frequencies);
        
        if (pi_flag)
            printf("pi:\t%lf\t", pi);
        if (ss_flag)
            printf("ss:\t%d\t", segsites);
        if (td_flag)
            printf("D:\t%lf\t", tajd(nsam, segsites, pi));
        if (tw_flag)
            printf("thetaW:\t%lf\t", theta_w(nsam, segsites));
        if (nh_flag)
            printf("num_haplotypes:\t%d\t", nh);
        if (ns_flag)
            printf("num_singletons:\t%d\t", num_singletons(nsam, hap_frequencies));
        if (ho_flag)
            printf("homozygosity:\t%lf\t", homozygosity(nsam, hap_frequencies));
        if (nss_flag)
            printf("nss:\t%d\t", num_singleton_sites(nsites, site_frequencies));
        if (r2_flag)
            printf("r2:\t%lf\t", R2(unic_frequencies, pi, nsam, segsites));
        if (fs_flag)
            printf("Fs:\t%lf\t", Fs(nsam, pi, nh));
        puts("");

        /* see if there's another replicate coming, check the number of samples
         * and sites advertised if there is, and reallocate space if the current
         * data structures aren't big enough to handle the next set */
        if (fgets(line, maxline, stdin) != NULL) {
       
            /* again, first line should be phylip-style " NSAM NSITES", 
             * so pull that data out */
            sscanf(line, " %d %d", &nextsam, &nextsites);

            /* if nextsam or nextsites are greater than the current,
             * we will have to update the sizes of our data structures */
            if (nextsam > nsam) {
                /* first free the old list */
                for (i=0; i<nsam; i++) 
                    free(list[i]);
                free(list);
                /* now create a new list */
                list = create_list(nextsam, nextsites+1);
                /* hap_frequencies needs to be expanded as well */
                hap_frequencies = (int *)realloc(hap_frequencies, 
                                                 nextsam*sizeof(int));
                /* expand unic_frequencies as well */
                unic_frequencies = (int *)realloc(unic_frequencies, 
                                                  nextsam*sizeof(int));
            }
            if (nextsites > nsites) {
                /* first, we'll need a bigger line buffer */
                maxline = nextsites + 100;
                line = (char *)realloc(line, (maxline+1)*sizeof(char));
                if (line == NULL)
                    perror("realloc error. couldn't make line bigger");

                /* free up and re-create site_frequencies structure */
                for (i=0; i<nsites; i++)
                    free(site_frequencies[i]);
                free(site_frequencies);
                /* now reallocate enough space for the number of sites in the data */
                site_frequencies = (int **)malloc(nextsites*sizeof(int*));
                /* initialize each site array */
                for (i=0; i<nextsites; i++)
                    site_frequencies[i] = (int *)malloc((4*sizeof(int)));
                
                /* finally, if nextsam is not bigger than the current nsam, 
                 * then we need to reallocate space in the list */
                biggerlist(nsam, nextsites, list);
            } 
            nsam = nextsam;
            nsites = nextsites;
        } else {
            /* we're at the end of the file, zero out nsam or segsites to
             * kill the while loop that's running this */
            nsam = 0;
            nsites = 0;
        }
    }
    
    exit(EXIT_SUCCESS);
}
