#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "simple_getopt.h"
#include "fs.h"
#include "r2.h"
#include "tajd.h"

#define PACKAGE "sample_stats2"
#define VERSION "0.0.1"

/* String containing name the program is called with. */
const char *program_name;

/* global count of the current maximum number of sites
 * used to allocate memory for the data array */
int maxsites = 1000 ;

/*  Calculates the frequency of an allele at a particular site
 *    across all chromosomes
 *
 *      allele      - the allele to count, either "1" or "0"
 *      site        - the position to count, in the range 0, number of sites - 1 
 *      nsam        - the number of samples in the dataset (number of rows in list)
 *      list        - a matrix containing the data, samples in rows, positions in columns
 *
 *  Returns an integer
 */
int frequency( char allele,int site,int nsam,  char **list) {
    int     i,                  /* iterator */
            count;              /* running count of allele */

    count = 0;

    for( i=0; i<nsam; i++) {
        count += ( list[i][site] == allele ? 1: 0 );
    }

    return( count);
} 

/*  Allocates space for a sample by positions list as
 *  a two dimensional matrix of chars
 *
 *      nsam        - the number of samples (rows)
 *      len         - the number of positions on the chromosome (columns)
 *
 *  Returns a two dimensional matrix of chars  
 */
char ** create_list(int nsam, int len) {
    int     i;                  /* iterator */
    char    **list;             /* pointer to the matrix we are creating here */

    /* Try to allocate enough memory for an array of char pointers
     * of length "nsam". Throw an error if it fails.  */
    if( ! ( list = (char **) malloc( (unsigned)( nsam*sizeof( char* )) ) ) ) {
        perror("alloc error in create_list");
    }

    /* Now walk through that intial array and try to allocate enough
     * memory for a char array of length "len" for each row (sample).
     * Throw an error if it fails. */
    for( i=0; i<nsam; i++) {
        if( ! ( list[i] = (char *) malloc( (unsigned) (len*sizeof( char )) ))) {
            perror("alloc error in cmatric. 2");
        }
    }
    
    return( list );
}

/*  Allocate more space for positions in each entry of the data list
 *
 *      nsam        - total number of samples in data list
 *      nmax        - new maximum number of positions
 *      list        - the data ( samples by positions matrix of chars )
 *
 *  Returns nothing
 */
void biggerlist(int nsam, unsigned nmax, char** list ) {
    int     i;                  /* iterator */

    /* change the global maxsites number to the new, larger value */
    maxsites = nmax;

    /* run through the list and attempt to expand the size of each 
     * positions array */
    for( i=0; i<nsam; i++) {
        list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
        if( list[i] == NULL ) {
            perror( "realloc error. bigger");
        }
    }
}                        

/*  Calculate pi (nucleotide diversity)
 *
 *      nsam            - total number of samples
 *      segsites        - total number of segregating sites 
 *                        ( length of data in site_freqs array )
 *      site_freqs      - array with counts of '1' per site 
 *
 *  Returns a double 
 */
double theta_pi( int nsam, int segsites, int *site_freqs) {
    int     s;                  /* iterator */
    double  pi,                 /* what we are calculating */
            p1,                 /* site frequency of derived allele */
            nd,                 /* nsam cast to double value */
            nnm1;               /* (N / (N Minus 1)) */

    pi = 0.0 ;

    /* create a double with the value of nsam */
    nd = nsam;

    nnm1 = nd/(nd-1.0);

    for( s = 0; s < segsites; s++) {
        /* calculate site frequency of '1' */
        p1 = site_freqs[s]/nd ;
        /* sum up the per site average distances */
        pi += 2.0*p1*(1.0 -p1)*nnm1 ;
    }
    
    return( pi ) ;
}

/*  Calculates Fay's theta_H
 *
 *      nsam            - total number of samples in data list
 *      segsites        - total number of segregating sites 
 *                        ( length of positions array )
 *      site_freqs      - array with counts of '1' per site 
 *
 *  Returns a double with the calculated value
 */
double theta_h( int nsam, int segsites, int *site_freqs) {
    int     s;                  /* iterator */
    double  pi,                 /* what we are calculating */
            p1,                 /* site count of derived allele */
            nd,                 /* nsam cast to double value */
            nnm1;               /* (N / (N Minus 1)) */

    pi = 0.0 ;

    /* create a double with the value of nsam */
    nd = nsam;

    nnm1 = nd/(nd-1.0) ;
    
    for( s = 0; s <segsites; s++) {
        p1 = site_freqs[s];
        pi += p1*p1 ; 
    }
    
    return( pi*2.0/(nd*(nd-1.0) ));
}

/*  Calculates Watterson's theta
 *
 *      nsam            - total number of samples in data list
 *      segsites        - total number of segregating sites 
 *                        ( length of positions array )
 *      site_freqs      - array with counts of '1' per site 
 *
 *  Returns a double with the calculated value
 */
double theta_w( int nsam, int segsites, int *site_freqs) {
    int     s;                  /* iterator */
    double  pi,                 /* what we are calculating */
            denom;              /* denominator of Watterson's theta */

    pi = 0.0 ;

    denom = 0.0;
    for (s=1; s < nsam; s++) {
        denom += 1.0/s;
    }

    for( s = 0; s <segsites; s++) {
        pi += 1.0/denom; 
    }
    
    return( pi ) ;
}

/*  Count up the haplotype frequencies in the data
 *
 *      nsam            - total number of samples in data list
 *      segsites        - total number of segregating sites 
 *                        ( length of positions array )
 *      list            - the data ( samples by positions matrix of chars )
 *      hap_freqs       - the (initialized) array of integers to fill (length nsam)
 *
 *  Returns nothing (fills in the array given)
 */
void count_haplotype_frequencies( int nsam, int segsites, char **list, int *hap_freqs ) {
    int     i, j;               /* iterators */

    /* first initialize all haplotype counts to 0 */
    for (i=0; i < nsam; i++) {
        hap_freqs[i] = 0;
    }

    /* If there are no segregating sites, then there is only one haplotype.
     * Set the haplotype count of the first haplotype to the number
     * of samples and all the rest to -9 */
    if ( segsites == 0 ) {
        hap_freqs[0] = nsam;
        for (i=1; i < nsam; i++) {
            hap_freqs[i] = -9;
        }
    } else {
        /* step through the data list, comparing each haplotype
         * to all those following it */
        for (i=0; i < nsam; i++) {
            /* if a haplotype has not been seen, it's count in the 
             * hap_freqs array will be 0 */
            if (hap_freqs[i] == 0) {
                /* we automatically have one of this haplotype */
                hap_freqs[i] = 1;
                /* compare the current haplotype to all following in the list */
                for (j=i+1; j < nsam; j++) {
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
int num_haplotypes(int nsam, int* hap_freqs) {
    int     i,                  /* iterator */
            count;              /* running total */

    count = 0;

    for (i=0; i < nsam; i++) {
        if ( hap_freqs[i] > 0 ) {
            count++;
        }
    }
    
    return (count);
}

/*  The size of the largest group of identical haplotypes
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
int max_identical_haplotypes(int nsam, int* hap_freqs) {
    int     i,                  /* iterator */
            max_num;            /* running largest value */

    max_num = 0;

    for (i=0; i < nsam; i++) {
        if ( hap_freqs[i] > max_num ) {
            max_num = hap_freqs[i];
        }
    }
    
    return (max_num);
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
int num_singletons(int nsam, int* hap_freqs) {
    int     i,                  /* iterator */
            count;              /* running total */

    count = 0;

    for (i=0; i < nsam; i++) {
        if ( hap_freqs[i] == 1 ) {
            count++;
        }
    }
    
    return (count);
}

/*  Count the total number of singleton sites in the data.
 *    This is the number of variant sites appearing once
 *    and only once in the entire dataset.
 *
 *      segsites        - total number of segregating sites 
 *                        ( length of positions array )
 *      site_freqs      - array with counts of '1' per site 
 *
 *  Returns an integer
 */
int num_singleton_sites(int segsites, int *site_freqs) {
    int     i,                  /* iterator */
            count;              /* running total */

    count = 0;

    for( i = 0; i < segsites; i++) {
        /* calculate site frequency of '1' */
        if (site_freqs[i] == 1) {
          count++;
        }
    }
    
    return (count);
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
double homozygosity(int nsam, int* hap_freqs) {
    int     i;                  /* iterator */
    double  total,              /* double value for the total number of samples */
            count,              /* temporary holder for each haplotype count */
            proportion,         /* temporary holder for each haplotype proportion */
            ho;                 /* calculated as we go (sum of per site homozygosities) */

    ho = 0.0;
    total = nsam;

    for (i=0; i < nsam; i++) {
        if ( hap_freqs[i] > 0 ) {
            count = hap_freqs[i];
            proportion = count/total;
            ho += proportion*proportion;
        }
    }
    
    return (ho);
}

/* Print help info. */
static void print_help (void) {
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
    -F        thetaH: Fay's H\n\
    -d        H:      difference between pi and Fay's H\n\
    -W        thetaW: Watterson's theta\n\
    -D        D:      Tajima's D\n\
    -H        ho:     Homozygosity\n\
    -n        nh:     number of haplotypes\n\
    -s        ns:     number of singletons\n\
    -N        nss:    number of singleton sites\n\
    -f        hf:     mean haplotype frequency\n\
    -i        ih:     max number identical haplotypes\n\
    -R        R2:     Ramos-Onsins & Rozas' R2\n\
    -U        Fs:     Fu's Fs\n", stdout);

  puts ("");
  fputs ("\
Calculate summary statistics from 'ms' output. With no options specified,\n\
the core set calculated by Hudson's original 'sample_stats' program will be\n\
written to the screen.  This original set is: pi, ss, D, thetaH, H\n\
\n\
This modified version of Hudson's original program has been optimized to\n\
calculate this original set about 40% faster.  In addition, it implements\n\
methods to calculate 4 more statistics: thetaW, ho, nh, and ns.\n\
\n\
Examples:\n\
      ms 10 1 -t 5 | sample_stats2\n\
      ms 10 1 -t 5 | sample_stats2 -SW\n\
      ms 10 1 -t 5 | sample_stats2 -pSFdWDHns", stdout);
  printf ("\n");
}

/* Print version and copyright information.  */
static void print_version (void) {
  printf ("(%s) version %s\n", PACKAGE, VERSION);
}

int main(int argc, char *argv[]) {
    int     i,                  /* iterator */
            nsam,               /* number of samples in the dataset */
            howmany;            /* number of replicates in the dataset */

    char    **list,             /* a matrix containing the data, 
                                 *   samples in rows, positions in columns*/
            line[1001],         /* temporary string to hold each line as it gets
                                 *   pulled from stdin */
            slashline[1001];    /* 'tbs' parameters are placed tab-delimited on a
                                 *   line beginning with "//". As the data are read in
                                 *   the text following any line starting with "//" is
                                 *   copied to the <slashline> and printed out following
                                 *   the summary statistics */

    int     segsites,           /* the number of sites for each replicate */
            count,              /* running tally of how many replicates have been processed */
            nh,                 /* the number of haplotypes */
            probflag;           /* 0 or 1, whether or not the input data includes 
                                 *   a "prob: ##" line*/

    double  pi,                 /* nucleotide diversity */
            th,                 /* Fay's theta H*/
            prob;               /* the prob value from the input */
    char    dum[20];            /* throwaway string, used when parsing first line of input file */
    
    int     *site_frequencies,  /* array holding the total count of '1's per site */
            *hap_frequencies,   /* array holding unique haplotypes counts */
            *unic_frequencies;  /* array holding count of unique sites per sequence */

    int     ss_flag,            /* 0 or 1; output the number of segregating sites */
            pi_flag,            /* 0 or 1; output nucleotide diversity */
            th_flag,            /* 0 or 1; output Fay's H (thetaH) */
            d_flag,             /* 0 or 1; output the difference between pi and Fay's H */
            tw_flag,            /* 0 or 1; output Watterson's theta (thetaW) */
            nh_flag,            /* 0 or 1; output the number of haplotypes */
            ns_flag,            /* 0 or 1; output the number of singleton haplotypes */
            ho_flag,            /* 0 or 1; output homozygosity */
            td_flag,            /* 0 or 1; output Tajima's D */
            nss_flag,           /* 0 or 1; output the number of singleton sites */
            hf_flag,            /* 0 or 1; output mean number of samples per haplotype */
            ih_flag,            /* 0 or 1; output max number of identical haplotypes */
            r2_flag,            /* 0 or 1; output Romas-Onsins & Rozas' R2 */
            fs_flag;            /* 0 or 1; output Fu's Fs */
    
    char    ch;                 /* current character iterator for getopt option parsing */
    int     remaining;          /* used to keep track of the number of remaining chars
                                 * when skipping through the 'positions' line */

    program_name = argv[0];

    /* by default, we print the same set as the original sample_stats */
    ss_flag = 1;
    pi_flag = 1;
    th_flag = 1;
     d_flag = 1;
    tw_flag = 0;
    nh_flag = 0;
    ns_flag = 0;
    ho_flag = 0;
    td_flag = 1;
    nss_flag = 0;
    hf_flag = 0;
    ih_flag = 0;
    r2_flag = 0;
    fs_flag = 0;

    /* However, if there are command line options, zero out all flags and
     * print only those specified in the options */
    if (argc > 1) {
      ss_flag = 0;
      pi_flag = 0;
      th_flag = 0;
       d_flag = 0;
      td_flag = 0;
    }

    /* Use getopt to parse the following flags:
     *      S - number of segregating sites
     *      p - pi
     *      F - Fay's H
     *      d - difference between pi and Fay's H (H in original sample_stats)
     *      W - Watterson's theta
     *      D - Tajima's D
     *      H - Homozygosity
     *      n - number of haplotypes
     *      s - number of singletons 
     *      N - number of singleton sites
     *      f - mean haplotype frequency
     *      i - max number identical haplotypes
     *      R - Ramos-Onsins & Rozas' R2
     *      U - Fu's Fs
     *      */
    while ((ch = getopt(argc, argv, "SpFdWDHnsNfiRUhv")) != -1) {
        switch (ch) {
	        case 'S':
		        ss_flag = 1;
		        break;
            case 'p':
                pi_flag = 1;
                break;
            case 'F':
                th_flag = 1;
                break;
            case 'd':
                d_flag = 1;
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
            case 'f':
                hf_flag = 1;
                break;
            case 'i':
                ih_flag = 1;
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

    /* read in first line of the ms output */
    fgets(line, 1000, stdin);
    /* the first line has the complete ms command that created this dataset
     * and is of the form "ms NSAMPLES NREPETITIONS [FLAGS]"
     * scan in the NSAMPLES and NREPETITIONS values, dropping "ms" via the
     * <dum> variable and ignoring everything else on this line */
    sscanf(line," %s  %d %d", dum,  &nsam, &howmany);

    /* pull off the second line (random number seeds) and throw it away */
    fgets(line, 1000, stdin);

    /* initialize the two dimensional char matrix <list> that will hold our data */
    list = create_list(nsam,maxsites+1);

    /* allocate enough space for the number of sites in the data */
    site_frequencies = (int *)malloc( maxsites*sizeof( int ));

    /* allocate enough space for the number of samples in the data */
    hap_frequencies = (int *)malloc( nsam*sizeof( int ));

    /* allocate enough space for the number of samples in the data */
    unic_frequencies = (int *)malloc( nsam*sizeof( int ));

    count=0;
    probflag=0;

    /* repeat as many times as we have replicates
     * (note that the post-increment of count happens after
     *  the sum has been calculated and the while test evaluated) */
    while( howmany - count++ ) {

        /* read in a sample */
        do {
            /* bail out if there's no data */
            if( fgets( line, 1000, stdin) == NULL ) {
                exit(0);
            }
            /* if this is the "//" line, then push the data into <slashline> */
            if( line[0] == '/' ) {
                strcpy(slashline,line+3);
            }
            /* otherwise, just read and throw away lines until we get to either a
             * "segsites: <...> " or a "prob: <...> " line */
        } while ( (line[0] != 's') && (line[0] != 'p' ) ) ;
 
        /* if we've hit the prob line, read it in. this line will only be present 
         * if both the "-s" and "-t" flags were used to generate the data in ms
         * (note that the ms documentation says that this will come after the
         *  "segsites: <...>" line, but in actuality it comes before).*/
        if( line[0] == 'p') {
            sscanf( line, "  prob: %lf", &prob );
            probflag = 1 ;
            /* bail out if the input ends */
            if( fgets( line, 1000, stdin) == NULL ) {
                exit(0);
            }
        }

        /* read in the number of segregating sites for this replicate */
        sscanf( line, "  segsites: %d", &segsites );

        /* increase the global maxsites variable if the current replicate has
         * more sites than we are currently prepared to deal with.  Also increase
         * the size of the data list to accomodate the larger amount of data about
         * to be read in. */
        if( segsites >= maxsites){
            maxsites = segsites + 10 ;
            site_frequencies = (int *)realloc( site_frequencies, maxsites*sizeof( int ));
            biggerlist(nsam,maxsites, list) ;
        }

        /* if this replicate has any segregating sites... */
        if( segsites > 0) {
            /* There is a line following segsites that looks like:
             *   positions: #.#### #.#### #.####
             * with as many numeric entries as there are segregating sites.
             * 
             * We're not using these data, so we want to pull this line off
             * and discard it.  However, as our line buffer is only 1000 characters,
             * anything over 142 segregating sites won't fit in one line read. 
             * So we need to make sure we get it all.
             * 
             * The 'remaining' variable is initialized with the total number
             * of characters on this line: 
             *   11 chars for 'position: '
             *   1  char  for the end of line '\0'
             *   7  chars for each segsite entry */
            remaining = segsites * 7 + 11 + 1;
            do {
                fgets( line, 1000, stdin);
                remaining = remaining - 1000;
            } while (remaining > 0);

            /* now pull off each line of data (each sample) into list */
            for( i=0; i<nsam;i++) fscanf(stdin," %s", list[i] );
        }

        /* fill in the site_frequencies array */
        for(i=0; i < segsites ; i++) {
            site_frequencies[i] = frequency('1', i, nsam, list);
        }

        /* fill in the unic_frequencies array if necessary */
        if (r2_flag) {
            count_binary_unic_frequencies(nsam, segsites, list, site_frequencies, unic_frequencies);
        }

        /* count up the haplotype frequencies if we are going to use them */
        if (nh_flag || ns_flag || ho_flag || hf_flag || fs_flag || r2_flag ) {
            count_haplotype_frequencies(nsam, segsites, list, hap_frequencies);
        }

        /* calculate pi if necessary */
        if ( pi_flag || d_flag || fs_flag || r2_flag )
            pi = theta_pi(nsam, segsites, site_frequencies);
    
        /* calculate Fay's H if necessary */
        if ( th_flag || d_flag )
            th = theta_h(nsam, segsites, site_frequencies);

        /* calculate the number of haplotypes if necessary */
        if (nh_flag || hf_flag || fs_flag)
            nh = num_haplotypes(nsam, hap_frequencies);
    
        if ( pi_flag )
            printf("pi:\t%lf\t", pi);
        if ( ss_flag )
            printf("ss:\t%d\t", segsites);
        if ( td_flag )
            printf("D:\t%lf\t", tajd(nsam,segsites,pi));
        if ( th_flag )
            printf("thetaH:\t%lf\t", th);
        if (  d_flag )
            printf("H:\t%lf\t", pi - th);
        if ( tw_flag )
            printf("thetaW:\t%lf\t", theta_w(nsam, segsites, site_frequencies));
        if ( nh_flag )
            printf("num_haplotypes:\t%d\t", nh);
        if ( ns_flag )
            printf("num_singletons:\t%d\t", num_singletons(nsam, hap_frequencies));
        if ( ho_flag )
            printf("homozygosity:\t%lf\t", homozygosity(nsam, hap_frequencies));
        if ( probflag )
            printf("prob:\t%g\t", prob);
        if ( nss_flag )
            printf("nss:\t%d\t", num_singleton_sites(segsites, site_frequencies));
        if ( hf_flag )
            printf("hf:\t%lf\t", (double)nsam/(double)nh);
        if ( ih_flag )
            printf("ih:\t%d\t", max_identical_haplotypes(nsam, hap_frequencies));
        if ( r2_flag )
            printf("r2:\t%lf\t", R2(unic_frequencies, pi, nsam, segsites));
        if ( fs_flag )
            printf("Fs:\t%lf\t", Fs(nsam, pi, nh));
        printf("%s\n", slashline);
       
    }
    
    exit (EXIT_SUCCESS);
}
