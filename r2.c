#include <stdlib.h>
#include <math.h>

/*  Count up the per sample unique site frequencies in the data;
 *    That is, how many sites are unique to each sequence.
 *
 *      nsam            - total number of samples in data list
 *      segsites        - total number of segregating sites 
 *                        ( length of positions array )
 *      list            - the data ( samples by positions matrix of chars )
 *      hap_freqs       - the (initialized) array of integers to fill (length nsam)
 *
 *  Returns nothing (fills in the array given)
 */
void count_binary_unic_frequencies(int nsam, int segsites, char **list, int *site_freqs, int *unic_freqs) 
{
    int     i, j;               /* iterators */

    /* first initialize all unic counts to 0 */
    for (i=0; i < nsam; i++) {
        unic_freqs[i] = 0;
    }

    /* If there are no segregating sites, then there are no unique sites. 
     * So only do the work if there are some segregating sites. */
    if (segsites > 0) {
        /* step through the site frequency spectrum, checking all unique sites */
        for (i=0; i<segsites; i++) {
            /* every time we find a site with frequency 1, find the sequence that is
             * the source of that site and up its' unic_freqs count. */
            if (site_freqs[i] == 1) {
                for (j=0; j<nsam; j++) {
                    if (list[j][i] == '1') {
                        unic_freqs[j] += 1;
                    }
                }
            }
        }
    }
}

/*  Count up the per sample unique site frequencies in the data;
 *    That is, how many sites are unique to each sequence.
 *
 *      nsam            - total number of samples in data list
 *      segsites        - total number of segregating sites 
 *                        ( length of positions array )
 *      list            - the data ( samples by positions matrix of chars )
 *      site_freqs      - array with counts of nucleotides per site 
 *      unic_freqs      - the (initialized) array of integers to fill (length nsam)
 *
 *  Returns nothing (fills in the array given)
 */
void count_agct_unic_frequencies(int nsam, int segsites, char **list, int site_freqs[][4], int *unic_freqs) 
{
    int     i, j, k;               /* iterators */
    char    *agct;

    agct = "AGCT\0";

    /* first initialize all unic counts to 0 */
    for (i=0; i<nsam; i++) {
        unic_freqs[i] = 0;
    }

    /* If there are no segregating sites, then there are no unique sites. 
     * So only do the work if there are some segregating sites. */
    if (segsites > 0) {
        /* step through the site frequency spectrum, checking all unique sites */
        for (i=0; i<segsites; i++) {
            for (j=0; j<4; j++) {
                /* every time we find a site with frequency 1, find the sequence that is
                 * the source of that site and up its' unic_freqs count. */
                if (site_freqs[i][j] == 1) {
                    for (k=0; k<nsam; k++) {
                        if (list[k][i] == agct[j]) {
                            unic_freqs[k] += 1;
                        }
                    }
                }
            }
        }
    }
}

/* Calculate Ramos-Onsins & Rozas' R2
 *
 *      pi              - average number of nucleotide differences
 *      segsites        - number of segregating sites 
 *      nsam            - number of samples
 *      unic_freqs      - array with number of unique sites per sample
 *
 *  Returns double R2 statistic
 */
double R2(int *unic_freqs, double pi, int nsam, int segsites)
{
    double  sm2;
    int     i;

    sm2 = 0.0;
    
    if(segsites == 0 || nsam == 0) 
        return(-10000);

    for (i=0; i<nsam; i++)
        sm2 += ((double)unic_freqs[i] - pi/2.0)*((double)unic_freqs[i] - pi/2.0);
    
    sm2 = sqrt(sm2/((double)nsam))/(double)segsites;
            
    if (sm2 < 1.0E-15)
        sm2 = 0.0;

    return (double)sm2;
}


