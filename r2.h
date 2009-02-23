void count_binary_unic_frequencies(int nsam, int segsites, char **list, int *site_freqs, int *unic_freqs);
void count_agct_unic_frequencies(int nsam, int segsites, char **list, int site_freqs[][4], int *unic_freqs);
double R2(int *unic_freqs, double pi, int nsam, int segsites);
