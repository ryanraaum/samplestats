#include <stdlib.h>
#include <assert.h>

#include "r2.h"

int main(int argc, char *argv[]) {
  char *list[2];
  int binary_site_freqs[10] = {0,1,0,0,1,0,1,0,0,0};
  int agct_site_freqs[4][4];

  int unic_freqs[2];

  list[0] = "0000001000\0";
  list[1] = "0100100000\0";
  count_binary_unic_frequencies(2, 10, list, binary_site_freqs, unic_freqs);
  assert(1 == unic_freqs[0]);
  assert(2 == unic_freqs[1]);

  list[0] = "AAAG\0";
  list[1] = "AGAA\0";

  agct_site_freqs[0][0] = 2;
  agct_site_freqs[0][1] = 0;
  agct_site_freqs[0][2] = 0;
  agct_site_freqs[0][3] = 0;
  agct_site_freqs[1][0] = 1;
  agct_site_freqs[1][1] = 1;
  agct_site_freqs[1][2] = 0;
  agct_site_freqs[1][3] = 0;
  agct_site_freqs[2][0] = 2;
  agct_site_freqs[2][1] = 0;
  agct_site_freqs[2][2] = 0;
  agct_site_freqs[2][3] = 0;
  agct_site_freqs[3][0] = 1;
  agct_site_freqs[3][1] = 1;
  agct_site_freqs[3][2] = 0;
  agct_site_freqs[3][3] = 0;

  count_agct_unic_frequencies(2, 4, list, agct_site_freqs, unic_freqs);

  assert(2 == unic_freqs[0]);
  assert(2 == unic_freqs[1]);

  exit(0);
}

