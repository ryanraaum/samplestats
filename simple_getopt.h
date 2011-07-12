#ifndef SIMPLE_GETOPT_H
#define SIMPLE_GETOPT_H

#include <stdio.h>
#include <string.h>

extern int optind, opterr;
extern char *optarg;
extern int optopt;
extern int opterr;

int getopt(int argc, char *argv[], char *optstring);

#endif //SIMPLE_GETOPT_H

