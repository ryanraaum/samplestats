#include "simple_getopt.h"

char        *optarg;            /* global pointer to option parameter */
int         optopt;             /* global pointer to unknown option char */
int         optind = 0;         /* global argv index */
int         opterr = 1;         /* global report errors flag */

int getopt(int argc, char *argv[], char *optstring)
{
  char        c;                /* current option char */
  char        *cp;              /* pointer to option char match in optstring */
  static char *next = NULL;     /* pointer to next char in current option */

  optarg = NULL;
  optopt = 0;

  /* If we're just starting (next == NULL) or 
   * we've reached the end of the current string (*next == '\0'),
   * then advance the global option index (optind)
   * and look at what we've got. 
   *
   * Return -1 if we find any problems, specifically:
   *  - no options to be found
   *  - option doesn't start with '-'
   *  - option consists only of '-'
   *
   * Stop scanning if we hit '--'
   *
   * Stop scanning and throw an error if we hit a long style 
   * option like '--something'
   *
   * If we don't find any problems, then:
   *  - point the 'next' pointer at the current argv entry and 
   *    skip it past the opening '-'
   *  - increment the global optind for the next round
   *  - proceed to process the current option. */
  if (next == NULL || *next == '\0') {
    /* if we're just starting, move into the first option */
    if (optind == 0)
      optind++;

    /* here we discover if:
     *  - our index is greater than the number of arguments (optind >= argc)
     *  - we've hit an improperly formated option:
     *     - that doesn't start with '-' (argv[optind][0] != '-')
     *     - that is _only_ '-' (argv[optind][1] == '\0') */
    if (optind >= argc || argv[optind][0] != '-' || argv[optind][1] == '\0') {
      optarg = NULL;
      if (optind < argc)
        optarg = argv[optind];
      return -1;
    }
    
    /* Stop scanning if we hit '--', increment optind so it points to the first
     * argv entry (if any) after '--' */
    if (strncmp(argv[optind], "--\0", 3) == 0) {
      optind++;
      return -1;
    }

    /* Stop scanning if we hit long-style option like '--something' 
     * If the user does not check for a non-null optarg containing
     * a string that starts with '--' then this functions like '--',
     * which is, stop scanning and point optind to the next argv entry. */
    if (strncmp(argv[optind], "--", 2) == 0) {
      optarg = argv[optind];
      optind++;
      return -1;
    }

    /* if we make it here, we've found no problems,
     * so point 'next' at the current option. */
    next = argv[optind];
    /* and skip past the opening '-' character */
    next++;

    optind++;
  }

  /* c gets next's value, next moves along */
  c = *next;
  next++; 

  /* does c match anything in the optstring? 
   * (if so, cp will point to the first match) */
  cp = strchr(optstring, c);

  /* if there's no match then return the current char. */
  if (cp == NULL) {
    optopt = c;
    if (opterr) {
      fprintf (stderr, "Unknown option `-%c'.\n", optopt);
    }
    return '?';
  }
  
  /* if c is ':', that implies that the option string starts
   * with '-:', which is an unusual option. Return c.
   *
   * (':' would be an unusual option because ':' is what is used in 
   *  the optstring to indicate that the previous option takes a parameter.) */
  if (c == ':')
    return c;

  /* Handle option parameters. Advance cp to see if the next character
   * in the optstring is ':'.  If so, try to get a parameter, which is:
   *  - the rest of the option "word" following the parameter-requiring
   *    option character (if any)
   *  - otherwise, the next "word" in the argv list
   *  - if no characters follow the parameter-requiring option,
   *    and there are no more option "words" in argv, simply
   *    default to returning the option character. However, optarg
   *    (which should hold the option parameter) will be NULL, and
   *    the option handling code in the main program should be
   *    prepared for that indication of a missing parameter. */
  cp++;
  if (*cp == ':') {
    if (*next != '\0') {
      optarg = next;
      next = NULL;
    } else if (optind < argc) {
      optarg = argv[optind];
      optind++;
    } else {
      return '?';
    }
  }

  return c;
}
