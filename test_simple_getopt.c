#include "test_simple_getopt.h"

/* String containing name the program is called with. */
const char *program_name;

static void print_help (void);
static void print_version (void);

int main (int argc, char *argv[])
{
	char ch;
	
	program_name = argv[0];

	if (argc == 1) {
		print_help();
	}

	while ((ch = getopt(argc, argv, "bhvf:")) != -1) {
		switch (ch) {
			case 'b':
				printf ("You specified the -b option.\n");
				break;
			case 'f':
				printf ("You specified the -f option ");
				if (optarg == NULL) {
				  printf("\nExiting: 'f' option requires a parameter. \nTry '%s -h' for help\n", program_name);
          exit (EXIT_FAILURE);
				}
				printf ("with the parameter \"%s\".\n", optarg);
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
	
	if (optarg != NULL && strncmp(optarg, "--", 2) == 0) {
    printf("Exiting. This program doesn't take long options like '%s'\nTry '%s -h' for help\n", optarg, program_name);
    exit (EXIT_FAILURE);
  }
  
  if (optind < argc) {
    for ( ; optind < argc; optind++)
      printf("extra option: %s\n", argv[optind]);
  }

	exit (EXIT_SUCCESS);
}

/* Print help info. */
static void print_help (void) {
  printf ("Usage: %s [OPTIONS] arguments ...\n", program_name);
  fputs ("Test the option parsing system.\n", stdout);

  puts ("");

  fputs ("\
  -h                  display this help and exit\n\
  -v                  display version information and exit\n", stdout);

  puts ("");

  fputs ("\
  -b                  capture a simple, no-parameter option 'b'\n\
  -f VALUE            capture a valued parameter\n", stdout);

  printf ("\n");
}



/* Print version and copyright information.  */

static void print_version (void) {
  printf ("(%s) version %s\n", PACKAGE, VERSION);
  puts ("");
}
