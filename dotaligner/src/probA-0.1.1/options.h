#ifndef _OPTIONS_H_
#include <getopt.h>
#define _OPTIONS_H_

extern int verbose_flag;


/*  the struct option structure has these fields: */
/*     const char *name                                                  */
/*     int has_arg    three legitimate values: no_argument(0)            */
/*  	                                       required_argument(1)      */
/*  	                                       optional_argument(2)      */
/*     int *flag      flag=0; => return val
                      otherwise => return 0; flag points to a variable which
		      is set to val if the option is found    */
/*     int val        value to return, or to load into  the  variable pointed
                      to by flag.*/

extern struct option long_options[];

extern int option_index;    /* getopt_long stores the option index here. */


extern char *shortopts;
extern sequ *input(int argc, char **argv);
extern void start(int argc, char **argv);

#endif












