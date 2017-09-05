/*
 * nw.h for program nw.
 *
 */

#include <iostream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>

#define DEBUG 0

using namespace std;


extern float  nw( 
                 string, string, 
                 string&, string&,
                 bool
              );

extern int  nw_align( 
                      int **, char **,
                      string, string, 
                      string&, string&,
                      int 
                    );

extern void  dpm_init        ( int **, char **, int, int, int );
extern void  print_al        ( string&, string& );
extern void  print_matrix    ( int ** const, string, string );
extern void  print_traceback ( char ** const, string, string );
extern int   max             ( int, int, int, char * );
extern int   readinput       ( istream & is, string & name, string & seq );
extern void  usage           ( char *program );

