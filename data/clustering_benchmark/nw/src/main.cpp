/*---------------------------------------------------------------
 *
 *   main.c++ for nw program.
 *
 *   Implements the Needleman-Wunsch algorithm
 *   for global nucleotide sequence alignment.
 *
 *   Rolf Muertter,  rolf@dslextreme.com
 *   9/5/2006
 *
 ---------------------------------------------------------------*/


#include "nw.h"
#include <cstdlib>
#include <fstream>

using namespace std;


int  main( int argc, char ** argv )
{
        char *  program = *argv ;
        bool    prm = false;

        if( argc < 3 )
		usage(program);
			        
        // Sequences to be aligned
        string  file_1   =  argv[ 1 ] ;
        string  file_2   =  argv[ 2 ] ;

        if( argc == 4 )
        {
                string  pr_arg  =  argv[ 3 ] ;
                if( pr_arg == "-p" )  prm = true;   // Print matrices
        }                       

   	/* read input file */
	int len_1, len_2;
	string seq_1, seq_2;
	string name1, name2;
        if( file_1.length() > 0 )
        {
            ifstream inputfile1;
            inputfile1.open(file_1.c_str());
            if( !inputfile1.is_open() ) {
                 cerr << "Could not open file " << file_1 << endl;
                 exit(EXIT_FAILURE);
            }
            len_1 = readinput(inputfile1, name1, seq_1);
            inputfile1.close();
        }
        else {
        	usage(program);
        }
        if( file_2.length() > 0 )
        {
            ifstream inputfile2;
            inputfile2.open(file_2.c_str());
            if( !inputfile2.is_open() ) {
                 cerr << "Could not open file " << file_2 << endl;
                 exit(EXIT_FAILURE);
            }
            len_2 = readinput(inputfile2, name2, seq_2);
            inputfile2.close();
        }
        else {
        	usage(program);
        }

        #if DEBUG
            cout << "seq_1: " << seq_1 << endl;
            cout << "seq_2: " << seq_2 << endl;
            cout << "-p: " << prm << endl;
        #endif

         // Aligned sequences
        string  seq_1_al;
        string  seq_2_al;

        // Get alignment
        float score = nw( seq_1, seq_2, seq_1_al, seq_2_al, prm ) ;   

        // Print alignment
	cout << "Score = " << score << endl;
        print_al( seq_1_al, seq_2_al ) ;        

        return  0 ;
}

