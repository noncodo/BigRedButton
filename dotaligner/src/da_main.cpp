/*
 * main_v2.cpp
 *
 *  Created on: Jul 03, 2014
 *      Author: Stefan Seemann, seemann@rth.dk
 *
 *  Implements a heuristic for aligning to base pair probability matrices (dotplots)
 *
 *  Based on Needleman-Wunsch algorithm implementation in http://www.rolfmuertter.com/code/nw.php
 *
 */

#include "nwdp.h"
#include <getopt.h>
#include <map>
#include <ctime>
#include <cstdlib>
#include <algorithm>
extern "C" {
  #include "probA-0.1.1/pfgoto.h"      // definition for C library functions
}

using namespace std;

/* global arguments */
double kappa = 0.5;
double open = 0.3;
double ext  = 0.1;
double theta = 0.5;
double zeta = 0.;
int mutation_rates_flag = 0;
int free_endgaps_flag = 0;

float INFINITE = -1000;


int main( int argc, char ** argv )
{
		char * program = *argv ;
        vector<double> probDbl_1;
        vector<double> probDbl_2;
        vector<double> probSgl_1;
        vector<double> probSgl_2;
        vector<double> probDblUp_1;
        vector<double> probDblUp_2;
        int len_1 = 0, len_2 = 0;
        string name1, name2, seq_1, seq_2;
        int precision = 4;
        double pnull = 0.0005;
        int subopt = 10;
        double temp = 1;

        /* arguments */
        string  filename1, filename2;

        /* boolean flags */
        static int help_flag = 0;

    	/* parsing long options */
    	while(1)
    	{
    		static struct option long_options[] =
    		{
    			{"dotplot", required_argument, 0, 'd'},			// input dot plots
    			{"kappa", required_argument, 0, 'k'},			// weight of sequence / unpaired probability similarity
    			{"open", required_argument, 0, 'o'},		 	// affine gap costs = open + k * ext (k gaps)
    			{"ext", required_argument, 0, 'e'},			    // affine gap costs = open + k * ext (k gaps)
    			{"precision", required_argument, 0, 'p'},       // number of digits considered of base pair reliabilities
    			{"pnull", required_argument, 0, 'n'},			// minimal probability
    			{"subopt", required_argument, 0, 's'},			// number of suboptimal sequence/unpaired probability alignments examined
    			{"theta", required_argument, 0, 't'},			// weight of sequence similarity in sequence compared to unpaired probability similarity
    			{"zeta", required_argument, 0, 'z'},            // weight of unpaired probability similarity compared to paired probabilities to upstream nucleotides
    			{"temp", required_argument, 0, 'T'},            // measure of our interest in suboptimal alignments; analogous to the thermodynamic temperature
    			{"mutation-rates", no_argument, &mutation_rates_flag, 1},  // whether probabilities are multiplied with mutation rates
    			{"free-endgaps", no_argument, &free_endgaps_flag, 1},  // whether end gaps are free
    			{"help", no_argument, &help_flag, 1},
    			{"verbose", no_argument, &verbose_flag, 1},     // verbose mode
    			{0, 0, 0, 0}
    		};

    		/* print help text */
    		if( help_flag )
    			usage_dotaligner(program);

    		/* getopt_long stores the option index here. */
    		int option_index = 0;

    		int cmd = getopt_long(argc, argv, "d:k:o:e:p:n:s:t:z:T:", long_options, &option_index);

    		/* Detect the end of the options. */
    		if (cmd == -1)
    			break;

    		switch(cmd)
    		{
    			case 0:   break;
    			case 'd': if (!filename1.length())
    						  filename1 = optarg;
    			          else if(!filename2.length())
    						  filename2 = optarg;
    			          else
    			        	  usage_dotaligner(program);
    			          break;
    			case 'k': kappa = atof(optarg); break;
    			case 'o': open = atof(optarg); break;
    			case 'e': ext = atof(optarg); break;
    			case 'p': precision = atoi(optarg); break;
    			case 'n': pnull = atof(optarg); break;
    			case 's': subopt = atoi(optarg); break;
    			case 't': theta = atof(optarg); break;
    			case 'z': zeta = atof(optarg); break;
    			case 'T': temp = atof(optarg); BETA = 1.0/temp; break;
    			default:  abort();
    			/* no break */
    		}
    	}

    	/* set probA specific parameters */
    	typ_flag = 1;
    	Egap_flag = ( free_endgaps_flag ) ? 1 : 0;
    	dnamt.p_n.o = open;
    	dnamt.p_n.e = ext;
    	ENDGAP = 0;

    	/* read input file */
    	sequ *seq_ar = new sequ[2];
        if( filename1.length() > 0 )
        {
            ifstream inputfile1;
            inputfile1.open(filename1.c_str());
            if( !inputfile1.is_open() ) {
                 cerr << "Could not open file " << filename1 << endl;
                 exit(EXIT_FAILURE);
            }
            len_1 = readinput(inputfile1, name1, seq_1, probDbl_1);
            inputfile1.close();

            seq_ar[0].name = strcpy((char*)malloc(name1.length()+1), name1.c_str());
            seq_ar[0].seq = strcpy((char*)malloc(seq_1.length()+1), seq_1.c_str());
        }
        else {
        	usage_dotaligner(program);
        }
        if( filename2.length() > 0 )
        {
            ifstream inputfile2;
            inputfile2.open(filename2.c_str());
            if( !inputfile2.is_open() ) {
                 cerr << "Could not open file " << filename2 << endl;
                 exit(EXIT_FAILURE);
            }
            len_2 = readinput(inputfile2, name2, seq_2, probDbl_2);
            inputfile2.close();

            seq_ar[1].name = strcpy((char*)malloc(name2.length()+1), name2.c_str());
            seq_ar[1].seq = strcpy((char*)malloc(seq_2.length()+1), seq_2.c_str());
        }
        else {
        	usage_dotaligner(program);
        }

        /****************************************************************************************/
        /* STEP 1: Fold envelopment -- basepair probability matrices using RNAfold or RNAplfold */
        /****************************************************************************************/

        /* get unpaired probabilities */
        getunpaired(probDbl_1, len_1, probSgl_1);
        getunpaired(probDbl_2, len_2, probSgl_2);
       	/* get percentage of paired probabilities to upstream nucleotides */
       	getpaired_up(probDbl_1, probSgl_1, len_1, probDblUp_1 );
       	getpaired_up(probDbl_2, probSgl_2, len_2, probDblUp_2 );

        /* log-odds scores */
        getlogoddsDbl( probDbl_1, seq_1, len_1, pnull );
        getlogoddsDbl( probDbl_2, seq_2, len_2, pnull );
       	getlogoddsSgl( probSgl_1, seq_1, len_1, pnull );
       	getlogoddsSgl( probSgl_2, seq_2, len_2, pnull );

        /* reduce depth of probability matrices */
        reducematrix(probDbl_1, len_1*len_1, precision);
        reducematrix(probDbl_2, len_2*len_2, precision);
        reducematrix(probSgl_1, len_1, precision);
        reducematrix(probSgl_2, len_2, precision);
        reducematrix(probDblUp_1, len_1, precision);
        reducematrix(probDblUp_2, len_2, precision);

        #if DEBUG
            for( int i=0; i<len_1; i++ ) { for( int j=0; j<len_1; j++ ) cerr << probDbl_1.at( i*len_1+j ) << "\t"; cerr << endl; }; cerr << endl;
            for( int i=0; i<len_2; i++ ) { for( int j=0; j<len_2; j++ ) cerr << probDbl_2.at( i*len_2+j ) << "\t"; cerr << endl; }; cerr << endl;
        	for( int i=0; i<len_1; i++ ) cerr << probSgl_1[ i ] << " "; cerr << endl;
        	for( int i=0; i<len_2; i++ ) cerr << probSgl_2[ i ] << " "; cerr << endl << endl;
            if( zeta )
            {
            	for( int i=0; i<len_1; i++ ) cerr << probDblUp_1[ i ] << " "; cerr << endl;
            	for( int i=0; i<len_2; i++ ) cerr << probDblUp_2[ i ] << " "; cerr << endl << endl;
            }
        #endif


        double score_seq, score_dp, score_gap, tempsim, maxsim;
        aligm maxa;
		int len_subaln, len_submatch;
		int *idx_1_subaln = new int[len_1];
		for( int i=0; i<len_1; i++ ) idx_1_subaln[ i ] = i;
		int *idx_2_subaln = new int[len_2];
		for( int i=0; i<len_2; i++ ) idx_2_subaln[ i ] = i;

        /**********************************************************************************************************/
        /* STEP 2.1: Alignment envelopment -- calculate optimal pairwise sequence alignment -- adapted from probA */
		/**********************************************************************************************************/

        aligm aligg;
        if( len_1 >= len_2 )
        {
        	aligg = align_da( seq_ar, probSgl_1, probDblUp_1, probSgl_2, probDblUp_2 ); /* here add simbp scoring scheme; synchronize Egap_flag for end_gaps with flag <free-endgaps> */
        	len_submatch = convert(aligg.a, idx_1_subaln, idx_2_subaln);
        }
        else
        {
        	aligg = align_da( seq_ar, probSgl_2, probDblUp_2, probSgl_1, probDblUp_1 ); /* here add simbp scoring scheme; synchronize Egap_flag for end_gaps with flag <free-endgaps> */
        	len_submatch = convert(aligg.a, idx_2_subaln, idx_1_subaln);
        }
        len_subaln = strlen(aligg.a);
        if( Egap_flag )
        	len_subaln -= count_endgaps(idx_1_subaln, len_1, idx_2_subaln, len_2, len_submatch);
	    //printf(" %s\n %s\n %s\n Score = %5.4f\n", aligg.s0.seq, aligg.a, aligg.s1.seq, aligg.score);
        //score_seq = simbp( seq_1, probSgl_1, probDblUp_1, len_1, seq_2, probSgl_2, probDblUp_2, len_2, idx_1_subaln, idx_2_subaln, len_submatch, aligg.a );
        score_seq = aligg.score / len_subaln;
		score_dp = simdp( seq_1, probDbl_1, len_1, seq_2, probDbl_2, len_2, idx_1_subaln, idx_2_subaln, len_submatch );
		tempsim = maxsim = kappa * score_seq + (1 - kappa) * score_dp;
		maxa = copy_aligm(aligg);

		/**************************************************************************************************************************************/
		/* STEP 2.2: Alignment envelopment -- calculate pairwise stochastic sequence alignments and match probabilities -- adapted from probA */
		/**************************************************************************************************************************************/

		real **Mpr;
        if( len_1 >= len_2 )
        	Mpr = partition_f_da( aligg, probSgl_1, probDblUp_1, probSgl_2, probDblUp_2 );
        else
        	Mpr = partition_f_da( aligg, probSgl_2, probDblUp_2, probSgl_1, probDblUp_1 );

		aligm sta;
        if( len_1 >= len_2 )
        	sta = stoch_backtr_da( aligg, probSgl_1, probDblUp_1, probSgl_2, probDblUp_2 );
        else
        	sta = stoch_backtr_da( aligg, probSgl_2, probDblUp_2, probSgl_1, probDblUp_1 );
        maxa.prob = sta.prob;

        double si;
        if( verbose_flag )
		{
    		/* calculate sequence identity of optimal sequence alignment */
    		si = seqidentity( sta.a );

    		//cout << "Optimal sequence alignment:" << endl;
		    printf(" seq_score\tprobability\tnorm_seq_score\tnorm_str_score\tcombined_score\tlength\tsequence_identity\talignment\n");
			printf(" %.4f\t\t%.2e\t%.4f\t\t%.4f\t\t%.8f\t%i\t%2.2f\t%s\n", sta.score, sta.prob, score_seq, score_dp, tempsim, len_subaln, si, sta.a);
		}

        while( subopt )
        {
            if( len_1 >= len_2 )
            	sta = stoch_backtr_da( aligg, probSgl_1, probDblUp_1, probSgl_2, probDblUp_2 );
            else
            	sta = stoch_backtr_da( aligg, probSgl_2, probDblUp_2, probSgl_1, probDblUp_1 );

        	/* STEP 3: Dotplot alignment score */
            if( len_1 >= len_2 )
            	len_submatch = convert(sta.a, idx_1_subaln, idx_2_subaln);
            else
            	len_submatch = convert(sta.a, idx_2_subaln, idx_1_subaln);
            len_subaln = strlen(sta.a);
            if( Egap_flag )
            	len_subaln -= count_endgaps(idx_1_subaln, len_1, idx_2_subaln, len_2, len_submatch);
            //score_seq = simbp( seq_1, probSgl_1, probDblUp_1, len_1, seq_2, probSgl_2, probDblUp_2, len_2, idx_1_subaln, idx_2_subaln, len_submatch, sta.a );
            score_seq = sta.score / len_subaln;
    		score_dp = simdp( seq_1, probDbl_1, len_1, seq_2, probDbl_2, len_2, idx_1_subaln, idx_2_subaln, len_submatch );
    		tempsim = kappa * score_seq + (1 - kappa) * score_dp;
    		if( verbose_flag )
    		{
        		/* calculate sequence identity of optimal sequence alignment */
        		si = seqidentity( sta.a );

    	        //cout << "Sub-optimal sequence alignments:" << endl;
    			printf(" %.4f\t\t%.2e\t%.4f\t\t%.4f\t\t%.8f\t%i\t%2.2f\t%s\n", sta.score, sta.prob, score_seq, score_dp, tempsim, len_subaln, si, sta.a);
    		}

    		if( tempsim > maxsim )
    		{
    			maxsim = tempsim;
    			maxa = copy_aligm(sta);
    		}

            free_stoch(sta);

            subopt--;
        }

        /********************************************/
        /* STEP 3: return optimal dotplot alignment */
        /********************************************/

		if( verbose_flag )
		{
			cout << "Best dotplot alignment:" << endl;
			printf("seq_score = %.4f\n", maxa.score);
			printf("probability = %.2e\n", maxa.prob);
			printf("combined_score = %.8f\n", maxsim);
			printf("length = %i\n", (int) strlen(maxa.a));
		}
		else
		{
			printf("%s\t%s\t%.4f\t%.2e\t%.8f\t%i\t%4.2f", filename1.c_str(), filename2.c_str(), maxa.score,  maxa.prob, maxsim, (int) strlen(maxa.a), seqidentity(maxa.a));
		}
		decode_alig_da(maxa);
        free_aligm(maxa);


        /* free memory */
        free_partition_f_da(Mpr, seq_ar);
        free_align_da(aligg);
        for (int i=0; i<2; i++) {
            free(seq_ar[i].name);
            free(seq_ar[i].seq);
        }
        delete[](seq_ar);
        delete[] idx_1_subaln;
	    delete[] idx_2_subaln;
	    probDbl_1.clear();
        probDbl_2.clear();
        probSgl_1.clear();
        probSgl_2.clear();
        probDblUp_1.clear();
        probDblUp_2.clear();

        return 0;
}
