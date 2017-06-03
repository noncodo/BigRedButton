/*
 * nwdp_v2.h
 *
 *  Created on: Jul 03, 2014
 *      Author: Stefan Seemann, seemann@rth.dk
 */

#ifndef NWDP_H_
#define NWDP_H_

#include <iostream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
//#include "boost/assign.hpp"
extern "C" {
  #include "probA-0.1.1/pfgoto.h"      // definition for C library functions
}

#define DEBUG 0

#define MAX(a,b) (( (a) > (b) ) ? (a) : (b) )
#define EQUAL(A,B) (fabs((A)-(B)) < 100*FLT_EPSILON)
/* encode the alignment: the longer sequence is always written above the
   shorter one, each state of the alignment (match, mismatch, indel) is
   symbolized a self-evident char.*/
#define M1 '|' /*  match */
#define M2 ':' /*  mismatch */
#define G3 '.' /*  gap in the lower sequence (the shorter one)*/
#define G4 '^' /*  gap in the upper sequence (the longer one)*/


using namespace std;
//using namespace boost::assign;

/* arguments */
extern double kappa;
extern double theta;
extern double zeta;
extern double open;
extern double ext;

extern int mutation_rates_flag;
extern int unpaired_up_down_flag;
extern int free_endgaps_flag;

/* adaptation from probA-file pfgoto.c for DotAligner */
static sequ *so; /* array with 2 members (longer seq + name is in so[0]) */
static iseq *si; /* like '*so' but sequences are encoded as strings of digits
		    see function code sequence */
static u_sc scmat; /* scoring matrix */
static matrix pmat,rmat; /* partition function and revers pf */
static real **Pr; /* match probability matrix */
static int sran=0; /* for initialisation of of the "Mersenne Twister"*/


class MapClass {
public:
	typedef map<char, int> CharIntMap;
	typedef map<string, int> StringIntMap;

	static CharIntMap init_base_index() {
		CharIntMap base_idx;
		base_idx['A'] = 0;
		base_idx['U'] = 1;
		base_idx['G'] = 2;
		base_idx['C'] = 3;
		return base_idx;
	}
	static StringIntMap init_pair_index() {
		StringIntMap pair_idx;
		pair_idx["AA"] = 0;
		pair_idx["AU"] = 1;
		pair_idx["AG"] = 2;
		pair_idx["AC"] = 3;
		pair_idx["UA"] = 4;
		pair_idx["UU"] = 5;
		pair_idx["UG"] = 6;
		pair_idx["UC"] = 7;
		pair_idx["GA"] = 8;
		pair_idx["GU"] = 9;
		pair_idx["GG"] = 10;
		pair_idx["GC"] = 11;
		pair_idx["CA"] = 12;
		pair_idx["CU"] = 13;
		pair_idx["CG"] = 14;
		pair_idx["CC"] = 15;
		return pair_idx;
	}
};

//MyClass::OpMap MyClass::bpindex = init_map();

/* base frequencies */
/* A, U, G, C */
const double basefreq[ 4 ] = { 0.364097, 0.273013, 0.211881, 0.151009 };

/* base pair frequencies */
/* A, U, G, C */
const double pairfreq[ 4 ][ 4 ] = { { 0.001167, 0.177977, 0.001058, 0.001806 },
									{ 0.177977, 0.002793, 0.049043, 0.000763 },
									{ 0.001058, 0.049043, 0.000406, 0.266974 },
									{ 0.001806, 0.000763, 0.266974, 0.000391 } };

/* base substitution matrix */
const int a = 1;   /* Match */
const int b = 0;   /* Mismatch */
/* A, U, G, C */
const int basesub[ 4 ][ 4 ] = { { a, b, b, b },
                                  { b, a, b, b },
                                  { b, b, a, b },
                                  { b, b, b, a } } ;
/*
 * Robert J Klein and Sean R Eddy: RSEARCH: Finding homologs of single structured RNA sequences. BMC Bioinformatics (2003)
 * When only a single query sequence is given, log-odds position independent substitution matrices are used to give the alignment scores.
 * Log-odds RIBOSUM matrix gives the log-odds ratio for observing a given substitution relative to background nucleotide frequencies.
 * Log-odds scores in RIBOSUM85-60 in Figure 3 are transformed to probabilities: p = 2^logodds / (1 + 2^logodds)
 */
const double basemutationrates[ 4 ][ 4 ] = { { 0.823288, 0.276183, 0.266590, 0.215979 },
											{ 0.276183, 0.758357, 0.230396, 0.325677 },
											{ 0.266590, 0.230396, 0.671272, 0.151999 },
											{ 0.215979, 0.325677, 0.151999, 0.690840 } };

/* base pair substitution matrix */
const int c = 1;   /* Match */
const int d = 0;   /* Mismatch */
/* AA, AU, AG, AC, UA, UU, UG, UC, GA, GU, GG, GC, CA, CU, CG, CC */
const int pairsub[ 16 ][ 16 ] = { { d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d },
									{ d, c, d, d, c, d, c, d, d, c, d, c, d, d, c, d },
									{ d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d },
									{ d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d },
									{ d, c, d, d, c, d, c, d, d, c, d, c, d, d, c, d },
									{ d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d },
									{ d, c, d, d, c, d, c, d, d, c, d, c, d, d, c, d },
									{ d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d },
									{ d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d },
									{ d, c, d, d, c, d, c, d, d, c, d, c, d, d, c, d },
									{ d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d },
									{ d, c, d, d, c, d, c, d, d, c, d, c, d, d, c, d },
									{ d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d },
									{ d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d },
									{ d, c, d, d, c, d, c, d, d, c, d, c, d, d, c, d },
									{ d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d } } ;
/*
 * Robert J Klein and Sean R Eddy: RSEARCH: Finding homologs of single structured RNA sequences. BMC Bioinformatics (2003)
 * When only a single query sequence is given, log-odds position independent substitution matrices are used to give the alignment scores.
 * Log-odds RIBOSUM matrix gives the log-odds ratio for observing a given substitution relative to background nucleotide frequencies.
 * Log-odds scores in RIBOSUM85-60 in Figure 3 are transformed to probabilities: p = 2^logodds / (1 + 2^logodds)
 */
//a = -8.4100
//printf("%.6f\n", 2^a / ( 1 + 2^a ))
const double pairmutationrates[ 16 ][ 16 ] = {
		{ 0.151108, 0.047680, 0.003297, 0.007542, 0.058441, 0.001883, 0.013792, 0.000391, 0.008535, 0.017158, 0.002972, 0.029698, 0.002177, 0.000157, 0.037546, 0.000047 },
		{ 0.047680, 0.957394, 0.027764, 0.195601, 0.753240, 0.112488, 0.412533, 0.034422, 0.024257, 0.600838, 0.020064, 0.866631, 0.020757, 0.027026, 0.760888, 0.009462 },
		{ 0.003297, 0.027764, 0.364817, 0.002103, 0.022670, 0.000465, 0.016027, 0.002133, 0.002553, 0.010203, 0.023451, 0.017996, 0.000734, 0.000885, 0.040400, 0.000042 },
		{ 0.007542, 0.195601, 0.002103, 0.188079, 0.024257, 0.004376, 0.008134, 0.002449, 0.001176, 0.036557, 0.000471, 0.066553, 0.001509, 0.000714, 0.016926, 0.001844 },
		{ 0.058441, 0.753240, 0.022670, 0.024257, 0.969080, 0.087084, 0.687872, 0.112488, 0.017042, 0.402492, 0.018242, 0.751949, 0.157438, 0.036557, 0.870587, 0.008419 },
		{ 0.001883, 0.112488, 0.000465, 0.004376, 0.087084, 0.503466, 0.048955, 0.015703, 0.000336, 0.043751, 0.000564, 0.062780, 0.002972, 0.026306, 0.072841, 0.022979 },
		{ 0.013792, 0.412533, 0.016027, 0.008134, 0.687872, 0.048955, 0.911248, 0.035592, 0.005308, 0.190205, 0.049278, 0.486141, 0.019793, 0.006050, 0.714013, 0.002931 },
		{ 0.000391, 0.034422, 0.002133, 0.002449, 0.112488, 0.015703, 0.035592, 0.097528, 0.009996, 0.024755, 0.000242, 0.042606, 0.007337, 0.065697, 0.032191, 0.005886 },
		{ 0.008535, 0.024257, 0.002553, 0.001176, 0.017042, 0.000336, 0.005308, 0.009996, 0.325677, 0.014369, 0.002449, 0.032845, 0.003945, 0.004753, 0.015385, 0.000181 },
		{ 0.017158, 0.600838, 0.010203, 0.036557, 0.402492, 0.043751, 0.190205, 0.024755, 0.014369, 0.917225, 0.017996, 0.698193, 0.004084, 0.020064, 0.453349, 0.004227 },
		{ 0.002972, 0.020064, 0.023451, 0.000471, 0.018242, 0.000564, 0.049278, 0.000242, 0.002449, 0.017996, 0.202227, 0.054029, 0.000380, 0.000076, 0.038050, 0.000163 },
		{ 0.029698, 0.866631, 0.017996, 0.066553, 0.751949, 0.062780, 0.486141, 0.042606, 0.032845, 0.698193, 0.054029, 0.980072, 0.015919, 0.017158, 0.811921, 0.071449 },
		{ 0.002177, 0.020757, 0.000734, 0.001509, 0.157438, 0.002972, 0.019793, 0.007337, 0.003945, 0.004084, 0.000380, 0.015919, 0.027764, 0.002774, 0.077663, 0.000714 },
		{ 0.000157, 0.027026, 0.000885, 0.000714, 0.036557, 0.026306, 0.006050, 0.065697, 0.004753, 0.020064, 0.000076, 0.017158, 0.002774, 0.170742, 0.031128, 0.017996 },
		{ 0.037546, 0.760888, 0.040400, 0.016926, 0.870587, 0.072841, 0.714013, 0.032191, 0.015385, 0.453349, 0.038050, 0.811921, 0.077663, 0.031128, 0.976230, 0.018746 },
		{ 0.000047, 0.009462, 0.000042, 0.001844, 0.008419, 0.022979, 0.002931, 0.005886, 0.000181, 0.004227, 0.000163, 0.071449, 0.000714, 0.017996, 0.018746, 0.076676 } };

extern float INFINITE;

typedef struct tagLocalHit {
	double similarity;		/* Similarity */
	int lstart_1;			/* Start position of LS_a */
	int lend_1;				/* End position of LS_a */
	int lstart_2;			/* Start position of LS_b */
	int lend_2;				/* End position of LS_b */
} LocalHit;

extern int readinput( istream & is, string & name, string & seq, vector<double> & prob );
extern void getunpaired( vector<double> & prob, int len, vector<double> & probSgl );
extern void getpaired_up( vector<double> & probDbl, vector<double> & probSgl, int len, vector<double> & probDblUp );
extern void getlogoddsDbl( vector<double> & probDbl, string seq, int len, double pnull );
extern void getlogoddsSgl( vector<double> & probSgl, string seq, int len, double pnull );
extern void reducematrix(vector<double> & prob, int len, int prec );
vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);

extern double nwdp( string seq_1, vector<double> & probSgl_1, int * idx_1_aln, int len_1, string seq_2, vector<double> & probSgl_2, int * idx_2_aln, int len_2, int & len_aln, double ** F, double ** Q, double ** P, char ** trF, char ** trQ, char ** trP );
extern double simbp( string seq_1, vector<double> & probSgl_1, vector<double> & probDblUp_1, int len_1, string seq_2, vector<double> & probSgl_2, vector<double> & probDblUp_2, int len_2, int * idx_1_aln, int * idx_2_aln, int len_match, char * aln_code );
extern void prob_backtracking( int * idx_1_subaln, int len_1, int * idx_2_subaln, int len_2, int & len_subaln, double ** F );
extern double simdp( string seq_1, vector<double> & probDbl_1, int len_1, string seq_2, vector<double> & probDbl_2, int len_2, int * idx_1_aln, int * idx_2_aln, int len_aln );
extern double gappenalty( char * aln_code );
extern int convert(char * aliggAlignment, int * idx_1_aln, int * idx_2_aln);
extern int count_endgaps(int * idx_1_aln, int len_1, int * idx_2_aln, int len_2, int len_match);
extern double seqidentity(char * aln_code );

extern void printalign(string & seq_1, int * idx_1_aln, string & seq_2, int * idx_2_aln, int len_aln );
extern void freeMatrix(double ** matrix, int leny);
extern void usage_dotaligner(char * program);

extern void nwdp_initF( double ** F, int L1, int L2 );
extern void nwdp_initF_affinegaps( double ** F, int L1, int L2, bool local );
extern void nwdp_initGap( double ** Q, int L1, int L2 );
extern double max3( double f1, double f2, double f3, char* ptr );
extern double max( double f1, double f2 );
template <typename T>
extern void print_matrixdp( T ** F, double * prob_1, int len_1, double * prob_2, int len_2 );
extern void nwdp_initTB( char ** traceback, int L1, int L2 );
extern void reverse( int * list, int len );
extern int statistical_sampling( unsigned int *sam_prob, int sam_len );

extern aligm align_da( sequ *seq_array, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short );
extern aligm align_2_da( sequ *os, iseq *is, u_sc smat, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short );
extern aligm trace_back_da( matrix m, sequ *os, iseq *is, u_sc smat, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short );
extern real **partition_f_da( aligm alig, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short );
extern matrix partf_da( aligm a, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short );
extern matrix revers_partf_da(sequ *os, iseq *is, u_sc m, aligm a, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short );
extern aligm stoch_backtr_da( aligm alig, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short );
extern aligm stoch_btr_da( matrix pmat, aligm alig, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short );
extern double calc_score_da( char *tr, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short );
extern void free_partition_f_da(real **m, sequ *s);
extern void free_align_da(aligm alig);
extern void decode_alig_da(aligm ai);


#endif /* NWDP_H_ */
