/*
 * nwdp.cpp
 *
 *  Created on: Mar 08, 2015
 *      Author: Stefan Seemann, seemann@rth.dk
 */

#include "nwdp.h"
#include <math.h>
extern "C" {
  #include "probA-0.1.1/pfgoto.h"      // definition for C library functions
}


using namespace std;
int stochswi=0;


void usage_dotaligner(char *program)
{
	cerr << "\nDotAligner v0.3";
	cerr << "\n===============";
	cerr << "\n   by Stefan E Seemann (seemann@rth.dk)\n";
	cerr << "\n   Usage:   " << program << " -d <file> -d <file> [ options ]";
    cerr << "\n   Semi-local pairwise alignment of two base pair probability matrices (dotplots).\n";
    cerr << "\n   -d --dotplot <file>    ... dotplot file (matrix of base pair probabilities)";
    cerr << "\n   Debug modes:";
    cerr << "\n   -k --kappa <double>    ... weight of string similarity compared to dotplot similarity (0..1; default = 0.5)";
    cerr << "\n   -t --theta <double>    ... weight of sequence/basepair similarity compared to similarity of unpaired/pairing probability (0..1; default = 0.5)";
    cerr << "\n                              only active if flag <mutation-rates> is not set";
    cerr << "\n   -o --gapopen <double>  ... gap-open penalty: affine gap costs = gapopen + (k-1) * gapext (k gaps; gapopen default = 0.3)";
    cerr << "\n                              set <gapopen> equal <gapext> to make the score becomes a linear function of gap length";
    cerr << "\n   -e --gapext <double>   ... gap-extension penalty: affine gap costs = gapopen + (k-1) * gapext (k gaps; gapext default = 0.1)";
    cerr << "\n                              set <gapext> equal 0 to keep the score similar regardless of gap length";
    cerr << "\n   -s --subopt <int>      ... number of examined suboptimal alignments (default = 10)";
    cerr << "\n   -p --precision <int>   ... number of digits considered of log-odds of base pair reliabilities (default = 4)";
    cerr << "\n   -n --pnull <prob>      ... minimal probability (default = 0.0005 )";
    cerr << "\n   --mutation-rates       ... usage of statistical substitution models RIBOSUM85-60";
    cerr << "\n   -z --zeta              ... weight of unpaired probability compared to upstream pairing probability (0..1; default = 0)";
    cerr << "\n                              only active if flag <mutation-rates> is set";
    cerr << "\n   -T --temp <double>     ... measure of our interest in suboptimal alignments; analogous to thermodynamic temperature (0..undef; default = 1)";
    cerr << "\n   --free-endgaps         ... whether start and end gaps are free";
    cerr << "\n   --verbose              ... verbose mode";
	cerr << "\n   --help                 ... this output\n";
    cerr << "\n   Output:   Similarity and alignment\n\n";

    exit( 1 ) ;
}


int readinput ( istream & is, string & name, string & seq, vector<double> & prob )
{
	// read probability matrix
	string line;
	int lnr = 0;
	int i = -1;

	while( !is.eof() ) {
		getline(is, line);

		if( i == -1 )
			name = line.replace(0,1,"");
		else
			if( i == 0 ) {
				seq = line;
				replace( seq.begin(), seq.end(), 'T', 'U' );
				lnr = line.length();
				prob.reserve( lnr * lnr );
			}
			else {
				vector<string> probv = split(line, ' ');
				for( unsigned int k = 0; k < probv.size(); k++ )
					prob.push_back( atof( probv[k].c_str() ) );
			}

		if( i == lnr )
			break;
		i++;
	}

	return lnr;
}


void getunpaired( vector<double> & prob, int len, vector<double> & probSgl )
{
	double psingle;

	probSgl.reserve( len );

	for( int i = 0; i < len; i++ )
	{
		psingle = 1;
		for( int j = 0; j < len; j++ )
			psingle -= prob.at( i*len + j );
		probSgl.push_back( psingle );
	}
}


void getpaired_up( vector<double> & probDbl, vector<double> & probSgl, int len, vector<double> & probDblUp )
{
	double pdouble;

	probDblUp.reserve( len );

	probDblUp.push_back( 0. );
	for( int i = 1; i < len; i++ )
	{
		pdouble = 0;
		for( int j = 0; j < i; j++ )
			pdouble += probDbl.at( i*len + j );
		if( probSgl[ i ] == 1 )
			probDblUp.push_back(0);
		else
			probDblUp.push_back( pdouble / ( 1 - probSgl[ i ] ) );
	}
}


void getunpaired_up( vector<double> & prob, int len, vector<double> & probSglUp )
{
	double psingle;

	probSglUp.reserve( len );

	probSglUp.push_back( 0. );
	for( int i = 1; i < len; i++ )
	{
		psingle = 1;
		for( int j = 0; j < i; j++ )
			psingle -= prob.at( i*len + j );
		probSglUp.push_back( psingle * ( i ) / ( len-1 ) );
	}
}


void getunpaired_down( vector<double> & prob, int len, vector<double> & probSglDown )
{
	double psingle;

	probSglDown.reserve( len );

	for( int i = 0; i < len-1; i++ )
	{
		psingle = 1;
		for( int j = i+1; j < len; j++ )
			psingle -= prob.at( i*len + j );
		probSglDown.push_back( psingle * ( len-i-1 ) / ( len-1 ) );
	}
	probSglDown.push_back( 0. );
}


vector<string> &split(const string &s, char delim, vector<string> &elems)
{
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}


vector<string> split(const string &s, char delim)
{
    vector<string> elems;
    return split(s, delim, elems);
}


void reducematrix(vector<double> & prob, int len, int precision)
{
	 char p[100];
     for( int i = 0; i < len; i++ ) {
    	 	sprintf(p, "%.*lf", precision, prob[i]);
    	 	prob.at(i) = atof(p);
     }
}


/**
   alignment of the pairing probabilities of two reference nucleotides (one line per dot plot)

   \param[seq_1]  sequence 1 ( S_a )
   \param[probSgl_1]  unpaired probabilities of S_a
   \param[idx_1_aln]  indices of nucleotides in the subsequence of S_a that is considered for alignment ( SS_a )
   \param[len_1] length of considered subsequence SS_a
   \param[seq_2]  sequence 2 ( S_b )
   \param[probSgl_2]  unpaired probabilities of S_b
   \param[idx_2_aln]  indices of nucleotides in the subsequence of S_b that is considered for alignment ( SS_b )
   \param[len_2]  length of considered subsequence SS_b
   \param[len_pair]  number of aligned bases

   \returns  Similarity score ( 0 .. 1 )
*/
double nwdp( string seq_1, vector<double> & probSgl_1, int * idx_1_aln, int L1, string seq_2, vector<double> & probSgl_2, int * idx_2_aln, int L2, int & Laln, double ** F, double ** Q, double ** P, char ** trF, char ** trQ, char ** trP )
{
	int global = 1;
    int i = 0, j = 0;
    double score, sim;

    /* opening gap penalty */
    double gapopen = open + ext;

    /* Initialize dynamic programming matrix */
    nwdp_initF_affinegaps( F, L1, L2, !global );

    /* Initialize Q and P matrices for cost of alignments that end with a gap */
    nwdp_initGap( Q, L1, L2 );
   	nwdp_initGap( P, L1, L2 );

    /* Initialize traceback matrices */
    nwdp_initTB( trF, L1, L2 );
    nwdp_initTB( trQ, L1, L2 );
    nwdp_initTB( trP, L1, L2 );

    /* base pair similarity of bases pairing with nuc_1 and bases pairing with nuc_2 */
    /* create alignment */
    char ptr;
    int x, y;
    for( j = 1; j <= L2; j++ )
    	for( i = 1; i <= L1; i++ )
        {
    	    // calculate P
   			P[ j ][ i ] = max3( P[ j-1 ][ i ] + ext, F[ j-1 ][ i ] + gapopen, INFINITE, &ptr );
            if( !global && P[ j ][ i ] < 0 )
            	P[ j ][ i ] = 0;
            trP[ j ][ i ] =  ptr;

    		// calculate Q
           	Q[ j ][ i ] = max3( INFINITE, F[ j ][ i-1 ] + gapopen, Q[ j ][ i-1 ] + ext, &ptr );
            if( !global && Q[ j ][ i ] < 0 )
            	Q[ j ][ i ] = 0;
            trQ[ j ][ i ] =  ptr;

            // calculate F
            switch( seq_1[ i-1 ] )
            {
            	case 'A':  x = 0 ;  break ;
                case 'C':  x = 1 ;  break ;
                case 'G':  x = 2 ;  break ;
                case 'T':  x = 3 ;  break ;
                case 'U':  x = 3 ;  break ;
            }
            switch( seq_2[ j-1 ] )
            {
            	case 'A':  y = 0 ;  break ;
                case 'C':  y = 1 ;  break ;
                case 'G':  y = 2 ;  break ;
                case 'T':  y = 3 ;  break ;
                case 'U':  y = 3 ;  break ;
            }

            score = theta * basesub[ x ][ y ] + (1 - theta) * ( 1 - abs( probSgl_1[ i-1 ] - probSgl_2[ j-1 ] ) );
            //cerr << i << " " << j << " " << score << endl;
            F[ j ][ i ] = max3( P[ j ][ i ], F[ j-1 ][ i-1 ] + score, Q[ j ][ i ], &ptr );
            if( !global && F[ j ][ i ] < 0 )
            	F[ j ][ i ] = 0;
            trF[ j ][ i ] =  ptr;
        }
    sim = F[ L2 ][ L1 ];
    i--; j--;

    /* backtracking */
    int k = 0;
    int w = 0;
    int maxmat = 0;
    while( i > 0 || j > 0 )
    {
    	if( !global )
    		if( (maxmat == 0 && F[ j ][ i ] == 0) || (maxmat == 1 && P[ j ][ i ] == 0) || (maxmat == 2 && Q[ j ][ i ] == 0) )
    			break;

    	switch( maxmat )
    	{
    		case 0 : switch( trF[ j ][ i ] )
    				 {
        			 	 case '|' : maxmat = 1;
        			 	 	 	 	break;
        			 	 case '\\': idx_1_aln[ k ] = i-1;
        			 	 	 	 	idx_2_aln[ k++ ] = j-1;
        			 	 	 	 	i--; j--;
        			 	 	 	 	w++;
        			 	 	 	 	break ;
        			 	 case '-' : maxmat = 2;
        			 	 	 	 	break;
    				 }
    				 break;
    		 case 1: switch( trP[ j ][ i ] )
    				 {
    		 	 	 	 case '|' : j--;
    		 	 	 	 	 	 	w++;
    		 	 	 	 	 	 	break;
    		 	 	 	 case '\\': maxmat = 0;
    		 	 	 	 	 	 	j--;
    		 	 	 	 	 	 	w++;
    		 	 	 	 	 	 	break;
    				 }
    				 break;
    		 case 2: switch( trQ[ j ][ i ] )
    		 	 	 {
 	 	 	 	 	 	 case '\\': maxmat = 0;
 	 	 	 	 	 	 	 	 	i--;
 	 	 	 	 	 	 	 	 	w++;
 	 	 	 	 	 	 	 	 	break;
 	 	 	 	 	 	 case '-' : i--;
 	 	 	 	 	 	 	 	 	w++;
 	 	 	 	 	 	 	 	 	break;
    		 	 	 }
    		 	 	 break;
    	}
    }

    Laln = k;
    reverse( idx_1_aln, Laln );
    reverse( idx_2_aln, Laln );

	#if DEBUG
        cout << "\nDynamic programming matrix: " << endl;
        double* tempidx_1 = new double[L1+1];
        double* tempidx_2 = new double[L2+1];
        for( int k = 0; k <= L1; k++ ) tempidx_1[ k ] = (double) k;
        for( int k = 0; k <= L2; k++ ) tempidx_2[ k ] = (double) k;
        cout << "F:" << endl;
    	print_matrixdp( F, tempidx_1, L1, tempidx_2, L2);
    	//print_matrixdp( trF, tempidx_1, L1, tempidx_2, L2);
        //cout << "P:" << endl;
    	//print_matrixdp( P, tempidx_1, L1, tempidx_2, L2);
    	//print_matrixdp( trP, tempidx_1, L1, tempidx_2, L2);
        //cout << "Q:" << endl;
    	//print_matrixdp( Q, tempidx_1, L1, tempidx_2, L2);
    	//print_matrixdp( trQ, tempidx_1, L1, tempidx_2, L2);
        delete[](tempidx_1);
        delete[](tempidx_2);
        //for(int i=0; i<Laln; i++){cerr << idx_1_aln[i] << " ";} cerr << endl;
        //for(int i=0; i<Laln; i++){cerr << idx_2_aln[i] << " ";} cerr << endl;
	#endif

 	return sim/Laln;
}


/*
 * initialize matrix of costs for alignment of prefixes (a_1 ... a_i; b_1 ... b_j)
 * for affine gap costs
 */
void nwdp_initF_affinegaps( double ** F, int L1, int L2, bool local )
{
	F[ 0 ][ 0 ] =  0. ;

	for( int j = 1; j <= L2; j++ )
		for( int i = 1; i <= L1; i++ )
			F[ j ][ i ] = INFINITE;

	for( int i = 1; i <= L1; i++ )
		F[ 0 ][ i ] = ( !local ) ? open + i * ext : 0;
	for( int j = 1; j <= L2; j++ )
		F[ j ][ 0 ] = ( !local ) ? open + j * ext : 0;
}


/*
 * initialize matrix of costs for alignment of prefixes (a_1 ... a_i; b_1 ... b_j)
 * that ends with a gap
 */
void nwdp_initGap( double ** Q, int L1, int L2 )
{
	Q[ 0 ][ 0 ] =  0. ;

	for( int j = 0; j <= L2; j++ )
		for( int i = 0; i <= L1; i++ )
			Q[ j ][ i ] = INFINITE;
}


void nwdp_initTB( char ** traceback, int L1, int L2 )
//void nwdp_initTB( vector< vector<char> > traceback, int L1, int L2 )
{
    traceback[ 0 ][ 0 ] = 'n' ;

    for( int i = 1; i <= L1; i++ )
    	traceback[ 0 ][ i ] =  '-' ;
    for( int j = 1; j <= L2 ; j++ )
    	traceback[ j ][ 0 ] =  '|' ;
}


double simdp( string seq_1, vector<double> & probDbl_1, int len_1, string seq_2, vector<double> & probDbl_2, int len_2, int * idx_1_aln, int * idx_2_aln, int len_aln )
{
	int offset_1, offset_2;
    int x, y;
	double sim = 0;
	MapClass::StringIntMap bpindex = MapClass::init_pair_index();

	/* sum scores of aligned base pair probability matrices */
	for( int i = 0; i < len_aln; i++ )
	{
		offset_1 = idx_1_aln[ i ] * len_1;
		offset_2 = idx_2_aln[ i ] * len_2;
		for( int j = i+1; j < len_aln; j++ )
		{
			x = bpindex[ string(&seq_1[ idx_1_aln[ i ] ], &seq_1[ idx_1_aln[ i ] ] + 1) + string(&seq_1[ idx_1_aln[ j ] ], &seq_1[ idx_1_aln[ j ] ] + 1) ];
			y = bpindex[ string(&seq_2[ idx_2_aln[ i ] ], &seq_2[ idx_2_aln[ i ] ] + 1) + string(&seq_2[ idx_2_aln[ j ] ], &seq_2[ idx_2_aln[ j ] ] + 1) ];
			/*cout << string(&seq_1[ idx_1_aln[ i ] ], &seq_1[ idx_1_aln[ i ] ]+1) + string(&seq_1[ idx_1_aln[ j ] ], &seq_1[ idx_1_aln[ j ] ]+1) << "\t"
					<< x << "\t" << string(&seq_2[ idx_2_aln[ i ] ], &seq_2[ idx_2_aln[ i ] ]+1) + string(&seq_2[ idx_2_aln[ j ] ], &seq_2[ idx_2_aln[ j ] ]+1) << "\t"
					<< y << "\t" << pairsub[x][y] << endl;*/
			//cout << i << "\t" << j << "\t" << probDbl_1[ offset_1 + idx_1_aln[j] ] << "\t" << probDbl_2[ offset_2 + idx_2_aln[ j ] ] << endl;

	    	if( mutation_rates_flag )
	    		sim += pairmutationrates[ x ][ y ] * ( 1 - abs( probDbl_1[ offset_1 + idx_1_aln[ j ] ] - probDbl_2[ offset_2 + idx_2_aln[ j ] ] ) );
	    	else
	    	{
	    		sim += theta * pairsub[ x ][ y ] + (1 - theta) * ( 1 - abs( probDbl_1[ offset_1 + idx_1_aln[ j ] ] - probDbl_2[ offset_2 + idx_2_aln[ j ] ] ) );
	    	}
		}
	}

	return 2*sim/len_aln/len_aln;
}


double max3( double f1, double f2, double f3, char* ptr )
{
	double  max = 0;

    if( f1 >= f2 && f1 >= f3 )
    {
    	max = f1;
        *ptr = '|';
    }
    else if( f2 > f3 )
    {
    	max = f2;
    	*ptr = '\\';
    }
    else
    {
    	max = f3;
    	*ptr = '-';
    }

    return max;
}


double max( double f1, double f2 )
{
    return ( f1 >= f2 ) ? f1 : f2;
}


template <typename T>
void print_matrixdp( T ** F, double * prob_1, int L1, double * prob_2, int L2 )
{
    cout << "        ";
    for( int i = 0; i < L1; i++ )
    	cout << prob_1[ i ] << "   ";
    cout << "\n  ";

    for( int j = 0; j <= L2; j++ )
    {
    	if( j > 0 )
    		cout << prob_2[ j-1 ] << " ";
    	for( int i = 0; i <= L1; i++ )
        {
    		cout.width( 3 );
            cout << F[ j ][ i ] << " ";
        }
        cout << endl;
    }
}


void freeMatrix(double ** matrix, int len)
{
	for(int i = 0; i < len; i++)
		delete[] matrix[i];
	delete[] matrix;
}


void reverse( int * list, int len )
{
	int temp;
	for( int i = 0; i < len/2; i++ )
	{
		temp = list[ i ];
	    list[ i ] = list[ len-i-1 ];
	    list[ len-i-1 ] = temp;
	}
}


/*
 * return normalized log odds of base pair probabilities weighted by minimal considered paired probabilities
 */
void getlogoddsDbl( vector<double> & probDbl, string seq, int len, double pnull )
{
	int offset;

    for( int i = 0; i < len; i++ )
    {
    	offset = i * len;
    	for( int j = 0; j < len; j++ )
    		probDbl[ offset + j ] = max( 0, log( probDbl[ offset + j ] / pnull ) / log( 1 / pnull ) );
    }
}


/*
 * return normalized log odds of unpaired probabilities weighted by minimal considered unpaired probability
 */
void getlogoddsSgl( vector<double> & probSgl, string seq, int len, double pnull )
{
    for( int i = 0; i < len; i++ )
    	probSgl[ i ] = max( 0, log( probSgl[ i ] / pnull ) / log( 1 / pnull ) );
}


void printalign( string & seq_1, int * idx_1_aln, string & seq_2, int * idx_2_aln, int len_aln )
{
	int len_1 = seq_1.length();
	int len_2 = seq_2.length();

	char *seqar_1 = new char[len_1+1];
	strcpy(seqar_1, seq_1.c_str());
	char *seqar_2 = new char[len_2+1];
	strcpy(seqar_2, seq_2.c_str());
	string seq_1_aln, seq_2_aln;

	int offset_1 = 0, offset_2 = 0;
	for (int k=0; k < len_aln; k++) {
		while( idx_1_aln[ k ] != offset_1 ) { seq_1_aln.push_back( seqar_1[ offset_1 ] ); seq_2_aln.push_back( '-' ); offset_1++; }
		while( idx_2_aln[ k ] != offset_2 ) { seq_2_aln.push_back( seqar_2[ offset_2 ] ); seq_1_aln.push_back( '-' ); offset_2++; }
		seq_1_aln.push_back( seqar_1[ idx_1_aln[ k ] ] ); offset_1++;
		seq_2_aln.push_back( seqar_2[ idx_2_aln[ k ] ] ); offset_2++;
	}
	while( offset_1 < len_1 ) { seq_1_aln.push_back( seqar_1[ offset_1 ] ); seq_2_aln.push_back( '-' ); offset_1++; }
	while( offset_2 < len_2 ) { seq_2_aln.push_back( seqar_2[ offset_2 ] ); seq_1_aln.push_back( '-' ); offset_2++; }

	cout << seq_1_aln << endl;
	cout << seq_2_aln << endl;

    delete[] seqar_1;
    delete[] seqar_2;
}


/*
 * generate a random number from 0 to sum of entries in sam_prob
 */
int statistical_sampling(unsigned int *sam_prob, int sam_len)
{
	long sum = 0;
	int i, sample;

	for( i=0; i<sam_len; i++ ) {
		sum += sam_prob[i];
	}

	if( sum>0 )
		sample = rand()%sum;
	else
		return 0;
	//printf("Sample: %i\n",sample);

	sum = 0;
	for( i=0; i<sam_len; i++ )
	{
		sum += sam_prob[i];
		if (sample<=sum) {
			return i;
		}
	}

	return 0;
}


/*
 *
 */
void prob_backtracking( int * idx_1_aln, int len_1, int * idx_2_aln, int len_2, int & Laln, double ** F )
{
	/* run through dynamic programming matrix using probabilistic sampling */
	int i = len_1;
	int j = len_2;
    int k = 0;
    int w = 0;
    int samprob;
    unsigned int btscores[3];

    while( i > 0 && j > 0 )
    {
   		btscores[0] = round(F[ j ][ i ]*1000);
   		btscores[1] = ( round(F[ j-1 ][ i ]) > 0 ) ? round(F[ j-1 ][ i ]*1000) : 0;
   		btscores[2] = ( round(F[ j ][ i-1 ]) > 0 ) ? round(F[ j ][ i-1 ]*1000) : 0;

   		samprob = statistical_sampling(btscores, 3);
   		//cerr << i << " " << j << " " << btscores[0] << " " << btscores[1] << " " << btscores[2] << " " << samprob << endl;

   		switch( samprob )
   		{
   			case 0:	/* base pair */
   					idx_1_aln[ k ] = i-1;
   					idx_2_aln[ k++ ] = j-1;
   					i--; j--;
   					w++;
   					break;
   			case 1: /* base i unmatched */
   					j--;
   					w++;
   					break;
   			case 2: /* base j unmatched */
   					i--;
   					w++;
   					break;
   			default:  abort();
   		}
    }

    Laln = k;
    reverse( idx_1_aln, Laln );
    reverse( idx_2_aln, Laln );
}


/*
 * calculate sequence / unpaired similarity score of suboptimal alignment
 */
double simbp( string seq_1, vector<double> & probSgl_1, vector<double> & probDblUp_1, int len_1, string seq_2, vector<double> & probSgl_2, vector<double> & probDblUp_2, int len_2, int * idx_1_aln, int * idx_2_aln, int len_match, char * aln_code )
{
	double sim = 0;
	double nzeta = 1-zeta;
	double ntheta = 1-theta;
	int len_aln = strlen(aln_code);

	/* sum scores of aligned bases */
	MapClass::CharIntMap baseindex = MapClass::init_base_index();
    int x, y;
    for( int i = 0; i < len_match; i++ )
    {
    	x = baseindex[ seq_1[ idx_1_aln[ i ] ] ];
    	y = baseindex[ seq_2[ idx_2_aln[ i ] ] ];

    	if( mutation_rates_flag )
    		sim += basemutationrates[ x ][ y ] *
    		       ( zeta * ( 1 - abs( probSgl_1[ idx_1_aln[ i ] ] - probSgl_2[ idx_2_aln[ i ] ] ) ) +
    		       nzeta * ( 1 - abs( probDblUp_1[ idx_1_aln[ i ] ] - probDblUp_2[ idx_2_aln[ i ] ] ) ) );
    	else
    		sim += theta * basesub[ x ][ y ] +
    			   ntheta * ( 1 - abs( probSgl_1[ idx_1_aln[ i ] ] - probSgl_2[ idx_2_aln[ i ] ] ) );
    }

    /* sum scores of gaps */
    sim += gappenalty( aln_code );

    return sim / len_aln;
}


/*
 * add scores of gaps
 */
double gappenalty( char * aln )
{
	double penalty = 0;

	int i = 0 , x = 0, y = 0, z = 0, m = 0, n = 0;
	while( i < strlen(aln) )
	{
		switch( aln[i] )
		{
			case '|': /* match */
				      n = m = 0;
				      x = 1;
					  break;
			case ':': /* mismatch */
				      n = m = 0;
				      x = 1;
					  break;
			case '^': /* gaps in longer sequence */
			          if( free_endgaps_flag && x == 0 && z == 0 && m == 0 )
			          {
			        	  n = 1;
			        	  y = 1;
			        	  break;
			          }
					  if( n == 0 )
						  penalty -= open;
					  else
						  penalty -= ext;
					  n = 1;
					  m = 0;
					  y = 1;
					  break;
			case '.': /* gaps in shorter sequence */
     		          if( free_endgaps_flag && x == 0 && y == 0 && n == 0 )
     		          {
     		        	  m = 1;
     		        	  z = 1;
     		        	  break;
     		          }
				      if( m == 0 )
				    	  penalty -= open;
					  else
						  penalty -= ext;
					  m = 1;
					  n = 0;
					  z = 1;
					  break;
			default:  abort();
		}
		i++;
	}

    if( free_endgaps_flag )
    {
    	i--;
    	m = 0, n = 0;
    	while( i >= 0 )
    	{
    		switch( aln[i] )
    		{
				case '^': /* gaps in longer sequence */
						  if( m == 1 )
							  return penalty;
						  if( n == 0 )
							  penalty += open;
						  else
							  penalty += ext;
						  n = 1;
						  break;
				case '.': /* gaps in shorter sequence */
					      if( n == 1 )
					    	  return penalty;
					      if( m == 0 )
					    	  penalty += open;
					      else
					    	  penalty += ext;
					      m = 1;
					      break;
				default:  return penalty;
    		}
    		i--;
    	}
    }

	return penalty;
}


/*
 * convert alignment string of probA to lists of indices
 */
int convert(char * aliggAlignment, int * idx_1_aln, int * idx_2_aln)
{
	int i = 0 , k = 0, m = 0, n = 0;
	while( i < strlen(aliggAlignment) )
	{
		switch( aliggAlignment[i] )
		{
			case '|': /* match */
					  idx_1_aln[ k ] = m++;
					  idx_2_aln[ k++ ] = n++;
					  break;
			case ':': /* mismatch */
				  	  idx_1_aln[ k ] = m++;
				  	  idx_2_aln[ k++ ] = n++;
					  break;
			case '^': /* gaps in longer sequence */
					  n++;
					  break;
			case '.': /* gaps in shorter sequence */
					  m++;
					  break;
			default:  abort();
		}
		i++;
	}

	return k;
}


/*
 * count the flanking gaps of the alignment
 */
int count_endgaps(int * idx_1_aln, int len_1, int * idx_2_aln, int len_2, int len_match)
{
	int fgaps = 0;

	/* 5' gaps */
	fgaps += idx_1_aln[0] + idx_2_aln[0];
	/* 3' gaps */
	fgaps += ( len_1-idx_1_aln[len_match-1]-1 ) + ( len_2-idx_2_aln[len_match-1]-1 );

	//cout << "COUNT_ENDGAPS = " << fgaps << endl;
	return fgaps;
}


/* adaptation of probA-function <align> for DotAligner */
aligm align_da(sequ *seq_array, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short )
{
    int i,l;
    int l0,l1;
    iseq *sti;
    u_sc stm;
    score_m **matser;
    score_m scoremat;
    aligm alig;

    so = sequences(seq_array);/* looks for the longer sequence and stores it in
			     so[0],the shorter one is stored in so[1],
			     both sequences are converted to uppercase
			     letters! */

    l0=strlen(so[0].seq);
    l1=strlen(so[1].seq);

    if(typ_flag==1)
    {
        MAT_SER[0]='D';MAT_SER[1]='N';MAT_SER[2]='A'; MAT_SER[3]='\0';
        scmat=convert_matrix(dnamt);
        si=code_seq(so,scmat);
    }
    else
    {
       fprintf(stderr,"Problem with score matrix selection\n");
        exit(11);
    }

    alig=align_2_da(so, si, scmat, probSgl_long, probDblUp_long, probSgl_short, probDblUp_short);

    alig=add_sequences(so,alig);

    /* allocated in init_matseries */
    free(matser);

    return(alig);
}


/* adaptation of probA-function <align_2> for DotAligner */
aligm align_2_da( sequ *os, iseq *is, u_sc smat, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short )
{
	int i,j;
	int l0;
	int l1;
	int l;
	int *s0, *s1;
	matrix m;
	aligm a;

	double sim;
	double nzeta = 1-zeta;
	double ntheta = 1-theta;

	MapClass::CharIntMap baseindex = MapClass::init_base_index();
    int x, y;

	s0=is[0].s;
	s1=is[1].s;

	l=strlen(smat.monomers);
	l0=strlen(os[0].seq);
	l1=strlen(os[1].seq);

	m.M=(real**) calloc((l1+1),sizeof(real*));
	for(i=0; i<l1+1; i++)
	    m.M[i]=(real*) calloc((l0+1),sizeof(real));

	m.E=(real**) calloc((l1+1),sizeof(real*));
	for(i=0; i<l1+1; i++)
	    m.E[i]=(real*) calloc((l0+1),sizeof(real));

	m.F=(real**) calloc((l1+1),sizeof(real*));
	for(i=0; i<l1+1; i++)
	    m.F[i]=(real*) calloc((l0+1),sizeof(real));

	/* recursion:
	   Mi,j=max{ Mi-1,j-1+match/mismatch,
	           { Ei-1,j-1+match/mismatch,
	           { Fi-1,j-1+match/mismatch;

	   Ei,j=max{ Mi,j-1+open,
	           { Ei,j-1+ext;

       Fi,j=max{ Mi-1,j+open,
	           { Ei-1,j+open,
	           { Fi-1,j+ext;   */

	/* initialisation: M    M(0,0) = E0,0 = F0,0 = 0;
	                        k = length of the gap;
	                        M(0,j>0) = M(i>0,0) = 4*open + 10*(k-1)*ext;(*)
	   (*) the values are selected in order to avoid the use of these entries (imagine a kind of "empty" alignment (e.g. seq1 aligned only to gap
		   characters), which you simply don't want to consider (you are looking for the best alignment in which i and j are paired!)

	                   E    E(0,j>0) = open + (k-1)*ext;
			                E(i>0,0) = 4*open + 10*(k-1)*ext;(*)

			           F    F(0,j>0) = 4*open + 10*(k-1)*ext;(*)
			                F(i>0,0) = open + (k-1)*ext; */


	/* M(0,0) = E(0,0) = F(0,0) = 0; */
	m.M[0][0]=m.E[0][0]=m.F[0][0]=0;
	/* initialise M,E,F*/

	i=0;
	j=1;
	m.M[i][j]=m.F[i][j]=4*smat.mat[l-2][s0[0]]; /* 4*Open in s0 */
	m.E[i][j]=smat.mat[l-2][s0[0]]; /* Open in s0 */
	for(j=2; j< l0+1; j++)
	{
	    /* 10*Extent in s0 */
	    m.M[i][j]=m.F[i][j]=m.M[i][j-1]+10*smat.mat[l-1][s0[j-1]];
	    m.E[i][j]=m.E[i][j-1]+smat.mat[l-1][s0[j-1]]; /* Extent in s0 */
	}

	i=1;
	j=0;
	m.M[i][j]= m.E[i][j]=4*smat.mat[l-2][s1[0]]; /* 4*Open in s1 */
	m.F[i][j]=smat.mat[l-2][s1[0]]; /* Open in s1 */
	for(i=2; i< l1+1; i++)
	{
	    /* 10*Extent in s1 */
	    m.M[i][j]=m.E[i][j]=m.M[i-1][j]+10*smat.mat[l-1][s1[i-1]];
	    m.F[i][j]=m.F[i-1][j]+smat.mat[l-1][s1[i-1]]; /* Extent in s1 */
	}

	if(Egap_flag)
	{
	    /* initialise E(0,1....l0+1 */
	    i=0;
	    /* value for endgaps */
	    for(j=1; j< l0+1; j++)
	    	m.E[i][j]=m.E[i][j-1]+smat.endgaps; /* value for endgaps */
	    /* initialise F(i....l1+1,0) */
	    j=0;
	    for(i=1; i< l1+1; i++)
		    m.F[i][j]=m.F[i-1][j]+smat.endgaps; /* value for endgaps */
	}

    for(i=1; i<l1+1; i++)
	{
	    for(j=1; j<l0+1; j++)
		{
	    	double max, e, f;

		    if(Egap_flag && i == l1 && j < l0)
		    {
		    	/* E */
		        m.E[i][j]=m.M[i][j-1]+smat.endgaps; /*Open in s1 */
		        e=m.E[i][j-1]+smat.endgaps; /* Extent in s1 */
		        m.E[i][j]=MAX(m.E[i][j],e);
		        e=0;

		        /* F */
		        m.F[i][j]=m.M[i-1][j]+smat.mat[l-2][s0[j-1]]; /* Open in s0 */
		        e=m.E[i-1][j]+smat.mat[l-2][s0[j-1]];
		        f=m.F[i-1][j]+smat.mat[l-1][s0[j-1]]; /* Extent in s0 */

		        max=MAX(m.F[i][j],e);
		        m.F[i][j]=MAX(max,f);
		        max=0;
		        e=0;
		        f=0;
		    }
		    else if(Egap_flag && j == l0 && i < l1)
		    {
		    	/* F */
		    	m.F[i][j]=m.M[i-1][j]+smat.endgaps; /*Open in s0 */
		    	e=m.E[i-1][j]+smat.endgaps;
		    	f=m.F[i-1][j]+smat.endgaps; /* Extent in s0 */

		    	max=MAX(m.F[i][j],e);
		    	m.F[i][j]=MAX(max,f);
		    	max=0;
		    	e=0;
		    	f=0;

		    	/* E */
		    	m.E[i][j]=m.M[i][j-1]+smat.mat[l-2][s1[i-1]]; /* Open in s1 */
		    	e=m.E[i][j-1]+smat.mat[l-1][s1[i-1]]; /* Extent in s1 */

		    	m.E[i][j]=MAX(m.E[i][j],e);
		    	e=0;
		    }
		    else if(Egap_flag && j == l0 && i == l1)
		    {
		    	/* F */
		    	m.F[i][j]=m.M[i-1][j]+smat.endgaps; /*Open in s0 */
		    	e=m.E[i-1][j]+smat.endgaps;
		    	f=m.F[i-1][j]+smat.endgaps; /* Extent in s0 */
		    	max=MAX(m.F[i][j],e);
		    	m.F[i][j]=MAX(max,f);
		    	max=0;
		    	e=0;
		    	f=0;

		    	/* E */
		    	m.E[i][j]=m.M[i][j-1]+smat.endgaps; /*Open in s1 */
		    	e=m.E[i][j-1]+smat.endgaps; /* Extent in s1 */
		    	m.E[i][j]=MAX(m.E[i][j],e);
		    	e=0;
		    }
		    else
		    {
		    	/* E */
		    	m.E[i][j]=m.M[i][j-1]+smat.mat[l-2][s1[i-1]]; /* Open in s1 */
		    	e=m.E[i][j-1]+smat.mat[l-1][s1[i-1]]; /* Extent in s1 */

		    	m.E[i][j]=MAX(m.E[i][j],e);
		    	e=0;

		    	/* F */
		    	m.F[i][j]=m.M[i-1][j]+smat.mat[l-2][s0[j-1]]; /* Open in s0 */
		    	e=m.E[i-1][j]+smat.mat[l-2][s0[j-1]];
		    	f=m.F[i-1][j]+smat.mat[l-1][s0[j-1]]; /* Extent in s0 */

		    	max=MAX(m.F[i][j],e);
		    	m.F[i][j]=MAX(max,f);
		    	max=0;
		    	e=0;
		    	f=0;
		    }

		    /* HERE ADD THE DOTALIGNER SINGLE STRANDED SCORING SCHEME */
		    /* replace smat.mat[s1[i-1]][s0[j-1]] */
	    	x = baseindex[ os[0].seq[ j-1 ] ];
	    	y = baseindex[ os[1].seq[ i-1 ] ];

		    if( mutation_rates_flag )
		    	sim = basemutationrates[ x ][ y ] *
	  		       	  ( zeta * ( 1 - abs( probSgl_long[ j-1 ] - probSgl_short[ i-1 ] ) ) +
	  		          nzeta * ( 1 - abs( probDblUp_long[ j-1 ] - probDblUp_short[ i-1 ] ) ) );
		    else
		    	sim = theta * basesub[ x ][ y ] +
	  			      ntheta * ( 1 - abs( probSgl_long[ j-1 ] - probSgl_short[ i-1 ] ) );

		    /* M */
		    m.M[i][j]=m.M[i-1][j-1]+sim;
		    /* M or MM */
		    e=m.E[i-1][j-1]+sim; /* M or MM */
		    f=m.F[i-1][j-1]+sim; /* M or MM */

		    max=MAX(m.M[i][j],e);
		    m.M[i][j]=MAX(max,f);

		    max=0;
		    e=0;
		    f=0;
		}
	}

    a=trace_back_da(m, os, is, smat, probSgl_long, probDblUp_long, probSgl_short, probDblUp_short);

    /* free the space allocated for the alignment matrices */
    for(i=0; i<l1+1; i++)
	{
    	free(m.M[i]);
	    free(m.E[i]);
	    free(m.F[i]);
	}
    free(m.M);
	free(m.E);
	free(m.F);

	return(a);
}


/* adaptation of probA-function <trace_back> for DotAligner */
aligm trace_back_da(matrix m, sequ *os, iseq *is, u_sc smat, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short)
{
	int i,j,k,l;
	int l0;
	int l1;
	int la;
	int *s0, *s1;
	real **M, **E, **F;
	aligm a;
	double max;
	enum Zust { Mm, Ee, Ff };
	enum Zust state;
	double open, ext;

	double sim;
	double nzeta = 1-zeta;
	double ntheta = 1-theta;

	MapClass::CharIntMap baseindex = MapClass::init_base_index();
    int x, y;

    s0=is[0].s;
	s1=is[1].s;
	l0=strlen(os[0].seq);
	l1= strlen(os[1].seq);

	l=strlen(smat.monomers);
	la=l0+l1;

	a.a=(char*) calloc((la),sizeof(char));
	M=m.M;
	E=m.E;
	F=m.F;

	i=l1;
	j=l0;
	la-=1;
	k=0;

	state = Mm; max = M[i][j];
	if (E[i][j]>max) {max = E[i][j]; state = Ee; }
	if (F[i][j]>max) {max = F[i][j]; state = Ff; }

	/* determine the score of an optimal alignment */
	a.score=max;

	while((i>0 && j>=0) || (i>=0 && j>0) )
	{
		switch (state)
		{
			case Mm:
				//printf("i %d  j %d  M[i][j] %6.4f\n",i,j,M[i][j]);
				a.a[la] = (s0[j-1] == s1[i-1]) ? M1 : M2;

			    /* replace smat.mat[s1[i-1]][s0[j-1]] */
		    	x = baseindex[ os[0].seq[ j-1 ] ];
		    	y = baseindex[ os[1].seq[ i-1 ] ];

			    if( mutation_rates_flag )
			    	sim = basemutationrates[ x ][ y ] *
		  		       	  ( zeta * ( 1 - abs( probSgl_long[ j-1 ] - probSgl_short[ i-1 ] ) ) +
		  		          nzeta * ( 1 - abs( probDblUp_long[ j-1 ] - probDblUp_short[ i-1 ] ) ) );
			    else
			    	sim = theta * basesub[ x ][ y ] +
		  			      ntheta * ( 1 - abs( probSgl_long[ j-1 ] - probSgl_short[ i-1 ] ) );

				if( EQUAL(M[i-1][j-1]+sim, M[i][j]) )
				{/* nothing to do */}
				else if( EQUAL(E[i-1][j-1]+sim, M[i][j]) )
				{ state = Ee; }
				else if ( EQUAL(F[i-1][j-1]+sim, M[i][j]) )
				{ state = Ff; }
				else {
					fprintf(stderr, "Error in backtracking case M\n"); exit(27);
				}
				i--; j--;
				break;

			case Ee:
				//printf("i %d  j %d  E[i][j] %6.4f\n",i,j,E[i][j]);
				/* (3) _ =  gaps in the shorter seq(s1,E) */
				a.a[la] = G3;
				if(i==0 && j==1)
				{
					--j;
					break;
				}
				if (Egap_flag && (i==l1 || i==0) )
					open = ext = smat.endgaps;
				else {
					open = smat.mat[l-2][s0[j]]; ext = smat.mat[l-1][s0[j]];
				}

				if ( EQUAL(E[i][j], E[i][j-1]+ext))
				{/* do nothing */}
				else if (EQUAL(E[i][j], M[i][j-1]+open))
					state = Mm;
				else {
					fprintf(stderr, "Error in backtracking case E\n"); exit(28);
				}
				j--;
				break;

			case Ff:
				//printf("i %d  j %d  F[i][j] %6.4f\n",i,j,F[i][j]);
				/* (4) _ =  gaps in the longer seq(s0,F) */
				a.a[la] = G4;
				if(i==1 && j==0)
				{
					--i;
					break;
				}
				if (Egap_flag && (j==l0 || j==0) )
					open = ext = smat.endgaps;
				else {
					open = smat.mat[l-2][s1[i]]; ext = smat.mat[l-1][s1[i]];
				}

				if (EQUAL(F[i][j], F[i-1][j]+ext))
				{/* do nothing */}
				else if (EQUAL(F[i][j], M[i-1][j]+open))
					state = Mm;
				else if (EQUAL(F[i][j], E[i-1][j]+open))
					state = Ee;
				else {
					fprintf(stderr, "Error in backtracking case F\n"); exit(29);
				}
				i--;
				break;
		}
		la--;
	}

	for(j=0, i=la+1; i<l0+l1; j++, i++)
		a.a[j] = a.a[i];
	a.a[j] = '\0';

	i=strlen(a.a);

	if((a.a = (char*) realloc(a.a,(i+1)*sizeof(char))) == 0)
	{
		fprintf(stderr,"trace_back:error in reallocation\n");
		exit(12);
	}

	return(a);
}


/* adaptation of probA-function <partition_f> for DotAligner */
real **partition_f_da(aligm alig, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short)
{
	pmat = partf_da( alig, probSgl_long, probDblUp_long, probSgl_short, probDblUp_short );

	rmat = revers_partf_da( so, si, scmat, alig, probSgl_long, probDblUp_long, probSgl_short, probDblUp_short );

	Pr = prop_alig( so, si, scmat, pmat.M, rmat.M, alig );

	return(Pr);
}


/* adaptation of probA-function <partf> for DotAligner */
matrix partf_da(aligm a, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short)
{
	int i,j;
	double opt;
	int l0,l1,l;
	int *s0, *s1;
	double Z;
	double s;
	double sk;     /* scaling factor */
	real **zM;
	real **zE;
	real **zF;
	matrix pf;
	u_sc m;

	double sim;
	double nzeta = 1-zeta;
	double ntheta = 1-theta;

	MapClass::CharIntMap baseindex = MapClass::init_base_index();
    int x, y;

    m=scmat;

	s0=si[0].s;
	s1=si[1].s;
	l0=strlen(so[0].seq);
	l1= strlen(so[1].seq);
	l=strlen(m.monomers);

	/* print some infos */
	if( verbose_flag )
	{
	    printf("Partition function over all sequence alignments of the two input sequences:\n");
		printf(" scoring matrix: %s;\n",m.name);
		printf(" gap penalties open:%.3f, extend:%.3f;\n",m.mat[l-2][s0[0]],m.mat[l-1][s0[0]]);
		printf(" temperature: %.2f\n", 1.0/BETA);
		if(Egap_flag)
			printf(" endgap penalty: %.3f\n",m.endgaps);
	}

	/* [i][j] => i runs along the rows, j runs along the columns !!!!*/

	/* zM: calculate only those alignments in which i and j are matched
	   zM(i,j)=exp(BETA*(M or MM))*(zM(i-1,j-1)+zE(i-1,j-1)+zF(i-1,j-1)); */
	zM=(real**) calloc((l1+1), sizeof(real*));
	for(i=0; i<l1+1; i++)
		zM[i]=(real*) calloc((l0+1), sizeof(real));

	/* the recurrence for zF: zF(i,j) = zM(i-1,j)*exp(BETA*(-OPEN)) +
	                                    + zE(i-1,j)*exp(BETA*(-OPEN)) +
	         	             		    + zF(i-1,j)*exp(BETA*(-EXT));
	   you need the previous row to calculate zFi,j;
	   in the 2D array the values of the previous row are stored in the
	   first row ( zF[0][...] ), the newly calculated values are stored
	   in the second one ( zF[1][...] );
	   after the calculation of the last column the values from zF[1][...]
	   are copied to zF[0][...]....*/

	zF=(real**) calloc((l1+1), sizeof(real*));
	for(i=0; i<l1+1; i++)
		zF[i]=(real*) calloc((l0+1), sizeof(real));

	/* the recurrence for zE is zE(i,j)= zM(i,j-1)*exp(BETA*(-OPEN) +
	                                     + zE(i,j-1)*exp(BETA*(-EXT);
	   one therefore never leaves the current row but one needs the values
	   of the previous column to calculate zF(i,j)  */

	zE=(real**) calloc((l1+1), sizeof(real*));
	for(i=0; i<l1+1; i++)
		zE[i]=(real*) calloc((l0+1), sizeof(real));

	/* opt: optimal score of the alignment of the two sequences */
	opt=a.score;
	s=2.0*opt/(l0+l1);
	sk=exp(BETA*s/2);

	zM[0][0]=1;

	/* zM[0][1...l0], zM[1...l0][0], zF[0][1...l0] and zE[1...l0][0]
	   are all 0, calloc initializes the allocated memory with 0 */

	/* zE[0][0]=0;
	   zE[0][j]=zM[0][j-1]*exp(BETA*(-OPEN)+zE[0][j-1]*exp(BETA*(-EXT) */

	i=0;
	j=1;
	zE[0][1]=zM[0][0]*exp(BETA*m.mat[l-2][s0[0]])/sk; /* Open in s0 */

	/* zM[0][j>0]=0 */
	for(j=2; j<=l0; j++)
	{
		/* Extent in s0 */
	    zE[0][j]=zE[0][j-1]*exp(BETA*m.mat[l-1][s0[j-1]])/sk;/* Extent in s0 */
	}

	/* initialize zF[0][0]=0;
	              zF[1][0]=zM[0][0]*exp(BETA*(-OPEN)); */
	i=1;
	j=0;
	zF[1][0]=zM[0][0]*exp(BETA*m.mat[l-2][s1[0]])/sk; /* Open in s1 */


	for(i=2; i<=l1; i++)
	{
		zF[i][0]=zF[i-1][0]*exp(BETA*m.mat[l-1][s1[i-1]])/sk; /* Extend in s1 */
	}

	if(Egap_flag)
	{
		i=0;
	    j=1;
	    zE[0][1]=zM[0][0]*exp(BETA*m.endgaps)/sk;
	    for(j=2; j<=l0; j++)
	    	zE[0][j]=zE[0][j-1]*exp(BETA*m.endgaps)/sk;
	    i=1;
	    j=0;
	    zF[1][0]=zM[0][0]*exp(BETA*m.endgaps)/sk;
	    for(i=2; i<=l1; i++)
	    	zF[i][0]=zF[i-1][0]*exp(BETA*m.endgaps)/sk;
	}
	for(i=1; i<=l1; i++)
	{
		for(j=1; j<=l0; j++)
		{
			double z;

			if(Egap_flag && i == l1 && j < l0)
		    {
				zE[i][j]=zE[i][j-1]*exp(BETA*m.endgaps);
				zE[i][j]+=zM[i][j-1]*exp(BETA*m.endgaps);
				zE[i][j]/=sk;

				/* Extent in s1 */
				zF[i][j]=zF[i-1][j]*exp(BETA*m.mat[l-1][s1[i-1]]);
				/* Open in s1 */
				zF[i][j]+=(zM[i-1][j]+zE[i-1][j])*exp(BETA*m.mat[l-2][s1[i-1]]);
				zF[i][j]/=sk;
		    }
			else if(Egap_flag && j == l0 && i < l1)
		    {
				zF[i][j]=zF[i-1][j]*exp(BETA*m.endgaps);
				zF[i][j]+=(zM[i-1][j]+zE[i-1][j])*exp(BETA*m.endgaps);
				zF[i][j]/=sk;

				/* Extent in s0 */
				zE[i][j]=zE[i][j-1]*exp(BETA*m.mat[l-1][s0[j-1]]);
				/* Open in s0 */
				zE[i][j]+=zM[i][j-1]*exp(BETA*m.mat[l-2][s0[j-1]]);
				zE[i][j]/=sk;
		    }
			else if(Egap_flag && j == l0 && i == l1)
		    {
				zE[i][j]=zE[i][j-1]*exp(BETA*m.endgaps);
				zE[i][j]+=zM[i][j-1]*exp(BETA*m.endgaps);
				zE[i][j]/=sk;

				zF[i][j]=zF[i-1][j]*exp(BETA*m.endgaps);
				zF[i][j]+=(zM[i-1][j]+zE[i-1][j])*exp(BETA*m.endgaps);
				zF[i][j]/=sk;
		    }
			else
		    {
				/* Extent in s0 */
				zE[i][j]=zE[i][j-1]*exp(BETA*m.mat[l-1][s0[j-1]]);
				/* Open in s0 */
				zE[i][j]+=zM[i][j-1]*exp(BETA*m.mat[l-2][s0[j-1]]);
				zE[i][j]/=sk;

				/* Extent in s1 */
				zF[i][j]=zF[i-1][j]*exp(BETA*m.mat[l-1][s1[i-1]]);
				/* Open in s1 */
				zF[i][j]+=(zM[i-1][j]+zE[i-1][j])*exp(BETA*m.mat[l-2][s1[i-1]]);
				zF[i][j]/=sk;
		    }

		    /* replace m.mat[s1[i-1]][s0[j-1]] */
	    	x = baseindex[ so[0].seq[ j-1 ] ];
	    	y = baseindex[ so[1].seq[ i-1 ] ];

		    if( mutation_rates_flag )
		    	sim = basemutationrates[ x ][ y ] *
	  		       	  ( zeta * ( 1 - abs( probSgl_long[ j-1 ] - probSgl_short[ i-1 ] ) ) +
	  		          nzeta * ( 1 - abs( probDblUp_long[ j-1 ] - probDblUp_short[ i-1 ] ) ) );
		    else
		    	sim = theta * basesub[ x ][ y ] +
	  			      ntheta * ( 1 - abs( probSgl_long[ j-1 ] - probSgl_short[ i-1 ] ) );
			zM[i][j] = exp(BETA*sim) * (zM[i-1][j-1]+zE[i-1][j-1]+zF[i-1][j-1]);
			zM[i][j]/=sk*sk;
		}
	}

	Z=zM[i-1][j-1]+zF[i-1][j-1]+zE[i-1][j-1];
	zM[0][0]=Z;

	pf.M=zM;
	pf.E=zE;
	pf.F=zF;

	return(pf);
}


/* adaptation of probA-function <revers_partf> for DotAligner */
matrix revers_partf_da(sequ *os, iseq *is, u_sc m, aligm a, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short )
{
	int i,j;
	double opt;
	int l0;
	int l1;
	int l;
	int *s0, *s1;
	double s;
	double sk; /* scaling factor */
	double Zr;
	real **rM;
	real **rE;
	real **rF;
	matrix pf;

	double sim;
	double nzeta = 1-zeta;
	double ntheta = 1-theta;

	MapClass::CharIntMap baseindex = MapClass::init_base_index();
    int x, y;

    l0=strlen(os[0].seq);
	l1=strlen(os[1].seq);
	l=strlen(m.monomers);
	s0=is[0].s;
	s1=is[1].s;

	/*i runs along the rows, j runs along the columns */

	/* Algorithm for the reverse partition function */
	/* rMi,j = exp(BETA * Match) * (rM i-1,j-1 + rE i-1,j-1 + rF i-1,j-1 );

       rEi,j = exp(BETA * Open) * (rMi,j-1 + rFi,j-1) +
               exp(BETA * Extend) * (rEi,j-1);

       rFi,j = exp(BETA * Open) * (rMi-1,j) +
             exp(BETA * Extend) * (rFi-1,j); */

	rM=(real**) calloc((l1+1),sizeof(real*));
	for(i=0; i<l1+1; i++)
		rM[i]=(real*) calloc((l0+1),sizeof(real));

	rF=(real**) calloc((2),sizeof(real*));
	for(i=0; i<2; i++)
		rF[i]=(real*) calloc((l0+1),sizeof(real));

	rE=(real**) calloc((2),sizeof(real*));
	for(i=0; i<2; i++)
		rE[i]=(real*) calloc((l0+1),sizeof(real));

	/* opt: optimal score of the alignment of the two sequences */
	opt=a.score;
	s=2.0*opt/(l0+l1);
	sk=exp(BETA*s/2);

	rM[l1][l0]=1;
	/* Open in s0, => the l0st member of s0 is s0[l0-1] */

	rE[0][l0-1]=rM[l1][l0]*exp(BETA*m.mat[l-2][s0[l0-1]])/sk;
	for(j=l0-2; j>=0; j--)
    {
		/* Extent in s0 */
		rE[0][j]=rE[0][j+1]*exp(BETA*m.mat[l-1][s0[j]])/sk;
    }

	if(Egap_flag)
    {
		i=l1;
		j=l0-1;
		rE[0][l0-1]=rM[l1][l0]*exp(BETA*m.endgaps)/sk;
		for(j=l0-2; j>=0; j--)
		{
			/* Extent in s0 */
			rE[0][j]=rE[0][j+1]*exp(BETA*m.endgaps)/sk;
		}
    }

	for(i=l1-1; i>=0; i--)
    {
		for(j=l0-1; j>=0; j--)
		{
			if(i==l1-1 && j==l0-1)
			{
				if(Egap_flag)
				{
					rF[1][l0]=rM[l1][l0]*exp(BETA*m.endgaps)/sk;
					/*  printf("%d,%d rF=%.4e  ",l1-1,l0,rF[1][l0]); */
				}
				else
				{
					/* Open in s1 */
					rF[1][l0]=rM[l1][l0]*exp(BETA*m.mat[l-2][s1[i]])/sk;
				}
			}
			if(i<l1-1 && j==l0-1)
			{
				if(Egap_flag)
				{
					rF[1][l0]=rF[0][l0]*exp(BETA*m.endgaps)/sk;
					/*  printf("%d,%d rF=%.4e  ",i,l0,rF[1][l0]); */
				}
				else
				{
					/* Extent in s1 */
					rF[1][l0]=rF[0][l0]*exp(BETA*m.mat[l-1][s1[i]])/sk;
				}
			}

			if(Egap_flag && j==0 && i>0)
			{
				rF[1][j]=rF[0][j]*exp(BETA*m.endgaps);
				rF[1][j]+=rM[i+1][j]*exp(BETA*m.endgaps);
				rF[1][j]/=sk;
				/*  printf("%d,%d rF=%.4e  ",i,j,rF[1][j]); */

				rE[1][j]=rE[1][j+1]*exp(BETA*m.mat[l-1][s0[j]]);/* ext in s0 */
				/* o in s0 */
				rE[1][j]+=(rM[i][j+1]+rF[1][j+1])*exp(BETA*m.mat[l-2][s0[j]]);
				rE[1][j]/=sk;
			}
			else if(Egap_flag && i==0 && j>0)
			{
				rE[1][j]=rE[1][j+1]*exp(BETA*m.endgaps);
				rE[1][j]+=(rM[i][j+1]+rF[1][j+1])*exp(BETA*m.endgaps);
				rE[1][j]/=sk;

				rF[1][j]=rF[0][j]*exp(BETA*m.mat[l-1][s1[i]]);/* ext in s1 */
				rF[1][j]+=rM[i+1][j]*exp(BETA*m.mat[l-2][s1[i]]);/* o in s1 */
				rF[1][j]/=sk;
			}
			else if(Egap_flag && i==0 && j==0)
			{
				rF[1][j]=rF[0][j]*exp(BETA*m.endgaps);
				rF[1][j]+=rM[i+1][j]*exp(BETA*m.endgaps);
				rF[1][j]/=sk;

				rE[1][j]=rE[1][j+1]*exp(BETA*m.endgaps);
				rE[1][j]+=(rM[i][j+1]+rF[1][j+1])*exp(BETA*m.endgaps);
				rE[1][j]/=sk;
			}
			else
			{
				/* rF(i,j) = ( rM(i+1,j)*exp(BETA*Open) +
	                           rF(i+1,j)*exp(BETA*Extend) ) / sk; */
				rF[1][j]=rF[0][j]*exp(BETA*m.mat[l-1][s1[i]]);/* ext in s1 */
				rF[1][j]+=rM[i+1][j]*exp(BETA*m.mat[l-2][s1[i]]);/* o in s1 */
				rF[1][j]/=sk;

				/* rE(i,j) =( (rM(i,j+1)+rF(i,j+1))*exp(BETA*Open) +
	                     	   rE(i,j+1)*exp(BETA*Extend) ) / sk; */
				rE[1][j]=rE[1][j+1]*exp(BETA*m.mat[l-1][s0[j]]);/* ext in s0 */
				/* o in s0 */
				rE[1][j]+=(rM[i][j+1]+rF[1][j+1])*exp(BETA*m.mat[l-2][s0[j]]);
				rE[1][j]/=sk;
			}

			/* rM(i,j)= ( exp(BETA*match or mismatch)*
	                     (rM(i+1,j+1)+rE(i+1,j+1)+rF(i+1,j+1) ) /sk*sk; */
			/* replace m.mat[s0[j]][s1[i]] */
	    	x = baseindex[ so[0].seq[ j ] ];
	    	y = baseindex[ so[1].seq[ i ] ];

	    	if( mutation_rates_flag )
			   	sim = basemutationrates[ x ][ y ] *
			       	  ( zeta * ( 1 - abs( probSgl_long[ j ] - probSgl_short[ i ] ) ) +
			          nzeta * ( 1 - abs( probDblUp_long[ j ] - probDblUp_short[ i ] ) ) );
			else
			   	sim = theta * basesub[ x ][ y ] +
				      ntheta * ( 1 - abs( probSgl_long[ j ] - probSgl_short[ i ] ) );
			rM[i][j] = exp(BETA*sim);
			rM[i][j] *= (rM[i+1][j+1]+rF[0][j+1]+rE[0][j+1]);
			rM[i][j] /= sk*sk;


			if(j == 0)
			{
				int k;
				for(k=l0; k>=0; k--)
				{
					rF[0][k]=rF[1][k];
					rF[1][k]=0;

					rE[0][k]=rE[1][k];
					rE[1][k]=0;
				}
			}
		}
    }

	Zr=rM[0][0]+rE[0][0]+rF[0][0];
	rM[l1][l0]=Zr;

	/*  printf("\n\n revers_partf:Zr;%e\n\n",Zr);  */
	/*    print_Zrevers(os,rM);  */

	/*  for(i=0; i<2; i++) */
	/*      { */
	/*        free(rE[i]); */
	/*        free(rF[i]); */
	/*      } */
	/*    free(rE); */
	/*    free(rF); */
	pf.M=rM;
	pf.E=rE;
	pf.F=rF;

	return(pf);
}


/* adaptation of probA-function <stoch_backtr> for DotAligner */
aligm stoch_backtr_da( aligm alig, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short )
{
	aligm st;

	st = stoch_btr_da( pmat, alig, probSgl_long, probDblUp_long, probSgl_short, probDblUp_short );

	return(st);
}


/* adaptation of probA-function <stoch_btr> for DotAligner */
aligm stoch_btr_da( matrix pz, aligm a, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short )
{
	int i,j,t;
	int l0,l1,l,Eg;
	int *s0, *s1;
	double s,opt,sk;
	double vrand, score;
	double x,y,z,Zx,bla;
	real pZ;
	real **M;
	real **E;
	real **F;
	aligm sta;
	char *tr;
	u_sc m;
	int opt_l;
	double s_p;

    m=scmat;

    opt_l=strlen(a.a)+20;

	s0=si[0].s;
	s1=si[1].s;
	l0=strlen(so[0].seq);
	l1= strlen(so[1].seq);
	l=strlen(m.monomers);

	/* scaling */
	opt=a.score;
	s=2.0*opt/(l0+l1);
	sk=exp(BETA*s/2);

	Eg = Egap_flag;

	pZ=pz.M[0][0]; /* partition function over all alignments */

	if(stochswi == 0)
	{
		double tmp;
	    double prob_o; /* probability of the (or better an) optimal alignment */

	    s_p = calc_score_da(a.a, probSgl_long, probDblUp_long, probSgl_short, probDblUp_short);

	    if(!(fabs(a.score - s_p) < 100*FLT_EPSILON) )
		{
	    	printf("Error in backtracking\n");
	    	exit(44);
		}

	    tmp = (BETA*a.score) - (l0+l1)*log(sk);
	    prob_o = exp(tmp) /pZ;

	    sta.a = a.a;
	    sta.score = a.score;
	    sta.prob = prob_o;

	    stochswi=1;

	    return sta;
	}

	M=pz.M;
	E=pz.E;
	F=pz.F;

	if(!sran)
	{
		//sgenrand((unsigned long) time(NULL));
		srand((unsigned long) time(NULL));
	    sran=1;
	}

	tr=(char*) calloc((l0+l1+1),sizeof(char));

	/* Backtracking */
	//rand=genrand();
	vrand=(double) rand() / RAND_MAX;

	Zx=M[l1][l0]+E[l1][l0]+F[l1][l0];

	x=M[l1][l0]/Zx;
	y=E[l1][l0]/Zx;
	z=F[l1][l0]/Zx;

	t=l0+l1;s0=si[0].s;
	s1=si[1].s;
	l0=strlen(so[0].seq);
	l1= strlen(so[1].seq);
	l=strlen(m.monomers);
	i=l1;
	j=l0;
	score=1.0;

	while((i>0 && j>=0) || (i>=0 && j>0))
    {
		//cout << vrand << "\t" << x << "\t" << y << "\t" << z << endl;
		if(vrand < x  )
		{
			/* M => match = 1 */
			if(s0[j-1] == s1[i-1])
		    {
				tr[t--]=M1; /* 1 match */
		    }
			else
		    {
				tr[t--]=M2; /* 2 mismatch */
		    }

			score*=x;

		    --i;
		    --j;

		    Zx=M[i][j]+E[i][j]+F[i][j];

		    x=M[i][j]/Zx;
		    y=E[i][j]/Zx;
		    z=F[i][j]/Zx;
		}
		else if(vrand <  x+y )
		{
			/* E => gap in s1(shorter seq) = 3 */
			/* m.mat[l-1][s0[k]];  Extend in s0 (3) */
			/* m.mat[l-2][s0[k]];  Open in s0 (3) */

			tr[t--]=G3; /* 3 gaps in the shorter seq(s1,E) */

			score*=y;

			Zx = E[i][j]*sk;
			j--;
			if(i==0)
		    {
				x=0.0;
				y=1.0;
				z=0.0;
		    }
			else
		    {
				if(Eg && i == l1)
				{
					x = M[i][j]*exp(BETA*m.endgaps) / Zx;
					y = E[i][j]*exp(BETA*m.endgaps) / Zx;
					z=0.0;
				}
				else
				{
					x = M[i][j]*exp(BETA*m.mat[l-2][s1[i-1]]) / Zx;/* Open */
					y = E[i][j]*exp(BETA*m.mat[l-1][s1[i-1]]) / Zx;/* Ext */
					z=0.0;
				}
		    }
		}
	    else
		{
	    	/* F => gap in s0(longer seq) = 3*/
	    	/* m.mat[l-2][s1[j]]; Open in s1 (4) */
	    	/* m.mat[l-1][s1[j]]; Extend in s1 (4) */

	    	tr[t--]=G4; /* gaps in s0(longer seq) = 4 */

	    	score*=z;

	    	Zx=F[i][j]*sk;
	    	--i;
	    	if(j==0)
	    	{
	    		x=0.0;
	    		y=0.0;
	    		z=1.0;
	    	}
	    	else
	    	{
	    		if(Eg && j == l0)
	    		{
	    			x=M[i][j]*exp(BETA*m.endgaps) / Zx;
	    			y=E[i][j]*exp(BETA*m.endgaps) / Zx;
	    			z=F[i][j]*exp(BETA*m.endgaps) / Zx;
	    		}
	    		else
	    		{
	    			x=M[i][j]*exp(BETA*m.mat[l-2][s0[j-1]]) / Zx;/* Open */
	    			y=E[i][j]*exp(BETA*m.mat[l-2][s0[j-1]]) / Zx;/* Open */
	    			z=F[i][j]*exp(BETA*m.mat[l-1][s0[j-1]]) / Zx;/* Extend */
	    		}
	    	}
		}
	    //rand=genrand();
		vrand=(double) rand() / RAND_MAX;
    }

	i=0;
    while(tr[i++]==0);
    i-=1;

    t=l0+l1+1-i;

	sta.a= (char*) calloc((t+1), sizeof(char));

	j=0;
	for(; i<=l0+l1; i++)
	{
		sta.a[j]=tr[i];
		++j;
	}

	sta.score=0.0;
	sta.prob=0.0;

    /* fill in structure aligm */
    s_p=calc_score_da(sta.a, probSgl_long, probDblUp_long, probSgl_short, probDblUp_short);
    sta.score=s_p;

    bla = (BETA*sta.score) - (l0+l1)*log(sk);
	sta.prob = exp(bla)/pZ;

	sta.s0.name=strdup(so[0].name);
	sta.s0.seq=strdup(so[0].seq);
	sta.s1.name=strdup(so[1].name);
	sta.s1.seq=strdup(so[1].seq);

    free(tr);

    return(sta);
}


/* adaptation of probA-function <calc_score> for DotAligner */
double calc_score_da( char *tr, vector<double> & probSgl_long, vector<double> & probDblUp_long, vector<double> & probSgl_short, vector<double> & probDblUp_short )
{
	/* matrix Z contains the match probabilities */
	int i,j,t,k,g3,g4;
	int *s0, *s1;
	int l0,l1,l,Eg;
	u_sc m;
	double s; /* score of the stoch_aln */
	double result;

	double sim;
	double nzeta = 1-zeta;
	double ntheta = 1-theta;

	MapClass::CharIntMap baseindex = MapClass::init_base_index();
    int x, y;

	m=scmat;

	s0=si[0].s;
	s1=si[1].s;

	l0=strlen(so[0].seq);
	l1= strlen(so[1].seq);
	l=strlen(m.monomers);

	Eg = Egap_flag;

	t=strlen(tr);

	i=0;
	k=0; /* index for s0 */
	j=0; /* index for s1 */
	s=0.0;
	g3=0;
	g4=0;
	while( i < t)
    {
		if((tr[i] == M1) || (tr[i] == M2) )
		{
		    /* replace m.mat[s0[k]][s1[j]] */
	    	x = baseindex[ so[0].seq[ k ] ];
	    	y = baseindex[ so[1].seq[ j ] ];

		    if( mutation_rates_flag )
		    	sim = basemutationrates[ x ][ y ] *
	  		       	  ( zeta * ( 1 - abs( probSgl_long[ k ] - probSgl_short[ j ] ) ) +
	  		          nzeta * ( 1 - abs( probDblUp_long[ k ] - probDblUp_short[ j ] ) ) );
		    else
		    	sim = theta * basesub[ x ][ y ] +
	  			      ntheta * ( 1 - abs( probSgl_long[ k ] - probSgl_short[ j ] ) );

			g3=0;
			g4=0;
			/* this function was originally part of fkt stoch_btr !!*/
			/* sta.al.s0[i]=so[0].seq[k]; sta.al.s1[i]=so[1].seq[j]; */
			s+=sim;
			#if DEBUG
				printf("%c  %c  %7.6f  %i  ", so[0].seq[ k ], so[1].seq[ j ], basemutationrates[ x ][ y ], basesub[ x ][ y ]);
				printf("%c  %.3f  s=%.3f\n", tr[i], sim, s);
			#endif
			++k;
			++j;
		}

		if(tr[i] == G3) /* 3 gaps in the shorter seq(s1,E)*/
		{
			g4=0;
			/* sta.al.s0[i]=so[0].seq[k]; sta.al.s1[i]=GAP; */

			if(Eg && (j==0 || j==l1))
			{
				s+=m.endgaps;
				#if DEBUG
					printf("%c  %.3f  s=%.3f\n", tr[i], m.endgaps, s);
				#endif
			}
			else if(g3 == 0)
			{
				s+=m.mat[l-2][s0[k]];/* Open in s0 (3) */
				#if DEBUG
					printf("%c  %.3f  s=%.3f\n", tr[i], m.mat[l-2][s0[k]], s);
				#endif
			}
			else
			{
				s+=m.mat[l-1][s0[k]];/* Extend in s0 (3) */
				#if DEBUG
					printf("%c  %.3f  s=%.3f\n", tr[i], m.mat[l-1][s0[k]], s);
				#endif
			}
			++g3;
			++k;
		}

		if(tr[i] == G4) /* 4 gaps in the longer seq (s0,F) */
		{
			g3=0;
			/* sta.al.s0[i]=GAP;  sta.al.s1[i]=so[1].seq[j]; */

			if(Eg && (k==0 || k==l0))
			{
				s+=m.endgaps;
				#if DEBUG
					printf("%c  %.3f  s=%.3f\n", tr[i], m.endgaps, s);
				#endif
			}
			else if(g4 == 0)
			{
				s+=m.mat[l-2][s1[j]];/* Open in s1 (4) */
				#if DEBUG
					printf("%c  %.3f  s=%.3f\n", tr[i], m.mat[l-2][s1[j]], s);
				#endif
			}
			else
			{
				s+=m.mat[l-1][s1[j]];/* Extend in s1 (4) */
				#if DEBUG
					printf("%c  %.3f  s=%.3f\n", tr[i], m.mat[l-1][s1[j]], s);
				#endif
			}
			++g4;
			++j;
		}
		++i;

    }

	result = s;

	return (result);
}


/* adaptation of probA-function <free_partition_f> for DotAligner */
void free_partition_f_da(real **m, sequ *s)
{
	int i, l0, l1;

	l0=length_l(s);
	l1=length_s(s);

	/* allocated in prop_alig */
	for(i=0; i<l0+1; i++)
	{
		free(m[i]);
	}
	free(m);

	/* allocated in partf_da */
	for(i=0; i<l1+1; i++)
    {
		free(pmat.M[i]);
		free(pmat.E[i]);
		free(pmat.F[i]);
    }
	free(pmat.M);
	free(pmat.E);
	free(pmat.F);

	/* allocated in revers_partf_da */
	for(i=0; i<l1+1; i++)
    {
		free(rmat.M[i]);
    }

	for(i=0; i<2; i++)
    {
		free(rmat.E[i]);
		free(rmat.F[i]);
    }
	free(rmat.E);
	free(rmat.F);
	free(rmat.M);
}


/* adaptation of probA-function <free_align> for DotAligner */
void free_align_da(aligm alig)
{
	int i,l;

	/* allocated in input */
	for (i=0; i<2; i++) {
		free(si[i].name);
		free(si[i].s);
	}
	free(si);

	/* allocated in sequences */
	for (i=0; i<2; i++)
    {
		free(so[i].name);
		free(so[i].seq);
    }
	free(so);

	/* allocated in convert_matrix */
	l=(strlen(scmat.monomers));

	for(i=0; i<l; i++)
    {
		free(scmat.mat[i]);
    }
	free(scmat.mat);
	free(scmat.name);
	free(scmat.monomers);

    /* allocated in trace_back */
	free(alig.a);
	free(alig.s0.seq);
	free(alig.s0.name);
	free(alig.s1.seq);
	free(alig.s1.name);
}


/* adaptation of probA-function <decode_alig> for DotAligner */
void decode_alig_da(aligm ai)
{
	int i,j;
	int k,l;
	int la;      /*length of the alignment*/
	al a;
	char *buf0, *buf1, *bufa;

	buf0 = (char*) calloc((61),sizeof(char));
	buf1 = (char*) calloc((61),sizeof(char));
	bufa = (char*) calloc((61),sizeof(char));

	a.score = ai.score;

	a.s0.name = strdup(ai.s0.name);
	a.s1.name = strdup(ai.s1.name);

	/* length of the trace (alignment encodes in 1,2,3,4) array */
	j=strlen(ai.a);
	la=j;
	a.s0.seq = (char*) calloc((j+1),sizeof(char));
	a.s0.seq[j] = '\0';
	a.s1.seq = (char*) calloc((j+1),sizeof(char));
	a.s1.seq[j] = '\0';

	i=0;
	k=0; /* index for s0 */
	l=0; /* index for s1 */

	/* decode an alignment but how : see begin of this file*/
	while(i < j)
    {
		if( (ai.a[i] == M1) || (ai.a[i] == M2))
		{
			a.s0.seq[i] = ai.s0.seq[k++];
			a.s1.seq[i] = ai.s1.seq[l++];
		}

		if(ai.a[i] == G3) /* (3) gaps in the shorter seq(s1,E) */
		{
			a.s0.seq[i] = ai.s0.seq[k++];
			a.s1.seq[i] = GAP;
		}

		if(ai.a[i] == G4) /* (4) gaps in the longer seq (s0,F) */
		{
			a.s0.seq[i] = GAP;
			a.s1.seq[i] = ai.s1.seq[l++];
		}

		++i;
    }

	/* determine the length of the longer sequence name */
	k=strlen(a.s0.name);
	l=strlen(a.s1.name);

	if(k >= l)
		l=k;

	i = la/60;
	k=0;

	if( verbose_flag )
	{
		for(j=0; j <= i; j++)
		{
			if(j < i)
			{
				strncpy(buf0, &a.s0.seq[k], 60);
				strncpy(buf1, &a.s1.seq[k], 60);
				strncpy(bufa, &ai.a[k], 60);
				k+=60;
				printf(" %s\n", buf0);
				printf(" %s\n", bufa);
				printf(" %s\n", buf1);
			}
			if(i == j)
			{
				printf("%s\n", &a.s0.seq[k]);
				printf("%s\n", &ai.a[k]);
				printf("%s\n", &a.s1.seq[k]);
			}
		}
	}
	else
	{
		printf("\t%s", &a.s0.seq[k]);
	    printf("\t%s", &ai.a[k]);
	    printf("\t%s", &a.s1.seq[k]);
	    printf("\n");
	}

	free(buf0);
	free(buf1);
	free(bufa);
	free_al(a);
}


/* calculate sequence identity from alignment code */
double seqidentity(char * aln )
{
	double seqid = 0;
	int i = 0, l = 0, m = 0, len1 = 0, len2 = 0;

	while( i < strlen(aln) )
	{
		switch( aln[i] )
		{
			case '|': /* match */
				      l++;
				      m++;
				      len1++;
				      len2++;
				      break;
			case ':': /* mismatch */
				      l++;
				      len1++;
				      len2++;
				      break;
			case '^': /* gaps in longer sequence */
		              l++;
				      len2++;
		              break;
			case '.': /* gaps in shorter sequence */
				      l++;
				      len1++;
				      break;
			default:  abort();
		}
		i++;
	}

	int minl = ( len1 <= len2 ) ? len1 : len2;
	seqid = (double) 100 * m / minl;
	//printf("%.8f\t%i\t%i\n", si, m, minl);

	return seqid;
}
