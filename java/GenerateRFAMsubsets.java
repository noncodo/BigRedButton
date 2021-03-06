//  Created by Martin A. Smith on Dec 2013.
// 	Distributed under GPL
//	martinalexandersmith[at]gmail[dot]com

/*  TO DO 
	-Build treeMap for individual families. If less than X sampled, discard
	-Extract structure mask from consensus
*/
import java.util.*; import java.io.*; import java.math.*;

public class GenerateRFAMsubsets {

	static String	OUT_PATH = 			"./rfam_subset.fa", 
					DB_PATH = 			"./input/";

	static File		OUT_FILE = 			new File(OUT_PATH);

	static boolean 	VERBOSE = 			false ,
					STRIP = 			false,
					ALIGN = 			true; 

	static int	MIN_PER_FAM = 			2,
				MAX_PER_FAM = 			20,
				TOTAL_SEQS = 			200,
				MAX_WIN = 				175,
				MIN_WIN = 				125,
				MIN_PI = 				0,
				MAX_PI = 				50;
//				LEN_DIFF =				10;     

/*******************************************************************************
*** 						 GATHER INPUT 
********************************************************************************/
	public static void main(String[] Args) throws IOException {
		String[] FileNames = new String [1] ; 
		File [] Files = new File [1] ; 
		//usage info
		if (Args.length == 0 ){
			System.out.println(
				"\n***************************************************************************************\n" +
				"***   GenerateRFAMsubsets\n"+
				"***   Selects random (sub)sequences from a reference database \n\n"+
			    "Usage:  java GenerateRFAMsubsets [options] -i [path to .fasta database folder]\n"+
			    "\t*** INPUT folder must contain only fasta files***\n"+
			    "Options:\n"+
   				"  -v            \tVerbose output\n"+
			    "  -min_win (int)\tminimum sequence/window length (default 125)\n"+
			    "  -max_win (int)\tmaximum sequence/window length (default 175)\n"+
				"  -min_pi (int) \tMinimum pairwise identity percentage (default 10)\n"+
				"  -max_pi (int) \tMaximum pairwise identity percentage (default 80)\n"+
				"  -max_f (int)  \tMaximum sequences per RFAM family/input file (default 20)\n"+
				"  -min_f (int)  \tMinimum sequences per RFAM family/input file (default 2)\n"+
				"  -t (int)      \tMaximum total sequences to include (default 250)\n"+
				"                \t**N.B.** There is no minimum value at the moment\n"+
				"  -l            \tUse Needleman-Wunsch with free end gaps to calculate pairwise identity.\n"+
				"                \tBy default, this is calculated from the intrinsic RFAM alignment.\n"+
				"  -strip        \tStrip gaps in output\n"+
//				"  -len   (int)\tMaximum percentage of length difference (default 10) \n"+
			    "  -o (string)   \twrite sampled sequences to file (default ./rfam_subset.fa) \n");
 			System.exit(0); 
		}
		// parse arguments
		flags: for (int i = 0 ; i != Args.length ; i++ ) {
			if ( Args[ i ].equals("-max_win") ) { 
				MAX_WIN = Integer.parseInt( Args[ i+1 ] ) ;
				i++ ; 
				continue flags; 
			}
			else if ( Args[ i ].equals("-min_win") ) { 
				MIN_WIN = Integer.parseInt( Args[ i+1 ] ) ;
				i++ ; 
				continue flags; 
			}
			else if ( Args[ i ].equals("-min_pi") ) {  
				MIN_PI = Integer.parseInt( Args[ i+1 ] ) ;
				i++ ; 
				continue flags; 
			}
			else if ( Args[ i ].equals("-max_pi") ) {  
				MAX_PI = Integer.parseInt( Args[ i+1 ] ) ;
				i++ ; 
				continue flags; 
			}
			else if ( Args[ i ].equals("-max_f") ) {  
				MAX_PER_FAM = Integer.parseInt( Args[ i+1 ] ) ;
				i++ ; 
				continue flags; 
			}
			else if ( Args[ i ].equals("-strip") ) {  
				STRIP = true ;
				continue flags; 
			}
			else if ( Args[ i ].equals("-min_f") ) {  
				MIN_PER_FAM = Integer.parseInt( Args[ i+1 ] ) ;
				i++ ; 
				continue flags; 
			}
			else if ( Args[ i ].equals("-t") ) {  
				TOTAL_SEQS = Integer.parseInt( Args[ i+1 ] ) ;
				i++ ; 
				continue flags; 
			}
			else if ( Args[ i ].equals("-l") ) {  
				ALIGN = true ; 
				// check to see if the NW binary was compiled and is present 
				// might not be necessary if bundled in a .jar 
				continue flags; 
			}
			else if ( Args[ i ].equals("-o") ) {  // out dir
				OUT_PATH = Args[ i+1 ] ;
				System.err.println(OUT_PATH) ; 
				OUT_FILE = new File(OUT_PATH) ; 
				if ( OUT_FILE.getParentFile() != null ) 
					OUT_FILE.getParentFile().mkdirs() ;	
				i++ ; 
				continue flags;
			}
			/*else if ( Args[ i ].equals("-len") ) {  // out dir
				LEN_DIFF = Integer.parseInt(Args[ i +1 ]); 
				i++;
				continue flags;
			}*/
			else if ( Args[ i ].equals("-v") ) {
				VERBOSE=true;
				continue flags ; 
			}
			else if ( Args[ i ].equals("-i") ) { 
				DB_PATH= Args[ i+1 ] +"/";
				File dir = new File(DB_PATH);
				Files = dir.listFiles(new FilenameFilter() {
    				@Override
    				public boolean accept(File dir, String name) {
        				return name.matches(".*\\.f.*a$");
    				}
				});
				FileNames = new String [ Files.length ] ; 
				// Files = dir.listFiles();
				if (Files == null) {
					System.out.println("Either dir does not exist or argument isin't a directory");
					System.exit(0); 
				}
				else {
					for (int j = 0 ; j < Files.length ; j++ ) {
						FileNames[j]=Files[j].toString();
					}
					i++ ;
				}	
			}
			else { 
				System.out.println("Couldn't parse argument: "+Args[ i ]+"\nPlease try again." );
				System.exit(0); 
			}

		}		
		/*******************************************************************************
		*** 						 Load data
		********************************************************************************/
		if (VERBOSE) // print out some basic information
			System.err.println(	"[ Checkpoint ] Loading data... "  ); 

		String [][] 	SeqNames = new String [ Files.length ][] ; 						// Array of sequence IDs
		String [][] 	SeqChars = new String [ Files.length ][] ;  					// Array of sequences 
		boolean [][] 	isWrongSize = new boolean [ Files.length ][  ] ,
						beenSampled = new boolean [ Files.length ][  ] ; 				// marker array

		int fams = 0 ; 
		for ( int i = 0 ; i != Files.length  ; i++  ) {
			if (VERBOSE)
				System.err.println("Parsing file: "+Files[i]);

			//read in individual files
			BufferedReader ReadFile = new BufferedReader(new FileReader( Files[i] )); 
			String Line = "" ; 
			Map RfamAln = new HashMap() ; 
			RfamAln = new TreeMap() ; 
			while ( (Line = ReadFile.readLine()) != null ) {
				RfamAln.put( Line , ReadFile.readLine().toUpperCase() ) ;  			
			}
			ReadFile.close(); 
			Object [] Names = RfamAln.keySet().toArray();								// returns an array of keys
			Object [] Sequences = RfamAln.values().toArray();							// returns an array of values

			SeqNames [ i ] = new String [ Names.length ] ; 								// Array of sequence IDs
			SeqChars [ i ] = new String [ Names.length ] ; 	 							// Array of sequences (as char[]) 

			//concatenate files into MULTIDIMENSIONAL arrays
			int ss = -1 ;
			for ( int j = 0 ; j != Names.length ; j++ ) {
				SeqNames[ i ][ j ] = Names[ j ].toString() ; 							// converts object to string
				if (SeqNames[ i ][ j ].matches(".*SS_cons") )
					ss = j ; 
				SeqChars[ i ][ j ] = Sequences[ j ].toString() ; 						
			}
			String Temp = SeqNames[ i ][ ss ];
			SeqNames[ i ][ ss ] = SeqNames[ i ][ 0 ] ;
			SeqNames[ i ][ 0 ] = Temp; 
			Temp = SeqChars[ i ][ ss ];
			SeqChars[ i ][ ss ] = SeqChars[ i ][ 0 ] ;
			SeqChars[ i ][ 0 ] = Temp;

			beenSampled[ i ] = new boolean [ Names.length ] ; 	
			isWrongSize[ i ]  = new boolean [ Names.length ] ; 							// implements 2nd dimension array
			RfamAln.clear() ; 	
		}
		/*******************************************************************************
		*** 						 START SAMPLING
		********************************************************************************/
		int	seqs = 0 ; 
		BufferedWriter SampleOutput = new BufferedWriter(new FileWriter( OUT_PATH )); 
		
		if (VERBOSE) // print out some basic information
			System.err.println(	"[ Checkpoint ] Sampling data... "  ); 
		// iterate through files non-randomly
		families: for ( int fam = 0 ; fam != Files.length ; fam++ ) {
			int fam_seqs = 0 ;
			Map SampledSeqs = new HashMap() ; 
			SampledSeqs = new TreeMap() ; 
			//RFAM family name
			String RFAM_ID = FileNames[ fam ].substring( FileNames[ fam ].lastIndexOf( "RF0" ), 
														 FileNames[ fam ].lastIndexOf("."));
			if (VERBOSE) // print out some basic information
				System.err.println(	"[ Checkpoint ] Sampling entry: "+ RFAM_ID  ); 

			//get a random sequence to initialise the sampling
			int random =  (int) (Math.random() * ((double) SeqNames[ fam ].length -1 ) + 1)    ; 
			int query_len = SeqChars[ fam ][ random ].replaceAll("-","").length();
			if ( query_len < MIN_WIN || query_len > MAX_WIN ) {
				int counter = 0 ; 
				while ( counter < 250 ) {
					random =  (int) (Math.random() * ((double) SeqNames[ fam ].length -1) + 1)   ; 
					query_len = SeqChars[ fam ][ random ].replaceAll("-","").length();
					if (  query_len >= MIN_WIN && query_len <= MAX_WIN )
						break; 
					counter++ ; 
				}
				System.err.println("[ \u001B[31m"+RFAM_ID+"\u001B[0m ] out of size range after 250 random samples");
						continue families; 
			}
			if (VERBOSE)
				System.err.println(	"[ Checkpoint ] Found size compatible random entry: "
									+ RFAM_ID +" "+SeqNames[ fam ][ random ] ); 
			// mark it as sampled
			beenSampled [ fam ][ random ] = true ; 
			beenSampled [ fam ][ 0 ] = true ; 

			fam_seqs++ ; 

			///////////////////////////////////////////////	
			// 	
			// 	Pseudo-random sampling 
			//	If a random choice has already been selected, 
			//  increments index until an unsampled sequence is found. 
			//
			int num_sampled = 1 ; // includes random 'seed' sequence
			scan : while ( 	num_sampled < SeqNames[ fam ].length && num_sampled <= MAX_PER_FAM  ) { 
				// get another random sample index
				int sample = (int) (Math.random() * ((double) SeqNames[ fam ].length -1 ) + 1);
				//make sure it's not the 'anchor' 
				if ( sample == random ) 
					continue scan; 
				if ( beenSampled [ fam ][ sample ] ) {
				//increment until an unsampled sequence is found
					if ( sample + 1 >= SeqNames[ fam ].length  ) {
						// Reaching end of array. Check if all members have been sampled.
						boolean isFull = false ; 
						for (int x = 0 ; x != beenSampled [ fam ].length ; x++ ) {
							if ( beenSampled [ fam ][ x ] == false ) {   		// found unsampled sequence
								isFull = false ; 
								sample = x ; 									// adjust sample index
								break; 
							}
						}
						if ( isFull ) {
							if (  VERBOSE ) 
								System.err.println("[ Warning ] Exhausted sequences. Moving to new entry"); 
							continue families ; 
						}
					}
					else 
						sample++ ; 
				}
				beenSampled [ fam ][ sample ] = true ;
				num_sampled++; 
				//length of sequence without gaps
				query_len = SeqChars[ fam ][ sample ].replaceAll("-","").length();

 				if ( query_len >= MIN_WIN && query_len <= MAX_WIN ) {
					double pid = 0; 
					double[] PIDtab = GetPID( SeqChars[ fam ][ random ], SeqChars[ fam ][ sample ], RFAM_ID) ;  
					pid = ( ALIGN )? PIDtab[1] * 100 : PIDtab[0] * 100 ;
					
					System.out.println( RFAM_ID +"\t"+ PIDtab[0]+"\t"+ PIDtab[1]); 

					// Ensure already sampled sequences don't clash (PIDs >> or << )
					if ( (int) pid > MAX_PI || (int) pid < MIN_PI ) {
					/*	if (  VERBOSE ) 
							System.err.println("[ Checkpoint ] PID out of range: rejecting sequence. "); 
					*/	continue scan ; 
					}
					if ( fam_seqs > 1 ) {
						if (VERBOSE) System.err.println(" -cross-checking sample #"+fam_seqs+" "+SeqNames[ fam ] [ sample ]);
				    	Collection C = SampledSeqs.values();  
    					Iterator<String> ITR = C.iterator();
    					while(ITR.hasNext()) {
    						PIDtab = GetPID( SeqChars[ fam ][ sample ], ITR.next() , RFAM_ID ) ; 
							pid = ( ALIGN )? PIDtab[1] * 100 : PIDtab[0] * 100 ;
      					 	if ( (int) pid > MAX_PI || (int) pid < MIN_PI )
								continue scan ; 
							
							System.out.println( RFAM_ID +"\t"+ PIDtab[0]+"\t"+ PIDtab[1]); 

						}
    				}
					if (VERBOSE)	 												
						// print out some basic information
						System.err.println(	"[ NOTE ] adding :\n\t"+ RFAM_ID +"\t"+ SeqNames[ fam ] [ sample ] 
								+"\t"+ (int) pid +"\t"+ query_len +"\t"+ fam_seqs ); 
					// Looks good, add it to output hashmap
					for (int x = 0 ; x != SeqNames[ fam ].length ; x++ )
					SampledSeqs.put( ">"+RFAM_ID+"_"
										+SeqNames[ fam ] [ sample ].substring(1, SeqNames[ fam ] [ sample ].lastIndexOf(".") )
										+"_"+(fam_seqs+1) +"_"+ (int) pid , SeqChars[ fam ][ sample ]) 	;	
					fam_seqs++ ; 
				}
				if ( fam_seqs == MAX_PER_FAM) {
					if (VERBOSE)
						System.err.println( "[ NOTE ] Reached maximum seqs per family: "+fam_seqs);
					break ; 
				}
			}
			if ( fam_seqs >1 ){
				if (VERBOSE) 												// print out some basic information
					System.err.println(	"[ NOTE ] added Random Initial selection : \n\t"+SeqNames[ fam ] [ random ]);
				SampledSeqs.put( ">"+RFAM_ID+"_"
								+SeqNames[ fam ] [ random ].substring(1, SeqNames[ fam ] [ random ].lastIndexOf(".") )
								+"_1", SeqChars[ fam ][ random ]);
			}


			// write it out
			if ( (fam_seqs >= MIN_PER_FAM || fam_seqs == MAX_PER_FAM) && seqs <= TOTAL_SEQS - 2 ) {
				fams++ ; 
				Object [] Names = SampledSeqs.keySet().toArray();								// returns an array of keys
				Object [] Sequences = SampledSeqs.values().toArray();							// returns an array of values

				if ( fam_seqs != Names.length )
					System.err.println( "[ \u001B[31mERROR\u001B[0m ] # sampled sequences != hash map size (fam_seqs: "
										+ fam_seqs +" Hash.length: "+ Names.length);
				// print out all sampled names and sequences 
				for ( int x = 0 ; x != Names.length ; x++ ) {
					if ( seqs == TOTAL_SEQS  ) {
						System.err.println(	"[ \u001B[32m"+RFAM_ID+"\u001B[0m ] Hit MAX_SEQS threshold. Sampled " 
											+ (seqs + x -1)  + " to file "+OUT_PATH ); 
						SampleOutput.close(); 
						System.exit(0) ;		
					}
					else 
						seqs++;
					// Write name out
					if (VERBOSE) 												// print out some basic information
								System.err.println(	"[ NOTE ] writing to file.... ");
					SampleOutput.write( Names[ x ].toString() + "\n");
					String Aligned = Sequences[ x ].toString() ; 
					
					// Strip gaps form consensus structure
					if (STRIP) {  // take it off
						String SeqStripped = new String() ;
						String StructureStrip = new String() ; 
						for ( int y = 0 ; y != Aligned.length() ; y ++ ) {
							 if ( Aligned.charAt( y ) != '-') {
							 	SeqStripped = SeqStripped + Aligned.charAt( y ); 
							 	StructureStrip = StructureStrip + SeqChars[ fam ][ 0 ].charAt( y ) ;			 						 
							 }
						}
						SampleOutput.write( SeqStripped +"\n" + 
							StructureStrip.replaceAll("<","(").replaceAll(">",")").replaceAll("[^<>()]","\\.") + "\n"  ); 
					}
					else {
						String Consensus = SeqChars[ fam ][ 0 ].toString().replaceAll("<","(").replaceAll(">",")").replaceAll("[^<>()]","\\.") ; 
						SampleOutput.write( Aligned +"\n" + Consensus  + "\n"  ); 
					}
				}
				System.err.println(	"[ \u001B[32m"+RFAM_ID+"\u001B[0m ] Added " + Names.length + " sequences from "+ RFAM_ID  ); 
				//seqs = seqs + Names.length ; 
			} 
			else if (VERBOSE)
				System.err.println( "[ WARNING ] Insufficient sequences (<"+fam_seqs+") sampled for "+ RFAM_ID );
			else 
				System.err.println(	"[ \u001B[31m"+RFAM_ID+"\u001B[0m ] less than "+MIN_PER_FAM+" sequences sampled "  ); 

			if ( seqs >= TOTAL_SEQS ) {
				System.err.println(	"[ \u001B[32m"+RFAM_ID+"\u001B[0m ] Hit max total sequences. Sampled " + seqs + " to file "+ OUT_PATH    ); 
				SampleOutput.close(); 
				break;
			}
			SampledSeqs.clear(); 
		}
		System.err.println(	"[ DONE ] Sampled " + seqs + " sequences from "+fams+" familes to file "+OUT_PATH    );
		SampleOutput.close(); 
	} 

	//####################################################################################################
	//     method to calculate pairwise indentity of 2 aligned sequences
	//####################################################################################################
	public static double[] GetPID ( String Seq1, String Seq2, String RFAM) throws IOException {
		int 	matches = 		0,
				seq1_len = 		0,
				seq2_len = 		0;
		double 	Pids [] = new double [2] ;

		if ( Seq1.length() != Seq2.length()) {
			System.err.println("[ \u001B[31mERROR\u001B[0m ] Sequences have uneven length! Check input:\n"+RFAM +"\n"+Seq1 +"\n"+Seq2); 
			System.exit(0); 
		}
		for (int col = 0 ; col != Seq1.length() ; col ++ ) {
			if ( Seq1.substring(col, col+1).matches("[a-zA-Z]") ) 
				seq1_len++ ; 
			if ( Seq2.substring(col, col+1).matches("[a-zA-Z]") ) 
				seq2_len++ ; 
			if ( Seq1.charAt( col ) == Seq2.charAt( col ) && Seq1.substring(col, col+1).matches("[a-zA-Z]") ) 
				matches++ ; 
		}
		Pids[0] = (double) matches  / Math.min( seq1_len, seq2_len ) ; 
		//System.out.println( RFAM +"\tRFAM\t"+ Pids[0] ); 
		if ( VERBOSE && matches > Math.min( seq1_len, seq2_len ) )
				System.err.println( "[ ERROR ] more RF matches than sequences\n"
									+ Seq1 + "\n"
									+ Seq2 + "\n"
									+ matches +" "
									+ Math.min( seq1_len, seq2_len )); 
			String StripSeq1 = Seq1.replaceAll("-","");
			String StripSeq2 = Seq2.replaceAll("-","");
			// assumes "align" binary is in appropriate relative path 
			String Command = "align -f "+ StripSeq1 + " " + StripSeq2 ; 
				/**** Typical output ******
				 ./align -l GCTTTTT CGTTTTT
				 Score = 10
				 Length = 8
				 Match = 6
				 SI = 75
				 GC-TTTTT
				 -CGTTTTT 
				**************************/	
			String Output = "" ; 
			Process NW ; 
			try {
				NW = Runtime.getRuntime().exec( Command );
				BufferedReader In = new BufferedReader(new InputStreamReader(NW.getInputStream()));
				int lineNo = 1 ; 
				while ( ( Output = In.readLine() ) != null ) {
					if ( lineNo == 3 ) {
				  		matches = Integer.parseInt( Output.split("\\s")[2] ) ;  
				  		break; 
				  	}
				  	lineNo++;
				}
			} catch (Exception OhNo) {
				System.err.println( "[ \u001B[31mERROR\u001B[0m ] Cannot execute pairwise alignment command:\n"+Command);
				OhNo.printStackTrace();
			}
			if ( VERBOSE && matches > Math.min( seq1_len, seq2_len ) )
				System.err.println( "[ ERROR ] more SW matches than sequences\n"
									+ Seq1 + "\n"
									+ Seq2 + "\n"
									+ matches +" "
									+ Math.min( seq1_len, seq2_len )); 
		Pids[1] = (double)matches / Math.min( seq1_len, seq2_len ) ; 
		// System.out.println( RFAM +"\tRFAM\t"+ Pids[1] ); 
		return Pids;
	}
			
}





