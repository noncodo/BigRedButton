//
//  StockholmToFasta.java
//  
//
//  Created by Martin Smith on 7/05/10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//
import java.util.*;import java.io.*;import java.math.*;

public class StockholmToFasta {
	public static void main(String[] Args) throws IOException {
		// * * * *  ADD something to ouput RFAM family ? ? ? 
		boolean isPublished = false ; 
		scanFiles : for ( int arg = 0 ; arg != Args.length ; arg++ ) {
			isPublished = false ; 
			System.out.println("---> New file being scanned.... "+Args[arg] ) ;  
			BufferedReader ReadFile = new BufferedReader(new FileReader( Args[ arg ] )); 
			// strip header
			String Line = "" ; 
			Map SeqHash = new HashMap() ;
			SeqHash = new TreeMap() ;
			scanlines: while ( (Line = ReadFile.readLine()) != null) { 
				String[] Entry = Line.split( "\\s+") ; // make exception for consensus tags
				//System.out.println( Line ); 
				if ( Line.length() < 10 || Entry.length < 2) 
					continue scanlines ; 
				//if ( Line.charAt(0)== '#' && !Line.matches( "^.*SS.*Published.*$") )  // only scan published RFAM entries
				//	continue scanFiles ;	
				if ( Line.charAt(0)== '#' && !Line.matches( "^#=GC.*SS_cons\\s.*" )  ) 
					continue scanlines ;
				else {
					// multiple lines per title
					String TempSeq = (String) SeqHash.get( Entry[0]) ; 
					if ( TempSeq != null ) {
						// merge sequences
						TempSeq = (Entry.length == 3)? 	TempSeq + Entry[2] :
														TempSeq + Entry[1] ; 
						SeqHash.put( Entry[0], TempSeq ) ; 
					}
					else if ( Entry.length == 3) 
						SeqHash.put( Entry[1], Entry[2] ) ; 
					else 
						SeqHash.put( Entry[0], Entry[1] ) ; 

				}
 			} 
			// check consensus
			/*if( !SeqHash.containsKey("#=GC RF") )  // if no consensus, flush it. 
				continue scanFiles ;
			else { // get rid of consensus < 120 ?????? this is quite stringent
				String Ungapped =  ((String)SeqHash.get("#=GC RF")).replaceAll("[~\\.]", "") ; 
				//System.out.println( SeqHash.get("#=GC RF") +"\n" + Ungapped);
				if ( Ungapped.length() < 100 ) 
					continue scanFiles ;
				else 
					System.out.println ( Ungapped.length() ) ; 
			//}
			*/// write file out
			BufferedWriter Out = new BufferedWriter(new FileWriter( Args[arg]+".fasta"));
			Iterator iterator = SeqHash.entrySet().iterator();
			while( iterator.hasNext() ) {
				Map.Entry keyEntry = (Map.Entry) iterator.next();
				Out.write( ">"+keyEntry.getKey()+"\n" );
				Out.write( SeqHash.get(keyEntry.getKey())+"\n" );
			}
			Out.close(); 
		} 
	}
		
}
