import java.io.*; 
import java.lang.Math ; 
public class ps2pp {
    public static void main(String[] Args) throws Exception  {
        for (int i = 0 ; i < Args.length ; i ++ ){
            BufferedReader ReadFile = new BufferedReader(new FileReader( Args[ i ] ));
            String Line = "" ; 
            String Sequence = "" ; 
            double [][] PP = new double [1][1]; 
            read : while ( (Line = ReadFile.readLine()) != null )  {
                
                if (Line.length() > 8 && Line.substring(0,9).equals("/sequence")) {
                    Sequence = ReadFile.readLine() ;
                    Sequence = Sequence.substring(0, Sequence.length()-1 ); 
                    while ( (Line = ReadFile.readLine()).charAt(0) != ')' ) {
                        Sequence = Sequence + Line ;
                        Sequence = Sequence.substring(0, Sequence.length()-1 ) ; 
                    }
                    System.out.println( Sequence.length() ); 
                    System.out.println( Sequence.toUpperCase() );
                    PP = new double[ Sequence.length() ][ Sequence.length() ] ; 
                }
                else if ( Line.length() != 0 && Line.matches(".*ubox") && Line.charAt(0) != '%' && Line.charAt(0) != '/'){
                    String [] Fields = Line.split(" "); 
                    //System.out.println( Line + Sequence.length() ) ; 
                    int x = Integer.parseInt( Fields[ 0 ]); 
                    int y = Integer.parseInt( Fields[ 1 ]); 
                    PP [ x-1 ][ y-1 ] = Math.pow( Double.parseDouble( Fields[ 2 ] ), 2) ;
                    PP [ y-1 ][ x-1 ] = Math.pow( Double.parseDouble( Fields[ 2 ] ), 2) ; 
                }
            }
            ReadFile.close(); 

            double [] UP = new double [ PP.length ]; 
            for ( int x = 0 ; x != PP.length ; x++ ) {
                for ( int y = 0 ; y != PP.length ; y++ ) {
                    System.out.print( PP[ x ][ y ] +" " );
                    UP[x] = UP[x]+PP[ x ][ y ] ; 
                }
                System.out.println();
            }
            System.out.println();
            for ( int x = 0 ; x != UP.length ; x++ ) 
                System.out.print( 1-UP[ x ] +" " );
            System.out.println() ; 

        }
    }
}

