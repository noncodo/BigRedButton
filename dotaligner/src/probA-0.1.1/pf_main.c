#include "pfgoto.h"


/*extern aligm aligg; holds the alignment, which should be easily accessable */

int main(int argc, char **argv)
{
  int i,n;
  sequ *s;    /* sequ structure: contains 2*(char *name + char *sequence) */
  aligm aligg,sta;
  real **Mpr;
  
  
      
  start(argc,argv);
  
  s=input(argc,argv);

  aligg=align(s);

  
  Mpr=partition_f(aligg);
  
  if(PS_flag)
    ps_plot(Mpr,"aln_sub.eps");

  n=0;
  while(n++ < Nr)
    {
      sta=stoch_backtr(aligg);
      free_stoch(sta);
      /*  sleep(2); */ 
    }
  
  /*  if(calcS_flag){ free(track); } */ 
  
  free_partition_f(Mpr,s);
  
  free_align(aligg);
      
  for (i=0; i<2; i++) {
    free(s[i].name);
    free(s[i].seq);
  }
  free(s);

  return(0); /*  if everything is okay => Exit 0  */
}



























