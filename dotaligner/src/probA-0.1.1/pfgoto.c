#include "pfgoto.h"

#define MAX(a,b) (( (a) > (b) ) ? (a) : (b) )
#define MIN(a,b) (( (a) < (b) ) ? (a) : (b) )
/* encode the alignment: the longer sequence is always written above the
   shorter one, each state of the alignment (match, mismatch, indel) is
   symbolized a self-evident char.*/
#define M1 '|' /*  match */
#define M2 ':' /*  mismatch */
#define G3 '.' /*  gap in the lower sequence (the shorter one)*/
#define G4 '^' /*  gap in the upper sequence (the longer one)*/

static sequ *so; /* array with 2 members (longer seq + name is in so[0]) */
static iseq *si; /* like '*so' but sequences are encoded as strings of digits
		    see function code sequence */
static u_sc scmat; /* scoring matrix */
static matrix pmat,rmat; /* partition function and revers pf */
static real **Pr; /* match probability matrix */
static int sran=0; /* for initialisation of of the "Mersenne Twister"*/
static int betaswi=0,stochswi=0;
static double beta=0;
static float T;
char *control_AA = "ABCDEFGHIKLMNPQRSTVWXYZOJ";
char *control_NA = "ABCDGHKMNRSTUVWYOJ";

void print_P(sequ *s, real **matr);
void print_Z(sequ *s, real **matr);
void print_Zrevers(sequ *s, real **matr);
void print_M(sequ *os, real **matr);
void print_smat(u_sc smat, score_m sm);
aligm trace_back(matrix m, sequ *os, iseq *is, u_sc smat);
aligm add_sequences(sequ *os, aligm a);
/*  void print_cod(char *tra, int ltra, Cod Caln); */
/*  int comp_score(const void *s1, const void *s2); */ 
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* determines the length of the longer sequence */
int length_l(sequ *seq)
{
  int l0,l1;

  l0 = strlen(seq[0].seq);
  l1 = strlen(seq[1].seq);

  if(l0 > l1)
    {
      return(l0);
    }
  else
    {
      return(l1);
    }
}
/*---------------------------------------------------------------------------*/
/* determines the length of the shorter sequence */
int length_s(sequ *seq)
{
  int l0,l1;

  l0 = strlen(seq[0].seq);
  l1 = strlen(seq[1].seq);

  if(l0 < l1)
    {
      return(l0);
    }
  else
    {
      return(l1);
    }
}


/*---------------------------------------------------------------------------*/
      /*------------------------ alignment functions -------------------*/
/*---------------------------------------------------------------------------*/
/* this function converts an alignment encoded as a character string to an
   alignment in ClustalW format*/
al decode_alig(aligm ai)
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

  /* decode an alignment but how : see beginn of this file*/
  

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

  for(j=0; j <= i; j++)
    {
     
      if(j<i)
	{
	  strncpy(buf0, &a.s0.seq[k], 60);
	  strncpy(buf1, &a.s1.seq[k], 60);
	  strncpy(bufa, &ai.a[k], 60);
	  k+=60;
	  /*  buf0[k]='\0'; */
/*  	  buf1[k]='\0'; */
          
	  /*printf("#%-10s     %s\n",a.s0.name,buf0);
	  printf("#               %s\n",bufa);
	  printf("#%-10s     %s\n",a.s1.name,buf1);
	  printf("#\n");
	  printf("#\n");*/
	  	  printf(" %s\n", buf0);
	  	  printf(" %s\n", bufa);
	  	  printf(" %s\n", buf1);
	}
      if(i == j)
	{
	  /*printf("#%-10s     %s\n",a.s0.name,&a.s0.seq[k]);
	  printf("#               %s\n",&ai.a[k]);
	  printf("#%-10s     %s\n",a.s1.name,&a.s1.seq[k]);
	  printf("#\n");
      printf("#\n");*/
    	  printf("%s\n", &a.s0.seq[k]);
    	  printf("%s\n", &ai.a[k]);
    	  printf("%s\n", &a.s1.seq[k]);
	}
    }

  free(buf0);
  free(buf1);

  return(a);

  
}
/*---------------------------------------------------------------------------*/
score_m scm_name(score_m **matx)
{
  static char *mat[] = {"gon", "pam", "blo"};
  score_m scores;
  

  /* blosum series: 30,50,62,80 */
  if(strncasecmp(MAT_SER,mat[2],3) == 0)
    {

      if(DISTANCE <= 40.0)
        {
          scores=*matx[0];
        }

      if(40.0 < DISTANCE && DISTANCE <= 56.0)
        {
          scores=*matx[1];
        }
      
      if(56.0 < DISTANCE && DISTANCE <= 70.0)
        {
          scores=*matx[2];
        }
      /*  printf("select_scm: DISTANCE=%f\n",DISTANCE); */
      if(DISTANCE > 70.0)
        {
          scores=*matx[3];
          /*  printf("score matrix:%s\n",matx[0].name); */
        }
    }

  /* pam series: 20,60,120,350 */
  else if(strncasecmp(MAT_SER,mat[1],3) == 0)
    {
      /*  printf("pam distance:%f\n",DISTANCE); */
      
      if(DISTANCE <= 40.0)
        scores=*matx[0];

      if(40.0 < DISTANCE && DISTANCE <= 90.0)
        scores=*matx[1];

      if(90.0 < DISTANCE && DISTANCE <= 235.0)
        scores=*matx[2];

      if(DISTANCE > 235.0)
        {
          scores=*matx[3];
          /*  printf("score matrix:%s\n",matx[3].name); */
        }
    }

  /* gonnet series: 40,80,120,160,250,300,350 */
  else if(strncasecmp(MAT_SER,mat[0],3) == 0)
    {
      /*  printf("pam distance:%f\n",DISTANCE); */

      if(DISTANCE <= 60.0)
        {
          scores=*matx[0];
          /*  printf("score matrix:%s\n",matx[0].name); */
        }
      
      if(60.0 < DISTANCE && DISTANCE <= 100.0)
        {
          scores=*matx[1];
          /*  printf("score matrix:%s\n",matx[1].name); */
        }
      
      if(100.0 < DISTANCE && DISTANCE <= 140.0)
        {
          scores=*matx[2];
          /*  printf("score matrix:%s\n",matx[2].name); */
        }
      
       if(140.0 < DISTANCE && DISTANCE <= 205.0)
         {
           scores=*matx[3];
           /*  printf("score matrix:%s\n",matx[3].name); */
         }
       if(205.0 < DISTANCE && DISTANCE <= 275.0)
         {
           scores=*matx[4];
           /*  printf("score matrix:%s\n",matx[4].name); */
         }
       
       if(275.0 < DISTANCE && DISTANCE <= 325.0)
         {
           scores=*matx[5];
           /*  printf("score matrix:%s\n",matx[5].name); */
         }

       if(DISTANCE > 325.0)
         {
           scores=*matx[6];
           /*  printf("score matrix:%s\n",matx[6].name); */
         }
    }
  else
  /*  if((strncasecmp(MAT_SER,mat[0],3)!=0) */
/*       &&(strncasecmp(MAT_SER,mat[1],3)!=0) */
/*       &&(strncasecmp(MAT_SER,mat[2],3)!=0)) */
    {
      fprintf(stderr,"select_score_mat:couldn't select a scoring matrix");
      usage(0);
      exit(3);
    }
  /*  printf("scm_name: matrix %s\n",scores.name); */

  return(scores);
    
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
void free_al(al a)
{
  /* allocated in decode_alig */
  free(a.s0.name);
  free(a.s0.seq);
  free(a.s1.name);
  free(a.s1.seq);
}
/*---------------------------------------------------------------------------*/
void free_align(aligm alig)
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

/*---------------------------------------------------------------------------*/
void free_aligm(aligm alig)
{
    free(alig.a);
    free(alig.s0.seq);
    free(alig.s0.name);
    free(alig.s1.seq);
    free(alig.s1.name);
}

/*---------------------------------------------------------------------------*/
aligm copy_aligm(aligm alig)
{
	aligm calig;

	calig.a = strdup(alig.a);
	calig.prob = alig.prob;
	calig.s0.name = strdup(alig.s0.name);
	calig.s0.seq = strdup(alig.s0.seq);
	calig.s1.name = strdup(alig.s1.name);
	calig.s1.seq = strdup(alig.s1.seq);
	calig.score = alig.score;

	return calig;
}
/*---------------------------------------------------------------------------*/
/* this function returns an optimal alignment between two sequences */
aligm align(sequ *seq_array)
{
  int i,l;
  int l0,l1;
  iseq *sti;
  u_sc stm;
  score_m **matser;
  score_m scoremat;
  aligm alig;
  
  so=sequences(seq_array);/* looks for the longer sequence and stores it in
			     so[0],the shorter one is stored in so[1],
			     both sequences are converted to uppercase
			     letters! */

  l0=strlen(so[0].seq);
  l1=strlen(so[1].seq);

  if(typ_flag==-1)
    {
      /* if the user didn't specify the kind of polymer used,
	 check_polymer sets the flag: Flag set by `-DNA'=1  or `-prot'=0 */
      typ_flag=check_polymer(so[0].seq,l0);
    }
  /* case protein typ_flag=0 */
  if((typ_flag==0) && (DISTANCE == -1))
    {
      /* for the initial aligment the score matrix gonnet_init is used,
	 MAT_SER decides wether a +/- or an all + matrix is used */
      stm=convert_matrix(gonnet_init);
      sti=code_seq(so,stm);
      matser=init_matseries();
      scoremat=select_scm(matser,so,sti,stm);
      
      /* free the initiation scoring matrix allocated in convert_matrix */ 
      l=(strlen(stm.monomers));
      for(i=0; i<l; i++)
	{
	  free(stm.mat[i]);
	}
      free(stm.mat);
      
      free(stm.name);
      free(stm.monomers);
      
      /* free the coded sequence allocated in code_seq */ 
      for (i=0; i<2; i++)
	{
	  free(sti[i].name);
	  free(sti[i].s);
	}
      free(sti);
      
      scmat=convert_matrix(scoremat);
      si=code_seq(so,scmat);
      }

  else if((typ_flag==0) && (DISTANCE >= 0))
    {
      matser=init_matseries();
      scoremat=scm_name(matser);
      scmat=convert_matrix(scoremat);
      si=code_seq(so,scmat);
    } 
  /* case nucleic acid typ_flag=1 */
  else if(typ_flag==1)
    {
      /*if(MAT_SER[0] != '\0')
	{
	  fprintf(stderr,"The input sequences are DNA therefor the program uses the default DNA scoring matrix \n");
	}*/
      MAT_SER[0]='D';MAT_SER[1]='N';MAT_SER[2]='A'; MAT_SER[3]='\0';
      scmat=convert_matrix(dnamt);
      si=code_seq(so,scmat);
    }
  else
    {
      fprintf(stderr,"Problem with score matrix selection\n");
      exit(11);
    }

  alig=align_2(so,si,scmat);

  alig=add_sequences(so,alig);

  /* allocated in init_matseries */ 
  free(matser);
   
  return(alig);
  
}
/*---------------------------------------------------------------------------*/
      /*------------------------ partition function -------------------*/
/*---------------------------------------------------------------------------*/
void free_partition_f(real **m, sequ *s)
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

   /* allocated in partf */
  
  for(i=0; i<l1+1; i++)
    {
      free(pmat.M[i]);
      free(pmat.E[i]);
      free(pmat.F[i]);
    }
  free(pmat.M);
  free(pmat.E);
  free(pmat.F);
    
  /* allocated in revers_partf */
  
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
/*---------------------------------------------------------------------------*/

/* partition_f calculates the partition function over all alignments of the
   two input sequences, the partition function is then used to determine the
   match probability for each possible match i,j between the two sequences. */

real **partition_f(aligm alig)
{
  /*  real **Pr; */
 
  pmat=partf(alig);
  
  rmat=revers_partf(so,si,scmat,alig);
  
  Pr=prop_alig(so,si,scmat,pmat.M,rmat.M,alig);
  
  return(Pr);

}

/*---------------------------------------------------------------------------*/
                 /*------------------dot plot-------------------*/
/*---------------------------------------------------------------------------*/
/* represents the match probabilities as a dot plot*/ 
void ps_plot(real **P, char *wastlfile)
{
    PS_plot(so, P, wastlfile);
}
  
/*---------------------------------------------------------------------------*/
      /*------------------ stochastic backtracking -------------------*/
/*---------------------------------------------------------------------------*/

aligm stoch_backtr(aligm alig2)
{
  aligm st;

  st=stoch_btr(pmat,alig2);

  return(st);

}
/*---------------------------------------------------------------------------*/
void free_stoch(aligm st)
{
   /* allocated in trace_back */

  free(st.a);

  free(st.s0.seq);
  free(st.s0.name);

  free(st.s1.seq);
  free(st.s1.name);

  
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* decides whether a all plus matrix or a +/- matrix is used */
int plus(void)
{
  int i, swi, lm;

  lm=strlen(MAT_SER);

  swi=0;
  
  for(i=0; i< lm; i++)
    {
      if(MAT_SER[i] == '_')
	{
	  if(MAT_SER[i+1] == 'p')
	    {
	      swi=1;  
	    }
	}
    }
  return(swi);
}

/*---------------------------------------------------------------------------*/
/* converts an N(N-1)/2 + N protein score matrix to an N*N score matrix */
u_sc convert_matrix(score_m sm)
{
  u_sc smat;
  int laa;
  int lname;
  int lsm,lm;
  int i,j;
  int sum;
  double min;
  int swi;    /* switch between +/- or a all + scoring matrix
		   swi=0 -> +/- matrix */
  swi=0;
  
  /* beta stores the original value of BETA */  
  if(!betaswi)
    {
      T=1/BETA;
      beta=BETA;
      ++betaswi;
    }
  
  laa=strlen(sm.monomers); /* lenght of the amino_acid_string */
  lname=strlen(sm.name);
  lm=strlen(MAT_SER);
  
  for(i=0; i<lm; i++)
    {
      tolower(MAT_SER[i]);
    }
  
  
  
  /* the gap penalties are added as the last but one (open) and last (extent)
     row of the score matrix array => laa + 2 */
  
  smat.name=(char*) calloc(20,sizeof(char));
  smat.monomers=(char*) calloc(30,sizeof(char));
  
  smat.mat=(double**) calloc((laa+2),sizeof(double*));
  for(i=0; i<laa+2; i++)
    {
      smat.mat[i]=(double*) calloc((laa+2),sizeof(double));
    }
 

  /* decide wether a +/- or a all + matrix is used */

  swi=plus();
  
  if(swi==1) /* all +  matrix */ 
    {      
 
      strcpy(smat.name, sm.name);
      
      /* add _p to the name of the matrix */ 
      smat.name[lname]='_';
      smat.name[lname+1]='p';
      smat.name[lname+2]='\0';
      
      smat.dist=sm.dist; 

      /* BETA contains the scoring matrix dependent partition modul (beta) and
	 a variable that influences the fidelity of the alignment (T):
	 
	 BETA = partition modul (beta) / T; */

      smat.beta=sm.beta; /* scoring matrix dependent partition modul (beta) */
      
      smat.beta*=beta;   /* BETA = 1/T */
      
      BETA=smat.beta;    /* => BETA = smat.beta/T */ 
      
      strcpy(smat.monomers, sm.monomers);
      
  /* add the gap penalties to the blosum 50 array: O for Open;
                                                   J for extent; */
      smat.monomers[laa]='O';
      smat.monomers[laa+1]='J';
      smat.monomers[laa+2]='\0';
      					   
  /* convert the  N(N-1)/2 + N protein score matrix to a N*N score matrix */

      lsm=laa*(laa+1)/2;

      min=sm.mat[0];
      for(i=1; i<lsm; i++)
	{
	  min=MIN(min,sm.mat[i]);
	}	  
      
      sum=0;
      i=0;
      while(i<laa)
	{
	  for(j=0; j<laa; j++)
	    {
	      if(i==j)
		{
		  smat.mat[i][j]=sm.mat[sum++]-min;
		  break;
		}
	      smat.mat[i][j]=smat.mat[j][i]=sm.mat[sum++]-min;
	    }
	  ++i;
	}

  /* add the gap penalties;  gap penalty = -open - (k-1) * extend;
     k length of the gap */

  /* open +/- */ 
      j=laa;
      for(i=0; i<laa; i++)
	{
	  smat.mat[i][j] = -sm.pos.o;
	}
      
      smat.mat[laa][j]=0.0;
      
      i=laa;
      for(j=0; j<laa; j++)
	{
	  smat.mat[i][j] = -sm.pos.o;
	}

  /* extend +/- */
  
      j=laa+1;
      for(i=0; i<laa+1; i++)
	{
	  smat.mat[i][j] = -sm.pos.e;
	}
      
      smat.mat[laa+1][j]=0.0;
      
      i=laa+1;
      for(j=0; j<laa+1; j++)
	{
	  smat.mat[i][j] = -sm.pos.e;
	}

      if(Egap_flag)
	{
	  smat.endgaps=ENDGAP;
	}
	  

     /*   print_smat(smat, sm); */   
    }
     /* matrix ready */

  
  if(swi==0) /* +/- matrix */ 
    {
      
      strcpy(smat.name, sm.name);
      
      smat.dist=sm.dist; 

      /* BETA contains the scoring matrix dependent partition modul (beta) and
	 a variable that influences the fidelity of the alignment (T):
	 
	 BETA = partition modul (beta) / T; */

      smat.beta=sm.beta; /* scoring matrix dependent partition modul (beta) */
      
      smat.beta*=beta;   /* BETA = 1/T */
      
      BETA=smat.beta;    /* => BETA = beta/T */ 

         
      strcpy(smat.monomers, sm.monomers);
  /* add the gap penalties to the blosum 50 array: O for Open;
                                                   J for extent; */
      smat.monomers[laa]='O';
      smat.monomers[laa+1]='J';
      smat.monomers[laa+2]='\0';
      
  /* convert the  N(N-1)/2 + N protein score matrix to a N*N score matrix */

      sum=0;
      i=0;
      while(i<laa)
	{
	  for(j=0; j<laa; j++)
	    {
	      if(i==j)
		{
		  smat.mat[i][j]=sm.mat[sum++];
		  break;
		}
	      smat.mat[i][j]=smat.mat[j][i]=sm.mat[sum++];
	    }
	  ++i;
	}

   
  /* add the gap penalties;  gap penalty = -open - (k-1) * extend;
     k length of the gap */

  /* open +/- */ 
      j=laa;
      for(i=0; i<laa; i++)
	{
	  smat.mat[i][j] = -sm.p_n.o;
	}
      
      smat.mat[laa][j]=0.0;
      
      i=laa;
      for(j=0; j<laa; j++)
	{
	  smat.mat[i][j] = -sm.p_n.o;
	}

  /* extend +/- */
  
      j=laa+1;
      for(i=0; i<laa+1; i++)
	{
	  smat.mat[i][j] = -sm.p_n.e;
	}
      
      smat.mat[laa+1][j]=0.0;
      
      i=laa+1;
      for(j=0; j<laa+1; j++)
	{
	  smat.mat[i][j] = -sm.p_n.e;
	}
      
      if(Egap_flag)
	{
	  smat.endgaps=ENDGAP;
	}

     /*   print_smat(smat, sm); */ 
      
    }
     /* matrix ready */

 
  
  
  
  /* build blosum 50 matrix as described in Henikoff, Curr Opin Struct Biol
     1996 Jun;6(3):353-60, Fig.3: add 3.4 to each entry and use these
     gap penalties: open -9.5; extent -0.6 */

  
  return(smat);

}

/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/* remember: Flag set by `-DNA'=1  or `-prot'=0;
   => ASCII Character set used !!!!!!!!!
   each character in the sequence is substituted by the index of this
   character in the corresponding monomer (amino acid or nucleic acid) array */
iseq *code_seq(sequ *os, u_sc s)
{
  int i,j;
  int l;
  int convert[127];
  int l0;
  int l1;
  iseq *sq;

  l=strlen(s.monomers);
  
  i=strlen(os[0].name);
  j=strlen(os[1].name);
  
  l0=strlen(os[0].seq);
  l1=strlen(os[1].seq);

  
  sq =(iseq*) calloc(2,sizeof(iseq));
  
 
  sq[0].name=(char*) calloc((i+1),sizeof(char));
  sq[0].s=(int*) calloc((l0+1),sizeof(int));
  sq[1].name=(char*) calloc((j+1),sizeof(char));
  sq[1].s=(int*) calloc((l1+1),sizeof(int));

  strcpy(sq[0].name, os[0].name);
  strcpy(sq[1].name, os[1].name);
  
  i=0;
  j=0;
  /* char array "convert": the position of the integer value of each character
     used in the monomer array contains the index of this character in the
     monomer array  e.g.:
     nucleic acid array "ABCDGHKMNRSTUVWYOJ"; A -> nucleic acid convert[0];
     ASCII;                                   A -> 65
                          => convert[65]=0;
  */
  

  for(i=0; i<127; i++)
    {
      convert[i]=-1;
    }
  for(j=0; j<l; j++)
    {
      for(i=0; i<127; i++)
	{
	  if(((int) s.monomers[j])==i)
	    {
	      convert[i]=j;
	      continue;
	    }
	}
    }

  /* O symbols gapopen, J symbols gapextent => both characters are not used
     for coding amino acids or nucleic acids =>
     convert[79 = O]=-1; convert[74 = J]-1; */

  convert[79]=-1;
  convert[74]=-1;
  
  /*  for(i=0; i<127; i++)
    {
      if(i%7 == 0)
	printf("\n");       
      printf("[%3d]=%3d ",i,convert[i]);
    } */
  /* convert the longer sequence */

  

  for(j=0; j<l0; j++)
    {
      for(i=0; i<127; i++)
	{
	  if(((int) os[0].seq[j])==i)
	    {
	      sq[0].s[j]=convert[i];
	    }
	}
    }

  
  /*printf("\n\n");
    printf("%s:",os[0].name); 
    for(i=0; i<l0; i++)
    {
    printf("%d ",sq[0].s[i]);
    }
    printf("\n%s",os[0].seq);
    printf("\n\n");*/ 

  /* convert the shorter sequence */

  for(j=0; j<l1; j++)
    {
      for(i=0; i<127; i++)
	{
	  if(((int) os[1].seq[j])==i)
	    {
	     sq[1].s[j]=convert[i];
	    }
	}
    }

  /*printf("%s:",os[1].name);
    for(i=0; i<l1; i++)
    {
    printf("%d ",sq[1].s[i]);
    }
    printf("\n%s",os[1].seq);
    printf("\n\n");*/

  /* check if there are characters not contained in monomers in the sequence */

  
  for(i=0; i<l0; i++)
    {
      if(sq[0].s[i] == -1)
	{
	 fprintf(stderr,"code_seq: unknown character in seq. %s, at position (%d)",os[0].name,(i+1));
	      exit(2); 
	}
    }

  for(i=0; i<l1; i++)
    {
      if(sq[1].s[i] == -1)
	{
	 fprintf(stderr,"code_seq: unknown character in seq. %s, at position (%d)",os[1].name,(i+1));
	      exit(2); 
	}
    }
 
  return(sq);

}

/*---------------------------------------------------------------------------*/
void print_smat(u_sc smat, score_m sm)
{
  int i,j;
  int sum;
  int l,ls;

  l=strlen(smat.monomers);
  ls=strlen(sm.monomers);

  ls+=2;
  if(l != ls)
    {
      fprintf(stderr,"print_smat:error in score matrix generation");
	      exit(2);
    }
  ls-=2;
  
  printf("%s\n",sm.name);
  printf("distance: %.1f\n",sm.dist);
  printf("beta: %.2f\n\n",sm.beta);
  printf("monomers: %s\n\n",sm.monomers);
  printf("all + matrix: open:%.2f, extend:%.2f  \n",sm.pos.o,sm.pos.e);
  printf("+/- matrix: open:%.2f, extend:%.2f  \n\n",sm.p_n.o,sm.p_n.e);
  
 
  sum=0;
  i=0; 
  while(i<ls)
    {
      for(j=0; j<ls; j++) 
	{
	  if(i==j)
	    {
	      printf("%5.1f,",sm.mat[sum++]);
	      break;
	    }
	  
	  printf("%5.1f,",sm.mat[sum++]);
	}
      printf("\n");    
	  ++i;
    }
  
  printf("\n\n");


  printf("scoring matrix: %s\n",smat.name);
  printf("distance: %.1f\n",smat.dist);
  printf("beta: %.2f\n\n",smat.beta);
  printf("monomers: %s\n\n",smat.monomers);
  
  for(i=0; i<l; i++)
    {
      for(j=0; j<l; j++)
	{
	  printf("%5.1f",smat.mat[i][j]);
	}
      printf("\n");
    }
  
 
  printf("\n\n");

}

/*---------------------------------------------------------------------------*/
 /* Select a scoring matrix : Flag set by `-DNA'=1  or `-prot'=0
    this function selects an appropriate scoring matrix and tinkers
    a structure containing the monomer order of the scoring matrix
    and the scoring matrix */

score_m select_scm(score_m **matx, sequ *os, iseq *is, u_sc smat)
{
  static char *mat[] = {"gon", "pam", "blo"};
  aligm a;
  float d, pam;
  score_m scores;
  
  /* align calculates an optimal alignment (using function trace_back
     for backtracking). align returns a structure of the typ aligm which
     contains an optimal alignment encoded in ints (see trace_back) and
     the score of this alignment */
  
  a=align_2(os,is,smat);

  

  /* the function pairwise_identity calculates the observed distance
     and the pam distance of the two sequences. its returns a structure
     of typ dist, which contains the observed and the pam distance */  


  d=observed_identity(a);

  /* a is not used anymore:free(a.a) allocated in trace_back ! */
  free(a.a);

  /* calculates the pam distance from the observed distance */ 
  pam=pam_distance(d);
  
  
  /* blosum series: 30,50,62,80 */
  if(strncasecmp(MAT_SER,mat[2],3) == 0)
    {
      double id;

      id=100*d; /* observed identity * 100 */ 


      if(id <= 40.0)
	{
	  scores=*matx[0];
	}

      if(40.0 < id && id <= 56.0)
	{
	  scores=*matx[1];
	}
      
      if(56.0 < id && id <= 70.0)
	{
	  scores=*matx[2];
	}
      
      if(id > 70.0)
	{
	  scores=*matx[3];
	}
    }

  /* pam series: 20,60,120,350 */
  else if(strncasecmp(MAT_SER,mat[1],3) == 0)
    {
      /*  printf("pam distance:%f\n",pam); */
      
      if(pam <= 40.0)
	scores=*matx[0];

      if(40.0 < pam && pam <= 90.0)
	scores=*matx[1];

      if(90.0 < pam && pam <= 235.0)
	scores=*matx[2];

      if(pam > 235.0)
	{
	  scores=*matx[3];
	  /*  printf("score matrix:%s\n",matx[3].name); */
	}
    }

  /* gonnet series: 40,80,120,160,250,300,350 */
  else if(strncasecmp(MAT_SER,mat[0],3) == 0)
    {
      /*  printf("pam distance:%f\n",pam); */

      if(pam <= 60.0)
	{
	  scores=*matx[0];
	  /*  printf("score matrix:%s\n",matx[0].name); */
	}
      
      if(60.0 < pam && pam <= 100.0)
	{
	  scores=*matx[1];
	  /*  printf("score matrix:%s\n",matx[1].name); */
	}
      
      if(100.0 < pam && pam <= 140.0)
	{
	  scores=*matx[2];
	  /*  printf("score matrix:%s\n",matx[2].name); */
	}
      
       if(140.0 < pam && pam <= 205.0)
	 {
	   scores=*matx[3];
	   /*  printf("score matrix:%s\n",matx[3].name); */
	 }
       if(205.0 < pam && pam <= 275.0)
	 {
	   scores=*matx[4];
	   /*  printf("score matrix:%s\n",matx[4].name); */
	 }
       
       if(275.0 < pam && pam <= 325.0)
	 {
	   scores=*matx[5];
	   /*  printf("score matrix:%s\n",matx[5].name); */
	 }

       if(pam > 325.0)
	 {
	   scores=*matx[6];
	   /*  printf("score matrix:%s\n",matx[6].name); */
	 }
    }
  else
  /*  if((strncasecmp(MAT_SER,mat[0],3)!=0) */
/*       &&(strncasecmp(MAT_SER,mat[1],3)!=0) */
/*       &&(strncasecmp(MAT_SER,mat[2],3)!=0)) */
    {
      fprintf(stderr,"select_score_mat:couldn't select a scoring matrix");
      usage(0);
      exit(3);
    }
  

  return(scores);

}

/*---------------------------------------------------------------------------*/
/* this function decides wether the input polymer is DNA or protein
   check if DNA or Protein
   The decision is based on counting all A,C,G,T,U 
   If >= 85% of all characters of the longer sequence are as above => DNA
   the idea for this function is taken from clustalw*/

int check_polymer(char *seq, int l)
{
  int i,j,c,residues;
  double ratio;
  char *dna_codes="ACGTU";
  int ldc;

  ldc=strlen(dna_codes);
  
  
  residues=0;
  
  for(i=0; i<l; i++)
    {
      c=seq[i];
      for(j=0; j<ldc; j++)
	{
	  if(dna_codes[j]==seq[i])
	    ++residues;
	}
    }
  
  ratio=(double)residues/(double)l;
  
  if(ratio >= 0.85)
    i=1;
  else
    i=0;
  
  return(i);
}
/*---------------------------------------------------------------------------*/
/* this function finds the longer one of the two input sequences, writes it
   in the first member of the sequence array (the shorter one is written in
   the secound member of sequence array) and converts both sequences to
   uppercase letters*/
sequ *sequences(sequ *seq_array)
{
  int i,lx,ly;
  int lnx,lny;
  char c;
  sequ *a, *b;
  sequ *s;


  a = (sequ*) calloc(1,sizeof(sequ));
  b = (sequ*) calloc(1,sizeof(sequ));

  s = (sequ*) calloc(2,sizeof(sequ));

  lnx=strlen(seq_array[0].name);
  lx=strlen(seq_array[0].seq);

  lny=strlen(seq_array[1].name);
  ly=strlen(seq_array[1].seq);

  if(lx >= ly)
    {
      a[0].name=(char*) calloc((lnx+1),sizeof(char));
      a[0].seq=(char*) calloc((lx+1),sizeof(char));
      b[0].name=(char*) calloc((lny+1),sizeof(char));
      b[0].seq=(char*) calloc((ly+1),sizeof(char));
      
      strcpy(a[0].name,seq_array[0].name);
      strcpy(a[0].seq,seq_array[0].seq);
      strcpy(b[0].name,seq_array[1].name);
      strcpy(b[0].seq,seq_array[1].seq);
    }
  else
    {
      a[0].name=(char*) calloc((lny+1),sizeof(char));
      a[0].seq=(char*) calloc((ly+1),sizeof(char));
      b[0].name=(char*) calloc((lnx+1),sizeof(char));
      b[0].seq=(char*) calloc((lx+1),sizeof(char));
      
      strcpy(a[0].name,seq_array[1].name);
      strcpy(a[0].seq,seq_array[1].seq);
      strcpy(b[0].name,seq_array[0].name);
      strcpy(b[0].seq,seq_array[0].seq);
      
    }
  lnx=0;
  lny=0;
  lx=0;
  ly=0;

  lnx=strlen(a[0].name);
  lx=strlen(a[0].seq);

  lny=strlen(b[0].name);
  ly=strlen(b[0].seq);
  
  s[0].name=(char*) calloc((lnx+1),sizeof(char));
  s[0].seq=(char*) calloc((lx+1),sizeof(char));
  s[1].name=(char*) calloc((lny+1),sizeof(char));
  s[1].seq=(char*) calloc((ly+1),sizeof(char));

  strcpy(s[0].name,a[0].name);
  strcpy(s[0].seq,a[0].seq);
  strcpy(s[1].name,b[0].name);
  strcpy(s[1].seq,b[0].seq);


  for(i=0; i<lx; i++)
    {
      c=toupper(s[0].seq[i]);
      s[0].seq[i]=c;
    }
  for(i=0; i<ly; i++)
    {
      c=toupper(s[1].seq[i]);
      s[1].seq[i]=c;
    }
  /*
    printf("sequences:\n%s\n%s\n\n%s\n%s\n\n",s[0].name,s[0].seq,s[1].name,s[1].seq);     */ 

  free(a[0].name);
  free(a[0].seq);
  free(a);
  

  free(b[0].name);
  free(b[0].seq);
  free(b);
  
  return(s);
}

/*---------------------------------------------------------------------------*/

/* calculate the score for an optimal alignment using the gotoh algorithm */
aligm align_2(sequ *os, iseq *is, u_sc smat)
{
  int i,j;
  int l0;
  int l1;
  int l;
  int *s0, *s1;
  matrix m;
  aligm a;

  
  s0=is[0].s;
  s1=is[1].s;
  
  l=strlen(smat.monomers);
  l0=strlen(os[0].seq);
  l1=strlen(os[1].seq);

  m.M=(real**) calloc((l1+1),sizeof(real*));
  for(i=0; i<l1+1; i++)
    {
     m.M[i]=(real*) calloc((l0+1),sizeof(real));
    }
  
  m.E=(real**) calloc((l1+1),sizeof(real*));
  for(i=0; i<l1+1; i++)
    {
      m.E[i]=(real*) calloc((l0+1),sizeof(real));
    }
  
  m.F=(real**) calloc((l1+1),sizeof(real*));
  for(i=0; i<l1+1; i++)
    {
       m.F[i]=(real*) calloc((l0+1),sizeof(real));
    }

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
                          k = lenght of the gap;
                          M(0,j>0) = M(i>0,0) = 4*open + 10*(k-1)*ext;(*)
			  (*) the values are selected in order to avoid the
			      use of these entries (imagine a kind of "empty"
			      alignment (e.g. seq1 aligned only to gap
			      characters), which you simply don't want to
			      consider (you are looking for the best
			      alignment in which i and j are paired!)

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
	{
	  m.E[i][j]=m.E[i][j-1]+smat.endgaps; /* value for endgaps */ 
	}
      /* initialise F(i....l1+1,0) */
      j=0;
      for(i=1; i< l1+1; i++)
	{
	  m.F[i][j]=m.F[i-1][j]+smat.endgaps; /* value for endgaps */ 
	}
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
	  /* M */
	  m.M[i][j]=m.M[i-1][j-1]+smat.mat[s1[i-1]][s0[j-1]];
	  /* M or MM */
	  e=m.E[i-1][j-1]+smat.mat[s1[i-1]][s0[j-1]]; /* M or MM */
	  f=m.F[i-1][j-1]+smat.mat[s1[i-1]][s0[j-1]]; /* M or MM */

	  max=MAX(m.M[i][j],e);
	  m.M[i][j]=MAX(max,f);

	  max=0;
	  e=0;
	  f=0;

	}      
    }

  /*  printf("align_2: score = %.4f\n",m.M[l1][l0]); */
  
  /*printf("M\n");
  print_M(os, m.M);
  printf("E\n");
  print_M(os, m.E);
  printf("F\n");
  print_M(os, m.F);*/
  
   
  a=trace_back(m, os, is, smat); 
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
/*---------------------------------------------------------------------------*/
void print_M(sequ *os, real **matr)
{
  int i,j;
  int l0;
  int l1;


  l0 = strlen(os[0].seq);
  l1 = strlen(os[1].seq);

  printf("      %c",GAP);
  for(j=0; j<l0; j++){
    printf("%8c",os[0].seq[j]);
  }

   printf("\n\n");
      printf("%c",GAP);

   for(j=0; j<=l0; j++){
       printf("  %6.1f",matr[0][j]);
   }
  
   printf("\n\n");
   for(i=1;i<=l1;i++){
   
     printf("%c",os[1].seq[i-1]);
    
     for(j=0;j<=l0;j++){
       printf("  %6.1f",matr[i][j]);
       
       if(j==l0){
         printf("\n\n");
       }
     }
   }   
      
   printf("\n\n");
}
/*---------------------------------------------------------------------------*/
/* calculate the pairwise identity and the pam distance of the two sequences */
float observed_identity(aligm a)
{
  int i;
  int match;
  int mismatch;
  int indel;
  float d;

  match=0;
  mismatch=0;
  indel=0;

  i=0;
  while(a.a[i])
    {
      if(a.a[i] == '|')
	{
	  ++match;
	}
      if(a.a[i] == ':')
	{
	  ++mismatch;
	}
      if((a.a[i] == '.') || (a.a[i] == '^'))
	{
	  ++indel;
	}
      ++i;
    }
  
  /*
    In Clustal V we used a simple formula to convert an observed distance
    to one that is corrected for multiple hits.  The observed distance is
    the mean number of differences per site in an alignment (ignoring sites
    with a gap) and is therefore always between 0.0 (for identical sequences)
    and 1.0 (no residues the same at any site).  These distances can be multi-
    plied by 100 to give percent difference values.  100 minus percent dif-
    ference gives percent identity.
    The formula we use to correct for multiple hits is from Motoo Kimura
    (Kimura, M. The neutral Theory of Molecular Evolution, Camb.Univ.Press,
    1983, page 75) and is:

    K = -Ln(1 - D - (D.D)/5)  where D is the observed distance and K is       
                              corrected distance.
			      
    This formula gives mean number of estimated substitutions per site and,
    in contrast to D (the observed number), can be greater than 1 i.e. more
    than one substitution per site, on average.
    This can also be expressed in PAM units by multiplying by 100 (mean
    number of substitutions per 100 residues).
    Dayhoff et al constructed an elaborate model of protein evolution based
    on observed frequencies of substitution between very closely related
    proteins.  Using this model, they derived a table relating observed
    distances to predicted PAM distances.
    Kimura's formula, above, is just a "curve fitting" approximation to this
    table. It is very accurate in the range 0.75 > D > 0.0 but becomes
    increasingly unaccurate at high D (>0.75) and fails completely at around
    D = 0.85.
*/


  /* observed identity => blosum:
     (a matrix derived from a data base of blocks in which sequence
     segments that are identical at >= 80% of aligned residues are
     clusted is referred to as BLOSUM 80, .. */
  if((match+mismatch) == 0)
    {
      d = 0.15; /* no matches or mismatches => high pam distance */ 
    }
  else
    {
      d = (double) match/(match + mismatch);
    }
  return(d);
}
/*---------------------------------------------------------------------------*/
float pam_distance(float d)
{
  float p_i,p_d;
  float k,dis;
  float pam;


  /* observed distance  => Clustal V*/ 
  dis =1 - d;

  /* percent difference */ 
  p_i=dis*100;

  /* pairwise identity */
  p_d=100 - (p_i);

  /* corrected distance */ 
  k = -log(1 - dis - (dis*dis)/5);

  /* pam distance */ 
  pam=k*100;

  if(fabs(k - 0.0 ) < 100*FLT_EPSILON)
    {
      k=0.0;
    }

  /* pam distance */ 
  pam=k*100;
  

  /* for observed distance > 0.751: a rough approximation of the pam
     distance is done */
  
  if(dis >= 0.751 && dis < 0.762 )
    pam=196;
  
  if(dis >= 0.762 && dis < 0.791 )
    pam=206;
   
  if(dis >= 0.791 && dis < 0.821)
    pam=236;

  if(dis >= 0.821 && dis < 0.848)
    pam=276;

  if(dis >= 0.848)
    pam=326;

  return(pam);

}
/*---------------------------------------------------------------------------*/
/* adds the sequences in the order used in the programm to the structure
   containing the alignment encoded as a string of digits -> This is
   important to decode the alignment */
aligm add_sequences(sequ *os, aligm a)
{
  int l0, l1, l, la;

  
  l0=strlen(os[0].seq);
  l1= strlen(os[1].seq);
  l=strlen(os[0].name);
  la=strlen(os[1].name);

  
  a.s0.name = (char*) calloc((l+1),sizeof(char));
  a.s1.name = (char*) calloc((la+1),sizeof(char));
  a.s0.seq = (char*) calloc((l0+1),sizeof(char));
  a.s1.seq = (char*) calloc((l1+1),sizeof(char));

  strcpy(a.s0.name, os[0].name);
  strcpy(a.s0.seq, os[0].seq);
  strcpy(a.s1.name, os[1].name);
  strcpy(a.s1.seq, os[1].seq);

  return(a);

}

/*---------------------------------------------------------------------------*/
/* the function trace_back finds an optimal alignment between the two
   sequences, the alignment is represented as a string of ints:
   alignment: encoded in ints
                          macro symbol meaning
	                1 M1    '|'    match
			2 M2    ':'    mismatch
	                3 G3    '.'    gaps in the shorter seq(s1,E)
	                4 G4    '^'    gaps in the longer seq(s0,F)
   see begin of this file - may be the symbols are not up to date!!! */
#define EQUAL(A,B) (fabs((A)-(B)) < 100*FLT_EPSILON)
aligm trace_back(matrix m, sequ *os, iseq *is, u_sc smat)
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

  /* determine the score of an optimal aligment */
  a.score=max;
  
  while((i>0 && j>=0) || (i>=0 && j>0) ) {

    switch (state)
      {
      case Mm:
	/*  printf("i %d  j %d\n",i,j); */ 
	a.a[la] = (s0[j-1] == s1[i-1]) ? M1 : M2;

	if( EQUAL(M[i-1][j-1]+smat.mat[s1[i-1]][s0[j-1]], M[i][j]) )
	  {/* nothing to do */}
	else if( EQUAL(E[i-1][j-1]+smat.mat[s1[i-1]][s0[j-1]],M[i][j]) )
	  { state = Ee; }
	else if (EQUAL(F[i-1][j-1]+smat.mat[s1[i-1]][s0[j-1]],M[i][j])) 
	  { state = Ff; }
	else {
	  fprintf(stderr, "Error in backtracking case M\n"); exit(27);
	}
	i--; j--;
	break;

      case Ee:
	/*  printf("i %d  j %d\n",i,j); */ 
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
	  /* smat.mat[l*][Hier richtig!] im faltungsteil noch aendern! */
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
	/*  printf("i %d  j %d\n",i,j); */ 
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
	  /* Hier richtig! im faltungsteil noch aendern! */
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
/*---------------------------------------------------------------------------*/

/* calculate the propability of each basepair i,j */
real **prop_alig(sequ *os, iseq *is, u_sc m,real **mat,real **revm,aligm a)
{
  int i,j;
  int k,l;
  int *s0, *s1;
  int l1;
  int l0;
  double opt;
  double s;
  real **P;
  real **X;

  s0=is[0].s;
  s1=is[1].s;
  
  l0=strlen(os[0].seq);
  l1=strlen(os[1].seq);

  /* l1 rows and l0 colums */ 
  P=(real**) calloc((l1+1),sizeof(real*));
  for(i=0; i<l1+1; i++)
    {
      P[i]=(real*) calloc((l0+1),sizeof(real));
    }
  /* dot plot l0 rows and l1 colums */ 
  X=(real**) calloc((l0+1),sizeof(real*));
  for(i=0;i<l0+1;i++){
  X[i] = (real*) calloc((l1+1),sizeof(real));  
  }

  /*  printf("prop_alig: %s\n",m.name); */

  P[0][0]=(mat[0][0]/revm[l1][l0])-1;

  P[0][0]=0;
 
  opt=a.score;

  /*  printf("opt:%d\n",opt);*/
  
  s=2.0*opt/(l0+l1); 
  s=exp(BETA*s);   /* Fehler beim Skalieren? */
  

  for(i=1; i<=l1; i++)
    {
      for(j=1; j<=l0; j++)
	{
	  double z;

	  z=(m.mat[s0[j-1]][s1[i-1]]);
	   /* Fehler beim Skalieren? */
	  P[i][j]=mat[i][j]*revm[i-1][j-1]*s;/* *s */
	  P[i][j]/=(exp(BETA*z)*mat[0][0]);
	  /*  P[i][j]*=exp(BETA*s*(i+j)/2);  */
	}
    }

  
  
  /*  printf("prop_alig: probs from partion function\n");     */
/*    print_P(os,P); */

  for(k=i=1; i&&k <=l1; i++,k++)
    {
      for(l=j=1; j&&l <=l0; j++,l++)
	{
	  X[l][k]=P[i][j];
	}
    }
  
  for(i=0; i<l1+1; i++)
    {
      free(P[i]);
    }
  free(P);
  
  return(X);
  
}

/*---------------------------------------------------------------------------*/

matrix revers_partf(sequ *os, iseq *is, u_sc m, aligm a)
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
 
  l0=strlen(os[0].seq);
  l1=strlen(os[1].seq);
  l=strlen(m.monomers);
  s0=is[0].s;
  s1=is[1].s;
  
  /*i runs along the rows, j runs along the columns */

  /* Algorithmus for the revers partion function */
  /* rMi,j = exp(BETA * Match) * (rM i-1,j-1 + rE i-1,j-1 + rF i-1,j-1 );

     rEi,j = exp(BETA * Open) * (rMi,j-1 + rFi,j-1) +
             exp(BETA * Extend) * (rEi,j-1);

     rFi,j = exp(BETA * Open) * (rMi-1,j) +
             exp(BETA * Extend) * (rFi-1,j); */
  
  rM=(real**) calloc((l1+1),sizeof(real*));
  for(i=0; i<l1+1; i++)
    {
      rM[i]=(real*) calloc((l0+1),sizeof(real));
    }

  rF=(real**) calloc((2),sizeof(real*));
  for(i=0; i<2; i++)
    {
      rF[i]=(real*) calloc((l0+1),sizeof(real));
    }
 
  rE=(real**) calloc((2),sizeof(real*));
  for(i=0; i<2; i++)
    {
      rE[i]=(real*) calloc((l0+1),sizeof(real));
    }

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
	  rM[i][j]=exp(BETA*m.mat[s0[j]][s1[i]]);
	  rM[i][j]*=(rM[i+1][j+1]+rF[0][j+1]+rE[0][j+1]);
	  rM[i][j]/=sk*sk;
	  

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
/*--------------------------------------------------------------------------*/
void print_Zrevers(sequ *s, real **matr)
{
  int i,j;
  int l0;
  int l1;
  
  l0=strlen(s[0].seq);
  l1=strlen(s[1].seq);
  
  for(i=0; i< l0; i++){
    printf("         %c",s[0].seq[i]);
  }

   printf("\n\n");
   for(i=0;i<=l1;i++){
     if(i<l1)
       printf("%c",s[1].seq[i]);
     if(i==l1)
        printf(" ");
     for(j=0;j<=l0;j++){
       printf("  %.2e",matr[i][j]);
       
       if(j==l0){
         printf("\n\n");
       }
     }
   }   
      
   printf("\n\n");
}

/*---------------------------------------------------------------------------*/
  /* calculate the partition function using a modification of the gotoh
     algorithm for sequence alignment */
matrix partf(aligm a)
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
	  printf(" Parameters:\n");
	  printf(" T = %.2f\n",T);
	  printf(" scoring matrix: %s;\n",m.name);
	  printf(" gap penalties open:%.3f, extend:%.3f;\n",m.mat[l-2][s0[0]],m.mat[l-1][s0[0]]);
	  if(Egap_flag)
	  {
		  printf(" endgap penalty: %.3f\n",m.endgaps);
	  }
  }

  /* [i][j] => i runs along the rows, j runs along the columns !!!!*/
 
  /* zM: calculate only those alignements in which i and j are matched
     zM(i,j)=exp(BETA*(M or MM))*(zM(i-1,j-1)+zE(i-1,j-1)+zF(i-1,j-1)); */
  zM=(real**) calloc((l1+1), sizeof(real*));
  for(i=0; i<l1+1; i++)
    {
      zM[i]=(real*) calloc((l0+1), sizeof(real));
    }
  
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
    {
      zF[i]=(real*) calloc((l0+1), sizeof(real));
    }
  
  /* the recurrence for zE is zE(i,j)= zM(i,j-1)*exp(BETA*(-OPEN) +
                                     + zE(i,j-1)*exp(BETA*(-EXT);
     one therefore never leaves the current row but one needs the values
     of the previous column to calculate zF(i,j)  */

  zE=(real**) calloc((l1+1), sizeof(real*));
  for(i=0; i<l1+1; i++)
    {
      zE[i]=(real*) calloc((l0+1), sizeof(real));
    }
  
  /* opt: optimale score of the alignement of the two sequences */
  opt=a.score;
  s=2.0*opt/(l0+l1);
  sk=exp(BETA*s/2);
  
  
  /* zM[0][0]=1 */
  zM[0][0]=1;
 
  /* zM[0][1...l0], zM[1...l0][0], zF[0][1...l0] and zE[1...l0][0]
     are all 0, calloc initializes the allocated memory with 0 */ 

  /* zE[0][0]=0;
     zE[0][j]=zM[0][j-1]*exp(BETA*(-OPEN)+zE[0][j-1]*exp(BETA*(-EXT) */

  i=0;
  j=1;
  /* zE[0][0]=0 */
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
	{
	  zE[0][j]=zE[0][j-1]*exp(BETA*m.endgaps)/sk;
	}
      i=1;
      j=0;
      zF[1][0]=zM[0][0]*exp(BETA*m.endgaps)/sk;
      for(i=2; i<=l1; i++)
	{
	  zF[i][0]=zF[i-1][0]*exp(BETA*m.endgaps)/sk; 
	}
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

	  /* HERE ADD THE DOTALIGNER SINGLE STRANDED SCORING SCHEME */
	  z=0;
	  z=(m.mat[s0[j-1]][s1[i-1]]);  
	  zM[i][j]=exp(BETA*z)*(zM[i-1][j-1]+zE[i-1][j-1]+zF[i-1][j-1]);
	  zM[i][j]/=sk*sk;
	}
    }

  Z=zM[i-1][j-1]+zF[i-1][j-1]+zE[i-1][j-1]; 
  zM[0][0]=Z;
  
 /*  printf("Zz=%f\n",Z); */
  
/*   printf("zM\n"); */
/*   print_Z(so,zM); */
 
/*   printf("zE\n"); */
/*   print_Z(so,zE); */

/*   printf("zF\n"); */
/*   print_Z(so,zF); */
  
  pf.M=zM;
  pf.E=zE;
  pf.F=zF;
  
  return(pf);
}

/*--------------------------------------------------------------------------*/
void print_Z(sequ *s, real **matr)
{
  int i,j;
  int la = strlen(s[0].seq);
  int lb = strlen(s[1].seq);
  
   printf("         %c",GAP);
  for(j=0; j<la; j++){
    printf("         %c",s[0].seq[j]);
  }

   printf("\n\n");
      printf("%c",GAP);

   for(j=0; j<=la; j++){
       printf("  %.2e",matr[0][j]);
   }
  
   printf("\n\n");
   for(i=1;i<=lb;i++){
   
       printf("%c",s[1].seq[i-1]);
    
     for(j=0;j<=la;j++){
       printf("  %.2e",matr[i][j]);
       
       if(j==la){
         printf("\n\n");
       }
     }
   }   
      
   printf("\n\n");
}
/*---------------------------------------------------------------------------*/
void print_P(sequ *s, real **matr)
{
  int i,j;
  int la = strlen(s[0].seq);
  int lb = strlen(s[1].seq);
  
   printf("\n");
   for(i=1;i<=la;i++)
     {
       printf("  %8c",s[0].seq[i-1]);
     }
   printf("\n");
   for(i=1;i<=lb;i++)
     {
       printf("%c",s[1].seq[i-1]);
       for(j=1;j<=la;j++)
	 {
	   printf("  %.2e",matr[i][j]);
       
	   if(j==la){
	     printf("\n\n");
	   }
	 }
     }   
      
   printf("\n\n");
}
/*-----------------------------------------------------------------------*/

aligm stoch_btr(matrix pz, aligm a)
{
  int i,j,t;
  int l0,l1,l,Eg;
  int *s0, *s1;
  double s,opt,sk;
  double rand, score;
  double x,y,z,Zx,bla;
  real pZ;
  real **M;
  real **E;
  real **F;
  /*real **Z; just for testing
    long **Test;*/
  aligm sta;
  char *tr;
  u_sc m;
  int opt_l;
  double s_p;
  /*  two_f P_P; */ 

  
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

  pZ=pz.M[0][0];/* partition function over all alignments */

  
  if(stochswi == 0)
    {
      double tmp;
      double prob_o; /*  probability of the (or better an) optimal alignment */

      /*  prob_o=exp(BETA*a.score)/(pZ*pow(sk,(l0+l1))); */
      /*tmp = (BETA*a.score) - (l0+l1)*log(sk);
      prob_o=exp(tmp) /pZ;

	  printf("pZ=%.3f, score=%.3f, trace=%s\n",pZ,a.score,a.a);
      
      printf("#number of stochastic alignments = %ld\n",Nr);
      printf("#score of the optimal alignment =  %.2f \n",a.score);
      printf("#probability of the optimal alignment =  %.3e \n",prob_o );
      
      decode_alig(a);
      
      printf("#\n");
      printf("#upper sequence\n");
      printf("# > %s\n# %s\n",so[0].name,so[0].seq);
      printf("#lower sequence\n");
      printf("# > %s\n# %s\n",so[1].name,so[1].seq);*/

      /*  if(calcS_flag)
	  {
	  printf("#\n# input alignment from file %s\n",TRACK);
	  s_p = tmp = prob_o = 0.0;
	  s_p = calc_score(track);
	  printf("# %s\t%.4f\t",track, s_p);
	  prob_o=exp(BETA*a.score)/(pZ*pow(sk,(l0+l1)));
	  tmp = (BETA*s_p) - (l0+l1)*log(sk);
	  prob_o=exp(tmp) /pZ;
	  printf("%.3e",prob_o);
	  P_P = comp_aln(track,track);
	  P_P.o %id with ref, P_P.e mean matchprob
	  printf("\t%.2f\t%.2e\n",P_P.o,P_P.e);
	  }  */
      //printf("#optimal alignment  \n");
      /* ask if the score of the optimal alignment calculated by the
	 Needleman-Wunsch algorithm is the same as the one calculated
	 by function calc_score(). */ 

      s_p=calc_score(a.a);
      if(!(fabs(a.score - s_p) < 100*FLT_EPSILON) )
	{
	  printf("Error in backtracking\n");
	  exit(44);
	}
      printf("%s\t",a.a);
     /*   printf("%.4f\t",calc_score(a.a)); */ 
      printf("%.4f\t",a.score);
      tmp = (BETA*a.score) - (l0+l1)*log(sk);
      prob_o=exp(tmp) /pZ;
      printf("%.2e",prob_o);
      //printf("\n");
  
      stochswi=1;

      return sta;
    }

  M=pz.M;
  E=pz.E;
  F=pz.F;
  
  if(!sran)
    {
      sgenrand((unsigned long) time(NULL));
      sran=1;
    }
  
  /*Test=(long* *) myspace((l1+1)*sizeof(long*));
    for(i=0; i<l1+1; i++)
    {
    Test[i]=(long*) myspace((l0+1)*sizeof(long));
    }
  
    Z=(real**) myspace((l1+1)*sizeof(real*));
    for(i=0; i<l1+1; i++)
    {
    Z[i]=(real*) myspace((l0+1)*sizeof(real));
    }*/
  
  tr=(char*) calloc((l0+l1+1),sizeof(char));
  
 
  /* Backtracking */
  rand=genrand();

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
      
      if(rand < x  )
	{
	  
	  /* M => match = 1 */
	  if(s0[j-1] == s1[i-1])
	    {
	      tr[t--]=M1;/* 1 match */
	    }
	  else
	    {
	      tr[t--]=M2;/* 2 mismatch */
	    }
	  score*=x;
	  
	  /*  ++Test[i][j]; */ 
	  
	  --i;
	  --j;
	  
	  Zx=M[i][j]+E[i][j]+F[i][j];
	      
	  x=M[i][j]/Zx;
	  y=E[i][j]/Zx;
	  z=F[i][j]/Zx;

	}
	  
      else if(rand <  x+y )
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
      rand=genrand();
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
      
      /*for(i=0; i<=l0+l1; i++)
	{
	  tr[i]=0;
	  }*/

      
      /* sta.al.s0=(char*) calloc((t+1),sizeof(char));
      sta.al.s0[t]='\0';
      sta.al.s1=(char*) calloc((t+1),sizeof(char));
      sta.al.s1[t]='\0';*/

      sta.score=0.0;
      sta.prob=0.0;
   

      /* fill in structure aligm */

      s_p=calc_score(sta.a);
      sta.score=s_p;
      
      
      /*  sta.prob=exp(BETA*sta.score)/(pZ*pow(sk,(l0+l1))); */

      bla = (BETA*sta.score) - (l0+l1)*log(sk);
      sta.prob=exp(bla) /pZ;
      
      sta.s0.name=strdup(so[0].name);
      sta.s0.seq=strdup(so[0].seq);
      sta.s1.name=strdup(so[1].name);
      sta.s1.seq=strdup(so[1].seq);

      /* output: aligment(string), Score(A), Prob(A) */

      /*  printf("%-*s\t",opt_l,sta.a); */
      printf("%s\t",sta.a);
      printf("%.4f\t",sta.score);
      printf("%.2e",sta.prob);
      //printf("\n");
      
      
      free(tr); 

    


      /*  test - test - test 
	  for(i=0; i<=l1; i++)   
	  {
	  for(j=0; j<=l0; j++)
	  {
	  Z[i][j]=0.0;
	  }
	  }
	  printf("n=%ld ",Nr);
	  for(i=1; i<l1+1; i++)
	  {
	  for(j=1; j<l0+1; j++)
	  {
	  Z[i][j]=(real) Test[i][j]/Nr;
	  }
	  }
	  printf("\n");
	  printf("\nTest\n");
	  for(i=1; i<l1+1; i++)
	  {
	  for(j=1; j<l0+1; j++)
	  {
	  printf("%12ld",Test[i][j]);
	  }
	  printf("\n");
	  }
	  printf("\n\n");
	  print_P(so,Z);*/

      
	
      return(sta);

}

/*-----------------------------------------------------------------------*/
void nrerror(char *message)
{
  fprintf(stderr, "\n%s\n", message);
  exit(-1);
}
/*-----------------------------------------------------------------------*/
/* calculate the score of an alignment given the alignment encoded as one
   string of symbols */ 
double calc_score(char *tr)
{
  /* matrix Z contains the match probabilities */ 
  int i,j,t,k,g3,g4;
  int *s0, *s1;
  int l0,l1,l,Eg;
  u_sc m;
  double s; /* score of the stoch_aln */
  int bla;
  double result; 

  bla=0; /* set bla to 1 to print the stuff in the while condition */
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
	  g3=0;
	  g4=0;
	  /* this function was originally part of fkt stoch_btr !!*/
	  /* sta.al.s0[i]=so[0].seq[k]; sta.al.s1[i]=so[1].seq[j]; */
	  s+=m.mat[s0[k]][s1[j]];
	  if(bla)
		{
		  printf("%.3f  s=%.3f\n",m.mat[s0[k]][s1[j]], s);
		}
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
	      if(bla)
		{
		  printf("%.3f  s=%.3f\n",m.endgaps, s);
		}
	    }
	  else if(g3 == 0)
	    {
	      s+=m.mat[l-2][s0[k]];/* Open in s0 (3) */
	      if(bla)
		{
		  printf("%.3f  s=%.3f\n",m.mat[l-2][s0[k]], s);
		}
	    }
	  else
	    {
	      s+=m.mat[l-1][s0[k]];/* Extend in s0 (3) */
	      if(bla)
		{
		  printf("%.3f  s=%.3f\n",m.mat[l-1][s0[k]], s);
		}
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
	      if(bla)
		{
		  printf("%.3f  s=%.3f\n",m.endgaps, s);
		}
	    }
	  else if(g4 == 0)
	    {
	      s+=m.mat[l-2][s1[j]];/* Open in s1 (4) */
	      if(bla)
		{
		  printf("%.3f  s=%.3f\n",m.mat[l-2][s1[j]], s);
		}
	    }
	  else
	    {
	      s+=m.mat[l-1][s1[j]];/* Extend in s1 (4) */
	      if(bla)
		{
		  printf("%.3f  s=%.3f\n",m.mat[l-1][s1[j]], s);
		}
	    }
	  ++g4;
	  ++j;
	}  
      ++i;
      
    }

  result = s;
  
  return (result);
  
}
/*-----------------------------------------------------------------------*/






























