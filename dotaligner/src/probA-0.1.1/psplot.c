/* Last changed Time-stamp: <2001-07-03 10:32:43 ulim> */
/* psplot.c */
#include "pfgoto.h"

#define PMIN 0.00001
#define  MAX(A,B)    (A)>(B)?(A):(B)
static char *time_stamp (void);
/*  int PS_plot(sequ so, real **P, char *wastlfile); */

/* produce PostScript dot plot from probabilities in pr[] array */
int PS_plot(sequ *so, real **P, char *wastlfile) {
  
  FILE *wastl;
  char name[31], *c;
  int i, j, la, lb;
  double tmp;
  char *A, *B;

  A=strdup(so[0].seq);
  B=strdup(so[1].seq);
  
  la = strlen(A);
  lb = strlen(B);
  if (wastlfile == NULL) wastl = stdout;
  if ( (wastl = fopen(wastlfile,"w")) == 0) {
    fprintf(stderr, "can't open %s for ps-plot\n", wastlfile);
    return 1;
  }
  
  strncpy(name, wastlfile, 30);
  if ((c=strrchr(name, '_'))!=0) *c='\0';
  fprintf(wastl,"%%!PS-Adobe-2.0 EPSF-1.2\n");
  fprintf(wastl,"%%%%Title: RNA DotPlot\n");
  fprintf(wastl,"%%%%Creator: PS_dot.c, ViennaRNA package\n");
  fprintf(wastl,"%%%%CreationDate: %s", time_stamp());
  fprintf(wastl,"%%%%BoundingBox: 20 2000 540 680\n");/* 255 - 2000 */
  fprintf(wastl,"%%%%DocumentFonts: Times-Roman\n");
  fprintf(wastl,"%%%%Pages: 1\n");
  fprintf(wastl,"%%%%EndComments\n\n");
  
  fprintf(wastl,"%%This file contains the square roots "
	  "of the base pair probabilities in the form\n");
  fprintf(wastl,"%% i  j  sqrt(p(i,j)) ubox\n");
  
  fprintf(wastl,"100 dict begin\n");  /* DSC says EPS should create a dict */
  fprintf(wastl,"\n/logscale false def\n\n");/*switch log (true) or not */
  fprintf(wastl,"%%delete next line to get rid of title\n"
	  "270 665 moveto /Times-Roman findfont 14 scalefont setfont "
	  "(%s) show\n\n", name);
  fprintf(wastl,"/lpmin {\n"
	  "   %g log  %% log(pmin) only probs>pmin will be shown\n"
	  "} bind def\n\n",PMIN);
  
  fprintf(wastl,"/ubox {\n"     /* upper triangle matrix */
	  "   logscale {\n"
	  "      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if\n"
	  "   } if\n"
	  "   3 1 roll\n"
	  "   exch lenA exch sub 1 add box\n"
	  "} bind def\n\n");

  fprintf(wastl,"/box { %%size x y box - draws box centered on x,y\n"
	  "   2 index 0.5 mul add            %% x += 0.5\n"
	  "   exch 2 index 0.5 mul add exch  %% x += 0.5\n"
	  "   newpath\n"
	  "   moveto\n"
	  "   dup neg   0 rlineto\n"
	  "   dup neg   0 exch rlineto\n"
	  "             0 rlineto\n"
	  "   closepath\n"
	  "   fill\n"
	  "} def\n\n");

  /* EPS should not contain lines >255 ( 255 - 2000 ) characters */

  fprintf(wastl,"/seqA { (\\\n");
  i=0;
  while (i<la) {
    fprintf(wastl, "%.2000s\\\n", A+i); /* 255 - 2000 */
    i+=2000;/* 255 - 2000 */
  }
  fprintf(wastl,") } def\n");

  fprintf(wastl,"/seqB { (\\\n");  
  i=0;
  while (i<lb) {
    fprintf(wastl, "%.2000s\\\n", B+i); /* 255 - 2000 */
    i+=2000;/* 255 - 2000 */
  }
  fprintf(wastl,") } def\n");
  
  fprintf(wastl,"/lenA { seqA length } def\n\n");
  fprintf(wastl,"/lenB { seqB length } def\n\n");
  
  fprintf(wastl,"72 216 translate\n");
  fprintf(wastl,"72 6 mul lenA 1 add div dup scale\n");
  fprintf(wastl,"/Times-Roman findfont 0.95 scalefont setfont\n\n");

  /* print sequence along all 4 sides */
  fprintf(wastl,"%% print seqA along all 4 sides\n");
  fprintf(wastl,"0 1 lenA 1 sub {\n");
  fprintf(wastl,"    dup\n");
  fprintf(wastl,"    0.7 add -0.3 lenA lenB sub add moveto\n");
  fprintf(wastl,"    seqA exch 1 getinterval\n");
  fprintf(wastl,"    show\n");
  fprintf(wastl,"} for\n");
  fprintf(wastl,"\n");
  fprintf(wastl,"0 1 lenA 1 sub {\n");
  fprintf(wastl,"    dup\n");
  fprintf(wastl,"    0.7 add 0.7 lenA add moveto\n");
  fprintf(wastl,"    seqA exch 1 getinterval\n");
  fprintf(wastl,"    show\n");
  fprintf(wastl,"} for\n\n");

  fprintf(wastl,"90  rotate\n");
  fprintf(wastl,"0 1 lenB 1 sub {\n");
  fprintf(wastl,"    dup\n");
  fprintf(wastl,"    lenA lenB sub add\n");  
  fprintf(wastl,"    0.7 add -0.2 moveto\n");
  fprintf(wastl,"    lenB 1 sub exch sub\n");
  fprintf(wastl,"    seqB exch 1 getinterval\n");
  fprintf(wastl,"    show\n");
  fprintf(wastl,"} for\n");
  fprintf(wastl,"270 rotate\n\n");

  fprintf(wastl,"270 rotate\n");
  fprintf(wastl,"0 1 lenB 1 sub {\n");
  fprintf(wastl,"    dup\n");
  fprintf(wastl,"    -0.3 add lenA sub  0.7 lenA add  moveto\n");
  fprintf(wastl,"    seqB exch 1 getinterval\n");
  fprintf(wastl,"    show\n");
  fprintf(wastl,"} for\n");
  fprintf(wastl,"90 rotate\n\n");

  /* do grid */
  fprintf(wastl,"0.5 dup translate\n");
  fprintf(wastl,"%%draw grid\n"
	  "0.01 setlinewidth\n"
	  "lenA log 0.9 sub cvi 10 exch exp  %% grid spacing\n"
	  "dup 1 gt {\n"
	  "   dup dup 20 div dup 2 array astore exch 40 div setdash\n"
	  "} { [0.3 0.7] 0.1 setdash } ifelse\n"
	  "0 exch lenA {\n"      /* for (i=0; i<=len; i++) */
	  "   dup dup\n"        
	  "   lenA lenB sub moveto\n"                     /* i 0 moveto   */
	  "   lenA lineto \n"                  /* i len lineto */
	  "   dup lenB gt {} {\n"
	  "    dup\n"
	  "    lenA exch sub 0 exch moveto\n"   /* 0 i moveto   */
	  "    lenA exch lenA exch sub lineto\n" /* len i lineto */
	  "   } ifelse\n"
	  "   stroke\n"
	  "} for\n"
	  "0.5 neg dup translate\n\n");

  /* print upper triangle */
  for (i=1; i<=la; i++)
    for (j=1; j<=lb; j++) {
      if (P[i][j]<PMIN) continue;
      tmp = sqrt(P[i][j]);
      fprintf(wastl,"%d %d %1.5f ubox\n", j, i, tmp);
    }

  fprintf(wastl,"showpage\n");
  /*  fprintf(wastl,"end\n"); */
  fprintf(wastl,"%%%%EOF\n");
  fclose(wastl);
  
  return 1; /* success */
}

static char *time_stamp (void) {
  time_t  cal_time;
  cal_time = time(NULL);
  return ( ctime(&cal_time) );
}




































