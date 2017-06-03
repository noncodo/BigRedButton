bin_PROGRAMS=probA

probA_SOURCES=pf_main.c 

lib_LIBRARIES = libprobA.a

libprobA_a_SOURCES = matrices.c mt.c gotoh_input.c options.c options.h pfgoto.c pfgoto.h utils.c getopt.c getopt1.c getopt.h psplot.h psplot.c utils.h

#  uncomment the following if proba requires the math library
probA_LDADD= libprobA.a -lm

EXTRA_DIST= probA.1

#  if you write a self-test script named `chk', uncomment the
#  following and add `chk' to the EXTRA_DIST list
#TESTS=chk

#  build and install the .info pages
info_TEXINFOS = probA_lib.texi
probA_TEXINFOS = gpl.texinfo
    
#  install the man pages
man_MANS=probA.1
