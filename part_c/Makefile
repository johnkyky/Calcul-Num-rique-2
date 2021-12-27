##########################################
# Makefile                               #
# Makefile for the code developed in TP1 #
#                                        #
# T. Dufaud                              #
##########################################
################################
# Variables for this makefile
################################
# 
CC=gcc

# 
# -- Compiler Option
#
OPTC=-O3 -fomit-frame-pointer -fPIC -mavx -DAdd_ -DF77_INTEGER=int -DStringSunStyle

#
# -- Directories
TPDIR=.
TPDIRSRC=$(TPDIR)/src

#
# -- librairies
LIBS=-llapacke -lblas -lm

# -- Include directories
INCLBLASLAPACK= -I /usr/include/x86_64-linux-gnu/

INCL= -I $(TPDIR)/include $(INCLBLASLAPACK) 
#
#################################################################
# makefile
############
#
OBJENV= tp_env.o
OBJTP2ITER= lib_poisson1D.o tp2_poisson1D_iter.o
OBJTP2DIRECT= lib_poisson1D.o tp2_poisson1D_direct.o
#

all: bin/tp_testenv bin/tp2poisson1D_iter bin/tp2poisson1D_direct bin/dgbmv bin/lu_trid bin/richardson_jacobi_gauss

testenv: bin/tp_testenv

tp2poisson1D_iter: bin/tp2poisson1D_iter

tp2poisson1D_direct: bin/tp2poisson1D_direct

tp_env.o: $(TPDIRSRC)/tp_env.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp_env.c 

lib_poisson1D.o: $(TPDIRSRC)/lib_poisson1D.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/lib_poisson1D.c 

tp2_poisson1D_iter.o: $(TPDIRSRC)/tp2_poisson1D_iter.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_poisson1D_iter.c  

tp2_poisson1D_direct.o: $(TPDIRSRC)/tp2_poisson1D_direct.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_poisson1D_direct.c  

dgbmv.o: $(TPDIRSRC)/dgbmv.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/dgbmv.c

lu_trid.o: $(TPDIRSRC)/lu_trid.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/lu_trid.c

richardson_jacobi_gauss.o: $(TPDIRSRC)/richardson_jacobi_gauss.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/richardson_jacobi_gauss.c

bin/tp_testenv: $(OBJENV) 
	$(CC) -o bin/tp_testenv $(OPTC) $(OBJENV) $(LIBS)

bin/tp2poisson1D_iter: $(OBJTP2ITER)
	$(CC) -o bin/tp2poisson1D_iter $(OPTC) $(OBJTP2ITER) $(LIBS)

bin/tp2poisson1D_direct: $(OBJTP2DIRECT)
	$(CC) -o bin/tp2poisson1D_direct $(OPTC) $(OBJTP2DIRECT) $(LIBS)

bin/dgbmv: lib_poisson1D.o dgbmv.o
	$(CC) -o bin/dgbmv $(OPTC) dgbmv.o lib_poisson1D.o $(LIBS)

bin/lu_trid: lib_poisson1D.o lu_trid.o
	$(CC) -o bin/lu_trid $(OPTC) lu_trid.o lib_poisson1D.o $(LIBS)

bin/richardson_jacobi_gauss: lib_poisson1D.o richardson_jacobi_gauss.o
	$(CC) -o bin/richardson_jacobi_gauss $(OPTC) richardson_jacobi_gauss.o lib_poisson1D.o $(LIBS)


run_testenv:
	bin/tp_testenv

run_tp2poisson1D_iter:
	bin/tp2poisson1D_iter

run_tp2poisson1D_direct:
	bin/tp2poisson1D_direct

clean:
	rm *.o bin/*
