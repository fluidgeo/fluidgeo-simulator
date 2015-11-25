#diretorios

DIRFONTES = ./fontes
DIRBIN    = ./bin
DIRINC    = ./include

#opcoes de compilacao
#PETSC     = withPetsc
OMP       = # withOMP
OPT       = -O3 -pg  #-fast
OPT       = -O3  
OPT       = -O3 -g -DwithGaussSkyline #-fast
#OPT       = -fast

TEMPOS    = mostrarTempos
DEBUG     = debug

#TIPOSOLUCAO = withLagrangeano
TIPOSOLUCAO = withRT

#compiladores
CC        = gcc
CC        = icc
FC        = ifort
#FC        = gfortran

ifeq ($(FC),ifort)
   FFLAGSI  = -module ${DIRINC}  
   FOMP     = -openmp ${OPT} -g  #-fast #-g #-warn  #para ifort
   LIBBLAS  = -mkl
endif

ifeq ($(FC),gfortran)
   #FFLAGSI  = -J ${DIRINC}  
   FFLAGSI  = -M ${DIRINC}  # ou
   FFLAGSI  = -J ${DIRINC}  
   FOMP     = -fopenmp  #-finit-local-zero  #-Wall #para gfortran
   LIBBLAS  = #-lblas -llapack
endif

ifeq ($(OMP),withOMP)
   FPPFLAGS  = -D${OMP} 
   FFLAGS    = ${FFLAGS0} ${FOMP} ${OPT}
   LFLAGS    = ${LFLAGS0} ${FOMP}  ${OPT}
else
   FFLAGS    = ${FFLAGS0}  ${OPT}
   LFLAGS    = ${LFLAGS0} ${OPT}
endif

#NOMEEXECUTAVELF=simuladorF.exe
#NOMEEXECUTAVELB=simuladorB.exe
NOMEEXECUTAVEL=simulador.exe

#NOMEEXECUTAVEL=escImiscivelSkyLine.exe

ifeq ($(LAPACK),withlapack)
   FPPFLAGSL = -D${LAPACK}   
   LFLAGS0   = $(LIBBLAS) 
#   NOMEEXECUTAVEL=escImiscivelLAPACK.exe
endif

LIB=$(LIBBLAS)

ifeq ($(UMFPACK),withumfpack)
   UMFPACK_DIR  = /usr/local/lib/
   UMFPACK_INC  = /local/include/
   UMFPACKCMP  = -DLP64 -I ${UMFPACK_DIR} 
   #UMFPACKLNK  = ${DIRBIN}/umf4_f77wrapper.o -L ${UMFPACK_DIR} -lumfpack -lm -lamd -lsuitesparseconfig -lrt
   UMFPACKLNK  = -L ${UMFPACK_DIR} -lumfpack -lm -lamd -lsuitesparseconfig -lrt
   FPPFLAGSU = -D${UMFPACK}   
#   NOMEEXECUTAVEL=escImiscivelUMF.exe
endif

ifeq ($(PARDISO),withpardiso)
   FPPFLAGSP = -D${PARDISO}
#   NOMEEXECUTAVEL=escImiscivelPardiso.exe 
endif

ifeq ($(TEMPOS),mostrarTempos)
   FPPFLAGSTime = -D${TEMPOS}
endif

ifeq ($(DEBUG),debug)
   FPPFLAGSDebug = -D${DEBUG}
endif

ifeq ($(UMFPACK),withumfpack)
   NSYMCRS=withcrs
   FPPFLAGSN = -D${NSYMCRS}
endif
ifeq ($(PARDISO),withpardiso)
   NSYMCRS=withcrs
   FPPFLAGSN = -D${NSYMCRS}
endif

ifeq ($(TIPOSOLUCAO),withLagrangeano)
FPPFLAGST  = -D${TIPOSOLUCAO}
else
ifeq ($(TIPOSOLUCAO),withRT)
FPPFLAGST  = -D${TIPOSOLUCAO}
endif
endif


FFLAGS0   = ${FPPFLAGS}  ${FPPFLAGSL}  ${FPPFLAGSU} ${FPPFLAGSP} ${FPPFLAGSN}  ${FFLAGSI} ${FPPFLAGST} ${FPPFLAGSTime} ${FPPFLAGSDebug}

NOMEMODULO01=mGlobais
NOMEMODULO02=algMatricial
NOMEMODULO03=malha
NOMEMODULO04=leituraEscrita
NOMEMODULO05=funcoesDeForma
NOMEMODULO06=utilitarios
NOMEMODULO07=parametros
NOMEMODULO08=coeficientes
NOMEMODULO09=blocoMacro
NOMEMODULO10=mInputReader
NOMEMODULO11=elasticidade
NOMEMODULO12=driverFluidgeo


OBJECTS0 = ${DIRBIN}/${NOMEMODULO01}.o ${DIRBIN}/${NOMEMODULO02}.o ${DIRBIN}/${NOMEMODULO03}.o \
	   ${DIRBIN}/${NOMEMODULO04}.o ${DIRBIN}/${NOMEMODULO05}.o ${DIRBIN}/${NOMEMODULO06}.o \
	   ${DIRBIN}/${NOMEMODULO07}.o ${DIRBIN}/${NOMEMODULO08}.o ${DIRBIN}/${NOMEMODULO09}.o \
	   ${DIRBIN}/${NOMEMODULO10}.o ${DIRBIN}/${NOMEMODULO11}.o ${DIRBIN}/${NOMEMODULO12}.o

SOURCES0 = ${DIRFONTES}/${NOMEMODULO01}.F90 ${DIRFONTES}/${NOMEMODULO02}.F90 ${DIRFONTES}/${NOMEMODULO03}.F90 \
	   ${DIRFONTES}/${NOMEMODULO04}.F90 ${DIRFONTES}/${NOMEMODULO05}.F90 ${DIRFONTES}/${NOMEMODULO06}.F90 \
	   ${DIRFONTES}/${NOMEMODULO07}.F90 ${DIRFONTES}/${NOMEMODULO08}.F90 ${DIRFONTES}/${NOMEMODULO09}.F90 \
	   ${DIRFONTES}/${NOMEMODULO10}.F90 ${DIRFONTES}/${NOMEMODULO11}.F90 ${DIRFONTES}/${NOMEMODULO12}.F90 

#OBJECTSF = ${OBJECTS0} ${DIRBIN}/${NOMEMODULO07F}.o ${DIRBIN}/${NOMEMODULO08F}.o
#SOURCESF = ${SOURCES0} ${DIRFONTES}/${NOMEMODULO07F}.F90 ${DIRFONTES}/${NOMEMODULO08F}.F90

#OBJECTSB = ${OBJECTS0} ${DIRBIN}/${NOMEMODULO07B}.o ${DIRBIN}/${NOMEMODULO08B}.o
#SOURCESB = ${SOURCES0} ${DIRFONTES}/${NOMEMODULO07B}.F90 ${DIRFONTES}/${NOMEMODULO08B}.F90


#all: ${DIRBIN}/${NOMEEXECUTAVELF}  ${DIRBIN}/${NOMEEXECUTAVELB} 
all: 	${DIRBIN}/${NOMEEXECUTAVEL} 
# 	clear; 
	  
# 	clear


ifeq ($(UMFPACK),withumfpack)
${DIRBIN}/umf4_f77wrapper.o: ${DIRFONTES}/umf4_f77wrapper.c
	@echo  0, compilando ${DIRFONTES}/umf4_f77wrapper.c 
	-$(CC) -c -I ${UMFPACK_DIR} ${DIRFONTES}/umf4_f77wrapper.c -o ${DIRBIN}/umf4_f77wrapper.o
endif

${DIRBIN}/${NOMEMODULO01}.o: ${DIRFONTES}/${NOMEMODULO01}.F90
	@echo  1, compilando ${NOMEMODULO01}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO01}.o  ${DIRFONTES}/${NOMEMODULO01}.F90

${DIRBIN}/${NOMEMODULO02}.o: ${DIRFONTES}/${NOMEMODULO02}.F90
	@echo  2, compilando ${NOMEMODULO02}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO02}.o  ${DIRFONTES}/${NOMEMODULO02}.F90

${DIRBIN}/${NOMEMODULO03}.o: ${DIRFONTES}/${NOMEMODULO03}.F90
	@echo  3, compilando ${NOMEMODULO03}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO03}.o  ${DIRFONTES}/${NOMEMODULO03}.F90

${DIRBIN}/${NOMEMODULO04}.o: ${DIRFONTES}/${NOMEMODULO04}.F90
	@echo  4, compilando ${NOMEMODULO04}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO04}.o  ${DIRFONTES}/${NOMEMODULO04}.F90

${DIRBIN}/${NOMEMODULO05}.o: ${DIRFONTES}/${NOMEMODULO05}.F90
	@echo  5, compilando ${NOMEMODULO05}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO05}.o  ${DIRFONTES}/${NOMEMODULO05}.F90

${DIRBIN}/${NOMEMODULO06}.o: ${DIRFONTES}/${NOMEMODULO06}.F90
	@echo  6, compilando ${NOMEMODULO06}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO06}.o  ${DIRFONTES}/${NOMEMODULO06}.F90

${DIRBIN}/${NOMEMODULO07}.o: ${DIRFONTES}/${NOMEMODULO07}.F90
	@echo  7, compilando ${NOMEMODULO07}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO07}.o  ${DIRFONTES}/${NOMEMODULO07}.F90

${DIRBIN}/${NOMEMODULO08}.o: ${DIRFONTES}/${NOMEMODULO08}.F90
	@echo  8, compilando ${NOMEMODULO08}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO08}.o  ${DIRFONTES}/${NOMEMODULO08}.F90
	
${DIRBIN}/${NOMEMODULO09}.o: ${DIRFONTES}/${NOMEMODULO09}.F90
	@echo  9, compilando ${NOMEMODULO09}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO09}.o  ${DIRFONTES}/${NOMEMODULO09}.F90

${DIRBIN}/${NOMEMODULO10}.o: ${DIRFONTES}/${NOMEMODULO10}.F90
	@echo 10, compilando ${NOMEMODULO10}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO10}.o  ${DIRFONTES}/${NOMEMODULO10}.F90

${DIRBIN}/${NOMEMODULO11}.o: ${DIRFONTES}/${NOMEMODULO11}.F90
	@echo 11, compilando ${NOMEMODULO11}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO11}.o  ${DIRFONTES}/${NOMEMODULO11}.F90
	
${DIRBIN}/${NOMEMODULO12}.o: ${DIRFONTES}/${NOMEMODULO12}.F90
	@echo 12, compilando ${NOMEMODULO12}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO12}.o  ${DIRFONTES}/${NOMEMODULO12}.F90
	
#fratura:
${DIRBIN}/${NOMEEXECUTAVELF}: ${OBJECTSF}
	@echo  gerando executavel: ${DIRBIN}/${NOMEEXECUTAVELF} 
	-${FC} ${LFLAGS} -o ${DIRBIN}/${NOMEEXECUTAVELF}  ${OBJECTSF} ${LIB} ${UMFPACKLNK} 

#bloco:
${DIRBIN}/${NOMEEXECUTAVELB}: ${OBJECTS0}
	@echo  gerando executavel: ${DIRBIN}/${NOMEEXECUTAVELB} 
	-${FC} ${LFLAGS} -o ${DIRBIN}/${NOMEEXECUTAVELB}  ${OBJECTS0} ${LIB} ${UMFPACKLNK} 

# 2 escalas:
${DIRBIN}/${NOMEEXECUTAVEL}: ${OBJECTS0}
	@echo  gerando executavel: ${DIRBIN}/${NOMEEXECUTAVEL} 
	-${FC} ${LFLAGS} -o ${DIRBIN}/${NOMEEXECUTAVEL}  ${OBJECTS0} ${LIB} ${UMFPACKLNK} 

#clean: cleanF cleanB

cleanF:
	@echo  FRATURA ... apagando:  ${DIRBIN}/${NOMEEXECUTAVELF} ${OBJECTSF} ${DIRINC}/*.mod
	rm -f ${DIRBIN}/${NOMEEXECUTAVELF} ${OBJECTSF} ${DIRINC}/*.mod

cleanB:
	@echo  BLOCO   ... apagando:  ${DIRBIN}/${NOMEEXECUTAVELB} ${OBJECTSB} ${DIRINC}/*.mod
	echo ${OBJECTSB} ${DIRINC}/*.mod
	rm -f ${DIRBIN}/${NOMEEXECUTAVELB} ${OBJECTSB} ${DIRINC}/*.mod

clean:
# 	clear
	@echo  2Escalas   ... apagando:  ${DIRBIN}/${NOMEEXECUTAVEL} ${OBJECTSO} ${DIRINC}/*.mod
	echo ${OBJECTS0} ${DIRINC}/*.mod
	rm -f ${DIRBIN}/${NOMEEXECUTAVEL} ${OBJECTS0} ${DIRINC}/*.mod
	
