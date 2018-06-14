#******************************************************************#
#                                                                  #
# Makefile for              P-MergerTree                           #
#                           Parallel Merger Tree tool for AHF	   #
#                                                            	   #
# Edoardo Carlesi (2018)                                           #
#                                                                  #
#******************************************************************#
include ../Makefile.config

.PHONY: clean

# flags for compiling all individual routines
# -------------------------------------------
CFLAGS	=	${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}

# flags for linking all parts together
# ------------------------------------
LFLAGS	=	${LNFLAGS}

# all object files
# ----------------
SRCS	=\
	ptree_io.c \
	ptree_vars.c \
	comm.c \
	common.c \
	ptree_funcs.c

OBJS	=\
	ptree_io.o \
	ptree_vars.o \
	common.o \
	comm.o \
	ptree_funcs.o

LIBS 	=\
	libio/libio.a \
	libsfc/libsfc.a \
	libutility/libutility.a \
	libahf/libahf.a \
	libamr_serial/libamr_serial.a

#========================
# Main binary executable
#========================

P-MergerTree	: P-MergerTree ${OBJS} ${LIBS}
		$(CC) $(CFLAGS) $(LFLAGS) -o P-MergerTree  ${LIBS} P-MergerTree.c  ${SRCS} ${LIBS} ${LNFLAGS}

# ----------------
# 	Clean
# ----------------

clean:
	rm -f *.o 
	cd libsfc ; ${MAKE} clean
	cd libio ; ${MAKE} clean
	cd libahf ; ${MAKE} clean
	cd libutility ; ${MAKE} clean
	cd libamr_serial ; ${MAKE} clean

# ------------------------------
#	 Make the libraries
# ------------------------------
libsfc/libsfc.a:
	cd libsfc ; ${MAKE} libsfc.a

libio/libio.a:
	cd libio ; ${MAKE} libio.a

libahf/libahf.a:
	cd libahf ; ${MAKE} libahf.a

libutility/libutility.a:
	cd libutility ; ${MAKE} libutility.a

libamr_serial/libamr_serial.a:
	cd libamr_serial ; ${MAKE} libamr_serial.a

#dependencies for .o files
# ------------------------
comm.o:	comm.c
	${CC} ${CFLAGS} -c comm.c

common.o:	common.c
	${CC} ${CFLAGS} -c common.c

ptree_io.o:	ptree_io.c
	${CC} ${CFLAGS} -c ptree_io.c

ptree_vars.o:	ptree_vars.c
	${CC} ${CFLAGS} -c ptree_vars.c

ptree_funcs.o:	ptree_funcs.c
	${CC} ${CFLAGS} -c ptree_funcs.c

P-MergerTree.o:	P-MergerTree.c
	${CC} ${CFLAGS} -c P-MergerTree.c
