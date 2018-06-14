##################################################################################
# WHATEVER YOU DO, PLEASE DO IT TO Makefile.config AND LEAVE THIS MAKEFILE ALONE #
##################################################################################
include Makefile.config


#------------------------------------------------------------------#
# the menu                                                         #
#------------------------------------------------------------------#
MASTER_DEFINEFLAGS	=	$(DEFINEFLAGS)

#------------------------------------------------------------------#
# make settings available to all other Makefiles                   #
#------------------------------------------------------------------#
export CC
export FC
export OPTIMIZE
export CCFLAGS
export LNFLAGS
export MASTER_DEFINEFLAGS
export MAKE

# everything in src/
#========================
all:	P-MergerTree


P-MergerTree:	
		cd src;\
		${MAKE} P-MergerTree;\
		mv -f P-MergerTree ../bin


#-------------------------------------------------------------------#
# "make clean" 
#-------------------------------------------------------------------#
clean:	
	cd src; ${MAKE} clean;\

dirs:
	mkdir -p bin
