# original file copied from /home/yoj/src/mech/BLD_2, 1/31/2008
# need to source scom command for compilers

  MODEL = CHEMMECH_ggb_SS4


FC = ifort
CC = icc
#F_FLAGS = -fixed -132 -O3 -fno-alias -mp1  -I . -g
#  LINK_FLAGS = $(myLINK_FLAG)
#C_FLAGS =  -O2  -DFLDMN=1
 f_FLAGS   = -fixed -132 -O0 -fpe0 -check uninit -check bounds -check format -g -traceback -fno-alias -mp1  -I .
 F_FLAGS   = $(f_FLAGS)
 f90_FLAGS = -free -132 -O0 -fpe0 -check uninit -check bounds -check format -g -traceback -fno-alias -mp1  -I .
 F90_FLAGS = $(f90_FLAGS)
 C_FLAGS   =  -O2  -DFLDMN=1

#LINK_FLAGS = -i-static
 LINK_FLAGS =


#FC = pgf90
#CC = pgcc

#F_FLAGS = -Mfixed -Mextend -O3 -I . -g -C
#C_FLAGS =  -O2  -DFLDMN=1
#LINK_FLAGS = -Bstatic  -Bstatic_pgi

# compiler options for subroutines

 
OBJECTS =\
    chemmech.o \
    eatcom.o \
    check_ss_spc.o \
    getchar.o \
    getlabel.o \
    get_operator.o \
    getprdct.o \
    getrate.o \
    getrctnt.o \
    getreal.o \
    getword.o \
    get_ss_data.o \
    junit.o \
    index1.o \
    lkupspec.o \
    rdline.o \
    size.o \
    wrexts.o \
    wrhdr1.o \
    wrbf6.o \
    wrbf16c.o \
    wrbf12d.o \
    wrbf12s.o \
    wrspecial_ext.o \
    wrss_ext.o  \
    vallabchr.o \
    valwrdchr.o \
    nameval.o
        
.SUFFIXES: .f .c

$(MODEL): $(OBJECTS)
	$(FC) $(LINK_FLAGS) $(OBJECTS) -o $@

.f.o:
	$(FC) $(F_FLAGS) -c $<

.c.o:
	$(CC) $(C_FLAGS) -c $<

clean:
	\rm -f *.o
