# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include Makefile.arch
-include ../MyConfig.mk

# Set the location of KNITRO.
KNDIR = /home/mendi/Programs/knitro-6.0.0-z
KNLIB_DYNAMIC = $(KNDIR)/lib/libknitro.so

#------------------------------------------------------------------------------

ifeq ($(shell $(ROOTCONFIG) --has-mathcore),yes)

MINLKHDO   = Min_Lkhd.$(ObjSuf) 
MINLKHDS   = Min_Lkhd.$(SrcSuf) 

PROBMAINO   = probmain.$(ObjSuf) 
PROBMAINS   = probmain.$(SrcSuf) 

PROBMAINMPIO   = probmain_mpi.$(ObjSuf) 
PROBMAINMPIS   = probmain_mpi.$(SrcSuf) 

PROBO   = prob.$(ObjSuf) 
PROBS   = prob.$(SrcSuf) 

MINLKHDLIBS = -lMathCore -lMinuit 


TMINLKHDSRC = TMinLkhdCint.$(SrcSuf)
TMINLKHDOBJ = TMinLkhdCint.$(ObjSuf)
TMINLKHDLIB = libTMinLkhdCint.$(DllSuf)

MINLKHD     = Min_Lkhd$(ExeSuf)
PROB        = prob$(ExeSuf)
PROBMPI     = prob_mpi$(ExeSuf)

endif

OBJS          = $(MINLKHDO) $(PROBO)

PROGRAMS      = $(MINLKHD) $(PROB) $(PROBMPI)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)
		@echo "$@ done"

$(TMINLKHDLIB):	$(TMINLKHDOBJ)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(LIBS) $(MINLKHDLIBS) $(KNLIB_DYNAMIC) $(OutPutOpt)$@
		$(MT_DLL)

$(MINLKHD):   $(MINLKHDO) $(TMINLKHDOBJ)
		$(LD) $(LDFLAGS) $(MINLKHDO) $(TMINLKHDOBJ) $(LIBS) $(MINLKHDLIBS)  $(OutPutOpt)$@
		@echo "$@ done"

$(PROB):      $(PROBO) $(PROBMAINO) $(TMINLKHDOBJ) 
		$(LD) $(LDFLAGS) $(PROBO) $(PROBMAINO) $(TMINLKHDOBJ) $(KNLIB_DYNAMIC) $(LIBS) $(MINLKHDLIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(PROBMPI):      $(PROBO) $(PROBMAINMPIO) $(TMINLKHDOBJ)
		$(LD) $(LDFLAGS) $(PROBO) $(PROBMAINMPIO) $(TMINLKHDOBJ) $(KNLIB_DYNAMIC) $(LIBS) $(MINLKHDLIBS)  $(OutPutOpt)$@
		@echo "$@ done"

clean:
		@rm -f $(OBJS) $(TMINLKHDSRC) core

distclean:      clean
		-@rm -f $(PROGRAMS) *Dict.* *.def *.exp *.o *Cint.* \
		   *.root *.ps *.so *.lib *.dll *.d *.log .def so_locations \
		   prob prob_mpi
		@rm -rf cxx_repository files

.SUFFIXES: .$(SrcSuf)

###


TMinLkhdCint.$(SrcSuf): TMinLkhd.cc TModel.cc gaushermite.C BitFunc.cc TMinLkhdLinkDef.h 
	@echo "Generating dictionary $@ using rootcint ..."
	@echo "rootcint -f $@ -c $^"	
	@rootcint -f $@ -c -I$(KNDIR)/include $^	

.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c -I$(KNDIR)/include $<
