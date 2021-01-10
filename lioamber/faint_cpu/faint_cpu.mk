################################################################################
INCLUDES :=
INCLUDES += subm_int1.o subm_int2.o subm_int3lu.o subm_int3mem.o subm_intsol.o \
            subm_int1G.o subm_intSG.o subm_intsolG.o subm_int2G.o subm_int3G.o \
            subm_intECPG.o subm_intECP.o subm_intECPG.o                        \
				subm_intdip.o subm_intfld.o

OBJECTS += $(INCLUDES)
$(OBJPATH)/subm_int3G.o: $(OBJPATH)/subm_int2G.o
$(OBJPATH)/subm_intECPG.o : $(OBJPATH)/subm_intECP.o
$(OBJPATH)/faint_cpu.o : $(INCLUDES:%.o=$(OBJPATH)/%.o) faint_cpu.mk
################################################################################
