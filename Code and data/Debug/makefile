################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init
RM := rm -rf

# to use "ifdef", add option "-x f95-cpp-input"
FFLAGS += -O0 -g -Wall -c -fmessage-length=0 -ffpe-trap=invalid,zero,overflow -finit-local-zero -x f95-cpp-input
CPPDEFS += -DFORTRANUNDERSCORE

USE_MPI = TRUE

ifeq ($(strip $(USE_MPI)),TRUE)
   COMP = /usr/bin/mpif90
   CPPDEFS += -DUSE_MPI
   MENDexe = MENDcn_mpi
else
   COMP = /usr/bin/gfortran
   MENDexe = MENDcn
endif

# All of the sources participating in the build are defined here
-include sources.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: MENDcn

# Tool invocations
MENDcn: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: MacOS X Fortran Linker'
#	$(COMP)  -o "MENDcn" $(OBJS) $(USER_OBJS) $(LIBS)
	$(COMP)  -o $(MENDexe) $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) $(MENDexe)
	-$(RM) ./*.mod
	-@echo ' '

.PHONY: all clean dependents

-include ../makefile.targets
