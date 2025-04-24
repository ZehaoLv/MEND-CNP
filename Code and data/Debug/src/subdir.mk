################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 

F90_SRCS += \
../src/MEND_READnml.f90 \
../src/MEND_main.f90 \
../src/MOD_MCMC.f90 \
../src/MOD_MEND.f90 \
../src/MOD_MEND_TYPE.f90 \
../src/MOD_OPT.f90 \
../src/MOD_OPT_TYPE.f90 \
../src/MOD_STRING.f90 \
../src/MOD_USRFS.f90 

OBJS += \
./src/MEND_READnml.o \
./src/MEND_main.o \
./src/MOD_MCMC.o \
./src/MOD_MEND.o \
./src/MOD_MEND_TYPE.o \
./src/MOD_OPT.o \
./src/MOD_OPT_TYPE.o \
./src/MOD_STRING.o \
./src/MOD_USRFS.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
#	$(COMP) -funderscoring -O0 -g -Wall -c -fmessage-length=0 -ffpe-trap=invalid,zero,overflow -finit-local-zero -o "$@" "$<"
	$(COMP) $(FFLAGS) $(CPPDEFS) -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/MEND_READnml.o: ../src/MEND_READnml.f90 src/MOD_MEND.o src/MOD_MEND_TYPE.o src/MOD_OPT_TYPE.o src/MOD_STRING.o src/MOD_USRFS.o

src/MEND_main.o: ../src/MEND_main.f90 src/MOD_MCMC.o src/MOD_MEND.o src/MOD_MEND_TYPE.o src/MOD_OPT.o src/MOD_OPT_TYPE.o src/MOD_STRING.o src/MOD_USRFS.o

src/MOD_MCMC.o: ../src/MOD_MCMC.f90 src/MOD_MEND.o src/MOD_OPT_TYPE.o

src/MOD_MEND.o: ../src/MOD_MEND.f90 src/MOD_MEND_TYPE.o src/MOD_USRFS.o

src/MOD_MEND_TYPE.o: ../src/MOD_MEND_TYPE.f90

src/MOD_OPT.o: ../src/MOD_OPT.f90 src/MOD_MEND.o src/MOD_OPT_TYPE.o src/MOD_USRFS.o

src/MOD_OPT_TYPE.o: ../src/MOD_OPT_TYPE.f90

src/MOD_STRING.o: ../src/MOD_STRING.f90

src/MOD_USRFS.o: ../src/MOD_USRFS.f90


