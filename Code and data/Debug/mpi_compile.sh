
mpif90 -O0 -g -Wall -c -fmessage-length=0 -ffpe-trap=invalid,zero,overflow -finit-local-zero -x f95-cpp-input -DFORTRANUNDERSCORE -DUSE_MPI -o "src/MOD_MEND_TYPE.o" "../src/MOD_MEND_TYPE.f90"
mpif90 -O0 -g -Wall -c -fmessage-length=0 -ffpe-trap=invalid,zero,overflow -finit-local-zero -x f95-cpp-input -DFORTRANUNDERSCORE -DUSE_MPI -o "src/MOD_USRFS.o" "../src/MOD_USRFS.f90"
mpif90 -O0 -g -Wall -c -fmessage-length=0 -ffpe-trap=invalid,zero,overflow -finit-local-zero -x f95-cpp-input -DFORTRANUNDERSCORE -DUSE_MPI -o "src/MOD_MEND.o" "../src/MOD_MEND.f90"
mpif90 -O0 -g -Wall -c -fmessage-length=0 -ffpe-trap=invalid,zero,overflow -finit-local-zero -x f95-cpp-input -DFORTRANUNDERSCORE -DUSE_MPI -o "src/MOD_OPT_TYPE.o" "../src/MOD_OPT_TYPE.f90"
mpif90 -O0 -g -Wall -c -fmessage-length=0 -ffpe-trap=invalid,zero,overflow -finit-local-zero -x f95-cpp-input -DFORTRANUNDERSCORE -DUSE_MPI -o "src/MOD_STRING.o" "../src/MOD_STRING.f90"
mpif90 -O0 -g -Wall -c -fmessage-length=0 -ffpe-trap=invalid,zero,overflow -finit-local-zero -x f95-cpp-input -DFORTRANUNDERSCORE -DUSE_MPI -o "src/MEND_IN.o" "../src/MEND_IN.f90"
mpif90 -O0 -g -Wall -c -fmessage-length=0 -ffpe-trap=invalid,zero,overflow -finit-local-zero -x f95-cpp-input -DFORTRANUNDERSCORE -DUSE_MPI -o "src/MOD_MCMC.o" "../src/MOD_MCMC.f90"
mpif90 -O0 -g -Wall -c -fmessage-length=0 -ffpe-trap=invalid,zero,overflow -finit-local-zero -x f95-cpp-input -DFORTRANUNDERSCORE -DUSE_MPI -o "src/MOD_OPT.o" "../src/MOD_OPT.f90"
mpif90 -O0 -g -Wall -c -fmessage-length=0 -ffpe-trap=invalid,zero,overflow -finit-local-zero -x f95-cpp-input -DFORTRANUNDERSCORE -DUSE_MPI -o "src/MEND_main.o" "../src/MEND_main.f90"

mpif90  -o "MENDcn"  ./src/MEND_IN.o ./src/MEND_main.o ./src/MOD_MCMC.o ./src/MOD_MEND.o ./src/MOD_MEND_TYPE.o ./src/MOD_OPT.o ./src/MOD_OPT_TYPE.o ./src/MOD_STRING.o ./src/MOD_USRFS.o   