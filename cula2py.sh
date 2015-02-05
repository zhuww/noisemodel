#echo "gcc -m64 -o $1.exe $1.c  -DNDEBUG -O3  -I$CULA_INC_PATH -L$CULA_LIB_PATH_64 -lcula_core -lcula_lapack -lcublas -lcudart"
#gcc -m64  -o $1.exe $1.c -DNDEBUG -O3 -I$CULA_INC_PATH -L$CULA_LIB_PATH_64  -lcula_core -lcula_lapack -lcublas -lcudart 
gccflags="-m64 -lm -DNDEBUG -I$CULA_INC_PATH -L$CULA_LIB_PATH_64  -lcula_core -lcula_lapack -lcublas -lcudart "
gcc -c -O3 -fPIC -shared $gccflags -o $1.o $1.c
gcc -o $1.so -shared $1.o $gccflags

