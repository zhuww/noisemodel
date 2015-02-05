echo "gcc -m64 -o $1.exe $1.c  -DNDEBUG -O3  -I$CULA_INC_PATH -L$CULA_LIB_PATH_64 -lcula_core -lcula_lapack -lcublas -lcudart"
gcc -m64  -o $1.exe $1.c -DNDEBUG -O3 -I$CULA_INC_PATH -L$CULA_LIB_PATH_64  -lcula_core -lcula_lapack -lcublas -lcudart 
