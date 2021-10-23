
main: matrixCalculation.o
	gcc ./lib/libMatrix.a main.c -o main;


matrixCalculation.o: basicCalculation.o
	gcc -I headers -c src/matrixCalculation.c -o lib/matrixCalculation.o
	ar rcs lib/libMatrix.a lib/matrixCalculation.o
	# nm lib/libMatrix.a
	

basicCalculation.o: 
	gcc -I headers -c src/basicCalculation.c -o lib/basicCalculation.o


clean:
	rm -d main lib/*.o lib/*.a
