
main: CMatLib.o
	gcc ./lib/libcml.a main.c -o main;


CMatLib.o: basic.o
	gcc -I headers -c src/CMatLib.c -o lib/CMatLib.o
	ar rcs lib/libcml.a lib/CMatLib.o
	# nm lib/libcml.a
	

basic.o: 
	gcc -I headers -c src/basic.c -o lib/basic.o


clean:
	rm -d main lib/*.o lib/*.a
