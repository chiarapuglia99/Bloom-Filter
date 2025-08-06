KmerCo.exe: KmerCo.o murmur.o prime.o initBF.o mask8.o
	gcc KmerCo.o murmur.o prime.o initBF.o mask8.o -o KmerCo.exe


KmerCo.o: KmerCo.c
	gcc -c KmerCo.c -o KmerCo.o


murmur.o: murmur.c
	gcc -c murmur.c -o murmur.o


prime.o: prime.c
	gcc -c prime.c -o prime.o

initBF.o: initBF.c
	gcc -c initBF.c -o initBF.o


mask8.o: mask8.c
	gcc -c mask8.c -o mask8.o


clean:
	del KmerCo.o murmur.o initBF.o mask8.o KmerCo.exe