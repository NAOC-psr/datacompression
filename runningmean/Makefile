.PHONY: all
all: runningmean.so test
runningmean.so:runningmean.c
	gcc -shared -I/usr/include/python2.7 -o runningmean.so -fPIC runningmean.c -w -lm -O3  -DHAVE_OPENMP  -std=c99
test: test.c
	gcc -shared -I/usr/include/python2.7 -o test -fPIC test.c -w -lm -O3  -DHAVE_OPENMP  -std=c99 -I/usr/include/cfitsio/  -lcfitsio

clean:
	rm *.so

#Winson: remove -lgomp (gnu openmp library), -fopenmp(openmp option in C) is removed to optimize the performance
#gcc -shared -I/usr/include/python2.7 -o runningmean.so -fPIC runningmean.c -w -lm -O3 -lgomp -fopenmp -DHAVE_OPENMP 

#test.so:test.c
	#gcc -shared -I/usr/include/python2.7 -o test.so -fPIC test.c 
#swap_bytes.o:swap_bytes.c
	#gcc -shared -fPIC -c swap_bytes.c
#send_stuff.o:send_stuff.c swap_bytes.o
	#gcc -shared -c -fPIC send_stuff.c swap_bytes.o
