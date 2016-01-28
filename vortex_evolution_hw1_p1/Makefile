cpp:
	g++ main.cpp -o main -std=c++17 -fopenmp
f2py:
	f2py -m fortutils --fcompiler=gnu95 --f90flags='-march=native -O3 -fopenmp' -lgomp -c fortutils.f90
