all:
	g++ solution_io.cpp dat2sol.cpp -O3 -Wall -o dat2solb
	g++ solution_io.cpp solb2dat.cpp -O3 -Wall -o solb2dat
debug:
	#g++ binary.cpp -g3 -Wall -fsanitize=address,bounds -Wfatal-errors -fuse-ld=gold -o dat2solb
	g++ solution_io.cpp dat2sol.cpp -g3 -Wall -fsanitize=address,bounds -Wfatal-errors -o dat2solb
	g++ solution_io.cpp solb2dat.cpp -g3 -Wall -fsanitize=address,bounds -Wfatal-errors -o solb2dat
