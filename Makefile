all:
	g++ binary.cpp -O3 -Wall -DNDEBUG -o dat2solb
debug:
	#g++ binary.cpp -g3 -Wall -fsanitize=address,bounds -Wfatal-errors -fuse-ld=gold -o dat2solb
	g++ solution_io.cpp dat2sol.cpp -g3 -Wall -fsanitize=address,bounds -Wfatal-errors -o dat2solb
