
all : ConstructMatricies

.PHONY : clean

clean :
	rm -f ConstructMatricies

ConstructMatricies : ConstructMatricies.cpp
	g++ -std=gnu++11 -O2 -Wall -Wextra -o ConstructMatricies ConstructMatricies.cpp
