# Compiler:
CC=g++

# General compiler flags:
C_FLAGS=-Wall -std=c++11

# Source files:
SRC=RK_2.cpp RK_4A.cpp RK_4B.cpp RK_45.cpp RK_5.cpp RK_10.cpp integrator.cpp main.cpp 

all:
	$(CC) $(SRC) $(C_FLAGS) -o main.out

