# Compiler:
CC=g++

# General compiler flags:
C_FLAGS=-Wall -std=c++11

# Source files:
SRC=integrator.cpp main.cpp 

all:
	$(CC) $(SRC) $(C_FLAGS) -o main.out

