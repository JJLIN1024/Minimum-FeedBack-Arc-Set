# CC and CFLAGS are varilables
CC = g++
CFLAGS = -O3 -std=c++11 -Wall
EXEC = cb

$(EXEC): src/cycleBreaking_new.cpp
	@$(CC) $< -o $@ $(CFLAGS)

clean:
	@rm -f $(EXEC)