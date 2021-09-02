#Makefile

CPPFLAGS =	-I/cygdrive/c/Users/#ballin/Documents/Libraries/boost_1_55_0/ -L /cygdrive/c/Users/#ballin/Documents/Libraries/boost_1_55_0/libs/ -Wall -O3 
CPP = g++

all:	lipidsdrugsacidandsolvent

lipidsdrugsacidandsolvent:
	${CPP} -o lipidsdrugsacidandsolvent_153 lipidsdrugsacidandsolvent_153.cpp ${CPPFLAGS}
	
clean:
	rm -f lipidsdrugsacidandsolvent_153