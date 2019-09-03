CC=g++
CPPFLAGS=-O3 -std=c++11
EXTER=myinclude/myutils_3.o

.PHONY: all
all: MyMath.o SAM_Base.o

SAM_Base.o: SAM_Base.hpp SAM_Base.cpp $(EXTER)
	$(CC) -c SAM_Base.cpp $(CPPFLAGS)
MyMath.o: MyMath.hpp MyMath.cpp
	$(CC) -c MyMath.cpp $(CPPFLAGS)

.PHONY: clean
clean:
	-rm *.o

