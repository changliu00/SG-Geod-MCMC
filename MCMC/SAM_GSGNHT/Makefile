CC=g++
CPPFLAGS=-O3 -std=c++11 -fopenmp
EXTER=../MyMath.o ../SAM_Base.o ../myinclude/myutils_3.o

samgsgnht: SAM_GSGNHT.o samgsgnht.cpp $(EXTER)
	$(CC) SAM_GSGNHT.o samgsgnht.cpp $(EXTER) -o samgsgnht $(CPPFLAGS)
SAM_GSGNHT.o: SAM_GSGNHT.cpp SAM_GSGNHT.hpp ../GMC.hpp ../GSGNHT.hpp $(EXTER)
	$(CC) -c SAM_GSGNHT.cpp $(CPPFLAGS)

.PHONY: clean
clean:
	-rm *.o samgsgnht

