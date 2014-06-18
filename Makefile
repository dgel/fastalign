CXX = g++
CXXFLAGS = -std=c++11 -O3 -flto -Werror -Wall

fast_align_ng: src/fast_align_ng.o src/model.o src/opts.o
	$(CXX) $(CXXFLAGS) -o fast_align_ng src/fast_align_ng.o src/model.o src/opts.o

src/fast_align_ng.o: src/*.h

src/model.o: src/model.h src/ttables.h src/corpus.h

src/opts.o: src/opts.h