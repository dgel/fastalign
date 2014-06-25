CXX = g++
CXXFLAGS = -std=c++11 -O3 -flto -Werror -Wall -g

fast_align_ng: src/fast_align_ng.o src/model.o src/opts.o src/reader.o src/threaded_output.o
	$(CXX) $(CXXFLAGS) -o fast_align_ng src/fast_align_ng.o src/model.o src/opts.o src/reader.o src/threaded_output.o -lpthread

src/fast_align_ng.o: src/fast_align_ng.cc src/*.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ -lpthread

src/model.o: src/model.cc src/model.h src/ttables.h src/corpus.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ -lpthread

src/opts.o: src/opts.cc src/opts.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ -lpthread

src/reader.o: src/reader.cc src/reader.h src/model.h src/threaded_output.h src/corpus.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ -lpthread

src/threaded_output.o: src/threaded_output.cc src/threaded_output.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ -lpthread

