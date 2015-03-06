all: hmm_tests hmm2d_tests

hmm_tests: hmm_tests.cpp hmm.hpp
	$(CXX) hmm_tests.cpp -g -o hmm_tests

hmm2d_tests: hmm2d_tests.cpp hmm2d.hpp hmm.hpp
	$(CXX) hmm2d_tests.cpp -g -o hmm2d_tests

