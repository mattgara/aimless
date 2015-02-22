all: sgm_stereo gibbs_stereo

sgm_stereo: sgm_stereo.cpp
	$(CXX) sgm_stereo.cpp -O3 -o sgm_stereo

gibbs_stereo: gibbs_stereo.cpp
	$(CXX) gibbs_stereo.cpp -O3 -fopenmp -o gibbs_stereo
