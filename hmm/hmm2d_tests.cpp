#include "hmm2d.hpp"
#include <iostream>
#include <cstdlib>

using aimless::HiddenMarkovModel2d;

typedef int energy_t; 

//Derive from HMM class to make custom HMM with no runtime penalty
template < class energy_t = energy_t>
class TestHmm: public HiddenMarkovModel2d<energy_t> {

    public:

        TestHmm(size_t nx, size_t ny, size_t nstate):
            HiddenMarkovModel2d<energy_t>(nx,ny,nstate) {}

        energy_t transition_cost(size_t x, size_t y, size_t s0,
                size_t s1, bool dirx) {
            return 0;
        }

        energy_t data_cost(size_t x, size_t y, size_t s0) {
            std::cout << " calling data cost " << std::endl;
            return 0;
        }

#if 0 
        template < class Iterator >
        energy_t get_configuration_energy( Iterator begin,
                Iterator end ) {
            energy_t energy = 0;
            Iterator cur, prev;
            size_t i = 0;

            cur = begin;
            prev = begin;
            energy += data_cost(i,*cur);
            cur++; i++;

            for ( ; cur != end; ++prev, ++cur, ++i ) {
                energy_t transitionterm, dataterm;
                dataterm = data_cost(i,*cur);
                transitionterm = transition_cost(i,*prev,*cur);
                energy += dataterm + transitionterm;
            }

            return energy;
        }

        template < class OutputIterator >
        energy_t brute_force_get_map_configuration( OutputIterator begin,
                OutputIterator end ) {

        }

        template < class Container, class Buffer > 
        void enumerate_powerset( size_t pos, Buffer &buf,
                Container &powerset ) {

            size_t nchain = this->nchain-1;
            size_t nstate = this->nstate;

             /* Handle base case */
            if ( pos == nchain ) {
                powerset.push_back(buf); /* Copies the buffer into powerset */
                return;
            }

            /* Recurse */
            for ( size_t istate = 0; istate < nstate; ++istate ) {
                buf[pos] = istate;
                this->enumerate_powerset(pos+1,buf,powerset);
            }

        }
#endif

};

int main( int argc, char * argv[] ) {

    size_t nx = 3;
    size_t ny = 3;
    size_t nstate = 2;


    TestHmm<> hmm2d(nx,ny,nstate);

    std::vector<size_t>


    return EXIT_SUCCESS;
}
