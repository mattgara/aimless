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
            return std::abs( (int) s1 - (int) s0 );
        }

        energy_t data_cost(size_t x, size_t y, size_t s0) {
            //std::cout << " calling data cost " << std::endl;
            return std::min( std::abs( (int) s0- (int) x),
                    std::abs( (int) s0- (int) y));
        }

        template < class Iterator >
        energy_t get_configuration_energy( Iterator begin,
                Iterator end ) {
            energy_t energy = 0;
            size_t i, j;
            size_t x,y;
            energy_t transitionterm, dataterm;

            size_t nx = this->nx;

            for ( size_t irow = 0; irow < this->ny; ++irow  ) {
                for ( size_t icol = 0; icol < this->nx; ++icol  ) {
                    i = irow * nx + icol;
                    dataterm = data_cost(icol,irow,*(begin+i));
                    transitionterm = 0;
                    if ( irow >= 1 ) {
                        j = (irow-1) * nx + icol;
                        transitionterm += transition_cost(icol,
                                irow,*(begin+j), *(begin+i),false);
                    }
                    if ( icol >= 1 ) {
                        j = (irow) * nx + icol-1;
                        transitionterm += transition_cost(icol,
                                irow,*(begin+j), *(begin+i),true);
                    }
                    energy += dataterm + transitionterm;
                }
            }


            return energy;
        }

        template < class Container, class Buffer > 
        void enumerate_powerset( size_t pos, Buffer &buf,
                Container &powerset ) {

            size_t nchain = this->nx * this->ny;
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

};

int main( int argc, char * argv[] ) {

    size_t nx = 4;
    size_t ny = 4;
    size_t nstate = 2;

    size_t nchain = nx * ny;

    TestHmm<> hmm2d(nx,ny,nstate);

    std::vector<size_t> buf(nchain);
    std::vector< std::vector<size_t> > powerset;
    hmm2d.enumerate_powerset(0,buf,powerset);

    std::cout << " # members in powerset: " << 
        powerset.size() << std::endl;

    energy_t minenergy =
        std::numeric_limits<energy_t>::max();
    for ( size_t j = 0; j < powerset.size(); ++j ) {
        std::cout << "      states:  ";
        for ( size_t i = 0; i < nchain; ++i ) {
            std::cout << powerset[j][i] << ", ";
        }
        energy_t energy = hmm2d.get_configuration_energy(powerset[j].begin(),
                powerset[j].end());
        std::cout << "  energy: " << energy << std::endl;
        if ( energy < minenergy ) {
            minenergy = energy;
        }
    }

    size_t xk = 6;
    size_t yk = 10;
    std::vector<size_t> mapstates(yk*nchain);
    std::vector<energy_t> mapenergies(yk);

    hmm2d.get_map_configuration(mapstates.begin(),
            mapstates.end(),
            mapenergies.begin(),
            mapenergies.end(),xk,yk);


    for ( size_t iyk = 0; iyk < yk; ++iyk ) {
        std::cout << " decoded MAP state of hmm,  " << std::endl;
        std::cout << "      energy:  " << mapenergies[iyk] << std::endl;
        std::cout << "      states:  ";

        for ( size_t i = 0; i < nchain; ++i ) {
            std::cout << mapstates[ iyk*nchain + i] << ", ";
        }
        std::cout << std::endl;
    }




    return EXIT_SUCCESS;

}
