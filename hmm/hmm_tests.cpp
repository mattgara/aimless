#include "hmm.hpp"
#include <iostream>
#include <cstdlib>

using aimless::HiddenMarkovModel;

typedef int energy_t; 

//Derive from HMM class to make custom HMM with no runtime penalty
template < class energy_t = energy_t>
class TestHmm: public HiddenMarkovModel<energy_t> {

    public:

        TestHmm(size_t nchain, size_t nstate):
            HiddenMarkovModel<energy_t>(nchain,nstate) {}

        /*  The cost of making transition from node i-1 to i with states s0 and
         *  s1 for respective nodes. */
        virtual energy_t transition_cost(size_t i, size_t s0, size_t s1) {
            return std::abs((int)s0 - (int)s1 - 1);
        }

        /* Cost of state s0 of node i. */
        virtual energy_t data_cost(size_t i, size_t s0) {
            return std::abs((int)s0 - (int)i);
        }

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

};

int main( int argc, char * argv[] ) {

    size_t nchain = 4;
    size_t nstate = 5;

    TestHmm<> hmm(nchain,nstate);

    std::vector<size_t> mapconfig(nchain);
    energy_t mapconfigenergy = 
        hmm.get_map_configuration(mapconfig.begin(),
                mapconfig.end());

    std::cout << " decoded MAP state of hmm,  " << std::endl;
    std::cout << "      energy:  " << mapconfigenergy << std::endl;
    std::cout << "      states:  ";

    for ( size_t i = 0; i < nchain; ++i ) {
        std::cout << mapconfig[i] << ", ";
    }
    std::cout << std::endl;

    std::cout << " independently computed energy of MAP config: " << 
        hmm.get_configuration_energy(mapconfig.begin(),
                mapconfig.end()) << std::endl;

    std::vector<size_t> buf(nchain);
    std::vector< std::vector<size_t> > powerset;
    hmm.enumerate_powerset(0,buf,powerset);

    std::cout << " # members in powerset: " << 
        powerset.size() << std::endl;

    energy_t minenergy =
        std::numeric_limits<energy_t>::max();
    for ( size_t j = 0; j < powerset.size(); ++j ) {
        std::cout << "      states:  ";
        for ( size_t i = 0; i < nchain; ++i ) {
            std::cout << powerset[j][i] << ", ";
        }
        energy_t energy = hmm.get_configuration_energy(powerset[j].begin(),
                powerset[j].end());
        std::cout << "  energy: " << energy << std::endl;
        if ( energy < minenergy ) {
            minenergy = energy;
        }
    }

    if ( minenergy == mapconfigenergy ) {
        std::cout << " ** SUCCESS ** Verified that MAP configuration energy is "
            "(one of) the minimal energy configurations." << std::endl;
    } else {
        std::cout << " ** FAILURE ** The MAP configuration energy is "
            " NOT one of the minimal energy configurations." << std::endl;
    }

    size_t k = 625;
    std::vector<energy_t> chainenergies(k);
    mapconfig.resize(k*nchain);
    hmm.get_map_configuration(mapconfig.begin(),
                mapconfig.end(),chainenergies.begin(),
                chainenergies.end(),k); /* Get top k chains */

    for ( size_t ik = 0; ik < k; ++ik ) {
        std::cout << " decoded (pseudo) MAP state " << ik
            << " of hmm,  " << std::endl;

        size_t offsetstart = ik*nchain;
        size_t offsetend = (ik+1)*nchain;

        energy_t energy = hmm.get_configuration_energy(mapconfig.begin()
                +offsetstart, mapconfig.begin()+offsetend);
        if ( energy != chainenergies[ik] ) {
            throw std::runtime_error("Energy calculations do not match.");
        }
        std::cout << "      energy: " << energy << std::endl;

        std::cout << "      states:  ";
        for ( size_t i = 0; i < nchain; ++i ) {
            std::cout << mapconfig[offsetstart+i] << ", ";
        }
        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}
