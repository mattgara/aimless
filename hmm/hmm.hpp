#ifndef AIMLESS_HMM_H
#define AIMLESS_HMM_H

#include <limits>
#include <cstdlib>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace aimless {

    template < class energy_t = int>
    class HiddenMarkovModel {

        protected:

            size_t nstate, nchain;
            std::vector<std::vector<energy_t> > hmm;
            std::vector<std::vector<size_t> > hmmpath;

        public:

            HiddenMarkovModel( size_t _nchain, size_t nstate  ): nstate(nstate), nchain(_nchain+1) {
                {
                    std::vector<energy_t> proto(nstate);
                    hmm.assign(_nchain,proto);
                }
                {
                    std::vector<energy_t> proto(1);
                    hmm.push_back(proto);
                }

                /* There is a difference in size of hmmpath and hmm by exactly
                 * 1, because hmmpath does not need to track the predecessor of
                 * the 0th node. */
                {
                    std::vector<size_t> proto(nstate);
                    hmmpath.assign(_nchain-1,proto);
                }
                {
                    std::vector<size_t> proto(1);
                    hmmpath.push_back(proto);
                }
            }

            /*  The cost of making transition from node i-1 to i with states s0
             *  and s1 for respective nodes. */
            virtual energy_t transition_cost(size_t i, size_t s0, size_t s1) { return 0; }

            /* Cost of state s0 of node i. */
            virtual energy_t data_cost(size_t i, size_t s0) { return 0; }


            template < class OutputStateIterator, class OutputEnergyIterator  >
                void get_map_configuration( OutputStateIterator begin,
                        OutputStateIterator end, OutputEnergyIterator ebegin,
                        OutputEnergyIterator eend, size_t k ) {

                    if ( k > std::pow(nstate,nchain-1) ) {
                        throw std::invalid_argument(
                                "Can not specify more chains to decode"
                                " than are available.");
                    }

                    if ( std::distance(ebegin,eend) < k ) {
                        throw std::invalid_argument(
                                "Energy iterator must have at least enough space to"
                                " store each hmm chain's energy.");
                    }

                    if ( std::distance(begin,end) < k*(nchain-1) ) {
                        throw std::invalid_argument(
                                "State iterator must have at least enough space to"
                                " store hmm chain.");
                    }

                    if ( k > 1 ) {

                        // Need to reallocate all datastructures to allow for
                        // enough space to work with
                        for ( size_t i = 1; i < nchain-1; ++i ) {
                            hmm[i].resize(k*nstate);
                        }
                        for ( size_t i = 0; i < nchain-2; ++i ) {
                            hmmpath[i].resize(k*nstate);
                        }

                        hmmpath[nchain-2].resize(k);
                        hmm[nchain-1].resize(k);

                    }


                    for ( size_t i = 0; i < nstate; ++i ) {
                        hmm[0][i] = data_cost(0,i); /* Initialize data terms
                                                       for zeroth node */
                    }

                    const size_t fullnstate = k*nstate;


                    for ( size_t ichain = 1; ichain < nchain; ++ichain ) {
                        bool islastnode = ichain == nchain-1;
                        bool isfirstnode = ichain == 1;
                        size_t _nstate = islastnode ?  1 : nstate;
                        for ( size_t istate = 0; istate < _nstate; ++istate ) {

                            size_t _fullnstate = isfirstnode ? nstate : fullnstate;

                            std::vector<bool> reservedidx( _fullnstate,
                                    false);

                            energy_t dataterm = islastnode ? 0 :
                                data_cost(ichain,istate);

                            for ( size_t ik = 0; ik < k; ++ik ) {

                                energy_t minenergy =
                                    std::numeric_limits<energy_t>::max();
                                size_t minenergyidx =
                                    std::numeric_limits<size_t>::max();


                                for ( size_t istateprev = 0; istateprev < _fullnstate;
                                        ++istateprev ) {


                                    if ( reservedidx[istateprev] ) {
                                        continue;
                                    }

                                    energy_t prevenergy =
                                        hmm[ichain-1][istateprev];

                                    if ( prevenergy ==
                                            std::numeric_limits<energy_t>::max() ){
                                        continue;
                                    }

                                    size_t _istateprev = isfirstnode ? istateprev :
                                        istateprev / k;

                                    energy_t transitionterm = islastnode ? 0 :
                                        transition_cost(ichain,_istateprev,istate);
                                    energy_t totalterm = dataterm + transitionterm
                                        + prevenergy;
                                    if ( totalterm < minenergy )  {
                                        minenergy = totalterm;
                                        minenergyidx = istateprev;
                                    }

                                }

                                if  ( minenergyidx ==
                                        std::numeric_limits<size_t>::max() ) {
                                    hmm[ichain][istate*k+ik] =
                                        std::numeric_limits<energy_t>::max();
                                } else {
                                    hmm[ichain][istate*k+ik] = minenergy;
                                    hmmpath[ichain-1][istate*k+ik] = minenergyidx;
                                    reservedidx[minenergyidx] = true;
                                }


                            }
                        }
                    }


                    OutputStateIterator lastend = begin+(nchain-1);
                    OutputEnergyIterator energyit = ebegin;

                    for ( size_t ik = 0; ik < k; ++ik ) {
                        //Decode
                        OutputStateIterator previt = lastend;
                        OutputStateIterator curit = --previt;
                        *(curit--) = hmmpath[nchain-2][ik];
                        for ( long long int i = nchain-3; i >= 0; i--, previt--, curit-- ) {
                            *curit = hmmpath[i][*previt];
                        }
                        curit++;
                        bool skip = true;
                        for ( OutputStateIterator it = curit; it != lastend; ++it ) {
                            if ( skip ) {
                                skip = false;
                                continue;
                            }
                            *it /= k;
                        }
                        lastend += nchain-1;
                        //Update energy
                        *(energyit++) = hmm[nchain-1][ik];
                    }


            }


            template < class OutputStateIterator>
                energy_t get_map_configuration( OutputStateIterator begin,
                        OutputStateIterator end ){

                    std::vector<energy_t> mapenergy(1);
                    get_map_configuration(begin,end,
                            mapenergy.begin(),
                            mapenergy.end(),1);

                    return mapenergy[0];

                }

    };

}

#endif
