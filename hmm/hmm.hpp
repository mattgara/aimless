#include <limits>
#include <cstdlib>
#include <vector>
#include <stdexcept>

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


            template < class OutputIterator >
                energy_t get_map_configuration( OutputIterator begin,
                        OutputIterator end ) {

                    if ( std::distance(begin,end) < nchain-1 ) {
                        throw std::invalid_argument(
                                "Iterator must have at least enough space to"
                                " store hmm chain.");
                    }

                    for ( size_t i = 0; i < nstate; ++i ) {
                        hmm[0][i] = data_cost(0,i); /* Initialize data terms
                                                       for zeroth node */
                    }
                    
                    for ( size_t ichain = 1; ichain < nchain; ++ichain ) {
                        bool islastnode = ichain == nchain-1;
                        size_t _nstate = islastnode ?  1 : nstate;
                        for ( size_t istate = 0; istate < _nstate; ++istate ) {
                            energy_t minenergy =
                                std::numeric_limits<energy_t>::max();
                            size_t minenergyidx =
                                std::numeric_limits<size_t>::max();
                            energy_t dataterm = islastnode ? 0 :
                                data_cost(ichain,istate);
                            for ( size_t istateprev = 0; istateprev < nstate;
                                    ++istateprev ) {
                                energy_t transitionterm = islastnode ? 0 :
                                    transition_cost(ichain,istateprev,istate);
                                energy_t totalterm = dataterm + transitionterm
                                    + hmm[ichain-1][istateprev];
                                if ( totalterm < minenergy )  {
                                    minenergy = totalterm;
                                    minenergyidx = istateprev;
                                }
                            }
                            hmm[ichain][istate] = minenergy;
                            hmmpath[ichain-1][istate] = minenergyidx;
                        }
                    }


                    //Decode
                    OutputIterator previt = end;
                    OutputIterator curit = --previt;
                    *(curit--) = hmmpath[nchain-2][0];
                    for ( long long int i = nchain-3; i >= 0; i--, previt--, curit-- ) {
                        *curit = hmmpath[i][*previt];
                    }

                    return hmm[nchain-1][0]; /* This is always lowest energy.*/



            }



    };

}
