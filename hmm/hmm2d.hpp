#ifndef AIMLESS_HMM2D_H
#define AIMLESS_HMM2D_H

#include "hmm.hpp"

namespace aimless {

    template < class energy_t = int>
    class HiddenMarkovModel2d {

        protected:

            size_t nx, ny, nstate;

        public:

            HiddenMarkovModel2d(size_t nx, size_t ny,
                    size_t nstate): nx(nx), ny(ny),
                    nstate(nstate) {
            }


            virtual energy_t transition_cost(size_t x, size_t y,
                    size_t s0, size_t s1, bool dirx) { return 0; }

            virtual energy_t data_cost(size_t x, size_t y,
                    size_t s0) { return 0; }

            class RowHmm: public HiddenMarkovModel<energy_t> {
                    private:
                        size_t row;
                        HiddenMarkovModel2d &parent;

                    public:

                        RowHmm(size_t nchain, size_t nstate, size_t row,
                                HiddenMarkovModel2d &parent):
                            HiddenMarkovModel<energy_t>(nchain,nstate),row(row),
                            parent(parent){}

                        energy_t transition_cost(size_t i, size_t s0, size_t s1) {
                            return parent.transition_cost(i,row,s0,s1,true);
                        }

                        energy_t data_cost(size_t i, size_t s0) {
                            return parent.data_cost(i,row,s0);
                        }

            };
            
            template < class OutputStateIterator, class OutputEnergyIterator  >
            void get_map_configuration( OutputStateIterator begin,
                        OutputStateIterator end, OutputEnergyIterator ebegin,
                        OutputEnergyIterator eend, size_t xk, size_t yk ) {

                {
                    size_t upperbound = std::pow(nstate,nx);
                    if ( xk > upperbound ) {
                        xk = upperbound;
                    }
                }

                std::vector< std::vector<size_t> > rowstates;
                {
                    std::vector<size_t> proto(nx);
                    rowstates.assign(ny,proto);
                }

                std::vector< std::vector<energy_t> > rowenergies;
                {
                    std::vector<energy_t> proto(nx);
                    rowenergies.assign(ny,proto);
                }

                for ( size_t iy = 0; iy < ny; ++iy ) {
                    RowHmm rhmm(nx,nstate,iy,*this);
                    rhmm.get_map_configuration(rowstates.begin(),
                            rowstates.end(), rowenergies.begin(),
                            rowenergies.end(), xk);
                }

            }

    };

}

#endif
