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

            class Hmm2d: public HiddenMarkovModel<energy_t> {

                    private:

                        HiddenMarkovModel2d &parent;
                        std::vector< std::vector<energy_t> > &rowenergies;
                        std::vector< std::vector<size_t> > &rowstates;

                    public:

                        Hmm2d(size_t nchain, size_t nstate,
                                std::vector< std::vector<energy_t> > &rowenergies,
                                std::vector< std::vector<size_t> > &rowstates,
                                HiddenMarkovModel2d &parent):
                            HiddenMarkovModel<energy_t>(nchain,nstate),
                            parent(parent),
                            rowenergies(rowenergies),
                            rowstates(rowstates){}

                        energy_t transition_cost(size_t i, size_t s0, size_t s1) {
                            energy_t totalenergy = 0;
                            size_t* row_s0 = &rowstates[i-1][s0*parent.nx];
                            size_t* row_s1 = &rowstates[i][s1*parent.nx];
                            for ( size_t icol = 0; icol < parent.nx; ++icol ){
                                totalenergy += 
                                    parent.transition_cost(icol,i,
                                            row_s0[icol],row_s1[icol],false);
                            }
                            return totalenergy;
                        }

                        energy_t data_cost(size_t i, size_t s0) {
                            return rowenergies[i][s0];
                        }

            };
            
            template < class OutputStateIterator, class OutputEnergyIterator  >
            void get_map_configuration( OutputStateIterator begin,
                        OutputStateIterator end, OutputEnergyIterator ebegin,
                        OutputEnergyIterator eend, size_t xk, size_t yk ) {

                if ( std::distance(begin,end) < yk*nx*ny ) {
                    throw std::invalid_argument(
                            "State iterator must have at least enough space to"
                            " store 2d hmm.");
                }

                {
                    size_t upperbound = std::pow(nstate,nx);
                    if ( xk > upperbound ) {
                        xk = upperbound;
                    }
                }

                std::vector< std::vector<size_t> > rowstates;
                {
                    std::vector<size_t> proto(nx*xk);
                    rowstates.assign(ny,proto);
                }

                std::vector< std::vector<energy_t> > rowenergies;
                {
                    std::vector<energy_t> proto(xk);
                    rowenergies.assign(ny,proto);
                }

                for ( size_t iy = 0; iy < ny; ++iy ) {
                    RowHmm rhmm(nx,nstate,iy,*this);
                    rhmm.get_map_configuration(rowstates[iy].begin(),
                            rowstates[iy].end(), rowenergies[iy].begin(),
                            rowenergies[iy].end(), xk);
                }

                {
                    size_t upperbound = std::pow(xk,ny);
                    if ( yk > upperbound ) {
                        throw std::invalid_argument(
                                "Can not choose yk lager than xk ^ ny.");
                    }
                }

                std::vector<size_t> hmm2dstates(ny*yk);
                Hmm2d hmm2d(ny,xk,rowenergies,rowstates,*this);
                hmm2d.get_map_configuration(hmm2dstates.begin(),
                        hmm2dstates.end(),ebegin,
                        eend,yk);

                // Decode the actual states in row-major order
                OutputStateIterator it = begin;
                for (size_t iy = 0; iy < yk; ++iy ) {
                    for ( size_t irow = 0; irow < ny; ++irow  ) {
                        for ( size_t icol = 0; icol < nx; ++icol  ) {
                            *(it++) = rowstates[irow][ nx * 
                                hmm2dstates[ iy*ny + irow ] + icol ];
                        }
                    }
                }


            }

    };

}

#endif
