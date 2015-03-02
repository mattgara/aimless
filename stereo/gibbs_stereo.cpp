#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <ctime>
#include <omp.h>
#include <limits>

void readpgm( const std::string &fn,
        unsigned char* &data,
        int &width,
        int &height ) {

    FILE*  stream = 0;

    stream = fopen(fn.c_str(),"r");
    if ( !stream ) {
        std::cerr << " Failed to read pgm." << std::endl;
        exit(-1);
    }


    int bsize;
    int read = fscanf(stream,"P5\n%d %d\n%d%*[\n]",&width,&height,&bsize);
    if ( read < 3 ) {
        std::cerr << " Failed to read pgm header." << std::endl;
        exit(-1);
    }

    data = new unsigned char[width*height];

    read = fread(data,1,width*height,stream);
    if ( read != width*height ) {
        std::cerr << " Failed to read pgm data." << std::endl;
        std::cerr << " Read " << read << " expected " << width*height
            <<  "." << std::endl;
        exit(-1);
    }

    fclose(stream);


}

void writepgm( const std::string &fn,
        unsigned char* &data,
        int &width,
        int &height ) {

    FILE*  stream = 0;

    stream = fopen(fn.c_str(),"w");
    if ( !stream ) {
        std::cerr << " Failed to write pgm." << std::endl;
        exit(-1);
    }


    int wrote = fprintf(stream,"P5\n%d %d\n%d\n",width,height,255);
    if ( wrote < 3 ) {
        std::cerr << " Failed to write pgm header." << std::endl;
        exit(-1);
    }

    wrote = fwrite(data,1,width*height,stream);
    if ( wrote != width*height ) {
        std::cerr << " Failed to write pgm data." << std::endl;
        std::cerr << " Wrote " << wrote << " expected " << width*height
            <<  "." << std::endl;
        exit(-1);
    }

    fclose(stream);


}

typedef struct _im {

    unsigned char *data;
    int width, height;

} im_t;

inline double smoothing_term( int *displut, unsigned char *disp, int width, int
        height, int y, int x, int thisdisp, double p1, double p2 ) {

    if ( y < 0 || y >= height || x < 0 || x >= width ) {
        return 0; /* TODO: Is this the correct behaviour? */
    }

    double smoothterm = 0;
    int absdisp = std::abs(displut[disp[y*width+x]] - thisdisp);
    if ( absdisp == 0 ) {
        smoothterm = 0;
    } else if (  absdisp == 1 ) {
        smoothterm = p1;
    } else { 
        smoothterm = p2;
    }
    return smoothterm;
}

inline int sample_distribution( int n, double *prob, double norm ) {

    double unirand = (double) std::rand() / (double) RAND_MAX;

    double cumsum = 0;
    double invnorm = 1/norm;

    int reti = -1;

    for ( int i = 0; i < n; ++i ) {
        cumsum += prob[i] * invnorm;
        if ( unirand <= cumsum ) { /* This must eventually be true of any prob.
                                      dist. */
            reti = i;
            break;
        }
    }

    if ( reti < 0  ){ /* This seems to happen only in the case when the
                         comparison above fails which is extremely unlikely. */
        printf(" ****warning corner case in sampling detected:\n"
                "   cumsum=%g\n   unirand=%g\n", cumsum, unirand);
        reti = n-1; /* It almost always seems to be the case that unirand == 1
                       && cumsum == 1 but it is a precision error that fails
                       the comparison above. */
    }

    assert( reti >= 0 );

    return reti;

}

double calculate_energy( int ndisp,
        int mindisp,
        double p1,
        double p2,
        double temperature,
        int width, /* Image width */
        int height, /* Image width */
        unsigned char *dataterms,
        unsigned char *ref,
        unsigned char *match,
        unsigned char *disp /* Disparity */) {

    int* displut = new int[ndisp];
    for ( int idisp = 0; idisp < ndisp; ++idisp ) {
        displut[idisp] = idisp + mindisp; 
    }

    const double invtemperature = 1/temperature;

    double energy = 0;
    for ( int irow = 0; irow < height; ++irow ) {
        for ( int icol = 0; icol < width; ++icol ) {
            int refidx = irow * width + icol;
            int idisp = disp[irow*width + icol];
            int thisdisp = displut[idisp];
            int matchx = icol + thisdisp;
            if ( matchx < 0 || matchx >= width ) {
                continue;
            }
            int matchidx = irow*width + matchx;
            double dataterm = dataterms[(irow*width+icol)*ndisp + idisp] * invtemperature;
            double smoothterm = 0;
            smoothterm += smoothing_term(displut,disp,width,height,irow+1,icol,thisdisp,p1,p2);
            smoothterm += smoothing_term(displut,disp,width,height,irow,icol+1,thisdisp,p1,p2);
            energy += dataterm + smoothterm;
        }
    }

    delete[] displut;

    return energy;

}


void block_gibbs_iter( bool smooth,
        int partition, /* 0 or 1 */
        int ndisp,
        int mindisp,
        double p1,
        double p2,
        double temperature,
        int width, /* Image width */
        int height, /* Image width */
        unsigned char *dataterms,
        unsigned char *ref,
        unsigned char *match,
        unsigned char *disp /* Disparity */) {


    int* displut = new int[ndisp];
    for ( int idisp = 0; idisp < ndisp; ++idisp ) {
        displut[idisp] = idisp + mindisp; 
    }

    const double invtemperature = 1/temperature;



    int nthread = 1;
#pragma omp parallel
    {
        nthread = omp_get_num_threads();
    }
    double **probs = new double*[nthread];
    for ( int i = 0; i < nthread; ++i ){
        probs[i] = new double[ndisp];
    }
#pragma omp parallel for
    for ( int irow = 0; irow < height; ++irow ) {
        int omptid = omp_get_thread_num();
        double* prob = probs[omptid];
        int shift = (partition + irow) % 2; /* 0 or 1 aternating between rows */
        for ( int icol = shift; icol < width; icol += 2 ) {
            int refidx = irow * width + icol; /* Every second element */

            //Need to compute the energy for every possible disparity:
            for ( int idisp = 0; idisp < ndisp; ++idisp ) {
                int thisdisp = displut[idisp];
                int matchx = icol + thisdisp;
                if ( matchx < 0 || matchx >= width ) {
                    prob[idisp] = -1;
                    continue;
                }
                int matchidx = irow*width + matchx;
                double dataterm = dataterms[(irow*width+icol)*ndisp + idisp] * invtemperature;


                //Smoothing term of the four neighbours.
                double smoothterm = 0;
                if ( smooth ) {
                    smoothterm += smoothing_term(displut,disp,width,height,irow+1,icol,thisdisp,p1,p2);
                    smoothterm += smoothing_term(displut,disp,width,height,irow-1,icol,thisdisp,p1,p2);
                    smoothterm += smoothing_term(displut,disp,width,height,irow,icol+1,thisdisp,p1,p2);
                    smoothterm += smoothing_term(displut,disp,width,height,irow,icol-1,thisdisp,p1,p2);
                }


                prob[idisp] = dataterm + smoothterm;

            }

            double partitionfunc = 0;
            for ( int idisp = 0; idisp < ndisp; ++idisp ) {
                if ( prob[idisp] >= 0 ) {
                    prob[idisp] = std::exp( - (prob[idisp]) );
                    partitionfunc += prob[idisp];
                } else {
                    prob[idisp] = 0;
                }
            }

            assert ( partitionfunc > 0 ); /* There must be some disparity with non-zero prob.*/

            disp[ refidx ] = sample_distribution(ndisp,prob,partitionfunc);

        }
    }

    delete[] displut;
    for ( int i = 0; i < nthread; ++i ) {
        delete[] probs[i];
    }
    delete[] probs;


}


void block_gibbs( int niter,
        int ndisp,
        int mindisp,
        double _p1,
        double _p2,
        double temperature,
        unsigned char *dataterms,
        im_t &im1,
        im_t &im2,
        im_t &disp) {

    if ( !(im1.width == im2.width && im1.height == im2.height  ) ) {
        std::cerr << " Images are not same dimensions." << std::endl;
        exit(-1);
    }

    int width = im1.width;
    int height = im1.height;

    double p1,p2;
    p1 = _p1 / temperature;
    p2 = _p2 / temperature;

    unsigned char *ref   = im1.data;
    unsigned char *match = im2.data;
    unsigned char *out   = disp.data;

    for ( int iter = 0; iter < niter;
            ++iter ) {
        for ( int k = 0; k < 2; ++k ) {
            int partition = k;
            block_gibbs_iter( true, partition, ndisp, mindisp, p1, p2, temperature, width,
                    height, dataterms, ref, match, out);
        }
        double energy = calculate_energy( ndisp, mindisp, p1, p2, temperature, width,
                height, dataterms, ref, match, out );

        if ( iter % 2  == 0 ) {
            writepgm("disp_inprogress.pgm",
                    disp.data,
                    disp.width,
                    disp.height);
        }

        std::cout << " done block gibbs iter " << iter+1 << " of " <<
            niter << " energy: " << energy << std::endl;

    }

}

void calculate_dataterms( int ndisp,
        int mindisp,
        int rady,
        int radx,
        im_t &im1,
        im_t &im2,
        unsigned char *dataterms ) {

    if ( !(im1.width == im2.width && im1.height == im2.height  ) ) {
        std::cerr << " Images are not same dimensions." << std::endl;
        exit(-1);
    }

    int width = im1.width;
    int height = im1.height;

    unsigned char *ref   = im1.data;
    unsigned char *match = im2.data;

    int* displut = new int[ndisp];
    for ( int idisp = 0; idisp < ndisp; ++idisp ) {
        displut[idisp] = idisp + mindisp; 
    }

#pragma omp parallel for
    for ( int irow = 0; irow < height; ++irow ) {

        int rsizey = 2*rady+1;
        int rsizex = 2*radx+1;
        double *refpx = new double[rsizey*rsizex];
        double *matchpx = new double[rsizey*rsizex];

        for ( int icol = 0; icol < width; ++icol ) {

            int refidx = irow * width + icol;

            for ( int idisp = 0; idisp < ndisp; ++idisp ) {
                int thisdisp = displut[idisp];
                int matchx = icol + thisdisp;
                if ( matchx < 0 || matchx >= width ) {
                    continue;
                }
                int matchidx = irow*width + matchx;

                refpx[0] = ref[refidx];
                matchpx[0] = match[matchidx];

                double refaggrterm = refpx[0];
                double matchaggrterm = matchpx[0];
                int countfound = 1;
                for ( int dy = -rady; dy <= rady; ++dy ) {
                    for ( int dx = -radx; dx <= radx; ++dx ) {
                        if ( dx == 0 && dy == 0 ) {
                            continue;
                        }
                        int x0, y0, x1, y1;
                        y0 = irow + dy;
                        x0 = matchx + dx;
                        y1 = y0;
                        x1 = icol + dx;
                        if ( x0 < 0 || x0 >= width || x1 < 0 || x1 >= width ||
                                y0 < 0 || y0 >= height ) {
                            continue;
                        }
                        refpx[countfound] = ref[y1*width+x1];
                        matchpx[countfound] = match[y0*width+x0];
                        refaggrterm += refpx[countfound];
                        matchaggrterm += matchpx[countfound];
                        countfound++;
                    }
                }

                //Compute means & std
                double refmean = refaggrterm / countfound;
                double matchmean = matchaggrterm / countfound;

                double refstd = 0;
                double matchstd = 0;

                for ( int ipx = 0; ipx < countfound; ++ipx ) {
                    double t;
                    t = (refpx[ipx] - refmean );
                    refstd += t*t;
                    t = (matchpx[ipx] - matchmean );
                    matchstd += t*t;
                }
                refstd = std::sqrt(refstd);
                matchstd = std::sqrt(matchstd);

                int stdthresh = 0;
                
                double norm = 1. / refstd / matchstd;
                double ncc = 0;
                for ( int ipx = 0; ipx < countfound; ++ipx ) {
                    ncc += (refpx[ipx] - refmean) * (matchpx[ipx] - matchmean ) * norm;
                }

                if ( matchstd > stdthresh && refstd >  stdthresh ) {
                    if ( !( ncc <= 1+1e-6 && ncc >= -1-1e-6) ) {
                        std::cout << " countfound: " << countfound << std::endl;
                        std::cout << " ncc: " << ncc << std::endl;
                        printf("nccp: %lf\n",ncc);
                        std::cout << " refstd: " << refstd << std::endl;
                        std::cout << " matchstd: " << matchstd << std::endl;
                        std::cout << " refmean: " << refmean << std::endl;
                        std::cout << " matchmean: " << matchmean << std::endl;
                        std::cout << " refpx: ";
                        for ( int ipx = 0; ipx < countfound; ++ipx ) {
                            std::cout << refpx[ipx] << ", ";
                        }
                        std::cout << std::endl;
                        std::cout << " matchpx: ";
                        for ( int ipx = 0; ipx < countfound; ++ipx ) {
                            std::cout << matchpx[ipx] << ", ";
                        }
                        std::cout << std::endl;
                    }
                    assert( ncc <= 1.+1e-6 && ncc >= -1.-1e-6);

                    double output = (1-ncc) * 128;
                    output = output < 0 ? 0 : output > 255 ? 255 : output;
                    dataterms[ (irow*width + icol)*ndisp + idisp ] = output;
                } else {
                    dataterms[ (irow*width + icol)*ndisp + idisp ] = 255;
                }

            }

        }

        delete[] refpx;
        delete[] matchpx;
    }

    delete[] displut;


}


int main( int argc, char *argv[] ) {

    std::srand(time(NULL));


    im_t im1, im2, disp;

    std::cout << " Reading im1 " << std::endl;
    readpgm("cones_left.pgm",
            im1.data,
            im1.width,
            im1.height);

    std::cout << " Reading im2 " << std::endl;
    readpgm("cones_right.pgm",
            im2.data,
            im2.width,
            im2.height);

    disp.width  = im1.width;
    disp.height = im1.height;
    disp.data   = new unsigned char[disp.width*disp.height];


    bool usecheckpoint = false;
    std::string checkpointfile = "";
    if ( usecheckpoint ) {
        readpgm(checkpointfile,
                disp.data,
                disp.width,
                disp.height);
    } else {
        std::fill(disp.data,disp.data+disp.width*disp.height,0);
    }

    int ndisp = 128;
    int mindisp = -(ndisp-1);
    double temperature;       /* Temperature is unfortunately critical to
                                 methods based on block gibbs sampling. Higher
                                 temperature induces more variability and
                                 possibly non convergence, lower temperatures
                                 cause numeric instabilities as well as poor
                                 solutions that refuse to shift from their
                                 locally optimal solution. */

    std::cout << " ndisp: " << ndisp << std::endl;
    std::cout << " mindisp: " << mindisp << std::endl;

    double p1, p2; /* The model that is minimized uses the same p1, p2
                      penalties as SGM and countless other stereo algorithms.
                      */

    /* Precompute the dataterms based on radius size */
    int drady = 1;
    int dradx = 1;
    unsigned char *dataterms = new unsigned char[im1.width*im1.height*ndisp];
    std::cout << " precomputing dataterms ... " << std::endl;
    calculate_dataterms(ndisp,mindisp,drady,dradx,im1,im2,dataterms);
    std::cout << " done precomputing dataterms " << std::endl;

    int niter;

    if ( !usecheckpoint ) {
        std::cout << " calculating prior ... " << std::endl;
        p1 = 0;
        p2 = 0;
        niter = 1;
        temperature = 1.;
        block_gibbs(niter,ndisp,mindisp,p1,p2,temperature,dataterms,im1,im2,disp);
        std::cout << " done calculating prior " << std::endl;
    }

    /* A very reasonable set of settings is:
     *
     * p1 = 10; p2 = 215; temperature = 80;
     *
     */

    std::cout << " calculating posterior ... " << std::endl;
    p1 = 10;
    p2 = 215;
    niter = 2000;
    temperature = 20.;
    block_gibbs(niter,ndisp,mindisp,p1,p2,temperature,dataterms,im1,im2,disp);

    writepgm("disp.pgm",
            disp.data,
            disp.width,
            disp.height);


    delete[] im1.data;
    delete[] im2.data;
    delete[] disp.data;
    delete[] dataterms;


    return 0;
}
