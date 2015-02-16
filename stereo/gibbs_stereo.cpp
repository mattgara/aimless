#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <vector>
#include <stdexcept>
#include <cassert>


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

inline int smoothing_term( int *displut, unsigned char *disp, int width, int
        height, int y, int x, int thisdisp, int p1, int p2 ) {

    if ( y < 0 || y >= height || x < 0 || x >= width ) {
        return 0; /* TODO: Is this the correct behaviour? */
    }

    int smoothterm = 0;
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

    //std::cout << " prob: ";
    //for ( int i = 0; i < n; ++i ) {
    //    double p = prob[i] * invnorm;
    //    std::cout << p << ", ";
    //}
    //std::cout << std::endl;

    for ( int i = 0; i < n; ++i ) {
        cumsum += prob[i] * invnorm;
        if ( unirand < cumsum ) { /* This must eventually be true of any prob. dist. */
            reti = i;
            break;
        }
    }

    assert( reti >= 0 );

    return reti;

}


void block_gibbs_iter( bool smooth,
        int partition, /* 0 or 1 */
        int ndisp,
        int mindisp,
        int p1,
        int p2,
        int width, /* Image width */
        int height, /* Image width */
        unsigned char *ref,
        unsigned char *match,
        unsigned char *disp /* Disparity */ ) {

    int* displut = new int[ndisp];
    for ( int idisp = 0; idisp < ndisp; ++idisp ) {
        displut[idisp] = idisp + mindisp; 
    }

    double* prob = new double[ndisp];

    double D = 1/ 20.; /* A normalizer */

    for ( int irow = 0; irow < height; ++irow ) {
        int shift = (partition + irow) % 2; /* 0 or 1 aternating between rows */
        for ( int icol = shift; icol < width; icol += 2 ) {
            int refidx = irow * width + icol; /* Every second element */
            int maxenergy, minenergy;
            int numrejected = 0;
            maxenergy = INT_MIN;
            minenergy = INT_MAX;
            //Need to compute the energy for every possible disparity:
            for ( int idisp = 0; idisp < ndisp; ++idisp ) {
                int thisdisp = displut[idisp];
                int matchx = icol + thisdisp;
                if ( matchx < 0 || matchx >= width ) {
                    prob[idisp] = -1;
                    numrejected++;
                    continue;
                }
                int matchidx = irow*width + matchx;
                int dataterm = std::abs((int)(ref[refidx]) - (int)(match[matchidx]));

                //Smoothing term of the four neighbours.
                int smoothterm = 0;
                if ( smooth ) {
                    smoothterm += smoothing_term(displut,disp,width,height,irow+1,icol,thisdisp,p1,p2);
                    smoothterm += smoothing_term(displut,disp,width,height,irow-1,icol,thisdisp,p1,p2);
                    smoothterm += smoothing_term(displut,disp,width,height,irow,icol+1,thisdisp,p1,p2);
                    smoothterm += smoothing_term(displut,disp,width,height,irow,icol-1,thisdisp,p1,p2);
                }

                int totalterm = dataterm + smoothterm;

                assert( totalterm >= 0 );

                if ( totalterm < minenergy ) {
                    minenergy = totalterm;
                }
                if ( totalterm > maxenergy ) {
                    maxenergy = totalterm;
                }
                prob[idisp] = totalterm * D;
            }

            //if ( smooth ) {
            //    std::cout << " min energy, max energy: " <<
            //        minenergy << ", " << maxenergy << std::endl;
            //}

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
    delete[] prob;


}


void block_gibbs( int niter1,
        int niter2,
        int ndisp,
        int mindisp,
        int p1,
        int p2,
        im_t &im1,
        im_t &im2,
        im_t &disp) {

    if ( !(im1.width == im2.width && im1.height == im2.height  ) ) {
        std::cerr << " Images are not same dimensions." << std::endl;
        exit(-1);
    }

    int width = im1.width;
    int height = im1.height;

    unsigned char *ref   = im1.data;
    unsigned char *match = im2.data;
    unsigned char *out   = disp.data;

    int _niter;

    _niter = 2*niter1;

    for ( int iter = 0; iter < _niter; ++iter ) {
        int partition = iter % 2;
        block_gibbs_iter( false, partition, ndisp, mindisp, p1, p2, width,
                height, ref, match, out);
        if ( partition == 0 ) {
            std::cout << " done non-smooth block gibbs iter " << (iter/2+1) << " of " <<
                niter1 << std::endl;
        }
    }

    _niter = 2*niter2;
    for ( int iter = 0; iter < _niter; ++iter ) {
        int partition = iter % 2;
        block_gibbs_iter( true, partition, ndisp, mindisp, p1, p2, width,
                height, ref, match, out);
        if ( partition == 0 ) {
            std::cout << " done smooth block gibbs iter " << (iter/2+1) << " of " <<
                niter2 << std::endl;
        }
    }

}


int main( int argc, char *argv[] ) {


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

    std::fill(disp.data,disp.data+disp.width*disp.height,0);

    int niter1 = 2;
    int niter2 = 20;
    int ndisp = 128;
    int mindisp = -(ndisp-1);
    int p1 = 10;
    int p2 = 120;
    std::cout << " ndisp: " << ndisp << std::endl;
    std::cout << " mindisp: " << mindisp << std::endl;
    block_gibbs(niter1,niter2,ndisp,mindisp,p1,p2,im1,im2,disp);

    writepgm("disp.pgm",
            disp.data,
            disp.width,
            disp.height);


    delete[] im1.data;
    delete[] im2.data;
    delete[] disp.data;


    return 0;
}
