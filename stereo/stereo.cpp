#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <climits>


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

#define CORR_ALGO_NCC 

void patchnxcorr( int rad, im_t &im1, im_t &im2, im_t &disp ) {

    if ( !(im1.width == im2.width && im1.height == im2.height  ) ) {
        std::cerr << " Images are not same dimensions." << std::endl;
        exit(-1);
    }

    int width = im1.width;
    int height = im1.height;

    int disparityrange = width;

    int sz = 2*rad+1;
    unsigned char *cache = new unsigned char[sz*sz];

    // Loop over all scan lines
    for ( int irow = 0; irow < height; ++irow ) { 

        for ( int icol = 0; icol < width; ++icol ) { 

            for ( int rrow = -rad; rrow <= rad; ++rrow ) {
                for ( int rcol = -rad; rcol <= rad; ++rcol ) {
                    int r = irow + rrow; 
                    int c = icol + rcol;
                    if ( r < 0 || r >= height || c < 0 || c >= width ) {
                        int pr = rrow + rad;
                        int pc = rcol + rad;
                        cache[pr*sz+pc] = 0;
                        continue;
                    }
                    int pr = rrow + rad;
                    int pc = rcol + rad;
                    cache[pr*sz+pc] = im2.data[r*width+c];
                }
            }

#ifdef CORR_ALGO_NCC
            double maxncc = 0;
#else 
            int minssd = 0;
#endif

            int minidx = -1;

            for ( int scol = std::max(0,icol - disparityrange);
                    scol < std::min(width,icol+disparityrange); ++scol ) { 

#ifdef CORR_ALGO_NCC
                double ncc = 0;
                double mean1 = 0;
                double mean2 = 0;
                double std1 = 0;
                double std2 = 0;
                double xcorr = 0;
#else 
                int ssd = 0;
#endif

#ifdef CORR_ALGO_NCC
                for ( int rrow = -rad; rrow <= rad; ++rrow ) {
                    for ( int rcol = -rad; rcol <= rad; ++rcol ) {
                        int r = irow + rrow; 
                        int c = scol + rcol;
                        int pr = rrow + rad;
                        int pc = rcol + rad;
                        mean2 += cache[pr*sz+pc];
                        if ( r < 0 || r >= height || c < 0 || c >= width ) {
                            continue;
                        }
                        mean1 += im1.data[r*width+c];
                    }
                }
                mean1 /= sz*sz;
                mean2 /= sz*sz;
                for ( int rrow = -rad; rrow <= rad; ++rrow ) {
                    for ( int rcol = -rad; rcol <= rad; ++rcol ) {
                        int r = irow + rrow; 
                        int c = scol + rcol;
                        int pr = rrow + rad;
                        int pc = rcol + rad;
                        double v2 = cache[pr*sz+pc] - mean2;
                        std2 += v2*v2;
                        double e;
                        if ( r < 0 || r >= height || c < 0 || c >= width ) {
                            e = 0;
                        } else {
                            e = im1.data[r*width+c];
                        }
                        double v1 =  e - mean1;
                        std1 += v1*v1;
                        xcorr += v1 * v2;
                    }
                }
                std1 = std::sqrt(std1);
                std2 = std::sqrt(std2);
                ncc = xcorr / std1 / std2;

                //std::cout << " ncc: " << ncc << std::endl;

                if ( minidx < 0 || ncc > maxncc ) {
                    maxncc = ncc;
                    minidx = scol;
                }

#else
                for ( int rrow = -rad; rrow <= rad; ++rrow ) {
                    for ( int rcol = -rad; rcol <= rad; ++rcol ) {
                        int r = irow + rrow; 
                        int c = scol + rcol;
                        int pr = rrow + rad;
                        int pc = rcol + rad;
                        int delta;
                        if ( r < 0 || r >= height || c < 0 || c >= width ) {
                            delta =  (int)cache[pr*sz+pc];
                            ssd   += delta * delta;
                            continue;
                        }
                        delta = (int)cache[pr*sz+pc] - (int)im1.data[r*width+c];
                        ssd += delta*delta;
                    }
                }

                if ( minidx < 0 || ssd < minssd ) {
                    minssd = ssd;
                    minidx = scol;
                }
#endif

            }

            //int disparity = icol - minidx;
            //std::cout << " disparity: " << disparity << std::endl;

            disp.data[irow*width + icol] = (icol-minidx) / 16 + 127.;


        }


    }

}

void sgmscan(int ndisp,
        int mindisp,
        int p1,
        int p2,
        int npix,
        unsigned char *ref,
        unsigned char *match,
        unsigned char *out ) {

    int* matrix = new int[ndisp*npix];
    std::fill(matrix,matrix+ndisp*npix,0); /* May be not required. */

    int* pathmatrix = new int[ndisp*npix];
    std::fill(pathmatrix,pathmatrix+ndisp*npix,0); /* May be not required. */

    int* displut = new int[ndisp];
    int* dispdifflut = new int[ndisp*ndisp];
    for ( int idisp = 0; idisp < ndisp; ++idisp ) {
        displut[idisp] = idisp + mindisp; 
    }
    for ( int idisp1 = 0; idisp1 < ndisp; ++idisp1 ) {
        for ( int idisp2 = 0; idisp2 < ndisp; ++idisp2 ) {
            dispdifflut[idisp1*ndisp+idisp2] =
                std::abs( displut[idisp1] - displut[idisp2] );
        }
    }

    // Initialize the 0th column first dependent on the data term only
    for ( int idisp = 0; idisp < ndisp; ++idisp ) {
        int matchpos = 0+displut[idisp];
        if ( matchpos >= 0 && matchpos < npix ) {
            matrix[idisp] = std::abs((int)(ref[0])-(int)(match[matchpos]));
        } else {
            matrix[idisp] = -1; /* Negative values will indicate invalid matching */
        }
        pathmatrix[idisp] = -1;
    }
    
    for ( int ipix = 1; ipix < npix; ++ipix ) {
        for ( int idisp = 0; idisp < ndisp; ++idisp ) {
                int matchpos = ipix+displut[idisp];
                int dataterm = INT_MAX;
                if ( matchpos >= 0 && matchpos < npix ) {
                    dataterm = std::abs((int)(ref[ipix])-(int)(match[matchpos]));
                } else {
                    matrix[ipix*ndisp + idisp] = -1; /* Negative values will indicate invalid matching */
                    pathmatrix[ipix*ndisp + idisp] = -1;
                    continue;
                }
                int mintotalterm = INT_MAX;
                int minidx = -1;
                for ( int idispprev = 0; idispprev < ndisp; ++idispprev ) {
                    int smoothterm = -1;
                    int v = matrix[ (ipix-1)*ndisp + idispprev ];
                    if ( v < 0 ) {
                        continue;
                    }

			
                    // mg: According to profiling, the line below is incredibly slow,
                    // thus recommend a lut to replace for speed.
                    //int absdisp = std::abs(displut[idispprev]-displut[idisp]);
                    int absdisp = dispdifflut[idisp*ndisp + idispprev];

                    if ( absdisp == 0 ) {
                        smoothterm = v;
                    } else if (  absdisp == 1 ) {
                        smoothterm = v + p1;
                    } else { 
                        smoothterm = v + p2;
                    }

                    int totalterm = smoothterm + dataterm;

                    if ( totalterm  < mintotalterm ) {
                        mintotalterm = totalterm;
                        minidx = idispprev;
                    }

                }
                matrix[ipix*ndisp + idisp ] = mintotalterm;
                pathmatrix[ipix*ndisp + idisp ] = minidx;
        }
    }

    // Decode from optimal
    int minenergy = INT_MAX;
    int minenergyidx = -1;
    for ( int idisp = 0; idisp < ndisp; ++idisp ) {
        int energy = matrix[(npix-1)*ndisp + idisp ];
        if ( energy < 0 ) {
            continue;
        }
        if ( energy < minenergy ) {
            minenergy = energy;
            minenergyidx = idisp;
        }
    }

    out[npix-1] = minenergyidx;

    for ( int ipix = npix-2; ipix >= 0; --ipix ) {
        out[ipix] = pathmatrix[ (ipix+1)*ndisp + out[ipix+1] ];
    }


    delete[] matrix;
    delete[] pathmatrix;
    delete[] displut;
    delete[] dispdifflut;


}


void sgm( int ndisp,
        int mindisp,
        int p1,
        int p2,
        im_t &im1,
        im_t &im2,
        im_t &disp ) {

    if ( !(im1.width == im2.width && im1.height == im2.height  ) ) {
        std::cerr << " Images are not same dimensions." << std::endl;
        exit(-1);
    }

    int width = im1.width;
    int height = im1.height;

    for ( int irow = 0; irow < height; ++irow ) { 
        int npix = width;
        unsigned char *ref   = im1.data + irow*width;
        unsigned char *match = im2.data + irow*width;
        unsigned char *out   = disp.data + irow*width;
        sgmscan(ndisp,mindisp,p1,p2,npix,ref,match,out);
        std::cout << "\r done row " << irow+1 << " of " << height << std::flush;
    }
    std::cout << std::endl;


}


int main( int argc, char *argv[] ) {


    //int rad = 7;
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

    //patchnxcorr(rad,im1,im2,disp);
    int ndisp = 128;
    int mindisp = -(ndisp-1);
    int p1 = 10;
    int p2 = 150;
    std::cout << " ndisp: " << ndisp << std::endl;
    std::cout << " mindisp: " << mindisp << std::endl;
    sgm(ndisp,mindisp,p1,p2,im1,im2,disp);

    writepgm("disp.pgm",
            disp.data,
            disp.width,
            disp.height);


    delete[] im1.data, im2.data, disp.data;


    return 0;
}
