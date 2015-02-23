#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <limits>

typedef  long long int energy_t; /* This datatype should be signed */


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

void generate_x_grid( int width, int height,
       std::vector< std::vector< int > > &rays,
       std::vector< std::vector< int > > &idx2ray,
       bool flip = false ) {

    for ( int irow = 0; irow < height; ++irow ) {
        size_t rayidx = rays.size();
        std::vector<int> ray;
        for ( int icol = 0; icol < width; ++icol ) {
            int y = irow;
            int x = flip ? width - icol - 1 : icol;
            int idx = y * width + x;
            idx2ray[idx].push_back(rayidx);
            idx2ray[idx].push_back(ray.size()/2); /* the current idx into ray */
            ray.push_back(y);
            ray.push_back(x);
        }
        rays.push_back(ray); //FIXME: Alot of copying
    }

}

void generate_y_grid( int width, int height,
       std::vector< std::vector< int > > &rays,
       std::vector< std::vector< int > > &idx2ray,
       bool flip = false ) {

    for ( int icol = 0; icol < width; ++icol ) {
        size_t rayidx = rays.size();
        std::vector<int> ray;
        for ( int irow = 0; irow < height; ++irow ) {
            int x = icol;
            int y = flip ? height - irow - 1 : irow;
            int idx = y * width + x;
            idx2ray[idx].push_back(rayidx);
            idx2ray[idx].push_back(ray.size()/2); /* the current idx into ray */
            ray.push_back(y);
            ray.push_back(x);
        }
        rays.push_back(ray); //FIXME: Alot of copying
    }

}

void generate_yx_grid( int width, int height,
       std::vector< std::vector< int > > &rays,
       std::vector< std::vector< int > > &idx2ray,
       int dy,
       int dx,
       bool flipy = false,
       bool flipx = false ) {

    if ( dx <= 0 || dy <= 0 ) {
        throw std::invalid_argument("dx and dy must be positive and greater than zero.");
    }

    for ( int icol = 0; icol < width; ++icol ) {
        size_t rayidx = rays.size();
        std::vector<int> ray;
        int posx, posy;
        for ( posx = icol, posy = 0;
                posx < width && posy < height; 
                posx += dx, posy += dy ) {
            int x = flipx ? width - posx - 1 : posx;
            int y = flipy ? height - posy  - 1 : posy;
            int idx = y * width + x;
            idx2ray[idx].push_back(rayidx);
            idx2ray[idx].push_back(ray.size()/2); /* the current idx into ray */
            ray.push_back(y);
            ray.push_back(x);
        }
        rays.push_back(ray); //FIXME: Alot of copying
    }

    for ( int irow = 1; irow < height; ++irow ) {
        size_t rayidx = rays.size();
        std::vector<int> ray;
        int posx, posy;
        for ( posx = 0, posy = irow;
                posx < width && posy < height; 
                posx += dx, posy += dy ) {
            int x = flipx ? width - posx - 1 : posx;
            int y = flipy ? height - posy  - 1 : posy;
            int idx = y * width + x;
            idx2ray[idx].push_back(rayidx);
            idx2ray[idx].push_back(ray.size()/2); /* the current idx into ray */
            ray.push_back(y);
            ray.push_back(x);
        }
        rays.push_back(ray); //FIXME: Alot of copying
    }

}



void sgm_scan(int ndisp,
        int mindisp,
        int p1,
        int p2,
        int width, /* Image width */
        int *unravel,
        int npix,
        unsigned char *ref,
        unsigned char *match,
        energy_t *energymatrix /* Preallocated */ ) {

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
        int refy = unravel[0];
        int refx = unravel[1];
        int matchposx = refx+displut[idisp];
        int matchposy = refy;
        if ( matchposx >= 0 && matchposx < width ) { /* epipolar line scans only */
            energymatrix[idisp] = std::abs((int)(ref[refy*width+refx])-
                    (int)(match[matchposy*width+matchposx]));
        } else {
            energymatrix[idisp] = -1; /* Negative values will indicate invalid matching */
        }
    }
    
    for ( int ipix = 1; ipix < npix; ++ipix ) {
        int rejectioncount = 0;
        for ( int idisp = 0; idisp < ndisp; ++idisp ) {
            int refy = unravel[2*ipix+0];
            int refx = unravel[2*ipix+1];
            int matchposx = refx+displut[idisp];
            int matchposy = refy;
            int dataterm = INT_MAX;
            if ( matchposx >= 0 && matchposx < width ) { /* epipolar line scans only */
                dataterm = std::abs((int)(ref[refy*width+refx])-
                        (int)(match[matchposy*width+matchposx]));
            } else {
                energymatrix[ipix*ndisp + idisp] = -1; /* Negative values will indicate invalid matching */
                rejectioncount++;
                continue;
            }
            energy_t mintotalterm = std::numeric_limits<energy_t>::max();
            for ( int idispprev = 0; idispprev < ndisp; ++idispprev ) {
                int smoothterm = -1;
                energy_t v = energymatrix[ (ipix-1)*ndisp + idispprev ];
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

                energy_t totalterm = smoothterm + dataterm;

                if ( totalterm  < mintotalterm ) {
                    mintotalterm = totalterm;
                }

            }
            energymatrix[ipix*ndisp + idisp ] = mintotalterm;
        }
        assert(rejectioncount < ndisp);
    }

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

    unsigned char *ref   = im1.data;
    unsigned char *match = im2.data;
    unsigned char *out   = disp.data;

    //Create a 16 ray supported SGM algorithm to test if it works
    //Currently commented out in favor of 4 ray version.
    std::vector< std::vector<int> > rays;
    std::vector< std::vector<int> > idx2ray(width*height);
    generate_x_grid(width,height,rays,idx2ray);
    generate_x_grid(width,height,rays,idx2ray,true);
    generate_y_grid(width,height,rays,idx2ray);
    generate_y_grid(width,height,rays,idx2ray,true);
    generate_yx_grid(width,height,rays,idx2ray,1,1,false,false);
    generate_yx_grid(width,height,rays,idx2ray,1,1,true,false);
    generate_yx_grid(width,height,rays,idx2ray,1,1,true,true);
    generate_yx_grid(width,height,rays,idx2ray,1,1,false,true);
    //generate_yx_grid(width,height,rays,idx2ray,2,1,false,false);
    //generate_yx_grid(width,height,rays,idx2ray,2,1,true,false);
    //generate_yx_grid(width,height,rays,idx2ray,2,1,true,true);
    //generate_yx_grid(width,height,rays,idx2ray,2,1,false,true);
    //generate_yx_grid(width,height,rays,idx2ray,1,2,false,false);
    //generate_yx_grid(width,height,rays,idx2ray,1,2,true,false);
    //generate_yx_grid(width,height,rays,idx2ray,1,2,true,true);
    //generate_yx_grid(width,height,rays,idx2ray,1,2,false,true);

    std::vector< std::vector<energy_t> > energies(rays.size());

#pragma omp parallel for
    for ( size_t i = 0; i < rays.size(); ++i ) {
        std::vector<int> &ray = rays[i];
        int *unravel = &rays[i][0];
        int npix = rays[i].size()/2;
        energies[i].resize(npix*ndisp); /* The ith energy, corresponding to this ray. */
        sgm_scan(ndisp,mindisp,p1,p2,width,unravel,npix,ref,match,&energies[i][0]);
        //if ( (i+1) % 15 == 0 ) {  /* Don't print too much */
        //    std::cout << "\r done scan " << i+1 << " of " << rays.size() << std::flush;
        //}

    }
    std::cout << std::endl;

    std::cout << " decoding image: " << std::endl;

    //Decode the image in 2D, this requires taking the min of the energy for 
    //each incident ray over all disparities
    for ( int irow = 0; irow < height; ++irow ) {
        for ( int icol = 0; icol < width; ++icol ) {

            size_t idx = irow * width + icol;
            size_t nrays = idx2ray[idx].size()/2;

            energy_t minenergy = std::numeric_limits<energy_t>::max();
            int mindisp  = -1;
            for ( int idisp = 0; idisp < ndisp; ++idisp ) {
                int energy = 0;
                bool dispinvalid = false;
                for ( size_t iray = 0; iray < nrays; ++iray ) {
                    int rayidx   = idx2ray[idx][2*iray+0]; /* Tells us which energy */
                    int idxinray = idx2ray[idx][2*iray+1]; /* Tells us where to look within energy */
                    {
                    int lookupidx =  ndisp*idxinray + idisp;
                    assert( lookupidx >= 0 && lookupidx < ndisp*rays[rayidx].size()/2 );
                    }
                    energy_t rayenergy = energies[rayidx][ndisp*idxinray + idisp];
                    if ( rayenergy < 0 ) {
                        dispinvalid = true;
                        break;
                    }
                    energy += rayenergy; /* Add the minimum energy of this ray */
                }
                if ( dispinvalid ) {
                    continue; /* If disparity is invalid do not allow. */
                }
                if ( energy < minenergy ) {
                    minenergy = energy;
                    mindisp = idisp;
                }
            }
            assert(mindisp >= 0);

            out[idx] = mindisp;

        }
        std::cout << "\r done row " << irow+1 << " of " << height << std::flush;
    }
    std::cout << std::endl;

    



}


int main( int argc, char *argv[] ) {


    im_t im1, im2, disp;

    std::cout << " Reading im1 " << std::endl;
    //readpgm("cones_left.pgm",
    readpgm("im0.pgm",
            im1.data,
            im1.width,
            im1.height);

    std::cout << " Reading im2 " << std::endl;
    //readpgm("cones_right.pgm",
    readpgm("im1.pgm",
            im2.data,
            im2.width,
            im2.height);

    disp.width  = im1.width;
    disp.height = im1.height;
    disp.data   = new unsigned char[disp.width*disp.height];

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


    delete[] im1.data;
    delete[] im2.data;
    delete[] disp.data;


    return 0;
}
