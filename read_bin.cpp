#include "read_bin.h"


#include <stdio.h>
#include <string.h>
#include "stdint.h"
#define __STDC_FORMAT_MACROS
#include "inttypes.h"

#ifndef PRId64
#define PRId64 "ld"
#endif

#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG)
#define LAMMPS_SMALLBIG
#endif

#if defined(LAMMPS_SMALLBIG)
typedef int tagint;
typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#elif defined(LAMMPS_SMALLSMALL)
typedef int tagint;
typedef int bigint;
#define BIGINT_FORMAT "%d"
#else /* LAMMPS_BIGBIG */
typedef int64_t tagint;
typedef int64_t bigint;
#define BIGINT_FORMAT "%" PRId64
#endif

void Bin::to_txt(std::string name)
{
    // open file to read
    FILE* fp;
    fopen_s(&fp, name.c_str(), "rb");
    if (!fp)
    {
        // ERROR: Could not open file
        return;
    }

    // open txt to write to
    FILE* fptxt;
    {
        std::string filetxt = name + ".txt";
        fopen_s(&fptxt, filetxt.c_str(), "w");
    }


    int i, j, k, m, n;
    bigint ntimestep, natoms;
    int size_one, nchunk, triclinic;
    double xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz;
    int boundary[3][2];
    char boundstr[9];

    int maxbuf = 0;
    double* buf = NULL;

    // loop until eof
    while (true)
    {
        fread(&ntimestep, sizeof(bigint), 1, fp);  // 8 byte, time-step

        // check eof
        if (feof(fp))
        {
            fclose(fp);
            fclose(fptxt);
            break;
        }

        fread(&natoms, sizeof(bigint), 1, fp);  // 8 bytes, N atoms
        fread(&triclinic, sizeof(int), 1, fp);  // 4 bytes, triclinic
        fread(&boundary[0][0], 6 * sizeof(int), 1, fp);  // 6*4 = 24 bytes boundary
        fread(&xlo, sizeof(double), 1, fp);  // 8 bytes each, boundaries
        fread(&xhi, sizeof(double), 1, fp);
        fread(&ylo, sizeof(double), 1, fp);
        fread(&yhi, sizeof(double), 1, fp);
        fread(&zlo, sizeof(double), 1, fp);
        fread(&zhi, sizeof(double), 1, fp);
        if (triclinic)
        {
            fread(&xy, sizeof(double), 1, fp);  // 8 bytes each, more detailed boundary information
            fread(&xz, sizeof(double), 1, fp);
            fread(&yz, sizeof(double), 1, fp);
        }
        fread(&size_one, sizeof(int), 1, fp);  // 4 bytes, size_one
        fread(&nchunk, sizeof(int), 1, fp);  // 4 bytes, N chunk


        // write
        fprintf(fptxt, "ITEM: TIMESTEP\n");
        fprintf(fptxt, BIGINT_FORMAT "\n", ntimestep);  // Line 2
        fprintf(fptxt, "ITEM: NUMBER OF ATOMS\n");
        fprintf(fptxt, BIGINT_FORMAT "\n", natoms);  // Line 4

        m = 0;
        for (int idim = 0; idim < 3; ++idim)
        {
            for (int iside = 0; iside < 2; ++iside)
            {
                if (boundary[idim][iside] == 0) boundstr[m++] = 'p';
                else if (boundary[idim][iside] == 1) boundstr[m++] = 'f';
                else if (boundary[idim][iside] == 2) boundstr[m++] = 's';
                else if (boundary[idim][iside] == 3) boundstr[m++] = 'm';
            }
            boundstr[m++] = ' ';
        }
        boundstr[8] = '\0';

        // Line 5:
        //     ITEM: BOX BOUNDS pp pp pp
        if (!triclinic)
        {
            fprintf(fptxt, "ITEM: BOX BOUNDS %s\n", boundstr);
            fprintf(fptxt, "%g %g\n", xlo, xhi);
            fprintf(fptxt, "%g %g\n", ylo, yhi);
            fprintf(fptxt, "%g %g\n", zlo, zhi);
        }
        else
        {
            fprintf(fptxt, "ITEM: BOX BOUNDS %s xy xz yz\n", boundstr);
            fprintf(fptxt, "%g %g %g\n", xlo, xhi, xy);
            fprintf(fptxt, "%g %g %g\n", ylo, yhi, xz);
            fprintf(fptxt, "%g %g %g\n", zlo, zhi, yz);
        }


        // start to deal with atom data
        fprintf(fptxt, "ITEM: ATOMS\n");

        for (i = 0; i < nchunk; ++i)
        {
            fread(&n, sizeof(int), 1, fp);  // 4 bytes, n (number of doubles in this chunk)

            // extend buffer to fit chunk size
            if (n > maxbuf)
            {
                if (buf) delete[] buf;
                buf = new double[n];  // `buf` -> double[] of size `n`
                maxbuf = n;
            }

            // read chunk and write as size_one values per line
            fread(buf, sizeof(double), n, fp);  // 8*n bytes of doubles (coordinates)

            // `size_one` probably means "size of one line", for water.bin, it should be 5?
            // so now `n` means `n` lines of data
            n /= size_one;
            m = 0;
            for (j = 0; j < n; ++j)  // for each line
            {
                for (k = 0; k < size_one; ++k)  // for each data
                    fprintf(fptxt, "%g ", buf[m++]);  // 

                fprintf(fptxt, "\n");
            }
        }

        // show which timestep has been written
        printf(" " BIGINT_FORMAT, ntimestep);
        fflush(stdout);
    }
    printf("\n");

    // clear buf
    if (buf) delete[] buf;
}

ParticleDataContainter Bin::load_chunk(std::string_view name)
{
    ParticleDataContainter res;




    return res;
}
