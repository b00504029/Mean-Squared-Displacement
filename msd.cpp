// msd.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//#define USE_KEYBOARD_INPUT
#define LARGE_FILE  // load the 1.6GB guy

#ifdef LARGE_FILE
#define HAS_MOL
#endif

//#define OUTPUT_MSD

#include "atom.h"
#include "atom_refactor.h"

constexpr size_t N_TEST = 3;

int main()
{
    for ( size_t i = 0; i < N_TEST; ++i ) {
        ATOM atom;
        double calc_time = atom.run();
        std::cout << "original: " << calc_time << " [s]" << std::endl;
    }

    for ( size_t i = 0; i < N_TEST; ++i ) {
        AtomR atom_r;
        double total_time = atom_r.run();
        std::cout << "refactored: " << total_time << " [s]" << std::endl;
    }


    return 0;
}

