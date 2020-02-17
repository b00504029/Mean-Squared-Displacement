// msd.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//#define USE_KEYBOARD_INPUT
#define LARGE_FILE  // load the 1.6GB guy

#ifdef LARGE_FILE
#define HAS_MOL
#endif

//#define OUTPUT_MSD

//#include "atom.h"
#include "atom_refactor.h"
//#include "read_bin.h"

constexpr size_t N_TEST = 3;

int main() {
    //for ( size_t i = 0; i < N_TEST; ++i ) {
    //    ATOM atom;
    //    double calc_time = atom.run();
    //    std::cout << "original: " << calc_time << " [s]" << std::endl;
    //}

    double bin_avg = 0.0;
    for (size_t i = 0; i < N_TEST; ++i) {
        AtomR atom_r;
        double total_time = atom_r.run_bin("data/NIMs_nvt_core_FRAMENUM100001.bin");
        bin_avg += total_time;
        std::cout << "bin: " << total_time << " [s]" << std::endl;
    }

    double refactored_avg = 0.0;
    for ( size_t i = 0; i < N_TEST; ++i ) {
        AtomR atom_r;
        double total_time = atom_r.run("data/NIMs_nvt_core_FRAMENUM100001.txt");
        refactored_avg += total_time;
        std::cout << "refactored: " << total_time << " [s]" << std::endl;
    }

    refactored_avg /= N_TEST;
    bin_avg /= N_TEST;

    std::cout << std::endl;

    std::cout << "Refactored init time avg: " << refactored_avg << " [s]" << std::endl;
    std::cout << "Binary init time avg:     " << bin_avg << " [s]" << std::endl;

    std::cout << std::endl;


    return 0;
}

