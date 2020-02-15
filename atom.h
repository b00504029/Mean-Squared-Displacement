#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct ATOM {
    double run();
    int init();
    int MSD_calc();

    // Skipping frames unwanted until wi
    void skip_lines(int, std::ifstream*);

    std::string str, name;                          // Store the instruction, the file name

    // parameters for sorting atoms
    int i = 0, j = 0, k = 0;                        // Counters
    int tot_part = 0, tot_frame = 0, tot_type = 0;  // # of total particles, frames, and atom types
    int istep = 0;                                  // ith step

    std::vector<int> itype;                         // Required atom types for MSD calculation
    std::vector<int> ntype;                         // # of atom types: {1,2,3,...,n}
    std::vector<int> npart;                         // # of total particles with different atom types: {125,250,....}
    
    // coords[time][type][part][xyz]
    std::vector<std::vector<std::vector<std::vector<double>>>> coords;  // Coordinates

    // Parameters for MSD
    double dt = 0.0;                                // Timestep (ps) between frames
    int wi = 0, wf = 0, ws = 0;                     // Initial frame, final frame ,stride of frames
    int taui = 0, tauf = 0, taus = 0, tau = 0;      // Minimum, final, stride of, current span of frames
    double MSD = 0.0;                               // Mean-Squared Displacement
};



double ATOM::run() {
#ifdef USE_KEYBOARD_INPUT
    // File name to read
    std::cout << "Enter the input file name: ";  // dump
    std::cin >> this->name;

    // Atom type
    int a;
    std::cout << "Enter the atom type: ";  // 1
    do {
        std::cin >> a;
        atom.itype.emplace_back(a);
    } while (std::cin.get() != '\n');

    // Timestep in ps
    std::cout << "Enter the timestep in ps: ";  // 1
    std::cin >> atom.dt;

    // Lag times t from taui to tauf by taus
    std::cout << "Enter the initial lag time, final lag time,and a stride of lag times: ";  // 0 10 1
    std::cin >> atom.taui >> atom.tauf >> atom.taus;

    // For analysis
    std::cout << "Enter the initial frame, final frame,and a span of frames: ";  // 0 10 1
    std::cin >> atom.wi >> atom.wf >> atom.ws;
#else
#ifndef LARGE_FILE
    name = "dump_FRAMENUM11";
    this->itype.emplace_back(2);
    this->dt = 1.0;
    this->taui = 0;
    this->tauf = 10;
    this->taus = 1;
    this->wi = 0;
    this->wf = 10;
    this->ws = 1;
#else
    name = "NOHMs_nvt_long_FRAMENUM100001.lammpstrj";
    this->itype.emplace_back(2);
    this->dt = 1.0;
    this->taui = 0;
    this->tauf = 100000;
    this->taus = 100;
    this->wi = 0;
    this->wf = 90000;
    this->ws = 100;
#endif // !LARGE_FILE
#endif

    // Initialization
    auto init_start_time = std::chrono::system_clock::now();
    this->init();
    std::chrono::duration<double> init_time = std::chrono::system_clock::now() - init_start_time;
    std::cout << "    init: " << init_time.count() << " [s]" << std::endl;

    // Mean-Squared Displacement (MSD) calculation
    auto calc_start_time = std::chrono::system_clock::now();
    this->MSD_calc();
    std::chrono::duration<double> calc_time = std::chrono::system_clock::now() - calc_start_time;
    std::cout << "    calc: " << calc_time.count() << " [s]" << std::endl;

    return (init_time + calc_time).count();
}


int ATOM::init() {
    // Determine the values of tot_part, tot_frame, and tot_type
    std::ifstream intrj;
    intrj.open(name, std::ios::in);

    int temptype = 0, maxtype = 0;
    double coord;  // Temporal coordinates
    std::vector<std::vector<std::vector<double>>> frame_data;

    // Skip the lines unwanted
    for (tot_frame = 0; tot_frame < wi; ++tot_frame) {
        if (tot_frame == 0) {
            getline(intrj, str);
            getline(intrj, str);
            getline(intrj, str);

            intrj >> tot_part;
            getline(intrj, str);

            skip_lines(5 + tot_part, &intrj);
        } else {
            skip_lines(9 + tot_part, &intrj);
        }
    }

    // Read the desired lines
    for ( ;!intrj.eof() && ( tot_frame <= wf ) ; ++tot_frame) {
        getline(intrj, str);
        getline(intrj, str);
        getline(intrj, str);

        intrj >> tot_part;
        getline(intrj, str);

        getline(intrj, str);
        getline(intrj, str);
        getline(intrj, str);
        getline(intrj, str);
        getline(intrj, str);

        if (intrj.eof()) break;

        for (j = 0; j < maxtype; ++j) {
            npart[j] = 0;
        }

        for (i = 0; i < tot_part; ++i) {
            intrj >> str;
#ifdef HAS_MOL
            intrj >> str;  // skip mol
#endif
            intrj >> temptype;

            int temptype_idx = temptype - 1;
            if (tot_frame == wi) {
                if (temptype > maxtype) {
                    for (j = maxtype; j < temptype; ++j) {
                        ntype.emplace_back(j + 1);
                        npart.emplace_back(0);
                    }

                    frame_data.resize(temptype);
                    maxtype = temptype;
                }

                ++npart[temptype_idx];
                frame_data[temptype_idx].resize(npart[temptype_idx]);
                for (j = 0; j < 3; j++) {
                    intrj >> coord;
                    frame_data[temptype_idx][npart[temptype_idx] - (int64_t)1].emplace_back(coord);
                }
            } else {
                for (j = 0; j < itype.size(); ++j) {
                    if (temptype == itype[j]) {
                        ++npart[temptype_idx];
                        for (k = 0; k < 3; k++) {
                            intrj >> frame_data[temptype_idx][npart[temptype_idx] - (int64_t)1][k];
                        }
                    }
                }
            }

            getline(intrj, str);
        }

        coords.emplace_back(frame_data);
    }
    tot_type = maxtype;
    intrj.close();

    return 0;
}

int ATOM::MSD_calc() {
#ifdef OUTPUT_MSD
    std::ofstream outf;
#ifndef LARGE_FILE
    outf.open("msd.txt", std::ios::out);
#else
    outf.open("msd_large.txt", std::ios::out);
#endif // !LARGE_FILE
    outf << "t(ps)" << " " << "MSD(angstroms)" << std::endl;
#endif

    int t0, n = 0;                               // Referecnce frame, # of total used particles
    double temp;                                 // Temporal average
    for (tau = taui; tau <= tauf; tau += taus) {
        MSD = 0;
        for (t0 = 0; t0 <= (wf - wi - tau - 1); t0 += ws) {
            temp = n = 0;
            for (i = 0; i < itype.size(); ++i) {
                int curr_type = itype[i] - 1;
                for (j = 0; j < npart[curr_type]; ++j) {
                    for (k = 0; k < 3; ++k) {
                        temp += std::pow(coords[(int64_t)t0 + tau][curr_type][j][k] - coords[t0][curr_type][j][k], 2);
                    }
                }

                n += npart[curr_type];
            }

            MSD += temp;
        }

        MSD *= ((ws/((double)wf - wi - tau))/n);
#ifdef OUTPUT_MSD
        outf << tau*dt << " " << MSD << std::endl;
#endif
    }

#ifdef OUTPUT_MSD
    outf.close();
#endif

    return 0;
}

void ATOM::skip_lines(int num_lines, std::ifstream* data) {
    for (i = 0; i < num_lines; ++i) {
        getline( *data, str );
    }
}

