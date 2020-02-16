#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>  // transform_reduce
#include <string>
#include <unordered_map>
#include <vector>

namespace {

std::string name; // file name

void skip_lines(size_t num_lines, std::ifstream& data) {
    std::string dummy_str;
    for (size_t i = 0; i < num_lines; ++i) {
        std::getline(data, dummy_str);
    }
}

/**
    parse frame number from file `name`
    the `name` format must be *FRAMENUM_____*, where _____ is a integer, and * can be anything
 */
size_t parse_frame_num(const std::string name) {
    constexpr std::string_view keyword = "FRAMENUM";
    if (size_t pos = name.find(keyword); pos != std::string::npos) {
        return std::stoull(name.c_str() + keyword.length() + pos);
    }

    return 0;
}

}

struct AtomR {
    using ParticleInfo = std::unordered_map<size_t, size_t>;  // particle type (ntype) -> particle quantity (npart)
    using Particles = std::vector<double>;                    // [x0, y0, z0, x1, y1, z1, ... ]

    double run();
    int init();
    void preread();
    int init_from_bin();
    void preread_from_bin(std::string name);
    int MSD_calc();

    // parameters for sorting atoms
    int64_t tot_part = 0, tot_frame = 0, tot_type = 0;  // # of total particles, frames, and atom types
    size_t istep = 0;                                  // ith step

    ParticleInfo particle_info;
    std::vector<size_t> itype;                         // Required atom types for MSD calculation
    std::unordered_map<size_t, size_t> itype_to_idx;   // if `itype` = [4, 2, 8], then `itype_to_idx` = [4 -> 0], [2 -> 1], [8 -> 2]
                                                       //     i.e. `type_idx` = itype_to_idx[type];
    // coords[time][type_idx]
    std::vector<std::vector<Particles>> coords;

    // Parameters for MSD
    double dt = 0.0;                                   // Timestep (ps) between frames
    size_t wi = 0, wf = 0, ws = 0;                     // Initial frame, final frame ,stride of frames
    size_t taui = 0, tauf = 0, taus = 0, tau = 0;      // Minimum, final, stride of, current span of frames
    double MSD = 0.0;                                  // Mean-Squared Displacement
};

double AtomR::run() {
    #ifdef USE_KEYBOARD_INPUT
    // File name to read
    std::cout << "Enter the input file name: ";  // dump
    std::cin >> name;

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
    this->itype.emplace_back(1);
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

int AtomR::init() {
    this->preread();

    // Determine the values of tot_part, tot_frame, and tot_type
    std::ifstream intrj;
    intrj.open(name, std::ios::in);

    std::string dummy_str;

    size_t curr_frame = 0;
    {   // jump to the `wi`-th frame
        size_t num_skip_lines = wi * (9 + tot_part);
        skip_lines(num_skip_lines, intrj);
        curr_frame = wi;
    }

    {   // init `itype_to_idx`
        itype_to_idx.reserve(tot_type + 1);
        for (size_t type_idx = 0; type_idx < itype.size(); ++type_idx) {
            size_t type = itype[type_idx];
            itype_to_idx[type] = type_idx;
        }
    }

    {   // pre-allocate of `coords`
        const size_t num_frames_to_read = wf - curr_frame + 1;
        coords.resize(num_frames_to_read);
        for (auto& particles_vec : coords) {
            particles_vec.resize(itype.size());  // required types

            for (auto& type : itype) {  // for particles of each type, allocate 3*num_of_particles (xyz coordinates)
                const size_t type_idx = itype_to_idx[type];
                const size_t num_particle_of_type = 3 * particle_info[type];
                particles_vec[type_idx].resize(num_particle_of_type);
            }
        }
    }

    // read frame data
    for (; !intrj.eof() && (curr_frame <= wf); ++curr_frame) {
        skip_lines(9, intrj);  // first 9 lines are not needed

        if (intrj.eof()) break;

        std::unordered_map<size_t, size_t> write_pos;  // record the position of `coords` to write
        write_pos.reserve(itype.size());

        // read particle data
        for (size_t i = 0; i < tot_part; ++i) {
            size_t type;
            intrj >> dummy_str;  // id
            #ifdef HAS_MOL
            intrj >> dummy_str;  // mol
            #endif
            intrj >> type;

            if (itype_to_idx.count(type) > 0) {
                // write to `coords[frame][type][x, y, z]`
                const size_t type_idx = itype_to_idx[type];

                size_t write_pos_of_curr_type = write_pos[type];
                double x, y, z;
                intrj >> x >> y >> z;
                coords[curr_frame][type_idx][write_pos_of_curr_type] = x;
                coords[curr_frame][type_idx][write_pos_of_curr_type + 1] = y;
                coords[curr_frame][type_idx][write_pos_of_curr_type + 2] = z;
                write_pos[type] += 3;
            }

            std::getline(intrj, dummy_str);  // consume the rest of the line
        }
    }

    intrj.close();

    return 0;
}

/**
    preread the first frame to get the following information
        `tot_part`
        `tot_frame`
        `tot_type`
        `particle_info`
 */
void AtomR::preread() {
    tot_frame = parse_frame_num(name);

    std::ifstream intrj;
    intrj.open(name, std::ios::in);
    std::string str;

    skip_lines(3, intrj);

    // total number of particles
    intrj >> tot_part;
    std::getline(intrj, str);

    skip_lines(5, intrj);

    // number of particles by type
    for (size_t i = 0; i < tot_part; ++i) {
        size_t type;

        intrj >> str;  // skip id
        #ifdef HAS_MOL
        intrj >> str;  // skip mol
        #endif
        intrj >> type;
        ++particle_info[type];

        std::getline(intrj, str);  // skip the rest
    }

    tot_type = particle_info.size();
}

int AtomR::init_from_bin() {
    return 0;
}

/**
    preread the first frame to get the following information
        `tot_part` done
        `tot_frame` done
        `tot_type`
        `particle_info`
 */
void AtomR::preread_from_bin(std::string name) {
    tot_frame = parse_frame_num(name);

    std::ifstream intrj(name, std::ios::in | std::ios::binary);

    int64_t timestep;
    intrj.read(reinterpret_cast<char*>(&timestep), sizeof(int64_t));  // timestep

    if (intrj.eof()) {
        intrj.close();
        return;
    }

    intrj.read(reinterpret_cast<char*>(&this->tot_part), sizeof(int64_t));  // natoms = 375

    int triclinic;
    intrj.read(reinterpret_cast<char*>(&triclinic), sizeof(int));  // triclinic = 0

    int boundary[3][2];
    intrj.read(reinterpret_cast<char*>(&boundary[0][0]), 6*sizeof(int));  // 6*4 = 24 bytes boundary

    std::vector<double> bbx(6, 0.0);  // 6 doubles to read
    intrj.read(reinterpret_cast<char*>(bbx.data()), 6*sizeof(double));  // xlo, xhi, ylo, yhi, zlo, zhi

    if (triclinic != 0) {
        std::vector<double> tri(3, 0.0);  // 3 doubles to read
        intrj.read(reinterpret_cast<char*>(tri.data()), 3*sizeof(double));  // xy, xz, yz
    }

    int size_of_one_line;
    intrj.read(reinterpret_cast<char*>(&size_of_one_line), sizeof(int));

    int nchunk;
    intrj.read(reinterpret_cast<char*>(&nchunk), sizeof(int));

    std::vector<double> buf;
    constexpr size_t type_offset = 1;  // if HAS_MOL, 2
    for (size_t i = 0; i < nchunk; ++i) {
        int nelems;
        intrj.read(reinterpret_cast<char*>(&nelems), sizeof(int));

        buf.resize(nelems, 0.0);
        intrj.read(reinterpret_cast<char*>(buf.data()), nelems*sizeof(double));

        for (size_t idx = 0; idx < buf.size(); idx += size_of_one_line) {
            // assume `size_of_one_line` = 5, i.e. id, type, x, y, z
            size_t type = static_cast<size_t>(buf[idx + type_offset]);
            ++this->particle_info[type];
        }
    }

    this->tot_part = this->particle_info.size();

    intrj.close();
}

int AtomR::MSD_calc() {
    #ifdef OUTPUT_MSD
    std::ofstream outf;
    #ifndef LARGE_FILE
    outf.open("msd_refactored.txt", std::ios::out);
    #else
    outf.open("msd_large_refactored.txt", std::ios::out);
    #endif // !LARGE_FILE
    outf << "t(ps)" << " " << "MSD(angstroms)" << std::endl;
    #endif

    size_t n = 0;  // # of total used particles
    for (size_t i = 0; i < itype.size(); ++i) {
        n += particle_info[itype[i]];
    }

    auto square_error = [] (const double a, const double b) -> const double {
        const double e = a - b;
        return e * e;
    };

    for (tau = taui; tau <= tauf; tau += taus) {
        MSD = 0;

        //  t0 <= (wf - wi - tau - 1)   has been rearranged into   (t0 + wi + tau + 1) <= wf
        // because `t0`, `wi`, `wf`, `tau` are size_t, and the minus operation might lead to overflow (< 0)
        //
        // of course, plus operation could also overflow,
        // but that happens only when all of the numbers are huge
        for (size_t t0 = 0; (t0 + wi + tau + 1) <= wf; t0 += ws) {  // Reference frame
            for (size_t i = 0; i < itype.size(); ++i) {
                const size_t type_idx = itype_to_idx[itype[i]];

                // transform reduce is the c++ implementation of Google MapReduce
                // which is designed for large scale parallelization
                const auto& vec_a = coords[t0 + tau][type_idx];
                const auto& vec_b = coords[t0][type_idx];
                MSD += std::transform_reduce(vec_a.cbegin(), vec_a.cend(), vec_b.cbegin(), 0.0, std::plus<>(), square_error);
            }
        }

        MSD *= ((ws / ((double)wf - wi - tau)) / n);
        #ifdef OUTPUT_MSD
        outf << tau * dt << " " << MSD << std::endl;
        #endif
    }

    #ifdef OUTPUT_MSD
    outf.close();
    #endif

    return 0;
}
