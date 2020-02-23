#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>  // transform_reduce
#include <string>
#include <unordered_map>
#include <vector>


/**
    binary file layout

    1st section
        -------------------------------- --------------------------------
        |           ntimestep 8          |            natoms 8            |
        -------------------------------- --------------------------------

    2nd section: triclinic and boundary
        -----------------
        |   triclinic 4   |
        --------------------------------------------------------------------------------------------------
        |                                        boundary 4 * 6                                            |
        --------------------------------------------------------------------------------------------------
        |                                        boundary data                                             |
        |                                            8 * 6                                                 |
        --------------------------------------------------------------------------------------------------
        |                         triclinic data 8 * 3 (exist only if triclinic == 1)                      |
        --------------------------------------------------------------------------------------------------

    3rd section
        ---------------------------------
        |  line size 4   |    nchunk 4    |
        ---------------- ----------------

        repeat the following block for `nchunk` times
        ----------------
        |     nelem 4    |
        --------------------------------------------------------------------------------------------------
        |                                 elements = data 8 * nelems ...                                  |
        --------------------------------------------------------------------------------------------------
 */
namespace BinReader {
    inline void ignore_int(std::istream& bin_data) {
        bin_data.ignore(sizeof(int));
    }

    inline void ignore_n_int(std::istream& bin_data, size_t n) {
        bin_data.ignore(n*sizeof(int));
    }

    inline void read_int(std::istream& bin_data, int& i) {
        bin_data.read(reinterpret_cast<char*>(&i), sizeof(int));
    }

    inline void ignore_int64(std::istream& bin_data) {
        bin_data.ignore(sizeof(int64_t));
    }

    inline void ignore_n_double(std::istream& bin_data, size_t n) {
        bin_data.ignore(n * sizeof(double));
    }

    inline void read_n_double(std::istream& bin_data, size_t n, void* buf) {
        bin_data.read(reinterpret_cast<char*>(buf), n * sizeof(double));
    }

    inline void ignore_timestep(std::istream& bin_data) {
        BinReader::ignore_int64(bin_data);
    }

    inline void ignore_natoms(std::istream& bin_data) {
        BinReader::ignore_int64(bin_data);
    }

    inline void read_natoms(std::istream& bin_data, int64_t& natoms) {
        bin_data.read(reinterpret_cast<char*>(&natoms), sizeof(int64_t));
    }

    inline void ignore_triclinic_and_boundary(std::istream& bin_data) {
        int triclinic;
        BinReader::read_int(bin_data, triclinic);

        BinReader::ignore_n_int(bin_data, 6);                                              // boundary
        BinReader::ignore_n_double(bin_data, 6);                                           // boundary data
        BinReader::ignore_n_double(bin_data, (triclinic > 0) ? (3 * sizeof(double)) : 0);  // triclinic data (optional)
    }

    inline void ignore_line_size(std::istream& bin_data) {
        BinReader::ignore_int(bin_data);
    }

    inline void read_line_size(std::istream& bin_data, int& line_size) {
        BinReader::read_int(bin_data, line_size);
    }

    inline void ignore_all_chunks(std::istream& bin_data) {
        BinReader::ignore_line_size(bin_data);

        int nchunk;
        BinReader::read_nchunk(bin_data, nchunk);
        for (size_t i = 0; i < nchunk; ++i) {
            int nelem;
            BinReader::read_nelem(bin_data, nelem);
            BinReader::ignore_n_double(bin_data, nelem);
        }
    }

    inline void read_nchunk(std::istream& bin_data, int& nchunk) {
        BinReader::read_int(bin_data, nchunk);
    }

    inline void read_nelem(std::istream& bin_data, int& nelem) {
        BinReader::read_int(bin_data, nelem);
    }

    inline void read_elements(std::istream& bin_data, int nelem, void* buf) {
        BinReader::read_n_double(bin_data, nelem, buf);
    }

    inline void read_info(std::istream& bin_data, int ninfo, void* buf) {
        BinReader::read_n_double(bin_data, ninfo, buf);
    }

    inline void read_xyz(std::istream& bin_data, void* buf) {
        BinReader::read_n_double(bin_data, 3, buf);
    }
}


namespace {

void skip_lines(size_t num_lines, std::ifstream& data) {
    std::string dummy_str;
    for (size_t i = 0; i < num_lines; ++i) {
        std::getline(data, dummy_str);
    }
}

void skip_frame(std::istream& bin_data, size_t nframes) {
    for (; nframes > 0; --nframes) {
        BinReader::ignore_timestep(bin_data);

        if (bin_data.eof()) {
            return;
        }

        BinReader::ignore_natoms(bin_data);

        BinReader::ignore_triclinic_and_boundary(bin_data);

        BinReader::ignore_all_chunks(bin_data);
    }
}

/**
    parse frame number from file `name`
    the `name` format must be *FRAMENUM_____*, where _____ is a integer, and * can be anything
 */
size_t parse_frame_num(const std::string name) {
    constexpr std::string_view keyword = "FRAMENUM";
    if (const size_t pos = name.find(keyword); pos != std::string::npos) {
        return std::stoull(name.c_str() + keyword.length() + pos);
    }

    return 0;
}

}

struct AtomR {
    using ParticleInfo = std::unordered_map<size_t, size_t>;  // particle type (ntype) -> particle quantity (npart)
    using Particles = std::vector<double>;                    // [x0, y0, z0, x1, y1, z1, ... ]

    double run(std::string name);
    double run_bin(std::string name);
    int init(std::string name);
    void preread(std::string name);
    void init_from_bin(std::string name);
    void preread_from_bin(std::string name);
    int MSD_calc();

    // input file columns
    std::vector<std::string> input_pattern;
    size_t type_offset = 0;

    // parameters for sorting atoms
    int64_t tot_part = 0, tot_frame = 0, tot_type = 0;  // # of total particles, frames, and atom types

    ParticleInfo particle_info;
    std::vector<size_t> itype;                          // Required atom types for MSD calculation
    std::unordered_map<size_t, size_t> itype_to_idx;    // if `itype` = [4, 2, 8], then `itype_to_idx` = [4 -> 0], [2 -> 1], [8 -> 2]
                                                        //     i.e. `type_idx` = itype_to_idx[type];
    // coords[time][type_idx]
    std::vector<std::vector<Particles>> coords;

    // Parameters for MSD
    double dt = 0.0;                                    // Timestep (ps) between frames
    size_t wi = 0, wf = 0, ws = 0;                      // Initial frame, final frame ,stride of frames
    size_t taui = 0, tauf = 0, taus = 0, tau = 0;       // Minimum, final, stride of, current span of frames
    double MSD = 0.0;                                   // Mean-Squared Displacement
};

double AtomR::run(std::string name) {
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
    this->itype.emplace_back(1);
    this->dt = 1.0;
    this->taui = 0;
    this->tauf = 10;
    this->taus = 1;
    this->wi = 0;
    this->wf = 100;
    this->ws = 1;
#else
    this->itype.emplace_back(2);
    this->dt = 1.0;
    this->taui = 0;
    this->tauf = 100000;
    this->taus = 100;
    this->wi = 0;
    this->wf = 100000;
    this->ws = 100;
#endif // !LARGE_FILE
#endif

    // Initialization
    auto init_start_time = std::chrono::system_clock::now();
    this->init(name);
    std::chrono::duration<double> init_time = std::chrono::system_clock::now() - init_start_time;
    std::cout << "    init: " << init_time.count() << " [s]" << std::endl;

    // Mean-Squared Displacement (MSD) calculation
    auto calc_start_time = std::chrono::system_clock::now();
    //this->MSD_calc();
    std::chrono::duration<double> calc_time = std::chrono::system_clock::now() - calc_start_time;
    std::cout << "    calc: " << calc_time.count() << " [s]" << std::endl;

    return (init_time + calc_time).count();
}

int AtomR::init(std::string name) {
    this->preread(name);

    std::ifstream intrj;
    intrj.open(name, std::ios::in);

    std::string dummy_str;

    size_t curr_frame = 0;
    {   // jump to the `wi`-th frame
        size_t num_skip_lines = this->wi * (9 + this->tot_part);
        skip_lines(num_skip_lines, intrj);
        curr_frame = this->wi;
    }

    {   // init `itype_to_idx`
        this->itype_to_idx.reserve(this->tot_type + 1);
        for (size_t type_idx = 0; type_idx < this->itype.size(); ++type_idx) {
            size_t type = this->itype[type_idx];
            this->itype_to_idx[type] = type_idx;
        }
    }

    {   // pre-allocate `coords`
        const size_t num_frames_to_read = this->wf - curr_frame + 1;
        this->coords.resize(num_frames_to_read);
        for (auto& particles_vec : this->coords) {
            particles_vec.resize(this->itype.size());  // required types

            for (auto& type : this->itype) {  // for particles of each type, allocate 3*num_of_particles (xyz coordinates)
                const size_t type_idx = this->itype_to_idx[type];
                const size_t num_particle_of_type = 3 * this->particle_info[type];
                particles_vec[type_idx].resize(num_particle_of_type);
            }
        }
    }

    // read frame data
    for (; !intrj.eof() && (curr_frame <= this->wf); ++curr_frame) {
        skip_lines(9, intrj);  // first 9 lines are not needed

        if (intrj.eof()) break;

        std::unordered_map<size_t, size_t> write_pos;  // record the position of `coords` to write
        write_pos.reserve(this->itype.size());

        // read particle data
        for (int i = 0; i < this->tot_part; ++i) {
            size_t type;
            intrj >> dummy_str;  // id
#ifdef HAS_MOL
            intrj >> dummy_str;  // mol
#endif
            intrj >> type;

            if (this->itype_to_idx.count(type) > 0) {
                // write to `coords[frame][type][x, y, z]`
                const size_t type_idx = this->itype_to_idx[type];

                size_t write_pos_of_curr_type = write_pos[type];
                double x, y, z;
                intrj >> x >> y >> z;
                this->coords[curr_frame][type_idx][write_pos_of_curr_type    ] = x;
                this->coords[curr_frame][type_idx][write_pos_of_curr_type + 1] = y;
                this->coords[curr_frame][type_idx][write_pos_of_curr_type + 2] = z;
                write_pos[type] += 3;
            }

            std::getline(intrj, dummy_str);  // consume the rest of the line
        }
    }

    intrj.close();

    return true;
}

/**
    preread the first frame to get the following information
        `tot_part`
        `tot_frame`
        `tot_type`
        `particle_info`
 */
void AtomR::preread(std::string name) {
    this->tot_frame = parse_frame_num(name);

    std::ifstream intrj;
    intrj.open(name, std::ios::in);
    std::string str;

    skip_lines(3, intrj);

    // total number of particles
    intrj >> this->tot_part;
    std::getline(intrj, str);

    skip_lines(5, intrj);

    // number of particles by type
    for (int i = 0; i < this->tot_part; ++i) {
        size_t type;

        intrj >> str;  // skip id
#ifdef HAS_MOL
        intrj >> str;  // skip mol
#endif
        intrj >> type;
        ++this->particle_info[type];

        std::getline(intrj, str);  // skip the rest
    }

    this->tot_type = this->particle_info.size();
}

double AtomR::run_bin(std::string name) {
#ifndef LARGE_FILE
    this->itype.emplace_back(1);
    this->dt = 1.0;
    this->taui = 0;
    this->tauf = 10;
    this->taus = 1;
    this->wi = 0;
    this->wf = 100;
    this->ws = 1;
#else
    this->input_pattern = {"id", "mol", "type", "x", "y", "z", "ax", "ay", "az"};
    this->itype.emplace_back(2);
    this->dt = 1.0;
    this->taui = 0;
    this->tauf = 100000;
    this->taus = 100;
    this->wi = 0;
    this->wf = 100000;
    this->ws = 100;
#endif

    // Initialization
    auto init_start_time = std::chrono::system_clock::now();
    this->init_from_bin(name);
    std::chrono::duration<double> init_time = std::chrono::system_clock::now() - init_start_time;
    std::cout << "    init: " << init_time.count() << " [s]" << std::endl;

    // Mean-Squared Displacement (MSD) calculation
    auto calc_start_time = std::chrono::system_clock::now();
    //this->MSD_calc();
    std::chrono::duration<double> calc_time = std::chrono::system_clock::now() - calc_start_time;
    std::cout << "    calc: " << calc_time.count() << " [s]" << std::endl;

    return (init_time + calc_time).count();
}

void AtomR::init_from_bin(std::string name) {
    this->preread_from_bin(name);

    std::ifstream bin_data(name, std::ios::in | std::ios::binary);

    // skip to frame `wi`
    skip_frame(bin_data, this->wi);
    size_t curr_frame = this->wi;

    {   // init `itype_to_idx`
        this->itype_to_idx.reserve(this->tot_type + 1);
        size_t type_idx = 0;
        for (auto type : this->itype) {
            this->itype_to_idx[type] = type_idx++;
        }
    }
    
    {   // pre-allocate of `coords`
        const size_t num_frames_to_read = this->wf - curr_frame + 1;
        this->coords.resize(num_frames_to_read);
        for (auto& particles_vec : this->coords) {
            particles_vec.resize(this->itype.size());  // required types

            for (auto& type : this->itype) {  // for particles of each type, allocate 3*num_of_particles (xyz coordinates)
                const size_t type_idx = this->itype_to_idx[type];
                const size_t num_particle_of_type = 3 * this->particle_info[type];
                particles_vec[type_idx].resize(num_particle_of_type);
            }
        }
    }

    // initialize variables used in the loops below
    std::unordered_map<size_t, size_t> write_pos;  // record the position of `coords` to write
    write_pos.reserve(this->itype.size());

    // whether current input binary file contains 'mol' data
    const bool is_has_mol =
        this->input_pattern.cend() != std::find_if(this->input_pattern.cbegin(), this->input_pattern.cend(),
                                                   [] (const std::string_view str) {
                                                       return str.compare("mol") == 0;
                                                   });
    const size_t info_len = is_has_mol ? 3 : 2;  // the first `info_len` elements are not coordinate data
    std::vector<double> info(info_len, 0.0);

    // read frame data
    for (; !bin_data.eof() && (curr_frame <= this->wf); ++curr_frame) {
        write_pos.clear();

        BinReader::ignore_timestep(bin_data);

        if (bin_data.eof())
            bin_data.close();
            return;

        BinReader::ignore_natoms(bin_data);

        BinReader::ignore_triclinic_and_boundary(bin_data);

        int size_of_one_line;
        BinReader::read_line_size(bin_data, size_of_one_line);

        int nchunk;
        BinReader::read_nchunk(bin_data, nchunk);

        for (; nchunk > 0; --nchunk) {
            int nelems;
            BinReader::read_nelem(bin_data, nelems);

            int nlines = nelems / size_of_one_line;
            for (; nlines > 0; --nlines) {
                int nelems_unread = size_of_one_line;

                BinReader::read_info(bin_data, info_len, info.data());
                nelems_unread -= static_cast<int>(info_len);

                size_t type = static_cast<size_t>(info[this->type_offset]);
                if (this->itype_to_idx.count(type) > 0) {  // required atom type, read x, y, z
                    const size_t type_idx = this->itype_to_idx[type];
                    const size_t write_pos_of_curr_type = write_pos[type];
                    BinReader::read_xyz(bin_data, &this->coords[curr_frame][type_idx][write_pos_of_curr_type]);

                    write_pos[type] += 3;
                } else {  // not required atom type, ignore
                    BinReader::ignore_n_double(bin_data, 3);
                }

                nelems_unread -= 3;

                BinReader::ignore_n_double(bin_data, nelems_unread);  // ignore trailing elements
            }
        }
    }

    bin_data.close();
}

/**
    preread the first frame to get the following information
        `type_offset`
        `tot_part`
        `tot_frame`
        `tot_type`
        `particle_info`
 */
void AtomR::preread_from_bin(std::string name) {
    // `tot_frame`
    this->tot_frame = parse_frame_num(name);

    std::ifstream bin_data(name, std::ios::in | std::ios::binary);

    BinReader::ignore_timestep(bin_data);

    if (bin_data.eof()) {
        bin_data.close();
        return;
    }

    // `tot_part`
    BinReader::read_natoms(bin_data, this->tot_part);

    BinReader::ignore_triclinic_and_boundary(bin_data);

    int size_of_one_line;
    BinReader::read_line_size(bin_data, size_of_one_line);
    assert(size_of_one_line == this->input_pattern.size());  // make `size_of_one_line` matches the `input_pattern`

    int nchunk;
    BinReader::read_nchunk(bin_data, nchunk);

    // `type_offset`
    auto iter = std::find_if(this->input_pattern.cbegin(), this->input_pattern.cend(), [] (const std::string_view str) {
        return str.compare("type") == 0;
    });
    assert(iter != this->input_pattern.cend());
    this->type_offset = iter - this->input_pattern.cbegin();

    // `particle_info`
    std::vector<double> buf;
    for(; nchunk > 0; --nchunk) {
        int nelem;
        BinReader::read_nelem(bin_data, nelem);

        buf.resize(nelem, 0.0);
        BinReader::read_elements(bin_data, nelem, buf.data());

        for (size_t idx = 0; idx < buf.size(); idx += size_of_one_line) {
            // assume `size_of_one_line` = 5, i.e. id, type, x, y, z
            size_t type = static_cast<size_t>(buf[idx + this->type_offset]);
            ++this->particle_info[type];
        }
    }

    // `tot_type`
    this->tot_type = this->particle_info.size();

    bin_data.close();
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

    // # of total used particles
    const size_t N = std::accumulate(this->itype.cbegin(), this->itype.cend(), 0, [&] (size_t tmp_sum, size_t type) {
        return tmp_sum + this->particle_info[type];
    });

    auto square_error = [](const double a, const double b) -> const double {
        const double e = a - b;
        return e * e;
    };

    constexpr double INIT_VAL = 0.0;
    for (this->tau = this->taui; this->tau <= this->tauf; this->tau += this->taus) {
        this->MSD = 0;

        //  t0 <= (wf - wi - tau - 1)   has been rearranged into   (t0 + wi + tau + 1) <= wf
        // because `t0`, `wi`, `wf`, `tau` are size_t, and the minus operation might lead to overflow (< 0)
        //
        // of course, plus operation could also overflow,
        // but that happens only when all of the numbers are huge
        for (size_t t0 = 0; (t0 + this->wi + this->tau + 1) <= this->wf; t0 += this->ws) {  // Reference frame
            for (auto& type: this->itype) {
                const size_t type_idx = this->itype_to_idx[type];

                // transform reduce is the c++ implementation of Google MapReduce
                // which is designed for large scale parallelization
                const auto& vec_a = this->coords[t0 + this->tau][type_idx];
                const auto& vec_b = this->coords[t0            ][type_idx];
                this->MSD += std::transform_reduce(vec_a.cbegin(), vec_a.cend(), vec_b.cbegin(), INIT_VAL, std::plus<>(), square_error);
            }
        }

        this->MSD *= ((this->ws / ((double)this->wf - this->wi - this->tau)) / N);
#ifdef OUTPUT_MSD
        outf << this->tau * this->dt << " " << this->MSD << std::endl;
#endif
    }

#ifdef OUTPUT_MSD
    outf.close();
#endif

    return 0;
}
