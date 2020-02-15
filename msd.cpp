// Read the input one time and simultaneously construct the array
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <vector>
#include <ctime>
#include <chrono>
using namespace std;

class ATOM
{
    public:
    string str, name;                         // Store the instruction, the file name
    double xl, xh, yl, yh, zl, zh;            // Boundary coordinates

    // parameters for sorting atoms
    int i, j ,k, l;                           // Counters
    int id;                                   // ID for all atoms read in order
    int tot_part, tot_frame, tot_type;        // # of total particles, frames, and atom types
    int istep;                                // ith step

    std::vector< int > itype;                 // Required atom types for MSD calculation
    std::vector< int > ntype;                 // # of atom types: {1,2,3,...,n}
    std::vector< int > npart;                 // # of total particles with different atom types: {125,250,....}
    std::vector< std::vector< std::vector< std::vector< double > > > > x;      // Coordinates

    // Parameters for MSD
    double dt;                                // Timestep (ps) between frames
    int ti, tf, ts;                           // Initial frame, final frame ,stride of frames
    int taui, tauf, taus, tau;                // Minimum, final, stride of, current span of frames
    double MSD;                               // Mean-Squared Displacement
    int MSD_calc();

    // Initialization
    int init();

    // Skipping frames unwanted until ti
    void skip_lines(int ,ifstream* );
};

int main(){
    // Construct a system
    ATOM atom;

    // File name to read
    cout<<"Enter the input file name: ";
    cin>>atom.name;

    // Atom type
    int a;
    cout<<"Enter the atom type: ";
    do {
        cin>>a;
        atom.itype.emplace_back(a);
    } while (cin.get() != '\n');

    // Timestep in ps
    cout<<"Enter the timestep in ps: ";
    cin>>atom.dt;

    // Lag times t from taui to tauf by taus
    cout<<"Enter the initial lag time, final lag time,and a stride of lag times: ";
    cin>>atom.taui>>atom.tauf>>atom.taus;

    // For analysis
    cout<<"Enter the initial frame, final frame,and a span of frames: ";
    cin>>atom.ti>>atom.tf>>atom.ts;

    // Calculate the wall-clock time
    auto wcts = std::chrono::system_clock::now();

    // Initialization
    atom.init();

    // Mean-Squared Displacement (MSD) calculation
    atom.MSD_calc();

    // End of calculation
    std::chrono::duration< double > wctduration = (std::chrono::system_clock::now() - wcts);
    std::cout<< "Finished in "<< wctduration.count() << " seconds [Wall Clock]" <<std::endl;

    return 0;
}

int ATOM::init()
{
    // Determine the values of tot_part, tot_frame, and tot_type
    ifstream intrj;
    intrj.open(name,ios::in);

    int temptype=0, maxtype=0;
    double coord;                     // Temporal coordinates
    std::vector < std::vector < std::vector < double > > > frame_data;

    // Skip the lines unwanted
    for (tot_frame = 0; tot_frame < ti; tot_frame++){
        if (tot_frame == 0){
            getline(intrj,str);
            getline(intrj,str);
            getline(intrj,str);

            intrj>>tot_part;       getline(intrj,str);

            skip_lines(5+tot_part,& intrj);
        }
        else {
            skip_lines(9+tot_part,& intrj);
        }
    }

    // Read the desired lines
    for (; !intrj.eof() && tot_frame <= tf; tot_frame++){
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);

        intrj>>tot_part;       getline(intrj,str);

        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);
        if (intrj.eof()) break;

        for (j = 0; j < maxtype; j++){
            npart[j]=0;
        }
        for (i = 0; i < tot_part; i++){
            intrj>>str;
            intrj>>temptype;

            if (tot_frame == ti){
                if (temptype > maxtype) {
                    for (j = maxtype; j < temptype; j++){
                        ntype.emplace_back(j+1);
                        npart.emplace_back(0);
                    }
                    frame_data.resize(temptype);
                    maxtype = temptype;
                }
                npart[temptype-1] ++;
                frame_data[temptype-1].resize(npart[temptype-1]);
                for (j = 0; j < 3; j++){
                    intrj>>coord;
                    frame_data[temptype-1][npart[temptype-1]-1].emplace_back(coord);
                }
            }
            else{
                for (j = 0; j < itype.size(); j++){
                    if (temptype == itype[j]) {
                        npart[temptype-1] ++;
                        for (k = 0; k < 3; k++){
                            intrj>>frame_data[temptype-1][npart[temptype-1]-1][k];
                        }
                    }
                }
            }
            getline(intrj,str);
        }
        x.emplace_back(frame_data);
    }
    tot_type = maxtype;
    intrj.close();
    return 0;
}

int ATOM::MSD_calc()
{
    ofstream outf;
    outf.open("msd.txt",ios::out);
    outf<<"t(ps)"<<" "<<"MSD(angstroms)"<<endl;

    int t0, n=0;                                 // Referecnce frame, # of total used particles
    double temp;                                 // Temporal average
    for (tau = taui; tau <= tauf; tau += taus){
        MSD = 0;
        for (t0 = 0; t0 <= (tf-ti-tau-1); t0+=ts){
            temp = n = 0;
            for (i = 0; i < itype.size(); i++){
                for (j = 0; j < npart[itype[i]-1]; j++){
                    for (k = 0; k < 3; k++){
                        temp += pow( x[t0+tau][itype[i]-1][j][k] - x[t0][itype[i]-1][j][k] ,2);
                    }
                }
                n += npart[itype[i]-1];
            }
            MSD += temp;
        }
        MSD *= ((ts/double(tf-ti+1-tau))/n);
        outf<<tau*dt<<" "<<MSD<<endl;
    }
    outf.close();

    return 0;
}

void ATOM::skip_lines(int num_lines, ifstream* data)
{
    for (i = 0; i < num_lines; i++){
        getline(*data,str);
    }
}
