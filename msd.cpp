// Read the input one time and simultaneously construct the array
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <vector>
#include <algorithm>
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
    int tot_part, tot_frame, tot_type;        // # of total particles, frames, and atom types
    int istep;                                // ith step

    std::vector< std::vector < int > > id;    // ID for all atoms read in order
    std::vector< int >::iterator id_it;       // Record the position of the ID in the id array
    std::vector< int >::iterator it;          // Record the position of the number in an array
    std::vector< int > itype;                 // Required atom types for MSD calculation
    std::vector< int > ntype;                 // # of atom types: {1,2,3,...,n}
    std::vector< int > npart;                 // # of total particles with different atom types: {125,250,....}
    std::vector< std::vector< std::vector< std::vector< double > > > > x;      // Coordinates of atoms

    bool subtract_com;                        // Flag for subtracting the drift
    std::vector< std::vector< double > > xcm;   // Coordinates of the center of mass

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

    // Decide if the drift should be removed
    cout<<"Enter 1/0 (Yes/No) for the drift removal: ";
    cin>>atom.subtract_com;

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
    int ID;
    id.resize(itype.size());
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
            intrj>>ID;
            intrj>>temptype;
            frame_data.resize(itype.size());
            it = std::find(itype.begin(), itype.end(), temptype);

            if (tot_frame == ti){
                if (temptype > maxtype) {
                    for (j = maxtype; j < temptype; j++){
                        ntype.emplace_back(j+1);
                        npart.emplace_back(0);
                    }
                    maxtype = temptype;
                }
                npart[temptype-1] ++;
                if (it != itype.end()){
                    id[distance(itype.begin(), it)].emplace_back(ID);
                    frame_data[distance(itype.begin(), it)].resize(npart[temptype-1]);
                    for (j = 0; j < 3; j++){
                        intrj>>coord;
                        frame_data[distance(itype.begin(), it)][npart[temptype-1]-1].emplace_back(coord);
                    }
                }
            }
            else{
                for (j = 0; j < itype.size(); j++){
                    if (temptype == itype[j]) {
                        npart[temptype-1] ++;
                        for (k = 0; k < 3; k++){
                            id_it = std::find(id[distance(itype.begin(), it)].begin(), id[distance(itype.begin(), it)].end(), ID);
                            intrj>>frame_data[distance(itype.begin(), it)][distance(id[distance(itype.begin(), it)].begin(), id_it)][k];
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

    // Initialization of center of mass
    xcm.resize(tf-ti);
    for (i = 0; i < ( tf-ti ); i++){
        for (j = 0; j < 3; j++){
            xcm[i].emplace_back(0);
        }
    }

    // If the drift should be cancelled
    int n;                                       // # of total used atoms
    if (subtract_com){
        for (i = 0; i < ( tf-ti ); i++){
            n = 0;
            for (j = 0; j < itype.size(); j++){
                for (k = 0; k < npart[itype[j]-1]; k++){
                    for (l = 0; l < 3; l++){
                        xcm[i][l] += x[i][j][k][l];
                    }
                }
                n += npart[itype[j]-1];
            }
            for (j = 0; j < 3; j++){
                xcm[i][j] *= 1/double(n);
            }
        }
    }

    // MSD calculation
    int t0;                                      // Referecnce frame
    double temp;                                 // Temporal average
    for (tau = taui; tau <= tauf; tau += taus){
        MSD = 0;
        for (t0 = 0; t0 <= (tf-ti-tau-1); t0+=ts){
            temp = n = 0;
            for (i = 0; i < itype.size(); i++){
                for (j = 0; j < npart[itype[i]-1]; j++){
                    for (k = 0; k < 3; k++){
                        temp += pow( ( x[t0+tau][i][j][k]-xcm[t0+tau][k] ) - ( x[t0][i][j][k]-xcm[t0][k] ) ,2);
                    }
                }
                n += npart[itype[i]-1];
            }
            MSD += temp;
        }
        MSD *= ((1.0/((tf-ti+1-tau)/ts)/n));
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
