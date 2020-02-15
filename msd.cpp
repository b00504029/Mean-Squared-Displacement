// Read the input one time and simultaneously construct the array
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
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
    int istep, itype;                         // ith step, atom type

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
};

int main(){
    // Construct a system
    ATOM atom;

    // File name to read
    cout<<"Enter the input file name: ";
    cin>>atom.name;

    // Atom type
    cout<<"Enter the atom type: ";
    cin>>atom.itype;

    // Timestep in ps
    cout<<"Enter the timestep in ps: ";
    cin>>atom.dt;

    // Lag times t from taui to tauf by taus
    cout<<"Enter the initial lag time, final lag time,and a stride of lag times: ";
    cin>>atom.taui>>atom.tauf>>atom.taus;

    // For analysis
    cout<<"Enter the initial frame, final frame,and a span of frames: ";
    cin>>atom.ti>>atom.tf>>atom.ts;

    // Initialization
    atom.init();

    // Mean-Squared Displacement (MSD) calculation
    atom.MSD_calc();

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
    for (tot_frame = 0; !intrj.eof() && tot_frame <= tf; tot_frame++){
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

            if (tot_frame == 0){
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
                npart[temptype-1] ++;
                for (j = 0; j < 3; j++){
                    intrj>>frame_data[temptype-1][npart[temptype-1]-1][j];
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

    int t0;                                      // Referecnce frame
    double temp;                                 // Temporal average
    for (tau = taui; tau <= tauf; tau += taus){
          MSD = temp = 0;
          for (t0 = ti; t0 <= (tf-tau-1); t0+=ts){
                temp = 0;
                for (j = 0; j < npart[itype-1]; j++){
                      for (k = 0; k < 3; k++){
                            temp += pow( x[t0+tau][itype-1][j][k] - x[t0][itype-1][j][k] ,2);
                      }
                }
                MSD += temp;
          }
          MSD *= ((ts/double(tf-ti+1-tau))/npart[itype-1]);
          outf<<tau<<" "<<MSD<<endl;
    }
    outf.close();

    return 0;
}
