// Read the input 3 times to construct the array
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

    std::vector< int > ipart;                 // # of particles with different atom types when counting
    std::vector< int > npart;                 // # of total particles with different atom types
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

    int temptype=0, nexttype=0;
    for (tot_frame = -1; !intrj.eof(); tot_frame++){
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);

        intrj>>tot_part;       getline(intrj,str);

        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);

        for (i = 0; i < tot_part; i++){
            intrj>>str;
            intrj>>temptype;
            if (temptype > nexttype) nexttype = temptype;
            getline(intrj,str);
        }
    }
    tot_type = nexttype;
    intrj.close();

    // Cconstruct the array for the system
    for (i = 0; i < tot_type; i++){
        ipart.emplace_back(0);
        npart.emplace_back(0);
    }

    intrj.open(name,ios::in);
    getline(intrj,str);
    getline(intrj,str);
    getline(intrj,str);
    getline(intrj,str);
    getline(intrj,str);
    getline(intrj,str);
    getline(intrj,str);
    getline(intrj,str);
    getline(intrj,str);

    for (i = 0; i < tot_part; i++){
        intrj>>id;
        intrj>>temptype;
        npart[temptype-1] ++;
        getline(intrj,str);
    }
    intrj.close();

    for (i = 0; i < tot_frame; i++) {
        std::vector< std::vector< std::vector< double > > > frame_data;
        for (j = 0; j < tot_type ; j++) {
            std::vector< std::vector< double > > type_data;
            for (k = 0; k < npart[j] ; k++){
                type_data.emplace_back(std::vector< double >{0, 0, 0});
            }
            frame_data.emplace_back(type_data);
        }
        x.emplace_back(frame_data);
    }

    // Read and sort coodinate for all particles
    intrj.open(name,ios::in);
    for (i=0 ; !intrj.eof(); i++){
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);
        getline(intrj,str);

        if (intrj.eof()) break;
        for (j = 0; j < tot_type; j++){
            ipart[j]=0;
        }
        for (k = 0; k < tot_part; k++){
            intrj>>id;
            intrj>>temptype;
            for (l = 0; l < 3; l++){
                intrj>>x[i][temptype-1][ipart[temptype-1]][l];
            }
            ipart[temptype-1] ++;
            getline(intrj,str);
        }
    }
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
