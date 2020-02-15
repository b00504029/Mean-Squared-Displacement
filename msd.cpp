// Read the input only when the frames need to be calculated
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
    string str, name;                         					// Store the instruction, the file name
    double xl, xh, yl, yh, zl, zh;            					// Boundary coordinates

    // Parameters for sorting particles
    int i, j ,k, l;                           					// Counters
    int tot_part, tot_frame, tot_type;        					// Total # of particles, frames, and atom types
    int n_chain;                              					// # of chains in one NIMs particle
    int n_bead;                               					// # of beads per chain
    int n_NIMs;                               					// # of NIMs in the system
	int n_site;													// # of interaction sites in one NIMs particle
	int tot_chain;												// # of all chains in the system
	
	// Useful vectors and variables
	bool subtract_com;                        					// Flag for subtracting the drift
	vector< double > mass;                    					// Masses stored in order of atom types
    vector< vector< double > > xcm;           					// Coordinates of the center of mass: [frame][coordinate]
	vector< vector< vector< vector< double > > > > x;      		// Coordinates of atoms: [2 frames][type][number][coordinate]
	vector< vector< vector< double > > > xcm_chain;  			// Coordinates of C.O.M. of polymer chains: [2 frames][number][coordinate]    
	vector< int > itype;                      					// Required atom types for MSD calculation
    vector< int > npart;                      					// # of total particles with different atom types required: {125,250,....}
    vector< vector < int > > id;             					// ID for all atoms read in order
	vector< int >::iterator it_type;       						// Record the position of the type in the itype array
	vector< int >::iterator it_id;       						// Record the position of the ID in the id array

    // Parameters for MSD and the standard deviation of SD
    double dt;                                					// Timestep (ps) between frames
    int ti, tf, ts;                           					// Initial frame, final frame ,stride of frames
    int taui, tauf, taus, tau;                					// Minimum, final, stride of, current span of frames
    double MSD;                               					// Mean-Square Displacement
    double SD;                                					// Square displacement
    double upper_STD;                         					// Upper STD of square displacement
    double lower_STD;                         					// Lower STD of square displacement
    int MSD_calc();

    // Initialization
    int init();

	// Input and output
	ifstream intrj;												// Load the trajectory
	ofstream outf;												// Write the MSD result
	ofstream outf_error;										// Write the standard deviation of square difference

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
    cout<<"Enter the atom types in the chain used for MSD calculation: ";
    do {
        cin>>a;
        atom.itype.emplace_back(a);
    } while (cin.get() != '\n');

    // Timestep in ps
    cout<<"Enter the timestep in ps between frames: ";
    cin>>atom.dt;

    // Lag times t from taui to tauf by taus
    cout<<"Enter the initial lag time, the final lag time,and a stride of lag times: ";
    cin>>atom.taui>>atom.tauf>>atom.taus;

    // For analysis
    cout<<"Enter the initial frame, the final frame,and a span of frames: ";
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
    // Assign the initial values to all parameters
    // Atom types: 1=B, 2=C, 3=N, 4=P
	tot_frame = tf-ti+1;
	tot_type = 4;	
	id.resize(tot_type);
	npart.resize(tot_type,0);
	
	// Masses
    mass.resize(tot_type);
    mass[0] = 44;
    mass[1] = 5500;
    mass[2] = 44;
    mass[3] = 44;
	
	// System parameters
	n_chain = 25;
	n_bead = 15;
	n_NIMs = 100;
	n_site = n_chain*n_bead+1;
	tot_chain = n_chain*n_NIMs;
    tot_part = n_site * n_NIMs;

	// Read the # of particles of all atom types	
	intrj.open(name,ios::in);
	
	int temp_id;										// Temporary ID	
    int temp_type;										// Temporary atom type	
	skip_lines(9,& intrj);
	for (i = 0; i < tot_part; i++){
		intrj >> temp_id;
		intrj >> temp_type;	
		it_type = find(itype.begin(), itype.end(), temp_type);
		
		if (it_type != itype.end()){
			npart[ temp_type-1 ]++;
			id[ temp_type-1 ].emplace_back(temp_id);
		}	
		getline(intrj,str);
	}
	
	intrj.close();
	
	// Construct the array x
	x.resize(2);										// 2 frames for the current frame and the next frame
	for (i = 0; i < x.size(); i++) {
		x[i].resize( tot_type );
		for (j = 0; j < tot_type; j++) {
			x[i][j].resize( npart[j] );
			for (k = 0; k < npart[j]; k++) {
				x[i][j][k].resize(3, 0);
			}
		}
	}
	
	// Construct the coordinates of COM of atoms used for MSD
    xcm.resize(2);										// 2 frames for the current frame and the next frame
    for (i = 0; i < xcm.size(); i++)
        for (j = 0; j < 3; j++)
            xcm[i].emplace_back(0);
	
	// Construct the array xcm_chain
	xcm_chain.resize(2);								// 2 frames for the current frame and the next frame
	for (i = 0; i < xcm_chain.size(); i++){
		xcm_chain[i].resize(tot_chain);
        for (j = 0; j < xcm_chain[i].size(); j++)
			for (k = 0; k < 3; k++)
				xcm_chain[i][j].emplace_back(0);
	}
	
    return 0;
}

int ATOM::MSD_calc()
{
	// Output the MSD result
    outf.open("msd_chain.txt",ios::out);
    outf<<"t(ps)"<<" "<<"MSD(angstrom^2)"<<endl;

/*     outf_error.open("std_chain.txt",ios::out);
    outf_error<<"t(ps)"<<" "<<"STD(angstrom^2)"<<endl; */

	// Useful variables for MSD calculation
	int n;												// Count the total particles used for calculation
	int t0;                                      		// Referecnce frame
	int temp_id;										// Temporary ID			
    int temp_type;										// Temporary atom type
	int i_NIMs;											// i-th NIMs particle
	int i_chain;										// i-th chain
	long long init_pos;									// Record the position of the initial frame		
	long long cur_pos;									// Record the position of the current frame
	long long next_pos;									// Record the position of the next frame
	double tot_mass;									// Total mass of used particles
	double temp_x[3];									// Temporary coordinates
	
	// Read the required data and perform MSD calculation
	intrj.open(name,ios::in);
    for (tau = taui; tau <= tauf; tau += taus){
		MSD = 0;
		// Skip the lines unwanted until the first frame
		if (tau == taui){
			for (i = 0; i < ti; i++){
				skip_lines(9+tot_part,& intrj);
			}
			init_pos = intrj.tellg();	
		}
		else {
			intrj.seekg(init_pos, intrj.beg);
		}
		cur_pos = intrj.tellg();
		
        for (t0 = 0; t0 <= (tf-ti-tau-1); t0 += ts){
			intrj.seekg(cur_pos, intrj.beg);
			if (t0 != 0){
				for (i = 0; i < ts; i++){
					skip_lines(9+tot_part,& intrj);
				}
			}
			cur_pos = intrj.tellg(); 
			skip_lines(9,& intrj);
			
			// Initialize the xcm
			for (i = 0; i < xcm.size(); i++)
				for (j = 0; j < xcm[i].size(); j++)
					xcm[i][j] = 0;
			
			// Initialize the xcm_chain
			for (i = 0; i < xcm_chain.size(); i++)
				for (j = 0; j < xcm_chain[i].size(); j++)
					for (k = 0; k < xcm_chain[i][j].size(); k++)
						xcm_chain[i][j][k] = 0;
			
			// Read the current frame (frame index = 0)
			tot_mass = 0;
			for (i = 0; i < tot_part; i++){
				intrj >> temp_id;
				intrj >> temp_type;
				it_type = find(itype.begin(), itype.end(), temp_type);
				it_id = find(id[ temp_type-1 ].begin(), id[ temp_type-1 ].end(), temp_id);
				
				if (it_type != itype.end()) {
					for (j = 0; j < 3; j++) {
						// Save the coordinates of the current frame
						intrj >> temp_x[j];
						x[0][ temp_type-1 ][ distance( id[ temp_type-1 ].begin(), it_id)][j] = temp_x[j];
						
						// Determine the chain index
						i_NIMs = (temp_id-1)/n_site + 1;
						i_chain = ceil(double((temp_id-(i_NIMs-1)*n_site)-1)/n_bead)+(i_NIMs-1)*n_chain;
						xcm_chain[0][ i_chain-1 ][j] += mass[ temp_type-1 ] * temp_x[j];
						
						// Calculate the C.O.M. of all atoms in the system
						if (subtract_com) xcm[0][j] += mass[ temp_type-1 ] * temp_x[j];
					}
					tot_mass += mass[ temp_type-1 ];
				}
				getline(intrj,str);
			}
			
			for (i = 0; i < xcm_chain[0].size(); i++)
				for (j = 0; j < 3; j++)
					xcm_chain[0][i][j] /= (tot_mass/tot_chain);

			for (j = 0; j < 3; j++)
				xcm[0][j] /= tot_mass;
			
			// Skip the lines between two frames
			if (tau == 0) intrj.seekg(init_pos, intrj.beg);
			if (t0 == 0){
				for (i = 0; i < tau-1; i++){
					 skip_lines(9+tot_part,& intrj);
				}
			}
			else {
				intrj.seekg(next_pos, intrj.beg);
				for (i = 0; i < ts; i++){
					skip_lines(9+tot_part,& intrj);
				}	
			}
			next_pos = intrj.tellg();
			skip_lines(9,& intrj);
			
			// Read the next frame (frame index = 1)
			for (i = 0; i < tot_part; i++){
				intrj >> temp_id;
				intrj >> temp_type;
				it_type = find(itype.begin(), itype.end(), temp_type);
				it_id = find(id[ temp_type-1 ].begin(), id[ temp_type-1 ].end(), temp_id);
				
				if (it_type != itype.end()) {
					for (j = 0; j < 3; j++) {
						// Save the coordinates of the next frame
						intrj >> temp_x[j];
						x[1][ temp_type-1 ][ distance( id[ temp_type-1 ].begin(), it_id)][j] = temp_x[j];
						
						// Determine the chain index
						i_NIMs = (temp_id-1)/n_site + 1;
						i_chain = ceil(double((temp_id-(i_NIMs-1)*n_site)-1)/n_bead)+(i_NIMs-1)*n_chain;
						xcm_chain[1][ i_chain-1 ][j] += mass[ temp_type-1 ] * temp_x[j];
						
						// Calculate the C.O.M. of all atoms in the system
						if (subtract_com) xcm[1][j] += mass[ temp_type-1 ] * temp_x[j];
					}
				}
				getline(intrj,str);
			}
			
			for (i = 0; i < xcm_chain[1].size(); i++)
				for (j = 0; j < 3; j++)
					xcm_chain[1][i][j] /= (tot_mass/tot_chain);
					
			for (j = 0; j < 3; j++)
				xcm[1][j] /= tot_mass;
			
			// MSD calculation
			for (i = 0; i < tot_chain; i++){
				for (j = 0; j < 3 ; j++){
					MSD += pow( ( xcm_chain[1][i][j]-xcm[1][j] ) - ( xcm_chain[0][i][j]-xcm[0][j] ), 2);
				}
			}
        }
        MSD *= ((1.0/((tf-ti+1-tau)/ts)/tot_chain));
        outf<<tau*dt<<" "<<MSD<<endl;

/*         // Useful variables for STD calculation
        int n_upper;                         			// Counts for SD > MSD
        int n_lower;                         			// Counts for SD < MSD
        double temp_upper;                   			// Temporary value for upper standard deviation calculation
        double temp_lower;                   			// Temporary value for lower standard deviation calculation
        upper_STD = lower_STD = n_upper = n_lower = 0;
		
		intrj.seekg(init_pos, intrj.beg);
		cur_pos = intrj.tellg();
		
        for (t0 = 0; t0 <= (tf-ti-tau-1); t0+=ts){
			intrj.seekg(cur_pos, intrj.beg);
			if (t0 != 0){
				for (i = 0; i < ts; i++){
					skip_lines(9+tot_part,& intrj);
				}
			}
			cur_pos = intrj.tellg(); 
			skip_lines(9,& intrj);
			
			// Initialize the xcm
			for (i = 0; i < xcm.size(); i++)
				for (j = 0; j < xcm[i].size(); j++)
					xcm[i][j] = 0;
				
			// Initialize the xcm_chain
			for (i = 0; i < xcm_chain.size(); i++)
				for (j = 0; j < xcm_chain[i].size(); j++)
					for (k = 0; k < xcm_chain[i][j].size(); k++)
						xcm_chain[i][j][k] = 0;
				
			// Read the current frame (frame index = 0)
			tot_mass = 0;
			for (i = 0; i < tot_part; i++){
				intrj >> temp_id;
				intrj >> temp_type;
				it_type = find(itype.begin(), itype.end(), temp_type);
				it_id = find(id[ temp_type-1 ].begin(), id[ temp_type-1 ].end(), temp_id);
				
				if (it_type != itype.end()) {
					for (j = 0; j < 3; j++) {
						// Save the coordinates of the current frame
						intrj >> temp_x[j];
						x[0][ temp_type-1 ][ distance( id[ temp_type-1 ].begin(), it_id) ][j] = temp_x[j];
						
						// Determine the chain index
						i_NIMs = (temp_id-1)/n_site + 1;
						i_chain = ceil(double((temp_id-(i_NIMs-1)*n_site)-1)/n_bead)+(i_NIMs-1)*n_chain;
						xcm_chain[0][ i_chain-1 ][j] += mass[ temp_type-1 ] * temp_x[j];
						
						// Calculate the C.O.M. of all atoms in the system
						if (subtract_com) xcm[0][j] += mass[ temp_type-1 ] * temp_x[j];
					}
					tot_mass += mass[ temp_type-1 ];
				}
				getline(intrj,str);
			}
			
			for (i = 0; i < xcm_chain[0].size(); i++)
				for (j = 0; j < 3; j++)
					xcm_chain[0][i][j] /= (tot_mass/tot_chain);
			
			for (j = 0; j < 3; j++) 
				xcm[0][j] /= tot_mass;
			
			// Skip the lines between two frames
			if (tau == 0) intrj.seekg(init_pos, intrj.beg);
			if (t0 == 0){
				for (i = 0; i < tau-1; i++){
					 skip_lines(9+tot_part,& intrj);
				}
			}
			else {
				intrj.seekg(next_pos, intrj.beg);
				for (i = 0; i < ts; i++){
					skip_lines(9+tot_part,& intrj);
				}
			}
			next_pos = intrj.tellg();
			skip_lines(9,& intrj);
			
			// Read the next frame (frame index = 1)
			for (i = 0; i < tot_part; i++){
				intrj >> temp_id;
				intrj >> temp_type;
				it_type = find(itype.begin(), itype.end(), temp_type);
				it_id = find(id[ temp_type-1 ].begin(), id[ temp_type-1 ].end(), temp_id);
				
				if (it_type != itype.end()) {
					for (j = 0; j < 3; j++) {
						// Save the coordinates of the next frame
						intrj >> temp_x[j];
						x[1][ temp_type-1 ][ distance( id[ temp_type-1 ].begin(), it_id) ][j] = temp_x[j];
						
						// Determine the chain index
						i_NIMs = (temp_id-1)/n_site + 1;
						i_chain = ceil(double((temp_id-(i_NIMs-1)*n_site)-1)/n_bead)+(i_NIMs-1)*n_chain;
						xcm_chain[1][ i_chain-1 ][j] += mass[ temp_type-1 ] * temp_x[j];
						
						// Calculate the C.O.M. of all atoms in the system
						if (subtract_com) xcm[1][j] += mass[ temp_type-1 ] * temp_x[j];
					}
				}
				getline(intrj,str);
			}
			for (i = 0; i < xcm_chain[1].size(); i++)
				for (j = 0; j < 3; j++)
					xcm_chain[1][i][j] /= (tot_mass/tot_chain);
			
			for (j = 0; j < 3; j++)
				xcm[1][j] /= tot_mass;
			
			// Standard deviation calculation
            temp_upper = temp_lower = 0;
            for (i = 0; i < tot_chain; i++){
				SD = 0;
				for (j = 0; j < 3; j++){
					SD += pow( ( xcm_chain[1][i][j]-xcm[1][j] ) - ( xcm_chain[0][i][j]-xcm[0][j] ) ,2);
				}
				if (SD > MSD) {
					temp_upper += pow( SD - MSD ,2);
					n_upper ++;
				}
				else {
					temp_lower += pow( SD - MSD ,2);
					n_lower ++;
				}
            }
            upper_STD += temp_upper;
            lower_STD += temp_lower;
        }
        upper_STD /= n_upper;
        upper_STD = sqrt(upper_STD);
        lower_STD /= n_lower;
        lower_STD = sqrt(lower_STD);
        outf_error<<tau*dt<<" "<<upper_STD<<" "<<lower_STD<<endl; */

    }

	intrj.close();
    outf.close();
//    outf_error.close();
    return 0;
}

void ATOM::skip_lines(int num_lines, ifstream* data)
{
	int i;
    for (i = 0; i < num_lines; i++){
        getline(*data,str);
    }
}
