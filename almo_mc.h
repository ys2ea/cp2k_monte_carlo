/*This code runs monte carlo by calling cp2k externally.


*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <utility>
#include <time.h>

#ifndef ALMO_MC_H
#define ALMO_MC_H

const double PI = 3.1415926536;
const double au_2_kelvin = 315735.;

class atom {
private: 
	char type_[2];
	double x_, y_, z_;

public:
	atom(char*, double, double, double);
	void coord(double, double, double);
	char* type();
	double x();
	double y();
	double z();
};

class almo_mc {
private:
	int Natom;
	std::vector<atom> atom_list;
	double energy_;
	double T_;
	double stepsize;
	
public:
	almo_mc(double T, double ss, const char* filename);
	void propose_update(const char*);
	void eval_update(double, const char*, const char*);
	std::vector<double> move;
	void set_energy(double);
	int Nupdate();
	int Nstep;
	int Naccept;
};

atom::atom(char* t, double ix, double iy, double iz) {
	type_[0] = t[0];
	type_[1] = t[1];
	x_ = ix;
	y_ = iy;
	z_ = iz;
}

void atom::coord(double newx, double newy, double newz) {
	x_ = newx;
	y_ = newy;
	z_ = newz;
}

char* atom::type() {
	return type_;
}

double atom::x() {
	return x_;
}

double atom::y() {
	return y_;
}

double atom::z() {
	return z_;
}


almo_mc::almo_mc(double T, double ss, const char* filename) {
	std::ifstream input;
	input.open(filename);     //file should only contain data lines
	std::string line;
	Natom = 0;
	Nstep = 0;
	Naccept = 0;
	T_ = T;
	stepsize = ss;
	double cx,cy,cz;
	char ct[4];
	while(getline(input, line)) {
		sscanf(line.c_str(), "%s %lf %lf %lf\n", ct, &cx, &cy, &cz);
		atom_list.push_back(atom(ct, cx, cy, cz));
		Natom += 1;
		std::cout << "Atom #: " << atom_list.size() << " (" << atom_list[Natom-1].type() << ", " << atom_list[Natom-1].x() << ", " << atom_list[Natom-1].y() << ", " << atom_list[Natom-1].z() << ")\n";
	}
	input.close();
}

void almo_mc::propose_update(const char* filename) {
	srand(time(NULL));
	double r1, r2, r3, r4,dx,dy,dz;
	FILE* output;
	output = fopen(filename, "w");
	for(int k = 0; k < Natom; k++) {
		//generate normal distribution using Box–Muller method 
		 r1 = double(rand())/double(RAND_MAX);
		 r2 = double(rand())/double(RAND_MAX);
		 r3 = double(rand())/double(RAND_MAX);
		 r4 = double(rand())/double(RAND_MAX);
		
		dx = stepsize*sqrt(-log(r1))*cos(2*PI*r2);
		dy = stepsize*sqrt(-log(r1))*sin(2*PI*r2);
		dz = stepsize*sqrt(-log(r3))*cos(2*PI*r4);
		move.push_back(dx);
		move.push_back(dy);
		move.push_back(dz);

		//std::cout << "( " << dx << ", " << dy << ", " << dz << ")\n";
				
		//temperary solution for bug
		if(atom_list[k].type()[0]=='S')
			fprintf(output, "%c%c  %lf  %lf  %lf\n", atom_list[k].type()[0], atom_list[k].type()[1], atom_list[k].x()+dx, atom_list[k].y()+dy, atom_list[k].z()+dz);
		else
			fprintf(output, "%c  %lf  %lf  %lf\n", atom_list[k].type()[0], atom_list[k].x()+dx, atom_list[k].y()+dy, atom_list[k].z()+dz);
		
	}
	fclose(output);
	Nstep ++;
}

void almo_mc::eval_update(double newe, const char* filename, const char* efname) {
	srand(time(NULL));
	double r = double(rand())/double(RAND_MAX);
	//std::cout << "coeff: " << (newe - energy_)*au_2_kelvin/T_ << " Exp: " << exp(-(newe - energy_)*au_2_kelvin/T_) << " r: " << r << "\n";
    //create energy file
    FILE* eout;
    eout = fopen(efname, "a");
    fprintf(eout, "Step: %d\t Old E: %lf\t New E: %lf\t exp: %lf\t rand: %lf Accepted: %d\n", Nstep, energy_, newe, exp(-(newe - energy_)*au_2_kelvin/T_), r, Naccept);
	fclose(eout);

	fflush(stdout);
	if(newe - energy_ < 0. || exp(-(newe - energy_)*au_2_kelvin/T_) > r) {
			
		if(move.size()!=3*Natom) 		std::cout << "move inconsistant!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		for(int k = 0; k < Natom; k ++) {
			double ox = atom_list[k].x();
			double oy = atom_list[k].y();
			double oz = atom_list[k].z();
			atom_list[k].coord(ox+move[3*k], oy+move[3*k+1], oz+move[3*k+2]);
			//if(k ==0) std::cout << "Move of a0: " << "(" << move[3*k] << ", " << move[3*k+1] << ", " << move[3*k+2] << ")\n";
		}
	
		Naccept ++;
		energy_ = newe;
		

		//create a trajectory input
		FILE* trajout;
		trajout = fopen(filename, "a");
		
		fprintf(trajout, "%d\n i = %d \t Energy = %lf \n", Natom, Naccept, energy_);
		std::cout << "Naccept: " << Naccept << "\n";
		for(int i = 0; i < Natom; i ++) {
			if(atom_list[i].type()[0]=='S')
			fprintf(trajout, "%c%c   %lf   %lf   %lf\n", atom_list[i].type()[0], atom_list[i].type()[1], atom_list[i].x(), atom_list[i].y(), atom_list[i].z());
			else
			fprintf(trajout, "%c   %lf   %lf   %lf\n", atom_list[i].type()[0], atom_list[i].x(), atom_list[i].y(), atom_list[i].z());
		}
		fclose(trajout);
	}

	move.clear();
}

void almo_mc::set_energy(double e) {
	energy_ = e;
}

int almo_mc::Nupdate() {
	return Naccept;
}
#endif
