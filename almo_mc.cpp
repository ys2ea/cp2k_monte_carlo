#include "almo_mc.h"

//const char e_command[200] = "mpirun -n 32 /home/ys2ea/dev/cp2k-mcgill/cp2k/exe/Graham-Linux-x86-64-intel-openmpi-2016/cp2k.popt";
//const char e_command[200] = "mpirun -n 32 /scratch/r/rzk/rzk/group/cp2k/cp2k-mcgill/cp2k/exe/SciNet-Linux-x86-64-intel-intelmpi/cp2k.popt";
//const char ifname[20] = "waters.inp";
//const char ofname[20] = "waters.out";
//const char trajname[20] = "pos-h2os.xyz";
//const char ename[20] = "e_waters.ener";
//const char initcoord[20] = "waters.xyz";
//const char runcoord[20] = "h2os.xyz";

const char e_command[200] = "mpirun -n 64 /home/ys2ea/dev/cp2k-mcgill/cp2k/exe/Graham-Linux-x86-64-gfortran-openmpi/cp2k.popt";
const char ifname[20] = "mc.inp";
const char ofname[20] = "mc.out";
const char trajname[20] = "pos-mc.xyz";
const char ename[20] = "e_mc.ener";
const char initcoord[20] = "coord_op.xyz";
const char runcoord[20] = "si.xyz";


//NEED TO MAKE SURE::
// the included coord file in the cp2k input matches runcoord, AND IT CAN'T HAVE BLANK LINES AT THE END
// If restarted, save the previous traj and ener files, or it will get removed
//
int main(int argc, char* argv[]) {
	if(argc!=3) {
		std::cout << "Usage: start try, start step\n";
	}

	int start_run, start_try;
	sscanf(argv[1], "%d", &start_try);
	sscanf(argv[2], "%d", &start_run);
	int nrun = 5000;
	int ntry = 0;
	almo_mc mcrun(500,0.013,initcoord);
	mcrun.Nstep = start_try;
	mcrun.Naccept = start_run;
	remove(trajname);
	remove(ename);
	char line0[200];
	sprintf(line0, "cp %s %s", initcoord, runcoord);
	system((char *)line0);
	//initial energy;
	char line1[200];
	char line2[200];
	sprintf(line1, "%s %s > %s0", e_command, ifname, ofname);
	std::cout << line1 << "\n";
	system((char *) line1);
	
	//grep energy from output
	sprintf(line2, "grep 'ENERGY| Total'  %s0 | sed 's/.* //g' > energy.dat", ofname);
	std::cout << line2 << "\n";
	system((char*) line2);

	//read in the energy
	FILE* ein;
	ein = fopen("energy.dat", "r");
	double e;
	fscanf(ein, "%lf", &e);
	mcrun.set_energy(e);
	std::cout << "Energy: " << e << "\n";

	while(mcrun.Nupdate()<nrun) {
	//while(ntry < 1) {
		mcrun.propose_update(runcoord);
		char linei[200];
		char linej[200];
		sprintf(linei, "%s %s > %s%d", e_command, ifname, ofname, ntry+1);
		//std::cout << linei << "\n";
		system((char *) linei);
		sprintf(linej, "grep 'ENERGY| Total'  %s%d | sed 's/.* //g' > energy.dat", ofname, ntry+1);
        //std::cout << linej << "\n";
        system((char *) linej);

		//char linek[20];
		//sprintf(linek, "%s%d", ofname, ntry+1);
		//remove(linek);

		FILE* iin;
		iin = fopen("energy.dat", "r");
		double e;
		fscanf(iin, "%lf", &e);
		fclose(iin);
		//std::cout << "Energy: " << e << "\n";
		mcrun.eval_update(e, trajname, ename);
		ntry ++;
	}  
	std::cout << "Total # of tries: " << ntry << "\n";
}
