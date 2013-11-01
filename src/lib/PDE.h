#ifndef PDE_H
#define PDE_H

#include "Residue.h"
#include <vector>

class PDE {
public:
	PDE(){
		;
	}

	void writeCycle0(vector<Residue> titGroupList, INI &ini){
		ofstream potfile;
		potfile.open("pka_cycle0.pde");
		potfile << "shared = /home/parallels/git/mfes/nglib/pointcharges" << endl;
		potfile << "shared = /home/parallels/git/mfes/nglib/energydiff" << endl;
		potfile << "shared = /home/parallels/git/mfes/nglib/writePotatAscii" << endl;
		potfile << endl;
		potfile << "define constant eps0 = 8.8541878e-22" << endl;
		potfile << "define constant q0 = 1.60217646e-19" << endl;
		potfile << endl;
		potfile << "define coefficient epsilon_solv" << endl;
		potfile << "(80*eps0),(4*eps0)" << endl;
		potfile << endl;
		potfile << "define coefficient epsilon_ref" << endl;
		potfile << "(4*eps0),(4*eps0)" << endl;
		potfile << endl;



		for (unsigned int i = 0; i < titGroupList.size(); i++){
			Residue currentTitGroup = titGroupList.at(i);
			string prefix = currentTitGroup.getIdentifier();
			unsigned int nrStates = currentTitGroup.getNrStates();
			string create = "";
			string stateIndex = "";

			potfile << "###################################################" << endl;
			potfile << "#" << endl;
			potfile << "#  Initialization: " << prefix << endl;
			potfile << "#" << endl;
			potfile << "###################################################" << endl;
			potfile << endl;
			potfile << endl;
			potfile << "mesh = " << prefix << ".vol" << endl;
			potfile << endl;
			potfile << "define fespace v -order=2 -type=h1ho -dirichlet=[1]" << endl;
			potfile << endl;
			potfile << "define gridfunction u_solv_"<<prefix<<" -fespace=v" << endl;
			potfile << "define gridfunction u_ref_"<<prefix<<" -fespace=v" << endl;
			potfile << endl;
			potfile << "define bilinearform a_solv_"<<prefix<<" -fespace=v -symmetric" << endl;
			potfile << "laplace epsilon_solv" << endl;
			potfile << endl;
			potfile << "define bilinearform a_ref_"<<prefix<<" -fespace=v -symmetric" << endl;
			potfile << "laplace epsilon_ref" << endl;
			potfile << endl;
			potfile << "define preconditioner c_solv_"<<prefix<<" -type=multigrid -bilinearform=a_solv_"<<prefix<<" -inverse=mumps" << endl;
			potfile << "define preconditioner c_ref_"<<prefix<<" -type=multigrid -bilinearform=a_ref_"<<prefix<<" -inverse=mumps" << endl;
			potfile << endl;
			potfile << "###################################################" << endl;
			potfile << "#" << endl;
			potfile << "#  Solvation calculations: " << prefix << endl;
			potfile << "#" << endl;
			potfile << "###################################################" << endl;
			potfile << "" << endl;
			potfile << "" << endl;
			for (unsigned int j = 1; j <= nrStates; j++){
				if (j == 1){
					create = "-file=create";
					stringstream ss;
					ss << nrStates;
					stateIndex = ss.str();
				} else {
					create = "";
					stateIndex = "-1";
				}
				potfile << "# state " << j << endl;
				potfile << "define linearform f_state_solv_" << j << "_"<<prefix<<" -fespace=v" << endl;
				potfile << "numproc pointcharges ps1_solv_"<<prefix<<" -linearform=f_state_solv_" << j << "_"<< prefix << " -pqrfile=" << prefix <<".state" << j <<".pqr   -interpolate " << endl;
				potfile << "numproc bvp np_solv_"<< j <<"_"<<prefix<<" -gridfunction=u_solv_" << prefix <<" -bilinearform=a_solv_" << prefix <<" -linearform=f_state_solv_"<< j <<"_"<<prefix<<" -preconditioner=c_solv_"<< prefix <<" -maxsteps=2" << endl;
				potfile << "numproc writepotat npeval_solv_"<<j<<"_"<<prefix<<" -gridfunction=u_solv_"<<prefix<<" -pqrfile=" << prefix << ".state" << j <<".pqr  -potatfile=cycle0."<<prefix<<".potat -statenr="<< stateIndex << " " << create << endl;
				potfile  << endl;

			}
			potfile << endl;
			potfile << "###################################################" << endl;
			potfile << "#" << endl;
			potfile << "#  Reference calculation: " << prefix << endl;
			potfile << "#" << endl;
			potfile << "###################################################" << endl;
			potfile << "" << endl;
			for (unsigned int j = 1; j <= nrStates; j++){
				if (j == 1){
					stringstream ss;
					ss << nrStates;
					stateIndex = ss.str();
				} else {
					create = "";
					stateIndex = "-1";
				}
				potfile << "# state " << j << endl;
				potfile << "define linearform f_state_ref_" << j <<"_"<<prefix<<" -fespace=v" << endl;
				potfile << "numproc pointcharges ps_ref_" << j <<"_"<<prefix<<" -linearform=f_state_ref_" << j << "_"<<prefix<<" -pqrfile=" << prefix <<".state" << j << ".pqr -interpolate" << endl;
				potfile << "numproc bvp np_ref_" << j << "_"<<prefix<<" -gridfunction=u_ref_"<<prefix<<" -bilinearform=a_ref_"<<prefix<<" -linearform=f_state_ref_" << j << "_"<<prefix<<" -preconditioner=c_ref_"<<prefix<<"  -maxsteps=2" << endl;
				potfile << "numproc writepotat npeval_ref_"<<j<<"_"<<prefix<<" -gridfunction=u_ref_"<<prefix<<" -pqrfile=" << prefix << ".state" << j <<".pqr -potatfile=cycle0." << prefix << ".potat -statenr=" << stateIndex << endl;

			}
			potfile << endl;
			potfile << endl;

		}
		potfile.close();

	}

	void writeCycle1(string molecule, vector<Residue> titGroupList, INI &ini){
		string mol = "protein";
		ofstream potfile;
		potfile.open("pka_cycle1.pde");
		potfile << "shared = /home/parallels/git/mfes/nglib/pointcharges" << endl;
		potfile << "shared = /home/parallels/git/mfes/nglib/energydiff" << endl;
		potfile << "shared = /home/parallels/git/mfes/nglib/writePotatAscii" << endl;
		potfile << endl;
		potfile << "define constant eps0 = 8.8541878e-22" << endl;
		potfile << "define constant q0 = 1.60217646e-19" << endl;
		potfile << endl;
		potfile << "define coefficient epsilon_solv" << endl;
		potfile << "(80*eps0),(4*eps0)" << endl;
		potfile << endl;
		potfile << "define coefficient epsilon_ref" << endl;
		potfile << "(4*eps0),(4*eps0)" << endl;
		potfile << endl;
		potfile << "###################################################" << endl;
		potfile << "#" << endl;
		potfile << "#  Initialization: " << mol << endl;
		potfile << "#" << endl;
		potfile << "###################################################" << endl;
		potfile << endl;
		potfile << endl;
		potfile << "mesh = "<<mol<<".vol" << endl;
		potfile << endl;
		potfile << "define fespace v -order=2 -type=h1ho -dirichlet=[1]" << endl;
		potfile << endl;
		potfile << "define gridfunction u_solv_"<<mol<<" -fespace=v" << endl;
		potfile << "define gridfunction u_ref_"<<mol<<" -fespace=v" << endl;
		potfile << endl;
		potfile << "define bilinearform a_solv_"<<mol<<" -fespace=v -symmetric" << endl;
		potfile << "laplace epsilon_solv" << endl;
		potfile << endl;
		potfile << "define bilinearform a_ref_"<<mol<<" -fespace=v -symmetric" << endl;
		potfile << "laplace epsilon_ref" << endl;
		potfile << endl;
		potfile << "define preconditioner c_solv_"<<mol<<" -type=multigrid -bilinearform=a_solv_"<<mol<<" -inverse=mumps" << endl;
		potfile << "define preconditioner c_ref_"<<mol<<" -type=multigrid -bilinearform=a_ref_"<<mol<<" -inverse=mumps" << endl;
		potfile << endl;



		for (unsigned int i = 0; i < titGroupList.size(); i++){
			Residue currentTitGroup = titGroupList.at(i);
			string prefix = currentTitGroup.getIdentifier();
			unsigned int nrStates = currentTitGroup.getNrStates();
			string create = "";
			string stateIndex = "";

			potfile << "###################################################" << endl;
			potfile << "#" << endl;
			potfile << "#  Solvation calculations: " << prefix << endl;
			potfile << "#" << endl;
			potfile << "###################################################" << endl;
			potfile << "" << endl;
			potfile << "" << endl;
			for (unsigned int j = 1; j <= nrStates; j++){
				if (j == 1){
					create = "-file=create";
					stringstream ss;
					ss << nrStates;
					stateIndex = ss.str();
				} else {
					create = "";
					stateIndex = "-1";
				}
				potfile << "# state " << j << endl;
				potfile << "define linearform f_state_solv_" << j << "_"<<prefix<<" -fespace=v" << endl;
				potfile << "numproc pointcharges ps1_solv_"<<prefix<<" -linearform=f_state_solv_" << j << "_"<< prefix << " -pqrfile=" << prefix <<".state" << j <<".pqr   -interpolate " << endl;
				potfile << "numproc bvp np_solv_"<< j <<"_"<<prefix<<" -gridfunction=u_solv_" << mol <<" -bilinearform=a_solv_" << mol <<" -linearform=f_state_solv_"<< j <<"_"<<prefix<<" -preconditioner=c_solv_"<< mol <<" -maxsteps=2" << endl;
				potfile << "numproc writepotat npeval_solv_"<<j<<"_"<<prefix<<" -gridfunction=u_solv_"<<mol<<" -pqrfile=" << molecule <<"  -potatfile=cycle1."<<prefix<<".potat -statenr="<< stateIndex << " " << create << endl;
				potfile  << endl;

			}
			potfile << endl;
			potfile << "###################################################" << endl;
			potfile << "#" << endl;
			potfile << "#  Reference calculation: " << prefix << endl;
			potfile << "#" << endl;
			potfile << "###################################################" << endl;
			potfile << "" << endl;
			for (unsigned int j = 1; j <= nrStates; j++){
				if (j == 1){
					stringstream ss;
					ss << nrStates;
					stateIndex = ss.str();
				} else {
					create = "";
					stateIndex = "-1";
				}
				potfile << "# state " << j << endl;
				potfile << "define linearform f_state_ref_" << j <<"_"<<prefix<<" -fespace=v" << endl;
				potfile << "numproc pointcharges ps_ref_" << j <<"_"<<prefix<<" -linearform=f_state_ref_" << j << "_"<<prefix<<" -pqrfile=" << prefix <<".state" << j << ".pqr -interpolate" << endl;
				potfile << "numproc bvp np_ref_" << j << "_"<<prefix<<" -gridfunction=u_ref_"<<mol<<" -bilinearform=a_ref_"<<mol<<" -linearform=f_state_ref_" << j << "_"<<prefix<<" -preconditioner=c_ref_"<<mol<<"  -maxsteps=2" << endl;
				potfile << "numproc writepotat npeval_ref_"<<j<<"_"<<prefix<<" -gridfunction=u_ref_"<<mol<<" -pqrfile=" << molecule << " -potatfile=cycle1." << prefix << ".potat -statenr=" << stateIndex << endl;

			}
			potfile << endl;
			potfile << endl;

		}
		potfile.close();

	}

};


#endif
