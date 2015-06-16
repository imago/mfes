/** @file PDE.h
 *  @brief Input file generation for NGSolve.
 *
 *  Depending on the mode mFES is, here PDE files are generated for potential
 *  energy difference compuations or for generation of PKINT files to compute
 *  pKA values and W-matrix elements to be used by Karlsberg2 program.
 *
 *  @author Ilkay Sakalli
 */

#ifndef PDE_H
#define PDE_H

#include "Residue.h"
#include <vector>
#include "Defs.h"

/**
 * @class PDE
 *
 * @brief This class holds the information to generate PDE files which can are
 *        used to compute potential energies. 
 *
 * Potential energies are computed using PDE files. These input files are written 
 * for NETGEN to use its interface to solve the linear Poisson-Boltzmann equation
 * via NGSolve library and other libraries provided by mFES.
 *
 * @author Ilkay Sakalli
 *
  */ 

class PDE {
public:

	PDE(){
		;
	}
	void writeEnergy(string fileName, INI &ini){
		string solOrder = ini.get<string>("solver.solution_order");
		string maxsteps = ini.get<string>("solver.maxsteps");
		string eps_in = ini.get<string>("experiment.eps_in");
		string eps_out = ini.get<string>("experiment.eps_out");

		float ionc = 0;
		boost::optional<string> ionc_opt = ini.get_optional<string>("experiment.ionc");

		if(ionc_opt){
		  ionc = ini.get<float>("experiment.ionc");
		}


		ofstream energyFile;
		energyFile.open("energy.pde");

		energyFile << "mesh = protein.vol" << endl;
		energyFile << endl;
		energyFile << "shared = pointcharges" << endl;
		energyFile << "shared = energydiff" << endl;
		energyFile << "shared = writePotatAscii" << endl;
		energyFile << endl;
		energyFile << "define constant eps0 = 8.8541878e-22" << endl;
		energyFile << "define constant q0 = 1.60217646e-19" << endl;
		energyFile << "define constant ionc = " << ionc << endl;
		energyFile << endl;
		energyFile << "define constant heapsize = " << HEAPSIZE << endl;
		energyFile << endl;
		energyFile << "define coefficient epsilon_solv" << endl;
		if (ionc == 0)
		  energyFile << "(" << eps_out << "*eps0),(" << eps_in << "*eps0),(" << eps_out << "*eps0)" << endl;
		else
		  energyFile << "(" << eps_out << "*eps0),(" << eps_out << "*eps0),(" << eps_in << "*eps0),("<< eps_out <<"*eps0)" << endl;

		energyFile << endl;
		energyFile << "define coefficient epsilon_ref" << endl;
		energyFile << "(" << eps_in << "*eps0),(" << eps_in << "*eps0),("<< eps_in <<"*eps0),("<< eps_in <<"*eps0)" << endl;
		energyFile << endl;
		energyFile << "define coefficient kappa" << endl;
		energyFile << "(("<< DL <<"*eps0)*ionc),0,0,0" << endl;
		energyFile << endl;
		energyFile << "define fespace v -order="<<solOrder<<" -type=h1ho -dirichlet=[1]" << endl;
		energyFile << endl;
		energyFile << "define gridfunction u_solv -fespace=v" << endl;
		energyFile << "define gridfunction u_ref -fespace=v" << endl;
		energyFile << endl;
		energyFile << "define bilinearform a_solv -fespace=v -symmetric" << endl;
		energyFile << "laplace epsilon_solv" << endl;
		energyFile << "mass kappa" << endl;
		energyFile << endl;
		energyFile << "define bilinearform a_ref -fespace=v -symmetric" << endl;
		energyFile << "laplace epsilon_ref" << endl;
		energyFile << "mass kappa" << endl;
		energyFile << endl;
		energyFile << "define linearform f -fespace=v" << endl;
		energyFile << endl;
		energyFile << "numproc pointcharges ps1 -linearform=f -pqrfile="<<fileName<<" -interpolate" << endl;
		energyFile << endl;
		energyFile << "define preconditioner c -type=multigrid -bilinearform=a_solv -inverse=mumps" << endl;
		energyFile << "define preconditioner c0 -type=multigrid -bilinearform=a_ref -inverse=mumps" << endl;
		energyFile << endl;
		energyFile << "numproc bvp np1 -gridfunction=u_solv -bilinearform=a_solv -linearform=f -preconditioner=c  -maxsteps=" << maxsteps << endl;
		energyFile << "numproc bvp np10 -gridfunction=u_ref -bilinearform=a_ref -linearform=f -preconditioner=c0  -maxsteps=" << maxsteps << endl;
		energyFile << endl;
		energyFile << "numproc energydiff npeval -gridfunction=u_solv -gridfunction0=u_ref -pqrfile=" << fileName << endl;
		energyFile.close();
	}

	void writeCycle0(vector<Residue> titGroupList, INI &ini){
		string solOrder = ini.get<string>("solver.solution_order");
		string maxsteps = ini.get<string>("solver.maxsteps");
		string eps_in = ini.get<string>("experiment.eps_in");
		string eps_out = ini.get<string>("experiment.eps_out");

		string extend_preconditioner = " -cylce=6 -smoother=block -coarsetype=direct -coarsesmoothingsteps=0 -notest";

		float ionc = 0;
		boost::optional<string> ionc_opt = ini.get_optional<string>("experiment.ionc");

		if(ionc_opt){
		  ionc = ini.get<float>("experiment.ionc");
		}


		for (unsigned int i = 0; i < titGroupList.size(); i++){
			Residue currentTitGroup = titGroupList.at(i);
			string prefix = currentTitGroup.getIdentifier();

			//	if ( !boost::filesystem::exists( "cycle0."+prefix+".potat" ) ){
			ofstream potfile;

				//				if (i > 0){
				//  potfile.open("pka_cycle0_"+prefix+".pde", ios::in | ios::out | ios::app);
				//}
				//else {		
			if(boost::filesystem::exists("pka_cycle0_"+prefix+".pde")){
			  boost::filesystem::remove("pka_cycle0_"+prefix+".pde");
			  cout << "removing pde file: " << "pka_cycle0_"+prefix+".pde" << endl;
			}				  		
			potfile.open("pka_cycle0_"+prefix+".pde", ios::out);
				  //				}
				potfile << "###################################################" << endl;
				potfile << "#" << endl;
				potfile << "#  Initialization: " << prefix << endl;
				potfile << "#" << endl;
				potfile << "###################################################" << endl;
				potfile << endl;
				potfile << "mesh = " << prefix << ".vol" << endl;
				potfile << endl;
				potfile << "shared = pointcharges" << endl;
				potfile << "shared = energydiff" << endl;
				potfile << "shared = writePotatAscii" << endl;
				potfile << "shared = precalculation" << endl;
				potfile << endl;
				potfile << "define constant eps0 = 8.8541878e-22" << endl;
				potfile << "define constant q0 = 1.60217646e-19" << endl;
				potfile << "define constant ionc = " << ionc << endl;
				potfile << endl;
				potfile << "define constant heapsize = " << HEAPSIZE << endl;
				potfile << endl;
				potfile << "define coefficient epsilon_solv" << endl;
				if (ionc == 0)
				  potfile << "(" << eps_out << "*eps0),(" << eps_in << "*eps0),(" << eps_out << "*eps0)" << endl;
				else
				  potfile << "(" << eps_out << "*eps0),(" << eps_out << "*eps0),(" << eps_in << "*eps0),("<< eps_out <<"*eps0)" << endl;
				potfile << endl;
				potfile << "define coefficient epsilon_ref" << endl;
				potfile << "(" << eps_in << "*eps0),(" << eps_in << "*eps0),("<< eps_in <<"*eps0),("<< eps_in <<"*eps0)" << endl;
				potfile << endl;
				potfile << "define coefficient kappa" << endl;
				potfile << "((" << DL << "*eps0)*ionc),0,0,0" << endl;
				potfile << endl;


				unsigned int nrStates = currentTitGroup.getNrStates();
				string create = "";
				string stateIndex = "";
				potfile << endl;
				potfile << "define fespace v -order="<< solOrder <<" -type=h1ho -dirichlet=[1]" << endl;
				potfile << endl;
				potfile << "define gridfunction u_solv_"<<prefix<<" -fespace=v" << endl;
				potfile << "define gridfunction u_ref_"<<prefix<<" -fespace=v" << endl;
				potfile << endl;
				potfile << "define bilinearform a_solv_"<<prefix<<" -fespace=v -symmetric" << endl;
				potfile << "laplace epsilon_solv" << endl;
				potfile << "mass kappa" << endl;
				potfile << endl;
				potfile << "define bilinearform a_ref_"<<prefix<<" -fespace=v -symmetric" << endl;
				potfile << "laplace epsilon_ref" << endl;
				//				potfile << "mass kappa" << endl;

				potfile << endl;
				potfile << "define preconditioner c_solv_"<<prefix<<" -type=multigrid -bilinearform=a_solv_"<<prefix<<" -inverse=mumps" << extend_preconditioner << endl;
				potfile << "define preconditioner c_ref_"<<prefix<<" -type=multigrid -bilinearform=a_ref_"<<prefix<<" -inverse=mumps" << extend_preconditioner << endl << endl;

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
					potfile << "numproc bvp np_solv_"<< j <<"_"<<prefix<<" -gridfunction=u_solv_" << prefix <<" -bilinearform=a_solv_" << prefix <<" -linearform=f_state_solv_"<< j <<"_"<<prefix<<" -preconditioner=c_solv_"<< prefix <<" -maxsteps=" << maxsteps << endl;
					potfile << "numproc writepotat npeval_solv_"<<j<<"_"<<prefix<<" -gridfunction=u_solv_"<<prefix<<" -pqrfile=" << prefix << ".state" << j <<".pqr  -potatfile=cycle0."<<prefix<<".potat -statenr="<< stateIndex << " " << create << endl;
					potfile  << endl;
					// Debug
					//			potfile << "numproc energydiff npeval_solv_diff_"<<j<<"_"<<prefix<<" -gridfunction=u_solv_"<<prefix<<" -gridfunction0=u_ref_"<<prefix<<" -pqrfile=" << prefix << ".state" << j <<".pqr" <<endl;


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
					potfile << "numproc bvp np_ref_" << j << "_"<<prefix<<" -gridfunction=u_ref_"<<prefix<<" -bilinearform=a_ref_"<<prefix<<" -linearform=f_state_ref_" << j << "_"<<prefix<<" -preconditioner=c_ref_"<<prefix<<"  -maxsteps=" << maxsteps << endl;
					potfile << "numproc writepotat npeval_ref_"<<j<<"_"<<prefix<<" -gridfunction=u_ref_"<<prefix<<" -pqrfile=" << prefix << ".state" << j <<".pqr -potatfile=cycle0." << prefix << ".potat -statenr=" << stateIndex << endl;
					
					// Debug
					//				potfile << "numproc energydiff npeval_ref_diff_"<<j<<"_"<<prefix<<" -gridfunction=u_solv_"<<prefix<<" -gridfunction0=u_ref_"<<prefix<<" -pqrfile=" << prefix << ".state" << j <<".pqr" <<endl;
				}
				potfile << endl;
				potfile << endl;
				potfile.close();

				//			}

		}

	}

	void writeCycle1(string molecule, vector<Residue> titGroupList, INI &ini){
		string mol = "protein";
		string solOrder = ini.get<string>("solver.solution_order");
		string maxsteps = ini.get<string>("solver.maxsteps");
		string eps_in = ini.get<string>("experiment.eps_in");
		string eps_out = ini.get<string>("experiment.eps_out");

		float ionc = 0;
		boost::optional<string> ionc_opt = ini.get_optional<string>("experiment.ionc");

		if(ionc_opt){
		  ionc = ini.get<float>("experiment.ionc");
		}


		string extend_preconditioner = " -cylce=6 -smoother=block -coarsetype=direct -coarsesmoothingsteps=0 -notest";

		bool init = true;

		if(boost::filesystem::exists("pka_cycle1.pde")){
		  boost::filesystem::remove("pka_cycle1.pde");
		  cout << "removing pde file: " << "pka_cycle1.pde" << endl;
		}

		for (unsigned int i = 0; i < titGroupList.size(); i++){
			Residue currentTitGroup = titGroupList.at(i);
			string prefix = currentTitGroup.getIdentifier();

			//			if ( !boost::filesystem::exists( "cycle1."+prefix+".potat" ) ){
			ofstream potfile;	
			potfile.open("pka_cycle1.pde", ios::in | ios::out | ios::app);
							
				if (init){

					potfile << "###################################################" << endl;
					potfile << "#" << endl;
					potfile << "#  Initialization: " << mol << endl;
					potfile << "#" << endl;
					potfile << "###################################################" << endl;
					potfile << endl;
					potfile << endl;
					potfile << "mesh = "<<mol<<".vol" << endl;
					potfile << endl;

					potfile << "shared = pointcharges" << endl;
					potfile << "shared = energydiff" << endl;
					potfile << "shared = writePotatAscii" << endl;
					potfile << "shared = precalculation" << endl;
					potfile << endl;
					potfile << "define constant eps0 = 8.8541878e-22" << endl;
					potfile << "define constant q0 = 1.60217646e-19" << endl;
					potfile << "define constant ionc = " << ionc << endl;

					potfile << endl;
					potfile << "define constant heapsize = " << HEAPSIZE << endl;
					potfile << endl;
					potfile << "define coefficient epsilon_solv" << endl;
					if (ionc == 0)
					  potfile << "(" << eps_out << "*eps0),(" << eps_in << "*eps0),(" << eps_out << "*eps0)" << endl;
					else
					  potfile << "(" << eps_out << "*eps0),(" << eps_out << "*eps0),(" << eps_in << "*eps0),("<< eps_out <<"*eps0)" << endl;
					potfile << endl;
					potfile << "define coefficient epsilon_ref" << endl;
					potfile << "(" << eps_in << "*eps0),(" << eps_in << "*eps0),("<< eps_in <<"*eps0),("<< eps_in <<"*eps0)" << endl;
					potfile << endl;
					potfile << "define coefficient kappa" << endl;
					potfile << "((" <<  DL << "*eps0)*ionc),0,0,0" << endl;
					potfile << endl;

					potfile << endl;
					potfile << "define fespace v -order="<< solOrder<<" -type=h1ho -dirichlet=[1]" << endl;
					potfile << endl;
					potfile << "define gridfunction u_solv_"<<mol<<" -fespace=v" << endl;
					potfile << "define gridfunction u_ref_"<<mol<<" -fespace=v" << endl;
					potfile << endl;
					potfile << "define bilinearform a_solv_"<<mol<<" -fespace=v -symmetric" << endl;
					potfile << "laplace epsilon_solv" << endl;
					potfile << "mass kappa" << endl;

					potfile << endl;
					potfile << "define bilinearform a_ref_"<<mol<<" -fespace=v -symmetric" << endl;
					potfile << "laplace epsilon_ref" << endl;
					//					potfile << "mass kappa" << endl;

					potfile << endl;
					potfile << "define preconditioner c_solv_"<<mol<<" -type=multigrid -bilinearform=a_solv_"<<mol<<" -inverse=mumps" << extend_preconditioner << endl;
					potfile << "define preconditioner c_ref_"<<mol<<" -type=multigrid -bilinearform=a_ref_"<<mol<<" -inverse=mumps" << extend_preconditioner << endl << endl;
					potfile << "numproc precalculation pre -pqrfile=" << molecule << endl;
					potfile << endl;
					init = false;
				}

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
					potfile << "numproc bvp np_solv_"<< j <<"_"<<prefix<<" -gridfunction=u_solv_" << mol <<" -bilinearform=a_solv_" << mol <<" -linearform=f_state_solv_"<< j <<"_"<<prefix<<" -preconditioner=c_solv_"<< mol <<" -maxsteps=" << maxsteps << endl;
					potfile << "numproc writepotat npeval_solv_"<<j<<"_"<<prefix<<" -gridfunction=u_solv_"<<mol<<" -pqrfile=" << molecule <<"  -potatfile=cycle1."<<prefix<<".potat -statenr="<< stateIndex << " " << create << " -predef" << endl;

					// Debug
					// potfile << "numproc energydiff npeval_diff_"<<j<<"_"<<prefix<<" -gridfunction=u_solv_"<<mol<<" -gridfunction0=u_ref_"<<mol<<" -pqrfile="<<molecule<<endl;
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
					potfile << "numproc bvp np_ref_" << j << "_"<<prefix<<" -gridfunction=u_ref_"<<mol<<" -bilinearform=a_ref_"<<mol<<" -linearform=f_state_ref_" << j << "_"<<prefix<<" -preconditioner=c_ref_"<<mol<<"  -maxsteps=" << maxsteps << endl;
					potfile << "numproc writepotat npeval_ref_"<<j<<"_"<<prefix<<" -gridfunction=u_ref_"<<mol<<" -pqrfile=" << molecule << " -potatfile=cycle1." << prefix << ".potat -statenr=" << stateIndex << " -predef" <<  endl;

					// Debug
					//					potfile << "numproc energydiff npeval_diff_"<<j<<"_"<<prefix<<" -gridfunction=u_solv_"<<mol<<" -gridfunction0=u_ref_"<<mol<<" -pqrfile="<<molecule<<endl;

				}
				potfile << endl;
				potfile << endl;
				potfile.close();
				// }

		}

	}

};


#endif
