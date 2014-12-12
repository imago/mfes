/*
 ============================================================================
 Name        : mfes.cpp
 Author      : Ilkay Sakalli
 Version     :
 Copyright   : GNU Library or Lesser General Public License version 2.0 (LGPLv2)
 Description : mFES - molecular Finite Element Solver
 ============================================================================
 */

#include "lib/boost/property_tree/ptree.hpp"
#include "lib/boost/property_tree/ini_parser.hpp"
#include "lib/boost/program_options.hpp"
#include "lib/boost/filesystem.hpp"
#include "lib/boost/algorithm/string.hpp"
#include "lib/boost/timer.hpp"

#include <iostream>
#include <vector>
#include <string>

#include "lib/Defs.h"
#include "lib/PQR.h"

using namespace std;
namespace po = boost::program_options;


void calcModelVolume(vector<PQR> &pqrList, INI &ini){
	for (unsigned int i = 0; i < pqrList.size(); i++){
		PQR currentPQR = pqrList.at(i);
		currentPQR.calcModel(ini);
	}
}


void calcDeltaG(vector<PQR> &pqrList, INI &ini){
  //	if ( !boost::filesystem::exists("result.out" ) ) {
	  for (unsigned int i = 0; i < pqrList.size(); i++){
	    PQR currentPQR = pqrList.at(i);
	    currentPQR.calcModel(ini);
	    currentPQR.writePDE(ini, "energy");
	    currentPQR.calcDeltaG();
	  }
  /*     	} else {
	  cout << "Calculations already done." << endl;
	}
  */


}

void calcpKa(vector<PQR> &pqrList, boost::property_tree::ptree &ini){
	string jobName = ini.get<string>("general.jobname");

	//	if ( !boost::filesystem::exists( "protein.vol" ) ) {
	  for (unsigned int i = 0; i < pqrList.size(); i++){
	    PQR currentPQR = pqrList.at(i);
	    currentPQR.addInfo(ini);
	    if (!currentPQR.STparsed)
	      currentPQR.parseSTFiles();
	    currentPQR.calcResidues(ini);
	    currentPQR.calcModel(ini);

	    //	    if (currentPQR.explicitModels)
	    currentPQR.calcExplicitModels(ini);

	    currentPQR.writePDE(ini);
	    currentPQR.calcPotat("cycle0");	      
	    currentPQR.calcPotat("cycle1");
	    
	    currentPQR.calcPkint("cycle0", jobName+".reference.pqr");
	    currentPQR.calcPkint("cycle1", jobName+".reference.pqr");
	    
	    cout << "cycle 0 results: " << endl;
	    currentPQR.writeOutPkint("cycle0", "cycle0.pkint",ini);
	    cout << endl;
	    cout << "cycle 1 results: " << endl;
	    currentPQR.writeOutPkint("cycle1", "cycle1.pkint",ini);
	    cout << endl;
	    
	    cout << "cycle 1 - cycle 0 results: " << endl;
	    currentPQR.writeOutPkint("diff", jobName+".pkint",ini);
	    currentPQR.calcW("diff");
	    currentPQR.writeOutW(jobName+".g");
	    cout << endl;
	  }
	  /*	} else {
	  cout << "Calculations already done." << endl;
	  }*/
}

void checkBoundary(INI &ini){
  string boundary = ini.get<std::string>("model.boundary");
  if ( !boost::filesystem::exists( boundary ) ) {
    std::cout << "Boundary surface given in config not found!" << std::endl;
    exit(1);
  }
}

void cleanFiles(INI &ini){
  vector<string> fileList = {"test.out", "times"};
  string surface_stl = ini.get<std::string>("model.surface_stl");
  fileList.push_back(surface_stl);
  string volume_vol = ini.get<std::string>("model.volume_vol");
  fileList.push_back(volume_vol);

  for (int i = 0; i < fileList.size(); i++){
    string file = fileList.at(i);
    if (boost::filesystem::exists(file)) {
      boost::filesystem::remove(file);
      cout << file << " removed." << endl;
    }
  }
}


int main(int argc, char* argv[]) {

	cout << "Hello World from mFES" << endl; /* prints Hello World from mFES */
	boost::timer t;

    try {

        po::options_description desc("Allowed options");
        desc.add_options()
          ("help,h", "produce help message")
	  ("ini,i", po::value<string>()->default_value("config.ini"), "set configuration file")
	  ("cleanfiles,c", "clean files before new run")
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }


        if (vm.count("ini")) {
	  if ( !boost::filesystem::exists( vm["ini"].as<string>() ) ) {
            cout << "Configuration file: " << vm["ini"].as<string>() << " could not be found!" << endl;
	    exit(1);
	  } 
	}else { 
	  cout << "Configuration file is set to "
	       << vm["ini"].as<string>() << ".\n";
        }

        boost::property_tree::ptree ini;
        boost::property_tree::ini_parser::read_ini(vm["ini"].as<string>(), ini);

	if (vm.count("cleanfiles")) {
	  cout << "Cleaning files ..." << endl;
	  cleanFiles(ini);
	  cout << "done." << endl;
	  exit(1);
	}

        cout << "Starting job \"" << ini.get<std::string>("general.jobname") << "\"" << endl;
        string location = ini.get<std::string>("general.molecule");
        string mode = ini.get<std::string>("general.mode");

        vector<PQR> pqrList;

        if ( boost::filesystem::is_directory( location.c_str() ) ){
        	boost::filesystem::path full_path( location );

        	for ( boost::filesystem::directory_iterator it( full_path );
        	        it != boost::filesystem::directory_iterator(); ++it )
        	  {
    	    	std::string ext = it->path().extension().string();
    	    	boost::to_lower(ext);

        		if ( boost::filesystem::is_regular_file( it->status() )
        			&& ( ext == ".pqr" || ext == ".pdb"))
        	    {
					PQR currentPQR(it->path().string());
					pqrList.push_back( currentPQR );
        	    }
        	  }

        } else { // is a file
        	PQR currentPQR( location );
        	pqrList.push_back( currentPQR );
        }

        cout << pqrList.size() << " molecule(s) found." << endl;
	if ( mode == "model" ) {
	  cout << "Calculation of model volume selected." << endl;
	  calcModelVolume( pqrList, ini );
	} else if ( mode == "energy") {
	  cout << "Calculation of energy difference selected." << endl;
	  checkBoundary(ini);
	  calcDeltaG( pqrList, ini );
        } else if (mode == "pka") { // pKa calculation
	  cout << "Calculation pKa values selected." << endl;
	  checkBoundary(ini);
	  calcpKa( pqrList, ini);
        } else {
	  cout << "Please set 'mode' to 'model', 'energy' or 'pka' in your configuration!" << endl;
        }

    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }

    float totalTime = t.elapsed();
    ofstream time;
    time.open ("times", ios::app );
    time << "total_time " << totalTime << " s" << endl;
    cout << "mFES total execution time " << totalTime << " seconds." <<endl;

    return 0;
}
