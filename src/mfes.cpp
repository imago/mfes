/*
 ============================================================================
 Name        : mfes.cpp
 Author      : Ilkay Sakalli
 Version     :
 Copyright   : GNU Library or Lesser General Public License version 2.0 (LGPLv2)
 Description : Hello World in C++,
 ============================================================================
 */

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/timer.hpp>

#include <iostream>
#include <vector>
#include <string>

#include "lib/Defs.h"
#include "lib/PQR.h"

using namespace std;
namespace po = boost::program_options;

void calcDeltaG(vector<PQR> &pqrList, INI &ini){
	for (unsigned int i = 0; i < pqrList.size(); i++){
		PQR currentPQR = pqrList.at(i);
		currentPQR.calcModel(ini);
		currentPQR.writePDE(ini, "energy");
		currentPQR.calcDeltaG();
	}


}

void calcpKa(vector<PQR> &pqrList, boost::property_tree::ptree &ini){
	for (unsigned int i = 0; i < pqrList.size(); i++){
		PQR currentPQR = pqrList.at(i);
		currentPQR.addInfo(ini);
		if (!currentPQR.STparsed)
			currentPQR.parseSTFiles();
		currentPQR.calcResidues(ini);
		currentPQR.calcModel(ini);
		 if (currentPQR.explicitModels)
			currentPQR.calcExplicitModels(ini);
		currentPQR.writePDE(ini);
		currentPQR.calcPotat("cycle0");
		currentPQR.calcPotat("cycle1");

		currentPQR.calcPkint("cycle0");
		currentPQR.calcPkint("cycle1");

//		cout << "cycle 0 results: " << endl;
//		currentPQR.writeOutPkint("cycle0");
//		cout << endl;
//		cout << "cycle 1 results: " << endl;
//		currentPQR.writeOutPkint("cycle1");
//		cout << endl;
		cout << "cycle 1 - cycle 0 results: " << endl;
		currentPQR.writeOutPkint("diff");
		cout << endl;
	}
}


int main(int argc, char* argv[]) {

	cout << "Hello World from mFES" << endl; /* prints Hello World from mFES */
	boost::timer t;

    try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("ini", po::value<string>(), "set configuration file")
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }

        if (vm.count("ini")) {
            cout << "Configuration file is set to "
                 << vm["ini"].as<string>() << ".\n";
        } else {
            cout << "Please set configuration file!\n";
        }

        boost::property_tree::ptree ini;
        boost::property_tree::ini_parser::read_ini(vm["ini"].as<string>(), ini);
        cout << "Starting job " << ini.get<std::string>("general.jobname") << endl;
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
        if ( mode == "energy") {
        	cout << "Calculation of energy difference selected." << endl;
            calcDeltaG( pqrList, ini );
        } else if (mode == "pka") { // pKa calculation
        	cout << "Calculation pKa values selected." << endl;
            calcpKa( pqrList, ini);
        } else {
        	cout << "Please set 'mode' to 'energy' or 'pka' in your configuration!" << endl;
        }

    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }

    cout << "mFES total execution time " << t.elapsed() << " seconds." <<endl;

    return 0;
}
