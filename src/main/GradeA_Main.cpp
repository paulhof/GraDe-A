/*
GraDe-A: Grain Detection Algorithm.
Copyright (C) 2016 Paul Hoffrogge

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "../ComputationManager.h"
#include "../io/GrainTimeEvolutionWriter.h"
#include "../io/CFGEditor.h"
#include "../io/FileEditor.h"
#include "../io/GrainCSVFileFormat.h"
#include "../StopWatch.h"

int writeProgramInfo(){
	std::vector<std::string> parameter {
		"\"-inputcfgfile(s) as wildcard-\"",
		"boundaryconditions "
		"latticeconstant ",
		"angularthreshold",
	};
	std::vector<std::string> optParameter{
		"startFileNum",
		"endFileNum",
		"\"restartcsvfile\"",
		"elementname",
		"printorientations (ON)"
	};
	std::string usage = "grade-A";
	for (int i = 0; i < parameter.size(); i++){
		usage += " " + parameter[i];
	}
	usage += " | optional:";
	for (int i = 0; i < optParameter.size(); i++){
		usage += " " + optParameter[i];
	}
	std::cout << "GraDe-A -- HELP: " << std::endl;
	std::cout << "Usage: " + usage << std::endl;
	std::cout << "boundaryconditions: p: periodic n: non periodic" << std::endl;
	std::cout << "latticeconstant: in Angstrom" << std::endl;
	std::cout << "angularthreshold: in degree" << std::endl;
	std::cout << "printorientations: ON: print, else: no print" << std::endl;
	std::cout << "Example: grade-A \"input*.cfg\" p 4.05 1.0" << std::endl;
	std::cout << "Example with restart-file: grade-A \"input*.cfg\" p 4.05 1.0 0 10 \"restart.csv\" " << std::endl;
	std::cout << "Example with orientation output: grade-A \"input*.cfg\" p 4.05 1.0 0 10 \"\" Al ON" << std::endl;
	//
	return -1;
}

void run(std::string& inputFileNamesWildCard,std::string restartFilename,  double a, double angularThreshold, bool periodic, std::string chemElementName, int startFileNum = 0, int endFileNum = INT_MAX, bool printOrientations = false){
	ComputationManager manager(periodic, a, angularThreshold, chemElementName, printOrientations);
	manager.run(inputFileNamesWildCard, restartFilename, startFileNum, endFileNum);
}

int main(int argc, char** argv){
	StopWatch watch;
	watch.trigger();
	std::cout << "This is GraDe-A " << version.getString() <<"\n"<< FANCYLINE << std::endl;
//----BEGIN
{
	int numParameter = argc-1;
	std::cout << "Input " << numParameter << " parameters, 4 necessary." <<  std::endl;
	if (numParameter < 4 || numParameter > 9){
		return writeProgramInfo();
	}
	//necessary parameter
	std::string inputWildCard = argv[1];
	bool isPeriodic = (0==strcmp(argv[2],"p"));
	if(!isPeriodic && (0!=strcmp(argv[2],"n"))){
		std::cerr << "Unknown boundary conditions \""<< argv[2] << "\" specified." << std::endl;
		std::cerr << "Supported boundary conditions are \"p\" (periodic) or \"n\" (non-periodic)." << std::endl;
		std::cerr << "Exited." << std::endl;
		return -1;
	}
	double latticeParameter = atof(argv[3]);
	if (latticeParameter <= 0. || latticeParameter > 100.){
		std::cerr << "Wrong value \""<< latticeParameter << " Angstrom\" given for lattice parameter."<< std::endl;
		std::cerr << "Supported range : 0-100 Angstrom " << std::endl;
		std::cerr << "Exited." << std::endl;
		return -1;
	}
	double angularThreshold = atof(argv[4]);
	if (angularThreshold <= 0. || angularThreshold > 90.){
		std::cerr << "Wrong value \""<< angularThreshold << " degree\" given for local criterion threshold."<< std::endl;
		std::cerr << "Supported range : 0-90 degree." << std::endl;
		std::cerr << "Exited." << std::endl;
		return -1;
	}
	angularThreshold /= RADTODEG;
	//optional parameter
	int startFileNum = 0;
	int endFileNum = INT_MAX;
	if(numParameter >= 5){
		startFileNum = atoi(argv[5]);
	}
	if(numParameter >= 6){
		endFileNum = atoi(argv[6]);
	}
	std::string initFileName ="";
	if(numParameter >= 7){
		initFileName = argv[7];
	}
	std::string chemElementName = "Al";
	if(numParameter >= 8){
		chemElementName = argv[8];
	}
	bool printOrientations = false;
	if(numParameter >= 9){
		printOrientations = 0 == strcmp(argv[9], "ON");
	}
	if(printOrientations){
		std::cout << "Atom-Orientation printing ON" << std::endl;
	}
	std::cout << "Grain-Volume is given in " << VOLUMEUNIT << std::endl;
	run(inputWildCard,initFileName, latticeParameter,angularThreshold, isPeriodic, chemElementName, startFileNum, endFileNum, printOrientations);
}
//----END
	watch.trigger();
	std::cout << "Computation took " << watch.getString() << std::endl;
}
