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

#ifndef COMPUTATIONMANAGER_H_
#define COMPUTATIONMANAGER_H_
#include "StopWatch.h"
#include "OrientatorFileQueue.h"
#include "ContainerData.h"
#include "AtomContainer.h"
#include "ContainerData.h"
#include "GradeA_Defs.h"
#include "GradeA_Version.h"
#include "GrainTracker.h"
#include "io/AtomIO.h"
#include "io/CFGImporter.h"
#include "io/GrainTimeEvolutionWriter.h"

#define PERIODIC_STRING "p"

#define MAXFRAGMENT 80


//! Class, which organizes a whole GraDe-A-computation.
class ComputationManager {
public:
	ComputationManager(bool inPeriodic = true, double inLatticeParameter = 4.05, double inAngularThreshold = 0.5/RADTODEG, std::string inChemElemName = "Al", bool inPrintOrientations = false);
	virtual ~ComputationManager();
	//! Executes a computation of multiple files identified by a wildcard-string.
	void run(std::string fileNameWildCard, std::string inInitGrainFileName = "", int startFileNum = 0, int endFileNum = INT_MAX);
	void runFile(int fileNum, OrientatorFileQueue & inQueue);
private:
	void parallelRun();
	bool initPrevTimeStepDataFromCSVFile(std::string initGrainFileName);
	//! Method, which runs a computation for a single file.
	//!\param[in] fileNum file-identifier for the underlaying filequeue.
	void runSingleFile(int fileNum);
	void initContainer();
	void writeCfgFile(std::string fileName);
	void writeCsvTableFile(std::string fileName);
	void resetPrevData();
	//IO
	OrientatorFileQueue queue;
	CSVTableWriter * csvTable = nullptr;
	//

	//Orientation calculation data
	AtomContainer * container = nullptr;
	double boxSize = 0.;
	bool periodic = true;
	bool printOrientations = false;
	//material
	const FccLattice * material = nullptr;
	//! nearest-neighbor search radii squared
	double NN_searchRadiusSqrMin, NN_searchRadiusSqrMax;
	//! threshold for grain finder
	double grainAngularThreshold;
	//
	Orientator orient;
	std::string initGrainFileName;
	//
};

#endif /* COMPUTATIONMANAGER_H_ */
