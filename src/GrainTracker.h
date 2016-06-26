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

#ifndef GRAINTRACKER_H_
#define GRAINTRACKER_H_
#include "OrientatorFileQueue.h"
#include "ContainerData.h"
#include "io/CSVTableWriter.h"
#include "GrainIDMapper.h"
#include "io/CFGEditor.h"
#include "io/CSVTableReader.h"
#include "io/CSVTableWriter.h"
//grain-grain mapping defaults
#define GG_MAX_MISORI 5 //degree
#define GG_MAX_A_DISTANCE 5 //times lattice parameter
#define GG_MAX_VOLFRAC 50.0 //percent

class GrainTracker {
public:
	GrainTracker(OrientatorFileQueue * const inQueue, const CubicLattice * inMaterial, std::string inInitGrainFileName = "", bool inEditCfgFiles = true);
	virtual ~GrainTracker();
	void run();
private:
	//!\brief Calculates and saves the orientation and center change for each grain to the initial state
	void calculateGrainDataChangeToInitial();
	//!\brief Saves data for prevData from a csv file
	bool initPrevTimeStepDataFromCSVFile(std::string inFileName);
	void addInitialGrainState(const GrainData * grain);
	//!\brief Deletes prevData and sets it to nullptr
	void resetPrevData();
	OrientatorFileQueue * queue;
	GrainIDMapper * mapping = nullptr;
	std::string initGrainFileName = "";
	ContainerData * curData = nullptr;
	ContainerData * prevData = nullptr;

	std::vector<GrainData> initialGrainStates;
	long grainNumberingBegin = 0;
	const CubicLattice * material;
	double maxCosHalfMisOri;
	double maxVolFrac;
	bool editCfgFiles = true;
};

#endif /* GRAINTRACKER_H_ */
