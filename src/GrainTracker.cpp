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

#include "GrainTracker.h"
#include "StopWatch.h"
GrainTracker::GrainTracker(OrientatorFileQueue * const inQueue, const CubicLattice * inMaterial, std::string inInitGrainFileName, bool inEditCfgFiles) {
	initGrainFileName = inInitGrainFileName;
	material = inMaterial;
	maxCosHalfMisOri = ori::cosHalfFromRad(GG_MAX_MISORI/RADTODEG);
	maxVolFrac = GG_MAX_VOLFRAC/100.0;
	queue = inQueue;
	editCfgFiles = inEditCfgFiles;
}

GrainTracker::~GrainTracker() {
	delete prevData;
	delete curData;
}

void GrainTracker::run() {
	std::cout << "Serial: Starting Grain Tracking." << std::endl;
	StopWatch watch;
	watch.trigger();
	std::string initCsvFileName;
	int initFileNum = 0;
	bool prevDataInitialized = false;
	if( initGrainFileName != "") {
		std::cout << "Initializing from grain data file \"" << initGrainFileName << "\" ." << std::endl;
		initCsvFileName = initGrainFileName;
		if (initPrevTimeStepDataFromCSVFile(initGrainFileName)){
			prevDataInitialized = true;
			grainNumberingBegin = prevData->maxGrainId() + 1 ;
		} else {
			std::cerr << "Failed to initialize from file \"" << initGrainFileName << "\" will continue with standard grain numbering."<< std::endl;
		}
	}

	if (!prevDataInitialized){
		if(initPrevTimeStepDataFromCSVFile(queue->outCsvFileName(0))){
			initCsvFileName = queue->outCsvFileName(0);
			prevDataInitialized = true;
			initFileNum = 1;
			grainNumberingBegin = prevData->getNumberOfGrains();
		} else {
			std::cerr << "Failed to initialize from file \"" << queue->outCsvFileName(0) << "\" will continue with standard grain numbering."<< std::endl;
			std::cerr << "Exited " << std::endl;
			return;
		}
	}

	std::cout << "Tracking applied from file \""<< initCsvFileName <<"\" for " << queue->numFiles()-initFileNum <<" files "<< std::endl;

	std::vector<GrainIDMapping> mappings;
	std::vector<int> mappingFileNums;
	//This for loop MUST be run serially!
	for(int iF = initFileNum; iF < queue->numFiles(); iF++){
		//read in the csv file
		CSVTableReader curReader  (queue->outCsvFileName(iF),csvFormat.getNumHeaderLines());
		curReader.parse();
		if(!csvFormat.isRightFormat(&curReader)){
			std::cout << "Exiting \"" << queue->outCsvFileName(iF) << "\" has wrong format" << std::endl;
			continue;
		}
		//init curData object
		curData = new ContainerData();
		if(!csvFormat.init(&curReader, curData)){
			delete curData;
			curData = nullptr;
			continue;
		}
		mapping = new GrainIDMapper(*material,prevData, curData);
		mapping->init(maxCosHalfMisOri,maxVolFrac, grainNumberingBegin);
		mapping->map();
		calculateGrainDataChangeToInitial();
		//first construct copies of that mapping stuff in order to start a task independently
		mappings.push_back(mapping);
		mappingFileNums.push_back(iF);
		grainNumberingBegin = mapping->getNewGrainId();
		delete mapping;
		mapping = nullptr;
		//for the next timestep(file) make curData available as prevData
		delete prevData;
		prevData = curData;
		//make curData able to be deleted again
		curData = nullptr;
	}
	watch.trigger();
	int numFiles = mappings.size();
	std::cout << "Tracking computation took " << watch.getString() << std::endl;
	std::cout << "Modifying both csv and cfg files" << std::endl;

//start a parallel team
#pragma omp parallel default(none) shared(std::cout, mappings, mappingFileNums) firstprivate(numFiles)
	{
		#pragma omp for nowait
		for(int i = 0; i < numFiles; i++){
			std::cout << "Writing csv file \"" << queue->outCsvFileName(mappingFileNums[i]) << "\"" << std::endl;
			CSVTableWriter csvWriter(queue->outCsvFileName(mappingFileNums[i]));
			mappings[i].print(&csvWriter);
			csvWriter.write();
		}


		if(editCfgFiles){

			#pragma omp for nowait
			for(int i = 0; i < numFiles; i++){
				std::cout << "Editing cfg file \"" << queue->outCfgFileName(mappingFileNums[i]) << "\"" << std::endl;
				CFGEditor editor(queue->outCfgFileName(mappingFileNums[i]));
				mappings[i].edit(&editor);
				editor.close();
			}

		}
	}
}

bool GrainTracker::initPrevTimeStepDataFromCSVFile(std::string inFileName) {
	//create local reader object for parsing input file.
	CSVTableReader startTable (inFileName, csvFormat.getNumHeaderLines());
	//parse the file
	if(!startTable.parse()){
		return false;
	}
	//check whether table has right format
	if(!csvFormat.isRightFormat(&startTable)){
		std::cerr << "File \"" << inFileName << "\" has wrong csv-grain data format." << std::endl;
		return false;
	}
	//setup data object for previous timestep
	prevData = new ContainerData();
	//init prevData
	if(!csvFormat.init(&startTable,prevData)){
		resetPrevData();
		return false;
	}

	for(int iG = 0; iG < prevData->getNumberOfGrains(); iG++){
		addInitialGrainState(prevData->getGrain(iG));
	}
	return true;
}

void GrainTracker::addInitialGrainState(const GrainData* grain) {
	GrainData emptyGrain;
	gID assignedId = grain->getAssignedId();
	if( assignedId < initialGrainStates.size()){
		initialGrainStates[assignedId] = *grain;
		return;
	}
	initialGrainStates.resize(assignedId,emptyGrain);
	initialGrainStates.push_back(*grain);
}

void GrainTracker::resetPrevData() {
	delete prevData;
	prevData = nullptr;
}

void GrainTracker::calculateGrainDataChangeToInitial() {
	if(curData == nullptr){
		return;
	}
	GrainData * grain;
	gID assignedId;
	for(int iG = 0; iG < curData->getNumberOfGrains(); iG++){
		grain = (GrainData *)curData->getGrain(iG);
		assignedId = grain->getAssignedId();
		if(assignedId >= initialGrainStates.size()){
			addInitialGrainState(grain);
			continue;
		}
		grain->setMisOrientation(initialGrainStates[assignedId].getOrientation());
		grain->setDistanceToInit(sqrt(curData->sqrDistance(grain->getCenter(), initialGrainStates[assignedId].getCenter())));
	}
}
