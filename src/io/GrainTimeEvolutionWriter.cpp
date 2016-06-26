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

#include "GrainTimeEvolutionWriter.h"

GrainTimeEvolutionWriter::GrainTimeEvolutionWriter(std::string inCsvFileNameWildCard) {
	inputFileWildCard = inCsvFileNameWildCard;
	csvQueue.initByWildcard(inputFileWildCard);
}


GrainTimeEvolutionWriter::~GrainTimeEvolutionWriter() {
	resetReader();
	resetWriters();
}

void GrainTimeEvolutionWriter::run() {
	std::string curFileId;
	csvQueue.autoFindFiles();
	if (csvQueue.numFiles() <= 0){
		return;
	}
	//init columns from first file
	initPropertiesFromFirstFile();
	//each property is written into another file
	//therefor setup a writer object for each
	initWritersForEachProperty();

	//now parse each file (= each timestep)
	for(int iF = 0 ; iF < csvQueue.numFiles(); iF++){
		//make the next file active and setup a reader for it
		reader = new CSVTableReader(csvQueue.fileName(iF),csvFormat.getNumHeaderLines());
		//try to parse it, else go to next file
		if(!reader->parse()){
			resetReader();
			continue;
		}
		//now put the data into writers
		addCurrentTimeStepDataToWriters();
		resetReader();
	}
	//now write all files
	write();
	//and delete the data
	resetWriters();
}

void GrainTimeEvolutionWriter::resetReader() {
	delete reader;
	reader = nullptr;
}

void GrainTimeEvolutionWriter::initPropertiesFromFirstFile() {
	reader = new CSVTableReader(csvQueue.curFileName(),csvFormat.getNumHeaderLines());
	if(!reader->parse()){
		resetReader();
		return;
	}
	for (int iC = 0; iC < reader->getNumColumns(); iC++){
		if(reader->getColumnName(iC) != ""){
			propertyNames.push_back(reader->getColumnName(iC));
		}
	}
	resetReader();
}

void GrainTimeEvolutionWriter::initWritersForEachProperty() {
	if (propertyNames.size() <= 0){
		return;
	}
	writers = new CSVTableWriter [propertyNames.size()];
	numWriters  = propertyNames.size();
	std::string subDir;
	subDir += ".";
	subDir += DIRCHAR;
	subDir += TIMEEVOSUBDIR;
	subDir += DIRCHAR;
	for(int iC = 0; iC < propertyNames.size(); iC++){
		writers[iC].setFileName( subDir + csvQueue.getFileNamePreFix() + "_TimeEvo_" + propertyNames[iC]+".csv");
		writers[iC].addNewColumn("TimeStep");
	}
}

void GrainTimeEvolutionWriter::addNewGrainId(std::string grainId) {
	grainIdOrder.push_back(grainId);
	for(int iW = 0; iW < numWriters; iW++){
		writers[iW].addNewColumn(grainId);
	}
}

int GrainTimeEvolutionWriter::grainIdPos(std::string grainId) {
	for (int iG = 0; iG < grainIdOrder.size(); iG++){
		if (grainId == grainIdOrder[iG]){
			return iG;
		}
	}
	return NOT_INSIDE;
}

void GrainTimeEvolutionWriter::addCurrentTimeStepDataToWriters() {
	if (reader == nullptr){
		return;
	}
	//add new line to each writer object
	//and set
	for (int iC = 0; iC < numWriters; iC ++){
		writers[iC].addNewLine();
		writers[iC].addEntryToCurrentLine(std::to_string(csvQueue.getCurFileNameNum()),0);
	}
	//parse each line (=each grain) of the input file
	long curGrainIdPos;
	for(long iL = 0; iL < reader->getNumLines(); iL ++){
		//get position of the current grain, to stick to the order
		curGrainIdPos = grainIdPos(reader->entry(iL,0));
		//if grain id is not known,...
		if( NOT_INSIDE == curGrainIdPos){
			//...add a new grain id
			addNewGrainId(reader->entry(iL,0));
			curGrainIdPos = grainIdOrder.size()-1;
		}
		//now add each property to the files
		for (int iC = 0; iC < numWriters; iC ++){
			writers[iC].addEntryToCurrentLine(reader->entry(iL,iC),curGrainIdPos+1);
		}
	}
}

void GrainTimeEvolutionWriter::write() {
	for (int iW = 0; iW < numWriters; iW ++){
		writers[iW].write();
	}
}


void GrainTimeEvolutionWriter::resetWriters() {
	delete [] writers;
	writers = nullptr;
	numWriters = 0;
}
