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

#include "GrainCSVFileFormat.h"
using namespace GrainCSVFile;

bool sortById(const Grain * g1, const Grain * g2){
	return g1->getAssignedId() < g2->getAssignedId();
}

bool sortByGrainId(const GrainData * g1, const GrainData * g2){
	return g1->getAssignedId() < g2->getAssignedId();
}

Format::Format() {}

Format::~Format() {}

void Format::prepareTable(CSVTableWriter* csvTable) const{
	if(csvTable == nullptr){
		return;
	}
	csvTable->addHeaderLine("GraDe-A v" + version.getString() + " grain data CSV file");
	for(int iC = 0; iC < getNumColumns(); iC++){
		csvTable->addNewColumn(colNames[iC]);
	}
	for(int iCO = 0; iCO < optionalColNames.size(); iCO++){
		csvTable->addNewColumn(optionalColNames[iCO]);
	}
}

int Format::getNumColumns() const {
	return colNames.size();
}

void Format::fillTable(CSVTableWriter * csvTable, const ContainerData * container, const CubicLattice & material) const{
	if(csvTable == nullptr){
			return;
		}
		if(container == nullptr){
			return;
		}

		prepareTable(csvTable);
			std::vector<const GrainData *> sortedGrainList;
			for(gID iG = 0; iG < container->getNumberOfGrains(); iG++){
				sortedGrainList.push_back(container->getGrain(iG));
			}
			std::sort(sortedGrainList.begin(), sortedGrainList.end(),sortByGrainId);
			const GrainData * grain;
			double bungeAngles[3];
			char iEntry = 0;
			//write header line
			std::string periodicString = container->isPeriodic() ? getPeriodicFlag() : "";
			csvTable->addHeaderLine(
					ori::to_string(container->getSize()[0])		+ csvTable->getSeparator() +
					ori::to_string(container->getSize()[1])		+ csvTable->getSeparator() +
					ori::to_string(container->getSize()[2])		+ csvTable->getSeparator() +
					periodicString								+ csvTable->getSeparator() +
					ori::to_string(container->getOrigin()[0])	+ csvTable->getSeparator() +
					ori::to_string(container->getOrigin()[1])	+ csvTable->getSeparator() +
					ori::to_string(container->getOrigin()[2])
					);
			//write grains line by line
			for(gID iG = 0; iG < sortedGrainList.size(); iG++){
				grain = sortedGrainList[iG];
				iEntry = 0;
				csvTable->addNewLine();
				//ID
				csvTable->addEntryToCurrentLine(std::to_string(grain->getAssignedId()),iEntry++);
				//NumAtoms
				csvTable->addEntryToCurrentLine(std::to_string(grain->getNumberOfAtoms()),iEntry++);
				//NumRegularAtoms
				csvTable->addEntryToCurrentLine(std::to_string(grain->getNumberOfRegularAtoms()),iEntry++);
				//NumOrphanAtoms
				csvTable->addEntryToCurrentLine(std::to_string(grain->getNumberOfOrphanAtoms()),iEntry++);
				//Position data
				double reducedPoint[DIM];
				container->reducedPoint(grain->getCenter(),reducedPoint);
				csvTable->addEntryToCurrentLine(ori::to_string(reducedPoint[0]),iEntry++);
				csvTable->addEntryToCurrentLine(ori::to_string(reducedPoint[1]),iEntry++);
				csvTable->addEntryToCurrentLine(ori::to_string(reducedPoint[2]),iEntry++);
				//Orientation in euler-angles
				grain->getOrientation()->obtainBungeAngles(bungeAngles);
				csvTable->addEntryToCurrentLine(ori::to_string(bungeAngles[0]*RADTODEG),iEntry++);
				csvTable->addEntryToCurrentLine(ori::to_string(bungeAngles[1]*RADTODEG),iEntry++);
				csvTable->addEntryToCurrentLine(ori::to_string(bungeAngles[2]*RADTODEG),iEntry++);
				//Orientation spread
				csvTable->addEntryToCurrentLine(ori::to_string(grain->getOrientationSpread()*RADTODEG),iEntry++);
				//Orientation as quaternion
				csvTable->addEntryToCurrentLine(ori::to_string(grain->getOrientation()->getQuaternion()[0]),iEntry++);
				csvTable->addEntryToCurrentLine(ori::to_string(grain->getOrientation()->getQuaternion()[1]),iEntry++);
				csvTable->addEntryToCurrentLine(ori::to_string(grain->getOrientation()->getQuaternion()[2]),iEntry++);
				csvTable->addEntryToCurrentLine(ori::to_string(grain->getOrientation()->getQuaternion()[3]),iEntry++);
				//Additional Data
				csvTable->addEntryToCurrentLine(ori::to_string(grain->getVolumeInLatticeUnit(material)),iEntry++);
				csvTable->addEntryToCurrentLine(ori::to_string(grain->getMisOrientationToInit()*RADTODEG),iEntry++);
				csvTable->addEntryToCurrentLine(ori::to_string(grain->getRedMisOrientationToInit()*RADTODEG),iEntry++);
				csvTable->addEntryToCurrentLine(ori::to_string(grain->getDistanceToInit()),iEntry++);
			}
}

void Format::fillTable(CSVTableWriter * csvTable, const AtomContainer * container, const CubicLattice & material) const{
	if(csvTable == nullptr){
		return;
	}
	if(container == nullptr){
		return;
	}

	prepareTable(csvTable);
		std::vector< const Grain *> sortedGrainList;
		for(gID iG = 0; iG < container->getNumGrains(); iG++){
			sortedGrainList.push_back(container->getGrain(iG));
		}
		std::sort(sortedGrainList.begin(), sortedGrainList.end(),sortById);
		const Grain * grain;
		double bungeAngles[3];
		char iEntry = 0;
		//write header line
		std::string periodicString = container->isPeriodic() ? getPeriodicFlag() : "";
		csvTable->addHeaderLine(
				ori::to_string(container->getSize()[0])		+ csvTable->getSeparator() +
				ori::to_string(container->getSize()[1])		+ csvTable->getSeparator() +
				ori::to_string(container->getSize()[2])		+ csvTable->getSeparator() +
				periodicString								+ csvTable->getSeparator() +
				ori::to_string(container->getOrigin()[0])	+ csvTable->getSeparator() +
				ori::to_string(container->getOrigin()[1])	+ csvTable->getSeparator() +
				ori::to_string(container->getOrigin()[2])
				);
		//write grains line by line
		for(gID iG = 0; iG < sortedGrainList.size(); iG++){
			grain = sortedGrainList[iG];
			iEntry = 0;
			csvTable->addNewLine();
			//ID
			csvTable->addEntryToCurrentLine(std::to_string(grain->getAssignedId()),iEntry++);
			//NumAtoms
			csvTable->addEntryToCurrentLine(std::to_string(grain->getNumberOfAtoms()),iEntry++);
			//NumRegularAtoms
			csvTable->addEntryToCurrentLine(std::to_string(grain->getNumberOfRegularAtoms()),iEntry++);
			//NumOrphanAtoms
			csvTable->addEntryToCurrentLine(std::to_string(grain->getNumberOfOrphanAtoms()),iEntry++);
			//Position data
			double reducedPoint[DIM];
			container->reducedPoint(grain->getPosition(),reducedPoint);
			csvTable->addEntryToCurrentLine(ori::to_string(reducedPoint[0]),iEntry++);
			csvTable->addEntryToCurrentLine(ori::to_string(reducedPoint[1]),iEntry++);
			csvTable->addEntryToCurrentLine(ori::to_string(reducedPoint[2]),iEntry++);
			//csvTable->addEntryToCurrentLine(ori::to_string(grain->getPosition()[0]),iEntry++);
			//csvTable->addEntryToCurrentLine(ori::to_string(grain->getPosition()[1]),iEntry++);
			//csvTable->addEntryToCurrentLine(ori::to_string(grain->getPosition()[2]),iEntry++);
			//Orientation in euler-angles
			grain->getOrientation()->obtainBungeAngles(bungeAngles);
			csvTable->addEntryToCurrentLine(ori::to_string(bungeAngles[0]*RADTODEG),iEntry++);
			csvTable->addEntryToCurrentLine(ori::to_string(bungeAngles[1]*RADTODEG),iEntry++);
			csvTable->addEntryToCurrentLine(ori::to_string(bungeAngles[2]*RADTODEG),iEntry++);
			//Orientation spread
			csvTable->addEntryToCurrentLine(ori::to_string(grain->orientationSpread()*RADTODEG),iEntry++);
			//Orientation as quaternion
			csvTable->addEntryToCurrentLine(ori::to_string(grain->getOrientation()->getQuaternion()[0]),iEntry++);
			csvTable->addEntryToCurrentLine(ori::to_string(grain->getOrientation()->getQuaternion()[1]),iEntry++);
			csvTable->addEntryToCurrentLine(ori::to_string(grain->getOrientation()->getQuaternion()[2]),iEntry++);
			csvTable->addEntryToCurrentLine(ori::to_string(grain->getOrientation()->getQuaternion()[3]),iEntry++);
			//
			csvTable->addEntryToCurrentLine(ori::to_string(grain->getVolumeInLatticeUnit(material)),iEntry++);
			csvTable->addEntryToCurrentLine(std::to_string(0.),iEntry++);
			csvTable->addEntryToCurrentLine(std::to_string(0.),iEntry++);
			csvTable->addEntryToCurrentLine(std::to_string(0.),iEntry++);
		}
}

bool Format::isRightFormat(CSVTableReader* csvTable) const{
	//header must be existing
	if(csvTable->getNumHeaderLines() != numHeaderLines) {
#ifdef DEBUGMODE
	std::cout << "Exited because num header lines is wrong" << std::endl;
#endif
		return false;
	}
	//check size
	double size[3];
	for(char dim = 0; dim < 2; dim ++){
		//sizes must be positive numbers
		size[dim] = atof(csvTable->headerEntry(containerSizeLineNum,containerSizeColNum + dim).c_str());
		if(  size[dim] <= 0.){
#ifdef DEBUGMODE
	std::cout << "Exited because size is wrong" << std::endl;
#endif
			return false;
		}
	}
	// periodic info
	std::string periodicEntry = csvTable->headerEntry(periodicLineNum, periodicColNum);
	if (periodicEntry != getPeriodicFlag() && periodicEntry != ""){
#ifdef DEBUGMODE
	std::cout << "Exited because periodic flag is wrong" << std::endl;
#endif
		return false;
	}
	//check origin
	for(char dim = 0; dim < 2; dim ++){
		//origin not  allowed to be larger than size
		if(  fabs(atof(csvTable->headerEntry(containerOriginLineNum,containerOriginColNum + dim).c_str())) > size[dim]){
#ifdef DEBUGMODE
	std::cout << "Exited because origin is wrong" << std::endl;
#endif
			return false;
		}
	}
	//check title line
	for(int iC = 0; iC < getNumColumns(); iC++){
		if(csvTable->getColumnName(iC) != colNames[iC]){
#ifdef DEBUGMODE
	std::cout << "Exited because colname " << iC << " is wrong" << std::endl;
#endif
			return false;
		}
	}
	return true;
}

bool Format::init(CSVTableReader* csvTable, ContainerData* data) const{
	if(data == nullptr){
		return false;
	}
	if(!isRightFormat(csvTable)){
		return false;
	}
	//parse size and periodic
	double size[3];
	double origin[3];
	for (char iC = 0; iC < 3; iC ++){
		size[iC] = atof(csvTable->headerEntry(containerSizeLineNum,containerSizeColNum + iC).c_str());
	}
	for (char iC = 0; iC < 3; iC ++){
		origin[iC] = atof(csvTable->headerEntry(containerOriginLineNum,containerOriginColNum + iC).c_str());
	}
	data->setSize(size);
	data->setPeriodic(csvTable->headerEntry(periodicLineNum, periodicColNum) == periodicFlag);
	//Grain center position
	gID maxAssignedId = 0;
	gID assignedId;
	long numAtoms;
	double pos[DIM];
	double qOri[4];
	double oriSpread;
	long numOrphanAtoms;
	long numRegularAtoms;

	//parse all grains
	for (int iG = 0; iG < csvTable->getNumLines(); iG++){
		//assigned Id
		assignedId = atol(csvTable->entry(iG, getColNum(Column::GrainID)).c_str());
		if(assignedId > maxAssignedId){
			maxAssignedId = assignedId;
		}
		//volume
		numAtoms = atol(csvTable->entry(iG, getColNum(Column::NumAtoms)).c_str());
		numRegularAtoms = atol(csvTable->entry(iG, getColNum(Column::NumRegularAtoms)).c_str());
		numOrphanAtoms = atol(csvTable->entry(iG, getColNum(Column::NumOrphanAtoms)).c_str());
		//pos
		for(char iEntry = getColNum(Column::PosX); iEntry <= getColNum(Column::PosZ); iEntry++){
			pos[iEntry-getColNum(Column::PosX)] = atof(csvTable->entry(iG, iEntry).c_str());
		}
		//oriSpread
		oriSpread = atof(csvTable->entry(iG, getColNum(Column::OrientationSpread)).c_str())/RADTODEG;
		//orientation
		for(char iEntry = getColNum(Column::quaternion_0); iEntry <= getColNum(Column::quaternion_3); iEntry++){
			qOri[iEntry-getColNum(Column::quaternion_0)] = atof(csvTable->entry(iG, iEntry).c_str());
		}

		//now add a grain to prevData
		GrainData grain(pos,qOri, numAtoms, assignedId, numRegularAtoms, numOrphanAtoms,oriSpread);
		data->addGrain(grain);
	}
	//grainNumberingBegin = maxAssignedId+1;
	return true;

}

const int Format::getColNum(Column col) const {
	return static_cast<int> (col);
}
