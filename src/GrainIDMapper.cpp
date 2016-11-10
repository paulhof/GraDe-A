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

#include "GrainIDMapper.h"
#define LENFACTOR 0.7839017541
GrainIDMapper::GrainIDMapper(const CubicLattice & material, ContainerData * inOldContainer, ContainerData * inCurContainer, long inNewGrainIdBegin){
	this->material = &material;
	prevContainer = inOldContainer;
	curContainer = inCurContainer;
	numCurGrains = curContainer->getNumberOfGrains();
	numPrevGrains = prevContainer->getNumberOfGrains();
	maxCosHalfMisOri = 1.;//cos = 1 -> no misorientation
	maxSqrDistance = 0.;
	maxVolFraction = 0.;
	//
	curMappedIds.resize(numCurGrains);
	curOldAssignedIds.resize(numCurGrains);
	prevGrainIsMapped.resize(numPrevGrains);
	newGrainId = inNewGrainIdBegin;
	newGrainIdBegin = inNewGrainIdBegin;
	for(long i = 0; i < numCurGrains; i++){
		curMappedIds[i] = NO_GRAIN;
	}
	for(long i = 0; i < numPrevGrains; i++){
		prevGrainIsMapped[i] = false;
	}
}

GrainIDMapper::~GrainIDMapper() {}

void GrainIDMapper::map() {
StopWatch mapWatch;
mapWatch.trigger();
GrainData * grain;
		for(gID iG = 0; iG < curContainer->getNumberOfGrains(); iG++){
				grain = (GrainData*) curContainer->getGrain(iG);
				curMappedIds[iG] = correspondingOldGrainId(grain);
				curOldAssignedIds[iG] = grain->getAssignedId();
				grain->setAssignedId(curMappedIds[iG]);
#ifdef DEBUGMODE
				#pragma omp critical
				{
					std::cout << "Thread "<< omp_get_thread_num() << ": For current grain " << iG << " with id " << curOldAssignedIds[iG] <<" found " << curMappedIds[iG] << std::endl;
				}
#endif
		}
mapWatch.trigger();
std::cout << "Mapping took " << mapWatch.getString() << std::endl;
}

void GrainIDMapper::init(double inMaxCosHalfMisOri,
		double inMaxVolFraction, long inNewGrainsIdBegin) {
	maxCosHalfMisOri = inMaxCosHalfMisOri;
	maxVolFraction = inMaxVolFraction;
	newGrainIdBegin = inNewGrainsIdBegin;
	newGrainId = newGrainIdBegin;
}

gID GrainIDMapper::getMappedId(long curGrainNum) const{
	return curMappedIds[curGrainNum];
}

gID GrainIDMapper::getMappedIdToOldAssignedId(gID oldAssignedId) const {
	for(int iG = 0; iG < curOldAssignedIds.size(); iG++){
		if(oldAssignedId == curOldAssignedIds[iG]){
			return getMappedId(iG);
		}
	}
	return NO_GRAIN;
}

long GrainIDMapper::getNumCurGrains() const{
	return numCurGrains;
}

long GrainIDMapper::getNumPrevGrains() const{
	return numPrevGrains;
}
void GrainIDMapper::print(CSVTableWriter * csvWriter) {
	if(csvWriter == nullptr){
		return;
	}
	csvFormat.fillTable(csvWriter,curContainer,*material);
}

void GrainIDMapper::edit(CFGEditor * editor) {
	int grainIDField = -1;
	//find out the column of the grain id
	for(int iAux = 0; iAux < editor->getNumAuxFields(); iAux++){
		if(editor->getAuxName(iAux) == GRAIN_ID_NAME){
			grainIDField = iAux;
		}
	}
	//if no column has been found, cancel
	if(grainIDField == -1) {
		return;
	}
	//now write data to file
	editor->readNextAtom();
	while(!editor->eof()){
		editor->setCurAux(grainIDField, std::to_string(getMappedIdToOldAssignedId(atol(editor->getCurAux(grainIDField).c_str()))));
		editor->readNextAtom();
	}
	editor->flush();
}

void GrainIDMapper::assign(AtomContainer* inCurContainer){
	if(inCurContainer->getNumGrains() < numCurGrains){
		return;
	}

	for(gID iG = 0; iG < numCurGrains; iG++){
		((Grain*)inCurContainer->getGrain(iG))->assign(getMappedId(iG));
	}
}

long GrainIDMapper::getNewGrainId() const {
	return newGrainId;
}

const ContainerData* GrainIDMapper::getCurContainer() const {
	return curContainer;
}

const std::vector<gID>& GrainIDMapper::getCurMappedIds() const {
	return curMappedIds;
}

const std::vector<gID>& GrainIDMapper::getCurOldAssignedIds() const {
	return curOldAssignedIds;
}

gID GrainIDMapper::correspondingOldGrainId(const GrainData * grain) {
	const GrainData * oldGrain;
	maxSqrDistance = SQR(LENFACTOR*pow(grain->getVolume(*material)/FOURTHIRDPI,ONETHIRD));
#ifdef DEBUGMODE
	std::cout << "Max Distance set to " << sqrt(maxSqrDistance) << " Angstrom " << std::endl;
#endif
	double smallestSqrDistance = maxSqrDistance;
	double curSqrDistance;
	gID closestCandidateId = NO_GRAIN;
	for (long iG = 0; iG < prevContainer->getNumberOfGrains(); iG++) {
		if(prevGrainIsMapped[iG]){
			continue;
		}
		oldGrain = prevContainer->getGrain(iG);

		//distance criterion
		curSqrDistance = curContainer->sqrDistance(grain->getCenter(),oldGrain->getCenter());
		if( curSqrDistance > maxSqrDistance){
			continue;
		}
		//orientation criterion
		if(ori::cosHalfMisOrientation(
				grain->getOrientation()->getQuaternion(),
				oldGrain->getOrientation()->getQuaternion())
				< maxCosHalfMisOri){
			continue;
		}
		//the current oldGrain is in the vicinity of grain AND has a similar orientation.
		//hence it is a possible candidate!
		//throughout this set find the one with the smallest distance to grain!
		if(curSqrDistance <= smallestSqrDistance){
			smallestSqrDistance = curSqrDistance;
			closestCandidateId = iG;
		}
	}
	//if we found a lucky one, return its assigned id
	if(closestCandidateId != NO_GRAIN){
		prevGrainIsMapped[closestCandidateId] = true;
		return prevContainer->getGrain(closestCandidateId)->getAssignedId();
	}
	//else this grain has no partner and gets a new one
	return newGrainId++;
}

GrainIDMapping::GrainIDMapping(const GrainIDMapper* mapper) {
	material = mapper->getMaterial();
	numCurGrains = mapper->getNumCurGrains();
	curContainer = *(mapper->getCurContainer());
	curMappedIds = mapper->getCurMappedIds();
	curOldAssignedIds = mapper->getCurOldAssignedIds();
}

void GrainIDMapping::print(CSVTableWriter * csvWriter) {
	if(csvWriter == nullptr){
		return;
	}
	csvFormat.fillTable(csvWriter,&curContainer,*material);
}

void GrainIDMapping::edit(CFGEditor * editor) {
	int grainIDField = -1;
	//find out the column of the grain id
	for(int iAux = 0; iAux < editor->getNumAuxFields(); iAux++){
		if(editor->getAuxName(iAux) == GRAIN_ID_NAME){
			grainIDField = iAux;
		}
	}
	//if no column has been found, cancel
	if(grainIDField == -1) {
		return;
	}
	//now write data to file
	editor->readNextAtom();
	long curOldId;
	long curMappedId;
	while(!editor->eof()){
		curOldId = atol(editor->getCurAux(grainIDField).c_str());
		curMappedId = getMappedIdToOldAssignedId(curOldId);
		editor->setCurAux(grainIDField, std::to_string(curMappedId));
		editor->readNextAtom();
	}
	editor->flush();
}

void GrainIDMapping::assign(AtomContainer* inCurContainer){
	if(inCurContainer->getNumGrains() < numCurGrains){
		return;
	}
	for(gID iG = 0; iG < numCurGrains; iG++){
		((Grain*)inCurContainer->getGrain(iG))->assign(getMappedId(iG));
	}
}

gID GrainIDMapping::getMappedId(long curGrainNum) const{
	if(curGrainNum < 0 || curGrainNum >= curMappedIds.size()){
		return -1;
	}
	return curMappedIds[curGrainNum];
}

gID GrainIDMapping::getMappedIdToOldAssignedId(gID oldAssignedId) const {
	for(int iG = 0; iG < curOldAssignedIds.size(); iG++){
		if(oldAssignedId == curOldAssignedIds[iG]){
			return getMappedId(iG);
		}
	}
	return NO_GRAIN;
}
