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

#include "ContainerData.h"

ContainerData::ContainerData(AtomContainer * con) {
	size[0] = con->getSize()[0];
	size[1] = con->getSize()[1];
	size[2] = con->getSize()[2];
	origin[0] = con->getOrigin()[0];
	origin[1] = con->getOrigin()[1];
	origin[2] = con->getOrigin()[2];
	for(long iG = 0; iG < con->getNumGrains(); iG++){
		grains.push_back(con->getGrain(iG));
	}
	periodic = con->isPeriodic();
	con->getAtomPropertyNames(propertyNames, false);
}

ContainerData::~ContainerData() {
}

void ContainerData::addGrain(GrainData& grain) {
	grains.push_back(grain);
}

double ContainerData::sqrDistance(const double* p1, const double* p2) const{
	if (!periodic) {
			return SQR(p1[0] - p2[0]) + SQR(p1[1] - p2[1]) + SQR(p1[2] - p2[2]);
	}
	//points in periodic images are moved inside the container
	//with periodic boundary conditions, reduced coordinates are inside the container (in [0,1]*size+origin).
	double redPos1[DIM];
	double redPos2[DIM];
	reducedPoint(p1, redPos1);
	reducedPoint(p2, redPos2);
	//now check all possible distances
	double halfSize[DIM] = {
		.5 * size[0],
		.5 * size[1],
		.5 * size[2]
	};
	double distance[DIM];
	// x-direction:
	distance[0] = fabs(redPos2[0] - redPos1[0]);
	if(distance[0] > halfSize[0]){
		distance[0] = size[0] - distance[0];
	}
	// y-direction
	distance[1] = fabs(redPos2[1] - redPos1[1]);
	if(distance[1] > halfSize[1]){
		distance[1] = size[1] - distance[1];
	}
	// z-direction
	distance[2] = fabs(redPos2[2] - redPos1[2]);
	if(distance[2] > halfSize[2]){
		distance[2] = size[2] - distance[2];
	}
	return SQR(distance[0]) + SQR(distance[1]) + SQR(distance[2]);
}

long ContainerData::getNumberOfGrains() const{
	return grains.size();
}

const GrainData* ContainerData::getGrain(gID grainId) const {
	if(grainId < grains.size()) {
		return grains.data() + grainId;
	}
	return nullptr;
}

void ContainerData::reducedPoint(const double* p, double* redP) const{
	redP[0] = reducedCoordinate(p[0],0);
	redP[1] = reducedCoordinate(p[1],1);
	redP[2] = reducedCoordinate(p[2],2);
}

double ContainerData::reducedCoordinate(double pos, unsigned char dimension) const{
	if (!periodic || dimension > 2){
		return pos;
	}
	double relPos = pos - origin[dimension];
	//rest = position in [-1,1]*size
	double rest = fmod(relPos, size[dimension]);
	//result is position in [0,1]*size + origin
	if (rest >= 0.){
		return origin[dimension] + rest;
	}
	return origin[dimension] + rest + size[dimension];
}

void ContainerData::setSize(double* inSize) {
	size[0] = inSize[0];
	size[1] = inSize[1];
	size[2] = inSize[2];
}

void ContainerData::setPeriodic(bool inPeriodic) {
	periodic = inPeriodic;
}

void ContainerData::setOrigin(double* inOrigin) {
	origin[0] = inOrigin[0];
	origin[1] = inOrigin[1];
	origin[2] = inOrigin[2];
}

const double* ContainerData::getSize() const {
	return size;
}

const double* ContainerData::getOrigin() const {
	return origin;
}

gID ContainerData::maxGrainId() const {
	gID maxGrainId = 0;
	for(gID iG = 0; iG < getNumberOfGrains(); iG++){
		if (grains[iG].getAssignedId() > maxGrainId){
			maxGrainId = grains[iG].getAssignedId();
		}
	}
	return maxGrainId;
}

void ContainerData::getAtomPropertyNames(std::vector<std::string> & outProperties) const{
 outProperties = propertyNames;
}

void ContainerData::setAtomPropertyNames(const std::vector<std::string>& inProperties) {
 propertyNames = inProperties;
}
