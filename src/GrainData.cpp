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

#include "GrainData.h"

GrainData::GrainData(double* inCenter, double * inQuaternion, long inNumAtoms,	gID inAssignedId ,
		 long inNumRegularAtoms, long inNumOrphanAtoms, double inOrientationSpread, const std::vector<double> & properties) {
	center[0] = inCenter[0];
	center[1] = inCenter[1];
	center[2] = inCenter[2];
	meanOrient.initbyQuaternion(inQuaternion);
	numAtoms = inNumAtoms;
	assignedId = inAssignedId;
	numOrphanAtoms = inNumOrphanAtoms;
	numRegularAtoms = inNumRegularAtoms;
	orientationSpread = inOrientationSpread;
	meanProperties = properties;
}

GrainData::GrainData(const Grain *inGrain)
{
	center[0] = inGrain->getPosition()[0];
	center[1] = inGrain->getPosition()[1];
	center[2] = inGrain->getPosition()[2];
	meanOrient.initbyQuaternion(inGrain->getOrientation()->getQuaternion());
	numAtoms = inGrain->getNumberOfAtoms();
	numRegularAtoms = inGrain->getNumberOfRegularAtoms();
	numOrphanAtoms = inGrain->getNumberOfOrphanAtoms();
	assignedId = inGrain->getAssignedId();
	orientationSpread = inGrain->orientationSpread();
	meanProperties = inGrain->getProperties();
}

GrainData::~GrainData() {
}

const double * GrainData::getCenter() const{
	return center;
}

const Orientation* GrainData::getOrientation() const{
	return &meanOrient;
}

double GrainData::getVolume(const CubicLattice & material) const{
	return getNumberOfAtoms()*material.getVolumePerAtom();
}

double GrainData::getVolumeInLatticeUnit(const CubicLattice & material) const{
	return getNumberOfAtoms()*material.getVolumePerAtomInVolumeUnit();
}

long GrainData::getNumberOfRegularAtoms() const {
	return numRegularAtoms;
}

long GrainData::getNumberOfOrphanAtoms() const {
	return numOrphanAtoms;
}

long GrainData::getNumberOfAtoms() const{
	return numAtoms;
}

gID GrainData::getAssignedId() const{
	return assignedId;
}

void GrainData::setAssignedId(gID inAssignedId) {
	assignedId = inAssignedId;
}

double GrainData::getOrientationSpread() const {
	return orientationSpread;
}

void GrainData::setMisOrientation(const Orientation* initOrientation) {
	misOrientationToInit = ori::misOrientation(meanOrient.getQuaternion(),initOrientation->getQuaternion());
	redMisOrientationToInit = ori::cubicMisOrientation(meanOrient.getQuaternion(),initOrientation->getQuaternion());
}

void GrainData::setDistanceToInit(double inDistance) {
	distanceToInit = inDistance;
}

const std::vector<double>& GrainData::getProperties() const {
	return meanProperties;
}
