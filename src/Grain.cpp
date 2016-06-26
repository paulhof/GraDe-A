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

#include "Grain.h"
Grain::Grain() {
	center[0] = 0.;
	center[1] = 0.;
	center[2] = 0.;
}

void Grain::addToCenter(const double *vec)
{
	center[0] += vec[0];
	center[1] += vec[1];
	center[2] += vec[2];
	nCenter++;
}

void Grain::add(Atom * atom, const Orientator * orient)
{
	atoms.push_back(atom);
	meanOrient.add(orient->getOrientation(atom->getOrientationId())->getQuaternion());
}

void Grain::addOrphan(Atom *oAtom)
{
	orphanAtoms.push_back(oAtom);
}

long Grain::getNumberOfAtoms() const{
	return getNumberOfRegularAtoms() + getNumberOfOrphanAtoms();
}

long Grain::getNumberOfRegularAtoms() const {
	return atoms.size();
}

long Grain::getNumberOfOrphanAtoms() const{
	return orphanAtoms.size();
}

const Orientation * Grain::getOrientation() const{
	return &meanOrient;
}

void Grain::recalculateMeanOrientation(){
	meanOrient.refresh();
}

void Grain::assign(gID grainId)
{
	assignedId = grainId;
	assignRegularAtoms(grainId);
	assignOrphanAtoms(grainId);
}

void Grain::assignRegularAtoms(gID grainId)
{
	for(long iA = 0; iA < atoms.size(); iA++){
		atoms[iA]->setGrainId(grainId);
	}
}

void Grain::assignOrphanAtoms(gID grainId)
{
	for(long iA = 0; iA < orphanAtoms.size(); iA++){
		orphanAtoms[iA]->setGrainId(grainId);
	}
}

void Grain::unAssign(){
	assignedId = NO_GRAIN;
	for(long iA = 0; iA < atoms.size(); iA++){
			atoms[iA]->setGrainId(assignedId);
	}
}

void Grain::calcAverageCenter()
{
	center[0] /= nCenter;
	center[1] /= nCenter;
	center[2] /= nCenter;
}

const double * Grain::getPosition() const{
	return center;
}

void Grain::calcTotalOrientationSpread(const Orientator * orient){
	calcRegularOrientationSpread(orient);
	calcOrphanOrientationSpread(orient);
	long  nTotalAtoms = atoms.size() + orphanAtoms.size();
	totalOriSpread = ( atoms.size() * oriSpread + orphanAtoms.size() * orphanOriSpread ) / nTotalAtoms;
	totalCosHalfOriSpread = ( atoms.size() * cosHalfOriSpread + orphanAtoms.size() * orphanCosHalfOriSpread ) / nTotalAtoms;
}

void Grain::calcRegularOrientationSpread(const Orientator * orient)
{
	oriSpread = 0.;
	cosHalfOriSpread = 0.;
	double cosHalfMisOri;
	oID oId;
	long n = 0;
	for (long i = 0; i < atoms.size(); i++){
		oId = atoms[i]->getOrientationId();
		if (oId == NO_ORIENTATION) continue;
		cosHalfMisOri = ori::cosHalfMisOrientation(
		orient->getOrientation(oId)->getQuaternion(),
		meanOrient.getQuaternion()
		);
		n++;
		cosHalfOriSpread += cosHalfMisOri;
		oriSpread += ori::radFromCosHalf(cosHalfMisOri);
	}
	if( n > 0){
		cosHalfOriSpread /= n;
		oriSpread /= n;
	}
}

void Grain::calcOrphanOrientationSpread(const Orientator * orient){
	orphanOriSpread = 0.;
	orphanCosHalfOriSpread = 0.;
	double cosHalfMisOri;
	oID oId;
	long n = 0;
	for (long i = 0; i < orphanAtoms.size(); i++){
		oId = orphanAtoms[i]->getOrientationId();
		if (oId == NO_ORIENTATION) continue;
		cosHalfMisOri = ori::cosHalfMisOrientation(
		orient->getOrientation(oId)->getQuaternion(),
		meanOrient.getQuaternion()
		);
		n++;
		orphanCosHalfOriSpread += cosHalfMisOri;
		orphanOriSpread += ori::radFromCosHalf(cosHalfMisOri);
	}
	if( n > 0){
		orphanCosHalfOriSpread /= n;
		orphanOriSpread /= n;
	}
}

double  Grain::orientationSpread() const {
	return oriSpread;
}

double Grain::cosHalfOrientationSpread() const {
	return cosHalfOriSpread;
}

Grain::~Grain() {
	unAssign();
}

gID Grain::getAssignedId() const {
	return assignedId;
}

double Grain::getVolume(const CubicLattice& material) const{
	return getNumberOfAtoms()*material.getVolumePerAtom();
}

double Grain::getVolumeInLatticeUnit(const CubicLattice& material) const{
	return getNumberOfAtoms()*material.getVolumePerAtomInVolumeUnit();
}
