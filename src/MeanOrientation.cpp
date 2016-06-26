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

#include "MeanOrientation.h"

MeanOrientation::MeanOrientation() {
	init();
}

void MeanOrientation::init(){
	qSum[0] = 0.;
	qSum[1] = 0.;
	qSum[2] = 0.;
	qSum[3] = 0.;
}

void MeanOrientation::add(const double *qIn)
{
	qSum[0] += qIn[0];
	qSum[1] += qIn[1];
	qSum[2] += qIn[2];
	qSum[3] += qIn[3];
}

MeanOrientation::~MeanOrientation() {
}
void MeanOrientation::refresh(){
	double reciprocalLength = 1./sqrt(SQR(qSum[0])+SQR(qSum[1])+SQR(qSum[2])+SQR(qSum[3]));
	//double reciprocalLength = 1./atomIds.size();
	q[0] = qSum[0]*reciprocalLength;
	q[1] = qSum[1]*reciprocalLength;
	q[2] = qSum[2]*reciprocalLength;
	q[3] = qSum[3]*reciprocalLength;
}

void MeanOrientationFromAtoms::addAtomByQuaternion(AtomID& inAtomId,double * qIn){
	//For reference see: Cho, Rollett, Oh: "Determination of a Mean Orientation in Electron Backscatter Diffraction Measurements"
	//Met & Mat Transactions A Vol. 36A, 2005
	atomIds.push_back(inAtomId);
	qSum[0] += qIn[0];
	qSum[1] += qIn[1];
	qSum[2] += qIn[2];
	qSum[3] += qIn[3];
	refresh();

}

void MeanOrientationFromAtoms::addAtomsByList(AtomID * atomIdList, long numAtoms, double * qIn){
	for (long i = 0; i < numAtoms; i++){
		atomIds.push_back(atomIdList[i]);
	}
	qSum[0] += numAtoms*qIn[0];
	qSum[1] += numAtoms*qIn[1];
	qSum[2] += numAtoms*qIn[2];
	qSum[3] += numAtoms*qIn[3];
	refresh();
}

long MeanOrientationFromAtoms::getNumAssignedAtoms(){
	return atomIds.size();
}

AtomID * MeanOrientationFromAtoms::getAtomIdList(){
	return &(atomIds[0]);
}

AtomID MeanOrientationFromAtoms::getAtomId(long iAtom){
	return atomIds[iAtom];
}
