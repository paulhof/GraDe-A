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

#include "Orientation.h"

Orientation::Orientation() {
}

void Orientation::initbyQuaternion(const double * qIn){
	q[0] = qIn[0];
	q[1] = qIn[1];
	q[2] = qIn[2];
	q[3] = qIn[3];
}

void Orientation::initbyBungeEuler(const double * b){
	ori::bunge2Quaternion(b,q);
}

const double * Orientation::getQuaternion() const{
	return q;
}

void Orientation::obtainBungeAngles(double * b) const{
	ori::quaternion2Bunge(q,b);
	if (SQR(b[1]) < SQR(0.2/RADTODEG)){
		b[0]+=b[2];
		b[2]=0.;
	}
}

Orientation::~Orientation() {
}
