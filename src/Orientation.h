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

#ifndef ORIENTATION_H_
#define ORIENTATION_H_
#include "GradeA_Defs.h"
//Structure to store a single orientation
//Only stores quaternion data, because internally all calculations are applied using unit quaternions
//additionally supports the Bunge-Euler representation by the definition of convertion routines
class Orientation {
public:
	Orientation();
	virtual ~Orientation();
	void initbyQuaternion(const double * qIn);
	void initbyBungeEuler(const double * b);
	const double * getQuaternion() const;
	void obtainBungeAngles(double * b) const;
protected:
	double q[4];
};

#endif /* ORIENTATION_H_ */
