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

#ifndef ATOM_H_
#define ATOM_H_
#include "GradeA_Defs.h"
#include "Orientation.h"

class Atom {
public:
	Atom();
	void init(double * pos, oID oriID, gID inGrainID);
	double * getPos();
	void setOrientationId(oID oriID);
	void setGrainId(gID inGrainID);
	gID getGrainId() const ;
	oID getOrientationId() const;
	virtual ~Atom();
private:
	oID orientID = 0;
	gID grainID = 0;
	double pos[DIM];
};

#endif /* ATOM_H_ */
