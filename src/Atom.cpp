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

#include "Atom.h"

Atom::Atom() {
	double zPos[DIM] = {0.,0.,0.};
	init(zPos, NO_ORIENTATION, NO_GRAIN);
}

double * Atom::getPos(){
	return pos;
}

void Atom::setOrientationId(long inOriID){
	orientID = inOriID;
};

oID Atom::getOrientationId() const{
	return orientID;
}

void Atom::init(double * inPos, oID inOriID, gID inGrainID){
	pos[0] = inPos[0];
	pos[1] = inPos[1];
	pos[2] = inPos[2];
	orientID = inOriID;
	grainID = inGrainID;
}

gID Atom::getGrainId() const {
	return grainID;
}

void Atom::setGrainId(gID inGrainId)
{
	grainID = inGrainId;
}

Atom::~Atom() {
}

