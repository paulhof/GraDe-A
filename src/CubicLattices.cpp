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
#include "GradeA_Defs.h"
#include "CubicLattices.h"

CubicLattice::CubicLattice(double latticeParameter, char atomsPerElementCell,std::string volumeUnit, std::string name) {
	this->latticeParameter = latticeParameter;
	this->atomsPerElementCell = atomsPerElementCell;
	this->volumeUnit = volumeUnit;
	this->name = name;
	volumePerAtom = CUBE(latticeParameter)/atomsPerElementCell;
	volumePerAtomInVolumeUnit = volumePerAtom;
	if(!fitUnit(this->volumeUnit)){
		std::cerr << "Wrong Unit given for CubicLattice -- Volume is given in A^3" << std::endl;
	}
}

FccLattice::FccLattice(double latticeParameter, std::string volumeUnit, std::string name):CubicLattice(latticeParameter,4,volumeUnit,name) {}

FccLattice::~FccLattice() {}

bool CubicLattice::fitUnit(std::string volumeUnit) {
	//Unit String should be of the from "[number]unit^3" [...]=optional
	//Supported units are:
	//"A": Angstrom
	//"nm": Nanometer
	//"um": Micrometer
	//"mm": Millimeter
	//"m": Meter
	std::size_t powPos = volumeUnit.rfind("^3");
	if(powPos == std::string::npos){
			return false;
	}
	std::string lengthUnit = volumeUnit.substr(0,powPos);
	std::size_t endNumberPos = lengthUnit.find_last_of("0123456789.");
	std::string lengthUnitName = lengthUnit.substr(endNumberPos+1,std::string::npos);
	std::string number = lengthUnit.substr(0,endNumberPos+1);
	double factor = 1.;
	if (number.size() > 0){
		factor = atof(number.c_str());
	}
	double volumeFactor;
	double A = 1e-10;
	if (lengthUnitName == "A"){	volumeFactor = 1.; }
	else if (lengthUnitName == "nm"){ volumeFactor = CUBE(A/1e-9); }
	else if (lengthUnitName == "um"){ volumeFactor = CUBE(A/1e-6); }
	else if (lengthUnitName == "mm"){ volumeFactor = CUBE(A/1e-3); }
	else if (lengthUnitName == "m"){ volumeFactor = CUBE(A/1.); }
	else{ return false; }
	volumePerAtomInVolumeUnit = volumeFactor/factor * volumePerAtom;
	return true;
}
CubicLattice::~CubicLattice() {
}

