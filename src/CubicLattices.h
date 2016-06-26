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

#ifndef SRC_CUBICLATTICES_H_
#define SRC_CUBICLATTICES_H_

class CubicLattice {
public:
	CubicLattice(double latticeParameter,char atomsPerElementCell,std::string volumeUnit = "A^3",std::string name ="Al");
	virtual ~CubicLattice();
	double getLatticeParameter() const {
		return latticeParameter;
	}

	double getVolumePerAtom() const {
		return volumePerAtom;
	}

	double getVolumePerAtomInVolumeUnit() const{
		return volumePerAtomInVolumeUnit;
	}

	const std::string& getVolumeUnit() const {
		return volumeUnit;
	}

	const std::string& getName() const {
		return name;
	}

private:
	bool fitUnit(std::string volumeUnit);
	double latticeParameter;
	double volumePerAtom;
	double volumePerAtomInVolumeUnit;
	char atomsPerElementCell;
	std::string volumeUnit;
	std::string name;
};

class FccLattice:public CubicLattice {
public:
	FccLattice(double latticeParameter,std::string volumeUnit = "A^3", std::string name = "Al");
	virtual ~FccLattice();
};

#endif /* SRC_CUBICLATTICES_H_ */
