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

#ifndef ATOMIO_H_
#define ATOMIO_H_
#define ALUMINUMELEMENTMASS 26.9815385
#include "../GradeA_Defs.h"
#include <string>
#define MAXLINELEN 200
class AtomIO {
public:
	AtomIO();
	bool openCfgFile(std::string &fileName);
	void writeCfgFileHeader(const double * boxSize, long nAtoms, long nAuxData, std::string * auxDataNames, const std::string &chemElement, double atomMass = ALUMINUMELEMENTMASS);
	void appendAtomToCfgFile(const double * pos, const double * data);
	void appendAtomToCfgFile(const double * pos, const std::vector<std::string>& properties);
	void closeCfgFile();
private:
	void writeCfgAtomPos(const double * pos, unsigned int precision = 15);
	void writeCfgAux(const std::string & aux);
	void writeCfgAux(double aux, unsigned int precision = 10);
	void writeCfgAux(long aux);
	void writeCfgLineBreak();
	bool fileOpen = false;
	char buffer[MAXLINELEN];
	double cfgBoxSize[3] = {0.,0.,0.};
	std::string fileName;
	FILE *fOutCfg = nullptr;
	long nCfgAuxData = 0;
};

#endif /* ATOMIO_H_ */
