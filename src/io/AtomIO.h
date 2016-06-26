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
	bool openCfgFile(std::string &fileName, const char mode);
	void writeCfgFileHeader(const double * boxSize, long nAtoms, long nAuxData, std::string * auxDataNames, const std::string &chemElement, double atomMass = ALUMINUMELEMENTMASS);
	void writeCfgAtomPos(double * pos);
	void writeCfgAux(long aux);
	void writeCfgLineBreak();
	void readNextIntFromCfgFile(int &value);
	void readNextDoubleFromCfgFile(double * value);
	void skipLinesCfg(long nLines);
	void appendAtomToCfgFile(double * pos, double * data);
	void closeInCfgFile();
	void closeOutCfgFile();
private:
	int readNextInt(FILE * f, int *ival);
	int readNextFloat(FILE * f, double *dval);
	void skipLines(FILE * f, long nLines);
	void readNextLine(FILE * f, char * line);
	bool fileOpen;
	char buffer[MAXLINELEN];
	double cfgBoxSize[3];
	std::string fileName;
	FILE *fInCfg;
	FILE *fOutCfg;
	FILE *fOutLammps;
	long nCfgAuxData;
};

#endif /* ATOMIO_H_ */
