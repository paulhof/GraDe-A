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

#include "AtomIO.h"

#include <cstdio>
#include <cstring>
AtomIO::AtomIO() {
	fOutCfg = nullptr;
}

bool AtomIO::openCfgFile(std::string &inFileName){
	FILE ** f;
	char modestr[2] = {'w','\0'};
	fileName = inFileName;
	fOutCfg = fopen(fileName.c_str(),modestr);
	if (fOutCfg == nullptr){
		std::cerr << "Error opening file \"" << fileName <<"\" for writing" << std::endl;
		return false;
	}
	return true;
}


void AtomIO::writeCfgFileHeader(const double * boxSize, long nAtoms, long nAuxData, std::string * auxDataNames, const std::string &chemElement, double atomMass){
	nCfgAuxData = nAuxData;
	cfgBoxSize[0] = boxSize[0];
	cfgBoxSize[1] = boxSize[1];
	cfgBoxSize[2] = boxSize[2];
	fprintf(fOutCfg,"Number of particles = %ld\n\n", nAtoms);
	fprintf(fOutCfg,"H0(1,1) = %lf\nH0(1,2) = 0\nH0(1,3) = 0\n\n", cfgBoxSize[0]);
	fprintf(fOutCfg,"H0(2,1) = 0\nH0(2,2) = %lf\nH0(2,3) = 0\n\n", cfgBoxSize[1]);
	fprintf(fOutCfg,"H0(3,1) = 0\nH0(3,2) = 0\nH0(3,3) = %lf\n\n", cfgBoxSize[2]);
	fprintf(fOutCfg,".NO_VELOCITY.\n");
	fprintf(fOutCfg,"entry_count = %ld\n",3 + nAuxData);
	for (long i = 0; i < nAuxData; i++){
		fprintf(fOutCfg,"auxiliary[%ld] = %s\n",i , auxDataNames[i].c_str());
	}
	fprintf(fOutCfg,"%lf\n%s\n", atomMass, chemElement.c_str());
}

void AtomIO::appendAtomToCfgFile(const double * pos, const double * data){
	fprintf(fOutCfg,"%0.16le %0.16le %0.16le", pos[0]/cfgBoxSize[0], pos[1]/cfgBoxSize[1], pos[2]/cfgBoxSize[2]);
	for(long i = 0; i < nCfgAuxData; i++){
		fprintf(fOutCfg," %lf ", data[i]);
	}
	fprintf(fOutCfg,"\n");
}
void AtomIO::appendAtomToCfgFile(const double * pos, const std::vector<std::string>& properties){
	writeCfgAtomPos(pos);
	for (int i = 0; i < properties.size(); i++){
		writeCfgAux(properties[i]);
	}
	writeCfgLineBreak();
}

void AtomIO::writeCfgAtomPos(const double * pos, unsigned int precision){
	const std::string formatString ("%0." + std::to_string(precision)+ "lg");
	const std::string posFields (formatString + " " + formatString + " " + formatString);
	fprintf(fOutCfg, posFields.c_str(), pos[0]/cfgBoxSize[0], pos[1]/cfgBoxSize[1], pos[2]/cfgBoxSize[2]);
}

void AtomIO::writeCfgAux(const std::string & aux){
	fprintf(fOutCfg," %s", aux.c_str());
}

void AtomIO::writeCfgAux(long aux){
	fprintf(fOutCfg," %ld", aux);
};

void AtomIO::writeCfgAux(double aux, unsigned int precision){
	const std::string formatString = " %0." + std::to_string(precision)+ "lg";
	fprintf(fOutCfg, formatString.c_str(), aux);
};

void AtomIO::writeCfgLineBreak(){
	fprintf(fOutCfg,"\n");
}

void AtomIO::closeCfgFile(){
	fclose(fOutCfg);
}
