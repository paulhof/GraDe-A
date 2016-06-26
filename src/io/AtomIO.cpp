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
	fInCfg = nullptr;
	fOutCfg = nullptr;
	fOutLammps = nullptr;
}

bool AtomIO::openCfgFile(std::string &inFileName , const char mode){
	FILE ** f;
	char modestr[2] = {mode,'\0'};
	if (mode == 'w'){
		f = &fOutCfg;
	}
	else{
		f = &fInCfg;
	}
	fileName = inFileName;
	*f = fopen(fileName.c_str(),modestr);
	if (*f == nullptr){
		std::cerr << "Error opening file \"" << fileName <<"\" for "<< mode << "-ing" << std::endl;
		return false;
	}
	return true;
}

void AtomIO::skipLinesCfg(long nLines){
	skipLines(fInCfg, nLines);
}

void AtomIO::readNextIntFromCfgFile(int &value){
	readNextInt(fInCfg, &value);
}

void AtomIO::readNextDoubleFromCfgFile(double * value){
	readNextFloat(fInCfg, value);
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

void AtomIO::writeCfgAtomPos(double * pos){
	fprintf(fOutCfg, "%0.14lf %0.14lf %0.14lf", pos[0]/cfgBoxSize[0], pos[1]/cfgBoxSize[1], pos[2]/cfgBoxSize[2]);
}

void AtomIO::writeCfgAux(long aux){
	fprintf(fOutCfg," %ld ", aux);
};

void AtomIO::writeCfgLineBreak(){
	fprintf(fOutCfg,"\n");
}

void AtomIO::appendAtomToCfgFile(double * pos, double * data){
	fprintf(fOutCfg,"%0.16le %0.16le %0.16le", pos[0]/cfgBoxSize[0], pos[1]/cfgBoxSize[1], pos[2]/cfgBoxSize[2]);
	for(long i = 0; i < nCfgAuxData; i++){
		fprintf(fOutCfg," %lf ", data[i]);
	}
	fprintf(fOutCfg,"\n");
}

void AtomIO::closeInCfgFile(){
	fclose(fInCfg);
}

void AtomIO::closeOutCfgFile(){
	fclose(fOutCfg);
}

void AtomIO::skipLines(FILE * f, long nLines){
	for (long i = 0; i < nLines; i++){
	fgets(buffer,MAXLINELEN,f);
	}
}

int AtomIO::readNextInt(FILE * f, int *ival){
  if(!fileOpen) return EOF;
  if(EOF==fscanf(f,"%[^0-9\n]",buffer))
  {
    return EOF;
  }else {
  if(nullptr==strchr(buffer,'#')){
      if (0==fscanf(f,"%d",ival)){
	  fscanf(f," %[\n]",buffer);
	  return readNextInt(f,ival);
	  }
	  return 1;
    }else{
      fgets(buffer,MAXLINELEN,f);
      return readNextInt(f,ival);
    }
  }
}

int AtomIO::readNextFloat(FILE * f, double *dval){
  if(!fileOpen) return EOF;
  fflush(stdin);
  if(EOF==fscanf(f,"%[^0-9\n-]",buffer))
  {
    return EOF;
  }else {
  if(nullptr==strchr(buffer,'#')){
      if (0==fscanf(f,"%lf",dval)){
	  fscanf(f," %[\n]",buffer);
	  return readNextFloat(f,dval);
	  }
	  return 1;
    }else{
      fgets(buffer,MAXLINELEN,f);
      return readNextFloat(f,dval);
    }
  }
}

