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
//This file has been taken from Ovito

#include "../io/CFGImporter.h"

bool isNewType(const std::string & line){
	//if a new type is introduced, the line contains only a single column
	//<number>< \t>*<\n>
	//        |<- trimPos
	size_t trimPos = line.find_first_of(" \t");
	return line.find_first_not_of(" \t",trimPos) == std::string::npos;
}

/******************************************************************************
* Checks if the given file has format that can be read by this importer.
******************************************************************************/
bool CFGImporter::checkFileFormat()
{
	// Open input file.
	TextReader stream(filename);

	// Read first line.
	stream.readLine();
	// CFG files start with the string "Number of particles".
	if(stream.lineStartsWith("Number of particles")){
		return true;
	}

	return false;
}



/******************************************************************************
* Parses the given input file and stores the data in the given container object.
******************************************************************************/
void CFGImporter::parseFile()
{
	reader = new TextReader(filename);
	header.parse(reader);
	//atomPositions = new double[header.getNumParticles()*DIM];
	//add all corresponding auxFields
	for(int i = 0; i < header.getNumAuxFields(); i++){
		data->addAtomProperty(header.getAuxField(i));
	}
	// Create particle mass and type properties.
	int currentAtomType = 0;
	double currentMass = 0;

		double trVec[3] = {
		 0.,0.,0.
		};
		//AffineTransformation H((header.transform * header.H0).transposed());
		transform = header.getTransform()->multiply(header.getH0()).transposed();
		//H.translation() = H * Vector3(-0.5f, -0.5f, -0.5f);
		transform.multiply(trVec, translate);
		//simulationCell().setMatrix(H);
		double v100[3] = {1.,0.,0.};
		double v010[3] = {0.,1.,0.};
		double v001[3] = {0.,0.,1.};
		transform.multiplyAndTranslate(v100,translate,v100);
		transform.multiplyAndTranslate(v010,translate,v010);
		transform.multiplyAndTranslate(v001,translate,v001);
		double max[3] = {v100[0], v100[1], v100[2]};
		double min[3] = {v100[0], v100[1], v100[2]};
		for (char i=0; i < 3; i++){
			if( v010[i] > max[i]) max[i] = v010[i];
			if( v001[i] > max[i]) max[i] = v001[i];
			if( v010[i] < min[i]) min[i] = v010[i];
			if( v001[i] < min[i]) min[i] = v001[i];
		}
		double size[3] ={max[0]-min[0], max[1]-min[1], max[2]-min[2]};
		double origin[3] = {0.,0.,0.};
		data->setSize(size);
		data->setOrigin(origin);
		data->generate(header.getNumParticles());
		// Read per-particle data.

	bool isFirstLine = true;
	double position [DIM];
	for(int particleIndex = 0; particleIndex < header.getNumParticles(); ) {
		if(!isFirstLine) {
			reader->readLine();
		} else {
			isFirstLine = false;
		}
		try {
			readAtom();
			particleIndex++;
		}
		catch(Exception& ex) {
			std::cerr << "Parsing error in line "  << reader->getLineNumber()  << " of CFG file."<< std::endl;
		}
		if(curAtomData.size() >= DIM){
			position[0] = atof(curAtomData[0].c_str());
			position[1] = atof(curAtomData[1].c_str());
			position[2] = atof(curAtomData[2].c_str());
			curAtomData.erase(curAtomData.begin(),curAtomData.begin()+DIM);
			transform.multiplyAndTranslate(position,translate,position);
			data->addAtom(position, curAtomData);
		}
	}
	delete reader;
}

void CFGImporter::readAtom() {
	curAtomData.clear();
	skipExtendedLines();
	const std::string & line  = reader->getLineStr();
	size_t fieldBegin = 0;
	size_t fieldEnd = line.find_first_of(" \t",fieldBegin);
	while (fieldEnd != std::string::npos){
		if(fieldBegin != fieldEnd){
			curAtomData.push_back(line.substr(fieldBegin,fieldEnd-fieldBegin));
		}
		fieldBegin = fieldEnd + 1;
		fieldEnd = line.find_first_of(" \t",fieldBegin);
	}
	//parse last field
	if(fieldBegin != line.size()){
		curAtomData.push_back(line.substr(fieldBegin));
	}
	//indicate an error if the number of entries is not valid
	if(curAtomData.size() != header.getEntryCount()){
		throw Exception("Data line in input file does not contain the right number of columns. Expected " +
				std::to_string(header.getEntryCount())+ " file columns, but found only " + std::to_string(curAtomData.size()));

	}
	curAtomNum++;
}

void CFGImporter::skipExtendedLines() {
	if(!header.isExtendedFormat()) {
		return;
	}
	if( isNewType(reader->getLineStr()) ){
		//new type means skipping two lines (mass and name info)
		reader->readLine();
		reader->readLine();
	}
}

void CFGImporter::readAtomPositions(int particleIndex, const char* s)
{
	int columnIndex = 0;
	while(columnIndex < DIM) {
		while(*s == ' ' || *s == '\t')
			++s;
		const char* token = s;
		while(*s > ' ')
			++s;
		if(s != token) {
			atomPositions[particleIndex*DIM+columnIndex] = parseField(particleIndex, columnIndex, token, s);
			columnIndex++;
		}
		if(*s == '\0') break;
		s++;
	}
	if(columnIndex < DIM){
		throw Exception("Data line in input file does not contain enough columns. Expected" + std::to_string(DIM)+ " file columns, but found only " + std::to_string(columnIndex));
	}
	transform.multiplyAndTranslate(atomPositions+particleIndex*DIM,translate,atomPositions+particleIndex*DIM);
}


double CFGImporter::parseField(int particleIndex, int columnIndex, const char* token, const char* token_end)
{
	const char *c;
	std::string fieldString;
	for(c = token; c != token_end; c++){
		fieldString += *c;
	}
	return atof(fieldString.c_str());
}
