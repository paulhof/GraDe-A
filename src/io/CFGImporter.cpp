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

struct CFGHeader {

	long numParticles;
	double unitMultiplier;
	Matrix3 H0;
	Matrix3 transform;
	//FloatType rateScale;
	double rateScale;
	bool isExtendedFormat;
	bool containsVelocities;
	//QStringList auxiliaryFields;
	void parse(TextReader * stream);
};

/******************************************************************************
* Checks if the given file has format that can be read by this importer.
******************************************************************************/
bool CFGImporter::checkFileFormat()
{
	// Open input file.
	TextReader stream(filename);

	// Read first line.
	stream.readLine(20);
	// CFG files start with the string "Number of particles".
	if(stream.lineStartsWith("Number of particles"))
		return true;

	return false;
}

/******************************************************************************
* Parses the header of a CFG file.
******************************************************************************/
void CFGHeader::parse(TextReader * stream)
{
	numParticles = -1;
	unitMultiplier = 1;
	H0.setIdentity();
	transform.setIdentity();
	rateScale = 1;
	isExtendedFormat = false;
	containsVelocities = true;
	int entry_count = 0;

	while(!stream->eof()) {
		std::string line(stream->readLine());
		// Ignore comments
		size_t commentChar = line.find('#');
		if(commentChar != std::string::npos) line.resize(commentChar);

		// Skip empty lines.
		size_t trimmedLine = line.find_first_not_of(" \t\n\r");
		if(trimmedLine == std::string::npos) continue;
		if(trimmedLine != 0) line = line.substr(trimmedLine);

		size_t splitChar = line.find('=');
		if(splitChar == std::string::npos) {
			if(stream->lineStartsWith(".NO_VELOCITY.")) {
				containsVelocities = false;
				continue;
			}
			break;
		}

		std::string key = line.substr(0, line.find_last_not_of(" \t\n\r", splitChar - 1) + 1);
		size_t valuestart = line.find_first_not_of(" \t\n\r", splitChar + 1);
		if(valuestart == std::string::npos) valuestart = splitChar+1;
		std::string value = line.substr(valuestart);

		if(key == "Number of particles") {
			numParticles = atoi(value.c_str());
			if(numParticles < 0 || numParticles > 1e9)
				throw Exception("CFG file parsing error. Invalid number of atoms (line " + std::to_string(stream->getLineNumber()) + "): " + std::to_string(numParticles) + ")");
		}
		else if(key == "A") unitMultiplier = atof(value.c_str());
		else if(key == "H0(1,1)") H0.set(0,0,atof(value.c_str()) * unitMultiplier);
		else if(key == "H0(1,2)") H0.set(0,1,atof(value.c_str()) * unitMultiplier);
		else if(key == "H0(1,3)") H0.set(0,2,atof(value.c_str()) * unitMultiplier);
		else if(key == "H0(2,1)") H0.set(1,0,atof(value.c_str()) * unitMultiplier);
		else if(key == "H0(2,2)") H0.set(1,1,atof(value.c_str()) * unitMultiplier);
		else if(key == "H0(2,3)") H0.set(1,2,atof(value.c_str()) * unitMultiplier);
		else if(key == "H0(3,1)") H0.set(2,0,atof(value.c_str()) * unitMultiplier);
		else if(key == "H0(3,2)") H0.set(2,1,atof(value.c_str()) * unitMultiplier);
		else if(key == "H0(3,3)") H0.set(2,2, atof(value.c_str()) * unitMultiplier);
		else if(key == "Transform(1,1)") transform.set(0,0,atof(value.c_str()));
		else if(key == "Transform(1,2)") transform.set(0,1,atof(value.c_str()));
		else if(key == "Transform(1,3)") transform.set(0,2,atof(value.c_str()));
		else if(key == "Transform(2,1)") transform.set(1,0,atof(value.c_str()));
		else if(key == "Transform(2,2)") transform.set(1,1,atof(value.c_str()));
		else if(key == "Transform(2,3)") transform.set(1,2,atof(value.c_str()));
		else if(key == "Transform(3,1)") transform.set(2,0,atof(value.c_str()));
		else if(key == "Transform(3,2)") transform.set(2,1,atof(value.c_str()));
		else if(key == "Transform(3,3)") transform.set(2,2,atof(value.c_str()));
		else if(key == "eta(1,1)") {}
		else if(key == "eta(1,2)") {}
		else if(key == "eta(1,3)") {}
		else if(key == "eta(2,2)") {}
		else if(key == "eta(2,3)") {}
		else if(key == "eta(3,3)") {}
		else if(key == "R") rateScale = atof(value.c_str());
		else if(key == "entry_count") {
			entry_count = atoi(value.c_str());
			isExtendedFormat = true;
		}
		else if(key.compare(0, 10, "auxiliary[") == 0) {
			isExtendedFormat = true;
			size_t endOfName = value.find_first_of(" \t");
		}
		else {
			throw Exception("Unknown key in CFG file header at line " + std::to_string(stream->getLineNumber()) + ": " + line + ")");
		}
	}
	if(numParticles < 0)
		throw Exception("Invalid file header. This is not a valid CFG file.");
}

bool isNewTypeOld(const char * inLine){
	for(const char* line = inLine; *line != '\0'; ++line) {
		if(*line <= ' ') {
			//if any character is something like a ' ' or '\t'
			for(; *line != '\0'; ++line) {
				//from this character on go to the end
				if(*line > ' ') {
					//if something is different to anything like ' ', then
					return false;
					break;
				}
			}
			break;
		}
	}
	return true;
}

bool isNewTypeNew(const char * inLine){
	std::string curLine = inLine;
	size_t trimPos = curLine.find_first_of(" \t");
	if(curLine.find_first_not_of(" \t",trimPos) != std::string::npos){
		return false;
	}
	return true;
}
/******************************************************************************
* Parses the given input file and stores the data in the given container object.
******************************************************************************/
void CFGImporter::parseFile()
{
	CFGHeader header;
	reader = new TextReader(filename);
	header.parse(reader);
	atomPositions = new double[header.numParticles*DIM];

	// Create particle mass and type properties.
	int currentAtomType = 0;
	double currentMass = 0;

		double trVec[3] = {
		 0.,0.,0.
		};
		//AffineTransformation H((header.transform * header.H0).transposed());
		transform = header.transform.multiply(&header.H0).transposed();
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
		data->generate(header.numParticles);
		// Read per-particle data.
	bool isFirstLine = true;
	for(int particleIndex = 0; particleIndex < header.numParticles; ) {
		if(!isFirstLine)
			reader->readLine();
		else
			isFirstLine = false;

		if(header.isExtendedFormat) {
			bool isNewType = isNewTypeNew(reader->getLine());
			//go for each character

			std::string particleTypeName;
			if(isNewType) {
				// Parse mass and atom type name.
				currentMass = atof(reader->getLine());
				const char* line = reader->readLine();
				while(*line != '\0' && *line <= ' ') ++line;
				const char* line_end = line;
				while(*line_end != '\0' && *line_end > ' '){
					particleTypeName += *line_end;
					++line_end;
				}
				continue;
			}
		}

		try {
			readAtomPositions(particleIndex, reader->getLine());
			particleIndex++;
		}
		catch(Exception& ex) {
			std::cerr << "Parsing error in line "  << reader->getLineNumber()  << " of CFG file."<< std::endl;
		}
	}
	delete reader;
	data->addAtoms(atomPositions, header.numParticles);
	delete[] atomPositions;
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
		std::cerr << "Dataline error" << std::endl;
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
