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

#include "CFGHeaderData.h"

CFGHeaderData::CFGHeaderData() {
	H0.setIdentity();
	transform.setIdentity();
	rateScale = 1.;
	extendedFormat = false;
	unitMultiplier = 1;
	//initialize with containsVelocities = true
	containsVelocities = true;
	auxFields.push_back("v_x");
	auxFields.push_back("v_y");
	auxFields.push_back("v_z");
}

void CFGHeaderData::parse(FileEditor* editor) {
	numParticles = -1;	//indicates an error if value is
						//not overwritten during header-parsing
	bool isDone = false;
	std::string line;
	while(!editor->eof() && !isDone ) {
		line = editor->readNextLine();
		parseLine(line, isDone);
	}
	if(numParticles < 0)
		throw Exception("Invalid file header. This is not a valid CFG file.");
}

void CFGHeaderData::parse(TextReader* reader) {
	numParticles = -1;	//indicates an error if value is
						//not overwritten during header-parsing
	bool isDone = false;
	std::string line;
	while(!reader->eof() && !isDone) {
		line = reader->readLine();
		parseLine(line, isDone);
	}
	if(numParticles < 0)
		throw Exception("Invalid file header. This is not a valid CFG file.");
}

int CFGHeaderData::getNumAuxFields() const{
	return auxFields.size();
}

bool CFGHeaderData::containsVelFields() const{
	return containsVelocities;
}

std::string CFGHeaderData::getAuxField(int fieldNum) const{
	if (fieldNum < 0 || fieldNum >= auxFields.size()){
		return "";
	}
	return auxFields[fieldNum];
}

void CFGHeaderData::parseLine(std::string & line, bool& isDone) {
	isDone = false;
	numLines++;
	// Ignore comments
	size_t commentChar = line.find('#');
	if(commentChar != std::string::npos) line.resize(commentChar);

	// Skip empty lines.
	size_t trimmedLine = line.find_first_not_of(" \t\n\r");
	if(trimmedLine == std::string::npos) return;
	if(trimmedLine != 0) line = line.substr(trimmedLine);

	size_t splitChar = line.find('=');
	if(splitChar == std::string::npos) {
		if(line.find(".NO_VELOCITY.") == 0) {
			containsVelocities = false;
			//delete v_* fields
			auxFields.erase(auxFields.begin(), auxFields.begin()+3);
		} else {
			isDone = true;
		}
		return;
	}

	std::string key = line.substr(0, line.find_last_not_of(" \t\n\r", splitChar - 1) + 1);
	size_t valuestart = line.find_first_not_of(" \t\n\r", splitChar + 1);
	if(valuestart == std::string::npos) valuestart = splitChar+1;
	std::string value = line.substr(valuestart);

	if(key == "Number of particles") {
		numParticles = atoi(value.c_str());
		if(numParticles < 0 || numParticles > 1e9)
			throw Exception("CFG file parsing error. Invalid number of atoms (line " + std::to_string(numLines) + "): " + std::to_string(numParticles) + ")");
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
		entryCount = atoi(value.c_str());
		extendedFormat = true;
	}
	else if(key.compare(0, 10, "auxiliary[") == 0) {
		extendedFormat = true;
		size_t endOfName = value.find_first_of(" \t");
		auxFields.push_back(value.substr(0, endOfName));
	}
	else {
		throw Exception("Unknown key in CFG file header at line " + std::to_string(numLines) + ": " + line + ")");
	}
}

const Matrix3* CFGHeaderData::getH0() const {
	return &H0;
}

const Matrix3* CFGHeaderData::getTransform() const {
	return &transform;
}
