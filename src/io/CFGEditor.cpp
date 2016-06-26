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

#include "CFGEditor.h"

CFGEditor::CFGEditor(std::string inFileName) {
	editor.setFileName(inFileName);
	if(!editor.open()){
		std::cerr << "Failed opening file \"" << inFileName << "\" for editing" << std::endl;
		return;
	}
	parseHeader();
}

CFGEditor::~CFGEditor() {
	close();
}

void CFGEditor::readNextAtom() {
	curAtomData.clear();
	std::string curLine = editor.readNextLine();
	skipExtendedLines();
	size_t fieldBegin = 0;
	size_t fieldEnd = curLine.find_first_of(" \t",fieldBegin);
	while (fieldEnd != std::string::npos){
		if(fieldBegin != fieldEnd){
			curAtomData.push_back(curLine.substr(fieldBegin,fieldEnd-fieldBegin));
		}
		fieldBegin = fieldEnd + 1;
		fieldEnd = curLine.find_first_of(" \t",fieldBegin);
	}
	curAtomNum++;
}

int CFGEditor::getNumAuxFields() const{
	return header.getNumAuxFields();
}

void CFGEditor::setCurAux(int auxNum, std::string value) {
	int fieldNum = auxNum + numPosFields + numVelFields;
	if(fieldNum < 0 || fieldNum >= curAtomData.size()) {
		return;
	}
	curAtomData[fieldNum] = value;
	editor.replaceCurLine(curFieldString());
}

long CFGEditor::getNumAtoms() const{
	return header.getNumParticles();
}

long CFGEditor::getCurAtomNum() const{
	return curAtomNum;
}

void CFGEditor::close() {
	editor.close();
}

void CFGEditor::parseHeader() {
	//ensure that header is only parsed once
	if(headerParsed){
		return;
	}
	header.parse(&editor);
	if (header.containsVelFields()){
		numVelFields = 3;
	}
	headerParsed  = true;
}

void CFGHeaderData::parse(FileEditor* editor) {
	numParticles = -1;
	unitMultiplier = 1;
	H0.setIdentity();
	transform.setIdentity();
	rateScale = 1.;
	extendedFormat = false;
	containsVelocities = true;

	int entry_count = 0;
	while(!editor->eof()) {
			numLines++;
			std::string line(editor->readNextLine());
			// Ignore comments
			size_t commentChar = line.find('#');
			if(commentChar != std::string::npos) line.resize(commentChar);

			// Skip empty lines.
			size_t trimmedLine = line.find_first_not_of(" \t\n\r");
			if(trimmedLine == std::string::npos) continue;
			if(trimmedLine != 0) line = line.substr(trimmedLine);

			size_t splitChar = line.find('=');
			if(splitChar == std::string::npos) {
				if(line.find(".NO_VELOCITY.") == 0) {
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
					throw Exception("CFG file parsing error. Invalid number of atoms (line " + std::to_string(editor->getLineNum()) + "): " + std::to_string(numParticles) + ")");
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
				extendedFormat = true;
			}
			else if(key.compare(0, 10, "auxiliary[") == 0) {
				extendedFormat = true;
				size_t endOfName = value.find_first_of(" \t");
				auxFields.push_back(value.substr(0, endOfName));
			}
			else {
				throw Exception("Unknown key in CFG file header at line " + std::to_string(editor->getLineNum()) + ": " + line + ")");
			}
		}
		if(numParticles < 0)
			throw Exception("Invalid file header. This is not a valid CFG file.");
}

int CFGHeaderData::getNumAuxFields() const{
	return auxFields.size();
}

std::string CFGHeaderData::getAuxField(int fieldNum) const{
	if (fieldNum < 0 || fieldNum >= auxFields.size()){
		return "";
	}
	return auxFields[fieldNum];
}

void CFGEditor::skipExtendedLines() {
	if(!header.isExtendedFormat()) {
		return;
	}
	std::string curLine = editor.getCurLine();
	size_t trimPos = curLine.find_first_not_of(" \t");
	if(curLine.find_first_not_of(" \t",trimPos) != std::string::npos){
		//no new type found means nothing to skip
		return;
	}
	//new type means skipping two lines (mass and name info)
	editor.readNextLine();
	editor.readNextLine();
}

std::string CFGEditor::curFieldString() {
	std::string ret;
	for(int i = 0; i < curAtomData.size(); i++){
		ret += curAtomData[i] + " ";
	}
	ret.pop_back(); //delete last whitespace
	return ret;
}

bool CFGHeaderData::containsVelFields() const{
	return containsVelocities;
}

std::string CFGEditor::getAuxName(int auxNum) const{
	if (auxNum < 0 || auxNum >= getNumAuxFields()){
			return "";
	}
	return header.getAuxField(auxNum);
}

bool CFGEditor::eof() const {
	return editor.eof();
}

void CFGEditor::flush() {
	editor.flush();
}

std::string CFGEditor::getCurAux(int auxNum) const{
	if(auxNum < 0 || auxNum >= getCurNumAuxFields()){
		return "";
	}
	return curAtomData[auxFieldNum(auxNum)];
}

int CFGEditor::getCurNumAuxFields() const {
	return curAtomData.size() - numPosFields - numVelFields;
}

int CFGEditor::auxFieldNum(int auxNum) const {
	return auxNum + numPosFields + numVelFields;
}
