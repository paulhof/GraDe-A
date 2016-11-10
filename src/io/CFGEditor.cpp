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
	skipExtendedLines();
	std::string curLine;
	if(curAtomNum == -1){
		curLine = editor.getCurLine();
	} else {
		curLine = editor.readNextLine();
	}
	size_t fieldBegin = 0;
	size_t fieldEnd = curLine.find_first_of(" \t",fieldBegin);
	while (fieldEnd != std::string::npos){
		if(fieldBegin != fieldEnd){
			curAtomData.push_back(curLine.substr(fieldBegin,fieldEnd-fieldBegin));
		}
		fieldBegin = fieldEnd + 1;
		fieldEnd = curLine.find_first_of(" \t",fieldBegin);
	}
	//parse last field
	if(fieldBegin != curLine.size()){
		curAtomData.push_back(curLine.substr(fieldBegin));
	}
	curAtomNum++;
}

int CFGEditor::getNumAuxFields() const{
	return header.getNumAuxFields();
}

void CFGEditor::setCurAux(int auxNum, std::string value) {
	int fieldNum = auxFieldNum(auxNum);
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

void CFGEditor::skipExtendedLines() {
	if(header.isExtendedFormat()) {
		std::string curLine = editor.getCurLine();
		//if a new type is introduced, the line contains only a single column
		//< \t>*<number><\n>
		//      |<- trimPos
		size_t trimPos = curLine.find_first_not_of(" \t");
		if(curLine.find_first_of(" \t",trimPos) == std::string::npos){
			//new type means skipping two lines (mass and name info)
			editor.readNextLine();
			editor.readNextLine();
		}
	}
}

std::string CFGEditor::curFieldString() {
	std::string ret;
	for(int i = 0; i < curAtomData.size(); i++){
		ret += curAtomData[i] + " ";
	}
	ret.pop_back(); //delete last whitespace
	return ret;
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
