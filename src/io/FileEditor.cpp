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

#include "../io/FileEditor.h"

FileEditor::FileEditor(std::string fileName) {
	setFileName(fileName);
}

FileEditor::~FileEditor() {
	close();
}

std::string FileEditor::readNextLine() {
	std::getline(inputReader,inputLine);
	if(inputReader.eof()){
		return "";
	}
	numLines++;
	editedLines.push_back(inputLine);
	curLineLength = inputLine.length();
	inputLine = "";
	return getCurLine();
}

std::string FileEditor::getCurLine() {
	if(editedLines.size() == 0){
		return "";
	}
	return editedLines.back();
}

bool FileEditor::eof() const {
	return inputReader.eof();
}

long FileEditor::getLineNum() {
	return numLines - 1;
}

void FileEditor::replaceCurLine(std::string lineText) {
	if(editedLines.size() == 0){
		return;
	}
	editedLines.back() = lineText;
	curLineLength = lineText.length();
}

bool FileEditor::open() {
	inputReader.open(inputFileName);
	if(!inputReader.is_open()){
		return false;
	}
	tempOutputStream.open(tempOutputFileName);
	if(!tempOutputStream.is_open()){
		inputReader.close();
		return false;
	}
	return true;
}

void FileEditor::flush() {
	for(long iL = 0; iL < editedLines.size(); iL++){
		tempOutputStream << editedLines[iL] + "\n";
	}
	tempOutputStream << std::flush;
	editedLines.resize(0);
}


bool FileEditor::close() {
	if(!tempOutputStream.is_open()){
		return false;
	}
	if(!inputReader.is_open()){
		return false;
	}
	//flush everything to the stream and close
	flush();
	std::string line;
	std::getline(inputReader,line);
	while(!inputReader.eof()){
		tempOutputStream << line << "\n";
		std::getline(inputReader,line);
	}
	tempOutputStream.flush();
	tempOutputStream.close();
	inputReader.close();
	if (tempOutputStream.bad()){
		return false;
	}
	if(inputReader.bad()){
		return false;
	}
	//replace input file with temporary file
	if(0 != remove(inputFileName.c_str())){
		return false;
	}
	return 0 == rename(tempOutputFileName.c_str(), inputFileName.c_str());
}

void FileEditor::setFileName(std::string inFileName) {
	inputFileName = inFileName;
	tempOutputFileName = inputFileName + TEMP_POSTFIX;
}
