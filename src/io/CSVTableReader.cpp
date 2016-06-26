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

#include "CSVTableReader.h"

CSVTableReader::CSVTableReader() {
	fileName = "";
	separator = ";";
	numHeaderLines = 0;
}

CSVTableReader::CSVTableReader(std::string inFilename, int inNumHeaderLines, std::string inSeparator ){
	fileName = inFilename;
	numHeaderLines = inNumHeaderLines;
	setCustomSeparator(inSeparator);
}

CSVTableReader::~CSVTableReader() {
}

bool CSVTableReader::parse() {
	if(!openFile()){
			std::cerr << "Cannot open file \"" << fileName << "\" for reading." << std::endl;
			return false;
	}

	std::string lineString;
	//first extract header-lines
	for(int iH = 0; iH < numHeaderLines; iH++){
		if(!filestream.good()){
			return false;
		}
		std::getline(filestream, lineString);
		addHeaderLine(lineString);
	}

	//parse title line into columnNames
	std::getline(filestream, lineString);
	parseTitleLine(lineString);

	//parse stream line by line
	while(filestream.good()){
		//get the current line string
		std::getline(filestream, lineString);
		//now parse the corresponding entries into a new line
		parseLine(lineString);
	}
	//delete last line when its empty and there are at least 2 columns
	if(lines.size() > 0 && columnNames.size() > 1 ){
		if(lines.back().size() == 1){
			std::cout << fileName <<" back = \"" << lines.back().back() << "\"" << std::endl;
			if(lines.back().back() == ""){
				lines.pop_back();
			}
		}
	}


	if(!closeFile()){
		std::cerr << "Cannot close file \"" << fileName << "\" after reading." << std::endl;
		return false;
	}

	return true;
}

void CSVTableReader::parseHeaderLine(std::string& inHeaderLine) {
	size_t prevPos = 0;
	size_t curPos = 0;
	std::string curName;
	//parse all separated entries of this line
	do {
		curPos = inHeaderLine.find(separator, prevPos);
		curName = inHeaderLine.substr(prevPos, curPos-prevPos);
		prevPos = curPos + separator.length();
		//now add the entry to an empty column
		addEntryToHeaderLine(curName);
	} while (curPos != std::string::npos);
}

void CSVTableReader::parseTitleLine(std::string& inTitleLine) {
	size_t prevPos = 0;
	size_t curPos = 0;
	std::string curName;
	//parse all separated entries of this line
	do {
		curPos = inTitleLine.find(separator, prevPos);
		curName = inTitleLine.substr(prevPos, curPos-prevPos);
		prevPos = curPos + separator.length();
		//now add the entry to an empty column
		addNewColumn(curName);
	} while (curPos != std::string::npos);
}

void CSVTableReader::parseLine(std::string& inLine) {
	size_t prevPos = 0;
	size_t curPos = 0;
	std::string curEntry;
	// add a new line to parse in
	addNewLine();
	long columnId = 0;
	do {
		curPos = inLine.find(separator, prevPos);
		curEntry = inLine.substr(prevPos, curPos-prevPos);
		prevPos = curPos + separator.length();
		//now add the entry to an empty column
		addEntryToCurrentLine(curEntry, columnId++);
	} while (curPos != std::string::npos);
}


bool CSVTableReader::addNewColumn(std::string columnName)
{
	columnNames.push_back(columnName);
	return true;
}

bool CSVTableReader::addEntryToHeaderLine(std::string entry) {
	if (headerLines.size() <= 0){
		headerLines.resize(1);
	}
	headerLines.back().push_back(entry);
	return true;
}
bool CSVTableReader::addEntryToCurrentLine(std::string entry, std::string inColumnName)
{
	for (int i = 0; i < columnNames.size(); i++){
		if(inColumnName == columnNames[i] ){
			safeAddToLine(entry,i);
			return true;
		}
	}
	return false;
}

bool CSVTableReader::addEntryToCurrentLine(std::string entry, long  columnId)
{
	if(columnId >= columnNames.size()){
		return false;
	}
	safeAddToLine(entry,columnId);
	return true;
}

void CSVTableReader::addNewLine()
{
	lines.resize(lines.size()+1);

}

void CSVTableReader::setFileName(std::string inFilename)
{
	fileName = inFilename;
}

void CSVTableReader::setCustomSeparator(std::string inSeparator)
{
	separator = inSeparator;
}

bool CSVTableReader::openFile()
{
	try{
		filestream.open(fileName.c_str());
	} catch (...){
		return false;
	}
	return filestream.is_open();
}

bool CSVTableReader::closeFile()
{
	try{
		filestream.close();
	} catch (...){
		return false;
	}
	return !filestream.is_open();
}



void CSVTableReader::addHeaderLine(std::string line)
{
	headerLines.resize(headerLines.size()+1);
	parseHeaderLine(line);
}

std::string CSVTableReader::entry(long lineNum, int colNum) const{
	if(lineNum >= lines.size() || lineNum < 0){
		return "";
	}
	if(colNum >= lines[lineNum].size() || colNum < 0){
		return "";
	}
	return lines[lineNum][colNum];
}

long CSVTableReader::getNumLines() const{
	return lines.size();
}

long CSVTableReader::getNumColumns() const{
	return columnNames.size();
}

int CSVTableReader::getNumHeaderLines() const {
	return headerLines.size();
}

std::string CSVTableReader::headerEntry(int lineNum, int colNum) const {
	if(lineNum >= headerLines.size() || lineNum < 0){
			return "";
	}
	if(colNum >= headerLines[lineNum].size() || colNum < 0){
		return "";
	}
	return headerLines[lineNum][colNum];
}

std::string CSVTableReader::headerLine(int lineNum) const {
	if (lineNum >= headerLines.size()){
		return "";
	}
	std::string lineString = "";
	for	(int iC = 0; iC < headerLines[lineNum].size(); iC++){
		lineString += headerEntry(lineNum,iC) + separator;
	}
	return lineString;
}

void CSVTableReader::safeAddToLine(std::string entry, long  columnId)
{
	if(lines.size() <= 0) {
		lines.resize(1);
	}
	if(lines.back().size() <= columnId){
		lines.back().resize(columnId+1, "");
	}
	lines.back()[columnId] = entry;
}

std::string CSVTableReader::getColumnName(int colNum) const{
	if(colNum < 0 || colNum >= columnNames.size()){
		return "";
	}
	return columnNames[colNum];
}
