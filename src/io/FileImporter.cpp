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

#include "../io/FileImporter.h"
FileImporter::FileImporter(std::string &inFilename, AtomContainer * inData, FileType inFileType){
	filename = inFilename;
	data = inData;
	fileType = inFileType;
}

FileImporter::~FileImporter() {
}

const char* TextReader::readLine(int maxSize)
{
	lineNumber++;
	line.clear();
	int size = getline(line,maxSize);
	if( eof() ) throw EndOfFile();
	return line.data();
}

const char* TextReader::readLineTrimLeft(int maxSize ) {
	const char* s = readLine(maxSize);
	while(*s != '\0' && (*s == ' ' || *s == '\t')) ++s;
	return s;
}

bool TextReader::eof()
{
 return fileStream->eof();
}

long TextReader::getLineNumber()
{
	return lineNumber;
}

TextReader::TextReader(std::string & inFileName)
{
	filename = inFileName;
	lineNumber = 0;
	fileStream = new std::ifstream;
	fileStream->open(filename.c_str());
	if(fileStream->fail())
			throw Exception("TextReader-FileStream failed to open \"" + filename + "\"");
}

TextReader::~TextReader()
{
	if(!fileStream->is_open()) std::cerr << "TextReader-FileStream has no file opened - will fail when closing" << std::endl;
	if(fileStream->fail())
			throw Exception("TextReader-FileStream failed before closing \"" + filename +"\"");
	fileStream->close();
	if(fileStream->fail())
		throw Exception("TextReader-FileStream failed to close \"" + filename +"\"");
	delete fileStream;
}

std::string TextReader::lineString()
{
	std::string str;
	for(int i = 0; i < line.size(); i++){
		str += line[i];
	}
	return str;
}

int TextReader::getline(std::vector<char> & outLine, int maxNumChars)
{
	std::string s;
	std::getline(*fileStream,s);
	int size;
	if(maxNumChars > 0){
		size = maxNumChars;
		if(s.size() < maxNumChars) size = s.size();
	} else {
		size = s.size();
	}
	outLine.resize(size+1);
	memcpy(outLine.data(), s.c_str(), size);
	outLine.push_back('\0');
	return size;
}







