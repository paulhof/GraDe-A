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

#ifndef FILEIMPORTER_H_
#define FILEIMPORTER_H_
#include "../GradeA_Defs.h"
#include "../AtomContainer.h"
#include <string>
#include <exception>
//PH: All currently supported filetypes (only cfg)
enum FileType{cfg};
//PH: replace Ovito-Matrix3 class with lightweight counterpart
class TextReader;
class Matrix3{
public:
	double get(char i,char j) const{
		if (i > 2) return 0.;
		if (i < 0) return 0.;
		if (j > 2) return 0.;
		if (j < 0) return 0.;
		return data[i*3+j];
	}
	void set(char i, char j, double value){
		if (i > 2) return;
		if (i < 0) return;
		if (j > 2) return;
		if (j < 0) return;
		data[i*3+j] = value;
	}
	void setIdentity()
	{
		data[0] = 1.; data[1] = 0.; data[2]=0.;
		data[3] = 0.; data[4] = 1.; data[5]=0.;
		data[6] = 0.; data[7] = 0.; data[8]=1.;
	}

	void multiply(double * inVector, double * outVector){
		double outV[3];
		outV[0] = data[0] * inVector[0] + data[1] * inVector[1] + data[2] * inVector[2];
		outV[1] = data[3] * inVector[0] + data[4] * inVector[1] + data[5] * inVector[2];
		outV[2] = data[6] * inVector[0] + data[7] * inVector[1] + data[8] * inVector[2];
		outVector[0] = outV[0];
		outVector[1] = outV[1];
		outVector[2] = outV[2];
	}
	void multiplyAndTranslate(double * inVector, double * translation, double * outVector){
		multiply(inVector,outVector);
		outVector[0] += translation[0];
		outVector[1] += translation[1];
		outVector[2] += translation[2];
	}

	Matrix3 multiply(const Matrix3 * mIn) const {
		double curVal = 0.;
		Matrix3 mOut;
		for(char i = 0; i < 3; i++){
			for(char j = 0; j < 3; j++){
			mOut.set(i,j,
			get(i,0)*mIn->get(0,j) +
			get(i,1)*mIn->get(1,j) +
			get(i,2)*mIn->get(2,j));
			}
		}
		return mOut;
	}
	Matrix3 transposed(){
		Matrix3 mT;
		for(char i = 0; i < 3; i++){
			for(char j = 0; j < 3; j++){
			mT.set(i,j,get(j,i));
			}
		}
		return mT;
	}
private:
	double data[9];
};
struct Exception : public std::exception
{
   std::string s;
   Exception(std::string ss) : s(ss) {}
   ~Exception() throw () {}
   const char* what() const throw() { return s.c_str(); }
};

struct EndOfFile : public std::exception
{
  const char * what () const throw ()
  {
    return "End Of File Reached";
  }
};


class FileImporter {
public:
	FileImporter(std::string &filename, AtomContainer * inData, FileType inFileType);
	virtual ~FileImporter();
protected:
	std::string filename;
	FileType fileType;
	AtomContainer * data;
	TextReader * reader = nullptr;
};

class TextReader{
public:
	TextReader(std::string & inFileName);
	~TextReader();
	bool eof();
	long getLineNumber();
	const char* getLine() const { return line_str.c_str(); }
	const std::string & getLineStr() const { return line_str;}
	bool lineStartsWith(const char* s) const {
		for(const char* l = getLine(); *s; ++s, ++l) {
			if(*l != *s) return false;
		}
		return true;
	}
	const char* readLine();
	//const char* readLineTrimLeft(int maxSize = 0);
	std::string lineString();
private:
	int getline();
	// The name of the input file (if known).
	std::string filename;
	/// Buffer holding the current text line.
	std::string line_str;
	///Some beautiful input stream which contains a huge amount of beer.
	std::ifstream * fileStream = nullptr;
	/// The current line number.
	long lineNumber;
	/// The current position in the uncompressed data stream.
	//long byteOffset;

};
#endif /* FILEIMPORTER_H_ */
