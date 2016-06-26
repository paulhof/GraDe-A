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

#include "GradeA_Defs.h"
#ifndef WINDOWS
std::string ori::to_string(double val, int precision, bool forceScientific){
		std::string formatSpecifier = "%.*g";
		if (forceScientific){
			formatSpecifier = "%.*e";
		}
		char test[2];
		int stringSize = snprintf(test,2,formatSpecifier.c_str(),val,precision);
		char * output = new char [stringSize+1];
		if(snprintf(output,stringSize+1,formatSpecifier.c_str(),val,precision) != stringSize){
			return "";
		}
		std::string outputstr = output;
		delete [] output;
		return outputstr;
}
#endif
//Resolve missing snprintf from C++11 in Cygwin
#ifdef WINDOWS
#include <iomanip>
#include <sstream>
std::string ori::to_string(double val, int precision, bool forceScientific){
	std::stringstream output;
		output.precision(precision);
		if (forceScientific){
			output << std::scientific;
		}
		output << val;
		return output.str();
}
#endif
