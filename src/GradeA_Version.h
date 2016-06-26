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

#ifndef ORIENTATORVERSION_H_
#define ORIENTATORVERSION_H_
#include <string>
#define MAJOR_VERSION 1
#define MINOR_VERSION 0
#define BUILD_NUMBER 0

//
// Version History:
//
//   1.0.0: First runnable version:
//     - atom wise orientation calculation
//     - grain identification
//     - orphan atom adoption
//     - multifile time-tracking
//     - parallelized with openmp
//     - time evolution of grain-properties written in "./TimeEvo/" folder

class OrientatorVersion {
public:
	OrientatorVersion();
	virtual ~OrientatorVersion();
	static std::string getString();
};

const OrientatorVersion version;
#endif /* ORIENTATORVERSION_H_ */
