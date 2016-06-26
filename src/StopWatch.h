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

#ifndef STOPWATCH_H_
#define STOPWATCH_H_
#include "GradeA_Defs.h"
class StopWatch {
public:
	StopWatch();
	virtual ~StopWatch();
	void trigger();

	unsigned int getHours() const {
		return hours;
	}

	unsigned int getMinutes() const {
		return minutes;
	}

	double getSeconds() const {
		return seconds;
	}

	double getDuration() const{
		return durationInSeconds;
	}

	std::string getString() const{
		return std::to_string(hours) + " hours " + std::to_string(minutes) + " minutes " + std::to_string(seconds)+  " seconds" ;
	}

private:
	bool running = false;
	double startTime = 0.;
	double durationInSeconds = 0.;
	unsigned int hours = 0;
	unsigned int minutes = 0;
	double seconds = 0.;
};

#endif /* STOPWATCH_H_ */
