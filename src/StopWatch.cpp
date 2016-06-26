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

#include "StopWatch.h"

StopWatch::StopWatch() {}

StopWatch::~StopWatch() {}

void StopWatch::trigger() {
	if(running){
		durationInSeconds = omp_get_wtime()- startTime;
		double rest = durationInSeconds;
		hours = rest/3600;
		rest -= 3600 * hours;
		minutes = rest/60;
		rest -= 60 * minutes;
		seconds = rest;
		running = false;
		return;
	}
	running = true;
	startTime = omp_get_wtime();
}
