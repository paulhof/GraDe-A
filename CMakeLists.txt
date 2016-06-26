# GraDe-A: Grain Detection Algorithm.
# Copyright (C) 2016 Paul Hoffrogge

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8.4)
project(GraDe-A)
IF(CMAKE_BUILD_TYPE)
	string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
ENDIF()
IF(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE_UPPER STREQUAL "RELEASE")
	set(BINARY_NAME "")
	set(CMAKE_BUILD_TYPE Release)
	message("-- Default Build Type set to ${CMAKE_BUILD_TYPE}")
	message("-- Release flags used = ${CMAKE_CXX_FLAGS_RELEASE}")
ELSE()
	add_definitions(-DDEBUGMODE)
	set(BINARY_NAME "_${CMAKE_BUILD_TYPE}")
	message("-- ${CMAKE_BUILD_TYPE} build type chosen")
	message("-- ${CMAKE_BUILD_TYPE} flags used = ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
ENDIF()

IF(CYGWIN)
message("-- Detected CYGWIN Environment")
ENDIF()
IF(WIN32)
message("-- Detected Win32 Environment")
ENDIF()
IF(NOT (CYGWIN OR WIN32))
message("-- Assuming Linux Environment")
ENDIF()

#Add C++11 support
set (C++11FLAG "-std=c++11")
IF(CMAKE_COMPILER_IS_GNUCXX)
	set (C++11FLAG "-std=gnu++11")
ENDIF()

set(MANUAL_C++11_FLAG OFF)
IF(CMAKE_VERSION VERSION_LESS "3.1.0")
	set(MANUAL_C++11_FLAG ON)
ENDIF()

IF(MANUAL_C++11_FLAG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${C++11FLAG}")
	message("-- Added ${C++11FLAG} compiler flag for C++11 support.")
ENDIF()

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

add_executable("grade-A${BINARY_NAME}"
	${CMAKE_SOURCE_DIR}/src/GradeA_Version.cpp
	${CMAKE_SOURCE_DIR}/src/main/GradeA_Main.cpp
	${CMAKE_SOURCE_DIR}/src/ComputationManager.cpp
	${CMAKE_SOURCE_DIR}/src/GlobalMethods.cpp
	${CMAKE_SOURCE_DIR}/src/OrientationMath.cpp
	${CMAKE_SOURCE_DIR}/src/StopWatch.cpp
	${CMAKE_SOURCE_DIR}/src/OrientatorFileQueue.cpp
	${CMAKE_SOURCE_DIR}/src/AtomContainer.cpp
	${CMAKE_SOURCE_DIR}/src/AtomBox.cpp
	${CMAKE_SOURCE_DIR}/src/Atom.cpp
	${CMAKE_SOURCE_DIR}/src/GrainTracker.cpp
	${CMAKE_SOURCE_DIR}/src/CubicLattices.cpp
	${CMAKE_SOURCE_DIR}/src/Orientator.cpp
	${CMAKE_SOURCE_DIR}/src/Orientation.cpp
	${CMAKE_SOURCE_DIR}/src/MeanOrientation.cpp
	${CMAKE_SOURCE_DIR}/src/GrainIdentificator.cpp
	${CMAKE_SOURCE_DIR}/src/Grain.cpp
	${CMAKE_SOURCE_DIR}/src/GrainIDMapper.cpp
	${CMAKE_SOURCE_DIR}/src/ContainerData.cpp
	${CMAKE_SOURCE_DIR}/src/GrainData.cpp
	${CMAKE_SOURCE_DIR}/src/io/FileImporter.cpp
	${CMAKE_SOURCE_DIR}/src/io/CFGImporter.cpp
	${CMAKE_SOURCE_DIR}/src/io/FileEditor.cpp
	${CMAKE_SOURCE_DIR}/src/io/CFGEditor.cpp
	${CMAKE_SOURCE_DIR}/src/io/GrainCSVFileFormat.cpp
	${CMAKE_SOURCE_DIR}/src/io/AtomIO.cpp
	${CMAKE_SOURCE_DIR}/src/io/CSVTableWriter.cpp	
	${CMAKE_SOURCE_DIR}/src/io/CSVTableReader.cpp
	${CMAKE_SOURCE_DIR}/src/io/GrainTimeEvolutionWriter.cpp	
)

IF(NOT MANUAL_C++11_FLAG)
	#automatically add support for C++11
	set_property(TARGET "grade-A${BINARY_NAME}" PROPERTY CXX_STANDARD 11)
	set_property(TARGET "grade-A${BINARY_NAME}" PROPERTY CXX_STANDARD_REQUIRED ON)
	message("-- Automatically added C++11 support")
ENDIF()
set(ARMA_LIB_ENDING ".so")
set(BLAS_LIBFILE "${CMAKE_SOURCE_DIR}/lib/libblas.dll")
set(LAPACK_LIBFILE "${CMAKE_SOURCE_DIR}/lib/liblapack.dll")

#find correct library files for windows
IF (EXISTS "${CMAKE_SOURCE_DIR}/lib/libblas.dll")
	set(BLAS_LIBFILE "${CMAKE_SOURCE_DIR}/lib/libblas.dll")
ELSEIF (EXISTS "${CMAKE_SOURCE_DIR}/lib/blas_win64_MT.dll")
	set(BLAS_LIBFILE "${CMAKE_SOURCE_DIR}/lib/blas_win64_MT.dll")
ENDIF()

IF (EXISTS "${CMAKE_SOURCE_DIR}/lib/liblapack.dll")
	set(LAPACK_LIBFILE "${CMAKE_SOURCE_DIR}/lib/liblapack.dll")
ELSEIF (EXISTS "${CMAKE_SOURCE_DIR}/lib/lapack_win64_MT.dll")
	set(LAPACK_LIBFILE "${CMAKE_SOURCE_DIR}/lib/lapack_win64_MT.dll")
ENDIF()
IF(CYGWIN)
add_definitions(-DCYGWIN)
ENDIF()
IF(WIN32 OR CYGWIN)
	add_definitions(-DWINDOWS)
	message("-- Linking BLAS library ${BLAS_LIBFILE}")
	target_link_libraries("grade-A${BINARY_NAME}" ${BLAS_LIBFILE})
	message("-- Linking LAPACK library ${LAPACK_LIBFILE}")
	target_link_libraries("grade-A${BINARY_NAME}" ${LAPACK_LIBFILE})
	set(ARMA_LIB_ENDING ".dll.a")
ELSE()
	#Link lapack
	FIND_PACKAGE(LAPACK REQUIRED)
	IF(LAPACK_FOUND)
		message("-- Found LAPACK")
		target_link_libraries("grade-A${BINARY_NAME}" ${LAPACK_LIBRARIES})
	ENDIF()

	#Link blas
	FIND_PACKAGE(BLAS)
	IF(BLAS_FOUND)
		message("-- Found BLAS")
		target_link_libraries("grade-A${BINARY_NAME}" ${BLAS_LIBRARIES})
	ENDIF()
ENDIF()
FIND_PACKAGE(Armadillo)
IF(ARMADILLO_FOUND)
	target_link_libraries("grade-A${BINARY_NAME}" ${ARMADILLO_LIBRARIES})
ELSE()
	set(ARMADILLO_DIR "${CMAKE_SOURCE_DIR}/armadillo")
	IF(EXISTS "${ARMADILLO_DIR}/libarmadillo${ARMA_LIB_ENDING}")
		message("-- Found Armadillo in ${ARMADILLO_DIR}")
		INCLUDE_DIRECTORIES("${ARMADILLO_DIR}/include")
		target_link_libraries("grade-A${BINARY_NAME}" ${ARMADILLO_DIR}/libarmadillo${ARMA_LIB_ENDING})
	ELSE()
	message(FATAL_ERROR "Armadillo not found in ${ARMADILLO_DIR}. Please download from http://arma.sourceforge.net/ and build.")
	ENDIF()
ENDIF()