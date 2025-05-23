# CMake project file for generated FMI code
#
# This file is part of FMICodeGenerator (https://github.com/ghorwin/FMICodeGenerator)
# 
# BSD 3-Clause License
#
# Copyright (c) 2018, Andreas Nicolai
# All rights reserved.
#
# see https://github.com/ghorwin/FMICodeGenerator/blob/master/LICENSE for details.

project( Pipe )

cmake_minimum_required(VERSION 3.5)

# add include directories
include_directories(
	${PROJECT_SOURCE_DIR}/../../src
)

add_library( ${PROJECT_NAME} SHARED
	${PROJECT_SOURCE_DIR}/../../src/${PROJECT_NAME}.cpp
	${PROJECT_SOURCE_DIR}/../../src/fmi2common/fmi2Functions.cpp
	${PROJECT_SOURCE_DIR}/../../src/fmi2common/InstanceData.cpp
)

# link against the dependent libraries
target_link_libraries( ${PROJECT_NAME}
	${APPLE_FRAMEWORKS}
)

# Library version strings set by code generator
set( SO_VERSION_MAJOR 1)
set( SO_VERSION_MINOR 0)
set( SO_VERSION_PATCH 0)
set( SO_VERSION_STRING
	${SO_VERSION_MAJOR}.${SO_VERSION_MINOR}.${SO_VERSION_PATCH}
)

set_target_properties( ${PROJECT_NAME} PROPERTIES
	VERSION ${SO_VERSION_STRING}
	SOVERSION ${SO_VERSION_MAJOR}
)
