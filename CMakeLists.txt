# Copyright (C) 2021 Marco Krummenacher
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8)
project(Fingerprint Fortran)
enable_language(Fortran)

set(CMAKE_CXX_STANDARD 14)

option(INTEL "compile with ifort and mkl" OFF)
option(OPENMP "compile with openmp support" OFF)
option(DEBUG "compile with debug flags" OFF)
option(CRYSTAL "crystal version of fingerprint" OFF)
option(CLUSTER "cluster version of fingerprint" OFF)
option(DISTANCE "distance version of fingerprint" OFF)

if(INTEL)
    set(CMAKE_Fortran_COMPILER ifort)
else()
    set(CMAKE_Fortran_COMPILER gfortran)
endif()

if(INTEL)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -no-wrap-margin") # -Ofast -ipo
    if(DEBUG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -mkl -xHost -fp-model strict -check bounds -g -traceback") # -Ofast -ipo
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost -mkl -Ofast") # -ipo
    endif()
    if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -parallel -parallel-source-info=1 ")# -parallel")
    endif()
else()
    # module dir
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -std=legacy")
    if(DEBUG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -Wall -Wconversion -g -Og -fbacktrace") # -ffpe-trap=zero,invalid") #-flto
    else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -ffast-math -march=native") #-flto
    endif()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    set(libraries ${libraries} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    if(OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
    endif()
endif()


if(CRYSTAL)
add_executable(FP_crystal.x
     src/fingerprint.f90
     src/fp_for_crystals.f90
)
target_link_libraries(FP_crystal.x ${libraries})
endif()
if(CLUSTER)
add_executable(FP_cluster.x
     src/fingerprint.f90
     src/fp_for_clusters.f90
)
target_link_libraries(FP_cluster.x ${libraries})
endif()
