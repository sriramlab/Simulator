cmake_minimum_required(VERSION 2.8.12)
project(EM_C++)

include_directories(include/)

set(CMAKE_CXX_FLAGS "-std=c++0x -O3 -lrt -lpthread ")

IF (NOT DEFINED SSE_SUPPORT)
    SET(SSE_SUPPORT 0)
ENDIF()
ADD_DEFINITIONS(-DSSE_SUPPORT=${SSE_SUPPORT})

IF (NOT DEFINED DEBUG)
    SET(DEBUG 0)
ENDIF()
ADD_DEFINITIONS(-DDEBUG=${DEBUG})




add_executable(Simulator_mafld src/simulator_maf_ld_arch.cpp src/genotype.cpp src/storage.cpp)
add_executable(Simulator_annot src/simulator.annot.cpp src/genotype.cpp src/storage.cpp)

add_executable(Simulator_gxe src/simulator_maf_ld_arch_gxe.cpp src/genotype.cpp src/storage.cpp)
add_executable(Simulator_gxe_with_corr src/simulator_maf_ld_arch_gxe_with_corr.cpp src/genotype.cpp src/storage.cpp)
add_executable(Simulator_gxe_same_causal src/simulator_maf_ld_arch_gxe_same_causal.cpp src/genotype.cpp src/storage.cpp)
add_executable(Simulator_gxe_causal_from_file src/simulator_maf_ld_arch_gxe_causal_from_file.cpp src/genotype.cpp src/storage.cpp)
add_executable(Simulator_gxe_heavy_tail src/simulator_maf_ld_arch_gxe_heavy_tail.cpp src/genotype.cpp src/storage.cpp)
