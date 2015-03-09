project(BaderCycle)

cmake_minimum_required(VERSION 2.8)

add_option(ENABLE_MPI "Enable MPI support?" ON)

if(ENABLE_MPI)
  find_package(MPI REQUIRED)

  # Legg til mpi til include path
  include_directories(${MPI_C_INCLUDE_PATH})

  # Sett MPI flagg
  add_definitions(${MPI_C_COMPILE_FLAGS})

  # Sett define for MPI-støtte
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DENABLE_MPI=1")
endif()

# Legg til c99
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

# Liste over kildekode for bader-applikasjonen
set(SOURCES bader.c
            GraphUtilities.c
            DistGather.c
            mpiUtils.c
            baderDetect.c)

# Legg til baderapplikasjon
add_executable(bader ${SOURCES})

# Lenk bader mot MPI
target_link_libraries(bader m)

if(ENABLE_MPI)
  target_link_libraries(bader ${MPI_C_LIBRARIES})
endif()
