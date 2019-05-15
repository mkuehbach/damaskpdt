cmake_minimum_required(VERSION 3.0)

################################################################################################################
##USER INTERACTION##############################################################################################
################################################################################################################
#please name your project accordingly
set(MYPROJECTNAME "damaskpdt_vorocomposer_debug")

#tell the top directory where this local DAMASKPDT version is stored
set(MYPROJECTDIR "/home/m.kuehbach/DAMASKPDT_VoroComposer_FINAL")
#set(MYPROJECTDIR "/talos/u/mkuehbac/DAMASKPDT_VoroComposer_FINAL")

#choose compiler
set(EMPLOY_INTELCOMPILER ON) #either or
set(EMPLOY_GNUCOMPILER OFF)

#choose optimization level
##-O0 nothing, debug purposes, -O1 moderate optimization, -O2 -O3 or production level up to aggressive optimization
if(EMPLOY_INTELCOMPILER)
	set(MYOPTLEVEL "-O3 -g -march=skylake") #-O3 -g -march=skylake") #"-O0 -g"
endif()
if(EMPLOY_GNUCOMPILER)
	set(MYOPTLEVEL "-O3 -g -march=native")
endif()

#necessary libraries and thirdparty dependencies
set(EMPLOY_MYHDF5 ON)
##set(EMPLOY_MYTETGEN OFF)
set(EMPLOY_MYIMKL ON)
set(EMPLOY_MYBOOST ON)
##set(EMPLOY_MYBSIMD OFF)
set(VERBOSE_VECTORIZATION OFF)

################################################################################################################################
#END OF INTERACTION FOR NON PRO USERS###########################################################################################
##in this section advanced users might want/need to make modifications if they use non default places for thirdparty libraries##
################################################################################################################################ 

#HDF5 BINARY CONTAINER IO LIBRARY
if(EMPLOY_MYHDF5)
	#HDF5 local installation for advanced I/O, collecting metadata and bundle analysis results together
	set(MYHDFPATH "/home/m.kuehbach/APT3DOIM/src/thirdparty/HDF5/CMake-hdf5-1.10.2/build/HDF5-1.10.2-Linux/HDF_Group/HDF5/1.10.2")
	##set(MYHDFPATH "/talos/u/mkuehbac/HDF5/CMake-hdf5-1.10.2/build/HDF5-1.10.2-Linux/HDF_Group/HDF5/1.10.2")
	include_directories("${MYHDFPATH}/include")
	link_directories("${MYHDFPATH}/lib")
	set(MYHDFLINKFLAGS "-L${MYHDFPATH}/lib/ ${MYHDFPATH}/lib/libhdf5_hl.a ${MYHDFPATH}/lib/libhdf5.a -lz -ldl")
else()
	set(MYHDFPATH "")
	set(MYHDFLINKFLAGS "")
	message([STATUS] "HDF5 not activated")
endif()

if(EMPLOY_TETGEN)
	## no TetGen
else()
	message([STATUS] "TetGen library not activated")
endif()

#INTEL MATH KERNEL LIBRARY threaded numerical library Intel MKL
if(EMPLOY_MYIMKL)
#	#MAWS30
#	#manual finding of IntelMKL where do we find the IntelMKL library
#	#check also https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#	set(MYMKLROOT "/opt/intel/compilers_and_libraries_2017/linux/mkl//")
#	set(MYMKLPATH "${MYMKLROOT}/lib/intel64_lin")
#	include_directories("${MYMKLROOT}/include")
#	set(MYMKL_COMPILEFLAGS "")
#	set(MYMKL_LINKFLAGS "-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl")

	#TALOS
	#INTEL MATH KERNEL LIBRARY threaded numerical library Intel MKL
	#manual finding of IntelMKL where do we find the IntelMKL library
	set(MYMKLROOT "/mpcdf/soft/SLE_15/packages/x86_64/intel_parallel_studio/2018.4/mkl")
	include_directories("${MYMKLROOT}/include")
	set(MYMKLCOMP "-I${MYMKLROOT}/include")
	message([STATUS] "${MYMKLCOMP}")
	message([STATUS] "${MYMKLROOT}")
	set(MYMKLLINK "-L${MYMKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -Wl,-rpath,${MYMKLROOT}/lib/intel64")
	message([STATUS] "${MYMKLLINK}")
else()
	message([FATAL_ERROR] "We need the IntelMKL library!")
endif()



#Vectorization using BSIMD location
if(EMPLOY_MYBSIMD)
	include_directories("${MYPROJECTDIR}/src/thirdparty/boost.simd/include")
	include_directories("${MYPROJECTDIR}/src/thirdparty/boost.simd/include/boost/simd")
	set(MYSIMDFLAGS "-msse4.2")
	message([STATUS] "Use shipped with boost.simd")
else()
	message([STATUS] "boost.simd not used")
	set(MYSIMDFLAGS "")
endif()

################################################################################################################
##AUTOMATIC SECTION#############################################################################################
################################################################################################################
#user input sanity checks
if(EMPLOY_INTELCOMPILER AND EMPLOY_GNUCOMPILER)
	message([FATAL_ERROR] "You cannot utilize two compiler at the same time!")
endif()


#automatically assign project name and compiler flags
project(${MYPROJECTNAME})
set(CMAKE_BUILD_DIR "build")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MYOPTLEVEL}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MYOPTLEVEL}")

#setting up compiler-specifics
#intel path
if(EMPLOY_INTELCOMPILER)
	message([STATUS] "Employing the Intel compiler!")
	if (VERBOSE_VECTORIZATION)
		set(MYVERBOSE "-qopt-report=5")
	else()
		set(MYVERBOSE "")
	endif()
	add_definitions("${MYOPTLEVEL} ${MYVERBOSE} ${MYSIMDFLAGS}")
	set(MYOMP "-qopenmp -lpthread")
	add_definitions(${MYOMP})
	set(MYLARGEFILESUPPORT "")
	add_definitions("-std=c++0x")
	add_definitions("-Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat 
		-Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull 
		-Wparentheses -Wpointer-sign -Wreorder -Wreturn-type -Wsign-compare 
		-Wsequence-point -Wtrigraphs -Wunused-variable")
	add_definitions("${MYLARGEFILESUPPORT}")
elseif(EMPLOY_GNUCOMPILER)
	message([STATUS] "Employing the GNU compiler!")
	if(VERBOSE_VECTORIZATION)
		set(MYVERBOSE "-fopt-info-all")
	endif()
	add_definitions("${MYOPTLEVEL} ${MYVERBOSE} ${MYSIMDFLAGS}")
	
	set(MYOMP "-fopenmp -lpthread")
	add_definitions(${MYOMP})
	set(MYLARGEFILESUPPORT "-D_FILE_OFFSET_BITS=64")
	#set(BSIMDPATH "-I=/home/m.kuehbach/DAMASK/damaskpdt/src/thirdparty/boost.simd/include -I=/home/m.kuehbach/DAMASK/damaskpdt/src/thirdparty/boost.simd/include/boost/simd")
	#add_definitions(${BSIMDPATH})
	add_definitions("-std=c++11")
	add_definitions("-Wall -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat 
		-Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull 
		-Wparentheses -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point 
		-Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable")
	add_definitions("${MYLARGEFILESUPPORT}")
else()
	message([FATAL_ERROR] "You have to utilize a compiler!")
endif()

message([STATUS] "Projectname is ${MYPROJECTNAME}")
message([STATUS] "We utilize optimization level ${MYOPTLEVEL}")

#parallelization - MPI process-level
#query location of MPI library
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(MYSRCPATH "${MYPROJECTDIR}/src/")
set(MYVOROSRCPATH "${MYSRCPATH}/thirdparty/VoroRycroft/voro++-0.4.6/src/")

add_executable(${MYPROJECTNAME}
	${MYSRCPATH}PDT_Profiler.cpp
	${MYSRCPATH}PDT_Settings.cpp
	${MYSRCPATH}PDT_Performance.cpp
	${MYSRCPATH}PDT_Datatypes.cpp
	${MYSRCPATH}PDT_DatatypesMPI.cpp
	${MYSRCPATH}PDT_TensorMath.cpp
	${MYSRCPATH}PDT_OriMath.cpp
	${MYSRCPATH}PDT_Math.cpp
	${MYSRCPATH}PDT_AABBTree.cpp
	${MYSRCPATH}PDT_HDF5.cpp
	${MYSRCPATH}PDT_XDMF.cpp
	${MYSRCPATH}PDT_VTKIO.cpp

	${MYVOROSRCPATH}c_loops.cc
	${MYVOROSRCPATH}cell.cc
	#${MYVOROSRCPATH}cmd_line.cc
	${MYVOROSRCPATH}common.cc
	${MYVOROSRCPATH}container.cc
	${MYVOROSRCPATH}container_prd.cc
	${MYVOROSRCPATH}pre_container.cc
	${MYVOROSRCPATH}unitcell.cc
	${MYVOROSRCPATH}v_base.cc
	#${MYVOROSRCPATH}v_base_wl.cc
	${MYVOROSRCPATH}v_compute.cc
	#${MYVOROSRCPATH}voro++.cc
	${MYVOROSRCPATH}wall.cc

	${MYSRCPATH}PDT_VoroComposer.cpp
	${MYSRCPATH}PDT_IntelMKL.cpp
	${MYSRCPATH}PDT_Crystallography.cpp
	${MYSRCPATH}PDT_GrainObject.cpp
	${MYSRCPATH}LOUVAIN_GraphBinary.cpp
	${MYSRCPATH}LOUVAIN_Core.cpp
	${MYSRCPATH}PDT_SpecOutReader.cpp
	${MYSRCPATH}PDT_OptionParser.cpp
	${MYSRCPATH}PDT_Main.cpp
)

#BOOST CXX HEADER LIBRARY
if(EMPLOY_MYBOOST)
	message([STATUS] "We need Boost as well...")
	include_directories(${Boost_INCLUDE_DIR})
	link_directories(${Boost_LIBRARY_DIR})
	find_package(Boost COMPONENTS system filesystem REQUIRED)
else()
	message([FATAL_ERROR] "We need the Boost Cpp library!")
endif()


#linking process
if(EMPLOY_INTELCOMPILER)
	target_link_libraries(${MYPROJECTNAME} ${MYOMP} ${MYLARGEFILESUPPORT} ${Boost_LIBRARIES} ${MYMKLLINK} ${MPI_LIBRARIES} ${MYHDFLINKFLAGS} )
elseif(EMPLOY_GNUCOMPILER)
	target_link_libraries(${MYPROJECTNAME} ${MYOMP} ${MYLARGEFILESUPPORT} ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${MYHDFLINKFLAGS}  )
else()
	message([FATAL_ERROR] "You have to utilize a compiler!")
endif()

#MPI compilation settings
if(MPI_COMPILE_FLAGS)
  set_target_properties(${MYPROJECTNAME} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${MYPROJECTNAME} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
