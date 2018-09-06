cmake_minimum_required(VERSION 3.0)

################################################################################################################
##USER INTERACTION##############################################################################################
################################################################################################################
#please name your project accordingly
set(MYPROJECTNAME "damaskpdt")

#choose compiler
set(EMPLOY_INTELCOMPILER ON)
set(EMPLOY_GNUCOMPILER OFF)
set(VERBOSE_VECTORIZATION OFF)

#optimization level -O0, -O1, -O2, or -O3
set(MYOPTLEVEL "-O3 -pg")

#BoostSIMD location
#this utilizes the copy of the boost simd shipped with damaskpdt utilize absolute addressing...
if(EMPLOY_INTELCOMPILER)
	include_directories("/nethome/m.kuehbach/damaskpdt/src/thirdparty/boost.simd/include")
	include_directories("/nethome/m.kuehbach/damaskpdt/src/thirdparty/boost.simd/include/boost/simd")
elseif(EMPLOY_GNUCOMPILER)
	include_directories("/home/m.kuehbach/DAMASK/damaskpdt/src/thirdparty/boost.simd/include")
	include_directories("/home/m.kuehbach/DAMASK/damaskpdt/src/thirdparty/boost.simd/include/boost/simd")
else()
	message([FATAL_ERROR] "You have to utilize a compiler!")
endif()
set(MYSIMDFLAGS "-msse4.2")
message([STATUS] "Use shipped with boost.simd")

#Intel MKL
#manual finding of IntelMKL where do we find the IntelMKL library
set(MYMKLROOT "/opt/intel/compilers_and_libraries_2017/linux/mkl//")
set(MYMKLPATH "${MYMKLROOT}/lib/intel64_lin")
include_directories("${MYMKLROOT}/include")
#set(MYMKL "-lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
#set(MYMKL "-DMKL_ILP64 -parallel -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl")
#set(MYMKL_COMPILEFLAGS "-mkl -DMKL_ILP64")
set(MYMKL_COMPILEFLAGS "")
set(MYMKL_LINKFLAGS "-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl")

#check also https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#works for now
#set(MYMKL_COMPILEFLAGS "")
#set(MYMKL_LINKFLAGS "-parallel -Wl,--start-group ${MYMKLPATH}/libmkl_intel_lp64.a ${MYMKLPATH}/libmkl_sequential.a ${MYMKLPATH}/libmkl_core.a -Wl,--end-group -lpthread")


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
		set(MYVERBOSE "-qopt-report-phase=vec")
	endif()
	set(MYLARGEFILESUPPORT "") #no particular action for intel but record size
	add_definitions("${MYOPTLEVEL} ${MYVERBOSE} ${MYSIMDFLAGS}")
	set(MYOMP "-qopenmp -lpthread")
	##BSIMD
	##BSIMD
	add_definitions(${MYOMP})
	add_definitions("-std=c++0x")
	add_definitions("-Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat 
		-Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull 
		-Wparentheses -Wpointer-sign -Wreorder -Wreturn-type -Wsign-compare 
		-Wsequence-point -Wtrigraphs -Wunused-variable")
	add_definitions("${MYLARGEFILESUPPORT}")
#gcc path
elseif(EMPLOY_GNUCOMPILER)
	message([STATUS] "Employing the GNU compiler!")
	if(VERBOSE_VECTORIZATION)
		set(MYVERBOSE "-fopt-info-all")
	endif()
	set(MYLARGEFILESUPPORT "-D_FILE_OFFSET_BITS=64")
	add_definitions("${MYOPTLEVEL} ${MYVERBOSE} ${MYSIMDFLAGS}")
	set(MYOMP "-fopenmp -lpthread")
	#set(BSIMDPATH "-I=/home/m.kuehbach/DAMASK/damaskpdt/src/thirdparty/boost.simd/include -I=/home/m.kuehbach/DAMASK/damaskpdt/src/thirdparty/boost.simd/include/boost/simd")
	#add_definitions(${BSIMDPATH})
	add_definitions(${MYOMP})
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

add_executable(${MYPROJECTNAME}
	../src/PDT_Profiler.cpp
	../src/PDT_Settings.cpp
	../src/PDT_Performance.cpp
	../src/PDT_Datatypes.cpp
	../src/PDT_DatatypesMPI.cpp
	../src/PDT_TensorMath.cpp
	../src/PDT_OriMath.cpp
	../src/PDT_VTKIO.cpp
	../src/PDT_IntelMKL.cpp
	../src/PDT_Crystallography.cpp
	../src/PDT_GrainObject.cpp
	../src/LOUVAIN_GraphBinary.cpp
	../src/LOUVAIN_Core.cpp
	../src/PDT_SpecOutReader.cpp
	../src/PDT_OptionParser.cpp
	../src/PDT_Main.cpp
)

#tell us where to find the BoostSIMD library
#message([STATUS] "We need the boost probable vector intrinsics library...")


##assure the Boost package is there to allow and ease the reading of files and folders from the filesystem
message([STATUS] "We need Boost as well...")
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
find_package(Boost COMPONENTS system filesystem REQUIRED)

#linking process
if(EMPLOY_INTELCOMPILER)
	target_link_libraries(${MYPROJECTNAME} ${MYOMP} ${MYLARGEFILESUPPORT} ${Boost_LIBRARIES} ${MYMKL_LINKFLAGS} ${MPI_LIBRARIES} )
elseif(EMPLOY_GNUCOMPILER)
	target_link_libraries(${MYPROJECTNAME} ${MYOMP} ${MYLARGEFILESUPPORT} ${Boost_LIBRARIES} ${MPI_LIBRARIES} )
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


