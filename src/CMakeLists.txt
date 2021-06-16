cmake_minimum_required(VERSION 3.10)
project (pharmit LANGUAGES C CXX)
# set the path to the library folder
# set the module path
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules )

set( SMINALIB "${SMINA_DIR}/build/libsmina.a" CACHE FILEPATH "Location of libsmina.a")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pthread" )
set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

set( MAIN_SOURCES
     basis.cpp MMappedRegion.h PharmerQuery.cpp QueryTripletFingerprint.cpp ShapeObj.h ThreePointData.h
     basis.h MolFilter.h PharmerQuery.h QueryTripletFingerprint.h ShapeResults.cpp Timer.h
     BitSetTree.h Corresponder.cpp MolProperties.cpp 
     PharmerServerCommands.h
     ReadMCMol.h ShapeResults.h
     Corresponder.h MolProperties.h 
     PharmerServer.cpp
     SimpleFingers.h TripleIndexer.h
     BoundingBox.cpp cors.h MTQueue.h PharmerServer.h RMSD.cpp SphereGrid.cpp Triplet.cpp
     BoundingBox.h dbloader.cpp params.h pharminfo.cpp RMSD.h SphereGrid.h TripletFingerprint.cpp
     BumpAllocator.h dbloader.h pharmarec.cpp pharminfo.h ShapeConstraints.cpp SpinLock.h TripletFingerprint.h
     cgi.cpp 
     FloatCoord.h pharmarec.h PMol.cpp ShapeConstraints.h SPSCQueue.h Triplet.h
     cgi.h
     pharmerdb.cpp PMol.h ThreadCounter.h tripletmatching.cpp
     main.cpp pharmerdb.h queryparsers.h ShapeObj.cpp ThreePointData.cpp tripletmatching.h
    tinyxml/tinystr.cpp 
    tinyxml/tinyxml.cpp 
    tinyxml/tinyxmlerror.cpp 
    tinyxml/tinyxmlparser.cpp 
    shapedb/GSSTreeCreator.cpp 
    shapedb/GSSTreeSearcher.cpp 
    shapedb/GSSTreeStructures.cpp 
    shapedb/KSamplePartitioner.cpp 
    shapedb/MGrid.cpp 
    shapedb/MappableOctTree.cpp 
    shapedb/MemMapped.cpp 
    shapedb/ShapeDistance.cpp 
    shapedb/WorkFile.cpp
    shapedb/molecules/MolSphere.cpp 
    shapedb/molecules/OBMoleculeAnalytic.cpp 
    shapedb/molecules/RDMoleculeAnalytic.cpp
    shapedb/packers/FullMergePacker.cpp 
    shapedb/packers/GreedyPacker.cpp 
    shapedb/packers/MatcherPacker.cpp 
    shapedb/packers/Packer.cpp 
    shapedb/packers/SpectralPacker.cpp
    CommandLine2/CommandLine.cpp
   )

if (NOT LEMON_DIR)
  set(LEMON_DIR "/usr/lib/cmake/") #this is where ubuntu puts the cmake config file
endif()

# add external libraries
find_package( Boost 1.53 COMPONENTS system regex filesystem program_options iostreams thread serialization date_time REQUIRED )
find_package( OpenBabel3 REQUIRED )
find_package( jsoncpp REQUIRED )
find_package( Eigen3 REQUIRED )
find_package( LEMON REQUIRED )
find_package( Gperftools REQUIRED )
find_package( Cgicc REQUIRED )
find_package( Fcgi COMPONENTS streambuf REQUIRED )
find_package( ANN REQUIRED )
find_package( CURL REQUIRED )

get_target_property(OpenBabel3_LIBRARY openbabel LOCATION)

include_directories("${PROJECT_BINARY_DIR}")
include_directories( ${Boost_INCLUDE_DIR} )
include_directories( ${OpenBabel3_INCLUDE_DIRS})
include_directories( $<TARGET_PROPERTY:jsoncpp_lib,INTERFACE_INCLUDE_DIRECTORIES> )
include_directories( ${jsoncpp_INCLUDE_DIRS} )
include_directories( "${EIGEN3_ROOT_DIR}/include" )
include_directories( ${LEMON_INCLUDE_DIRS} )
include_directories( ${bm_SOURCE_DIR} )
include_directories ("${PROJECT_SOURCE_DIR}/CommandLine2")
include_directories( ${GPERFTOOLS_INCLUDE_DIR} )
include_directories( ${CGICC_INCLUDE_DIR} )
include_directories( ${ANN_INCLUDE_DIR} )
include_directories( ${CURL_INCLUDE_DIR} )
include_directories( ${FCGI_INCLUDE_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR})
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/shapedb)
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/shapedb/molecules)
include_directories( ${SMINA_DIR}/src/lib)

add_definitions(-DEIGEN_STACK_ALLOCATION_LIMIT=0)


# add the executable
add_executable(pharmit ${MAIN_SOURCES})
target_link_libraries(pharmit jsoncpp ${SMINALIB} ${GPERFTOOLS_LIBRARIES} ${CGICC_LIBRARY} ${Boost_LIBRARIES} 
                      ${OpenBabel3_LIBRARY} openbabel ${jsoncpp_LIBRARIES} ${EIGEN3_LIBRARY} ${LEMON_LIBRARY} ${ANN_LIBRARY} ${CURL_LIBRARY} ${FCGI_LIBRARIES})

# add the install targets
install (TARGETS pharmit DESTINATION bin)
