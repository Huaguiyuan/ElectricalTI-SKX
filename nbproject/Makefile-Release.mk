#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/Calculations.o \
	${OBJECTDIR}/ElectronSite.o \
	${OBJECTDIR}/ElectronSystem.o \
	${OBJECTDIR}/FiniteTemperature.o \
	${OBJECTDIR}/InputOutput.o \
	${OBJECTDIR}/LLG_Equation.o \
	${OBJECTDIR}/MovieWindow.o \
	${OBJECTDIR}/NeighborList.o \
	${OBJECTDIR}/Nodes.o \
	${OBJECTDIR}/OpenBoundary.o \
	${OBJECTDIR}/SpinSystem.o \
	${OBJECTDIR}/TopologicalCharge.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/mersenne.o \
	${OBJECTDIR}/userintf.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/usr/local/dislin -L/usr/local/armadillo6.1/lib64 -Wl,-rpath,/usr/local/dislin/lib

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/berrycurvaturejbttvar

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/berrycurvaturejbttvar: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/berrycurvaturejbttvar ${OBJECTFILES} ${LDLIBSOPTIONS} -larmadillo -ldiscpp

${OBJECTDIR}/Calculations.o: Calculations.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Calculations.o Calculations.cpp

${OBJECTDIR}/ElectronSite.o: ElectronSite.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ElectronSite.o ElectronSite.cpp

${OBJECTDIR}/ElectronSystem.o: ElectronSystem.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ElectronSystem.o ElectronSystem.cpp

${OBJECTDIR}/FiniteTemperature.o: FiniteTemperature.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FiniteTemperature.o FiniteTemperature.cpp

${OBJECTDIR}/InputOutput.o: InputOutput.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/InputOutput.o InputOutput.cpp

${OBJECTDIR}/LLG_Equation.o: LLG_Equation.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/LLG_Equation.o LLG_Equation.cpp

${OBJECTDIR}/MovieWindow.o: MovieWindow.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/MovieWindow.o MovieWindow.cpp

${OBJECTDIR}/NeighborList.o: NeighborList.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/NeighborList.o NeighborList.cpp

${OBJECTDIR}/Nodes.o: Nodes.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Nodes.o Nodes.cpp

${OBJECTDIR}/OpenBoundary.o: OpenBoundary.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OpenBoundary.o OpenBoundary.cpp

${OBJECTDIR}/SpinSystem.o: SpinSystem.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SpinSystem.o SpinSystem.cpp

${OBJECTDIR}/TopologicalCharge.o: TopologicalCharge.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/TopologicalCharge.o TopologicalCharge.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/mersenne.o: mersenne.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/mersenne.o mersenne.cpp

${OBJECTDIR}/userintf.o: userintf.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/local/dislin -I/usr/local/armadillo6.1/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/userintf.o userintf.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/berrycurvaturejbttvar

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
