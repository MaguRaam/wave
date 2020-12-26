include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CXXFLAGS += -pedantic -std=c++14

wave: wave.o
	-${CXXLINKER} -o wave wave.o  ${PETSC_LIB}
	${RM} wave.o
