# enumerate

# Files to compile
SET ( FILES
  treeconstraint
  treelp
  treetraversal
  typetrie
  )

# Prepend folder name
FOREACH ( SOURCE_FILE ${FILES} )
  SET ( SOURCES ${SOURCES} enumerate/${SOURCE_FILE})
ENDFOREACH(SOURCE_FILE)

SET(SOURCES ${SOURCES} PARENT_SCOPE)

if (${REGINA_INSTALL_DEV})
  INSTALL(FILES
    enumconstraints.h
    doubledescription.h
    doubledescription-impl.h
    hilbertcd.h
    hilbertcd-impl.h
    hilbertdual.h
    hilbertdual-impl.h
    hilbertprimal.h
    hilbertprimal-impl.h
    maxadmissible.h
    maxadmissible-impl.h
    ordering.h
    treeconstraint.h
    treelp.h
    treelp-impl.h
    treetraversal.h
    treetraversal-impl.h
    typetrie.h
    DESTINATION ${INCLUDEDIR}/enumerate COMPONENT Development)
endif (${REGINA_INSTALL_DEV})
