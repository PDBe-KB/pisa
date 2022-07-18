cmake_minimum_required(VERSION 2.8.3)
project(pisa C CXX)

option(USE_QT "Build QtPisa" ON)
option(USE_STATIC "Use static libraries when available" OFF)
option(INSTALL_MAC_BUNDLE "On OSX create QtPisa bundle" ON)

if(USE_STATIC)
  # dynamic libraries are allowed if static are not found,
  # because some Qt dependencies, like libX11, are normally dynamic.
  # On some systems librt.a can cause problems. Workaround is to remove it.
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

find_package(CCP4 REQUIRED mmdb2 ssm ccp4srs)
include_directories("${CMAKE_SOURCE_DIR}" "${CCP4_INCLUDE_DIRS}")

if(USE_QT)
  find_package(Qt4 REQUIRED QtCore QtGui QtNetwork QtWebKit QtXml QtXmlPatterns)
  include(${QT_USE_FILE})

  find_path(QTLIBS_INCLUDE qtx/qtx_dnlmanager.h
            HINTS ${CMAKE_INSTALL_PREFIX}/include)
  find_library(QTX_LIBRARY NAMES qtx
               HINTS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib)
  find_library(FGL_LIBRARY NAMES fgl
               HINTS ${LIB_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX}/lib)

  include_directories(${QTLIBS_INCLUDE})
  QT4_ADD_RESOURCES(RCS qtpisa/qtpisa.qrc)
endif()

set(PISALIB_SRC
	pisalib/chem_equilibrium.cpp
	pisalib/pisa_data.cpp
	pisalib/pisa_engunits.cpp
	pisalib/pisa_rcsbdat.cpp
	pisalib/pisa_asmstock.cpp
	pisalib/pisa_datastats.cpp
	pisalib/pisa_interface.cpp
	pisalib/pisa_symnumber.cpp
	pisalib/pisa_asmunit.cpp
	pisalib/pisa_detail.cpp
	pisalib/pisa_lists.cpp
	pisalib/pisa_types.cpp
	pisalib/pisa_assemblies.cpp
	pisalib/pisa_domain.cpp
	pisalib/pisa_molref.cpp
	pisalib/pisa_assembly.cpp
	pisalib/pisa_domains.cpp
	pisalib/pisa_prosurf.cpp
	pisalib/pisa_crystsplit.cpp
	pisalib/pisa_engine.cpp
	pisalib/pisa_view.cpp
	pisalib/pisa_query.cpp)

add_library(pisalib STATIC ${PISALIB_SRC})
#target_link_libraries(pisalib ${CCP4_LIBRARIES})

set(PISA_SRC
    pisa/pisa_analyse.cpp
    pisa/pisa.cpp)

add_executable(pisa ${PISA_SRC})
target_link_libraries(pisa pisalib ${CCP4_LIBRARIES})


set(JSPISA_SRC
    jspisa/jspisa.cpp
    jspisa/jspisa_init.cpp
    jspisa/jspisa_monomers.cpp
    jspisa/jspisa_splits.cpp
    jspisa/jspisa_base.cpp
    jspisa/jspisa_interfaces.cpp
    jspisa/jspisa_process.cpp
    jspisa/jspisa_stock.cpp)

add_executable(jspisa ${JSPISA_SRC})
target_link_libraries(jspisa pisalib rvapi ${CCP4_LIBRARIES})

if(USE_QT)

set(QTPISA_SRC
	qtpisa/qtpisa.cpp
	qtpisa/qtpisa_helppage.cpp
	qtpisa/qtpisa_montab.cpp
	qtpisa/qtpisa_viewer.cpp
	qtpisa/qtpisa_asmdtab.cpp
	qtpisa/qtpisa_intfd1tab.cpp
	qtpisa/qtpisa_msgtab.cpp
	qtpisa/qtpisa_viewsignal.cpp
	qtpisa/qtpisa_asmtab.cpp
	qtpisa/qtpisa_intfd2tab.cpp
	qtpisa/qtpisa_prefdialog.cpp
	qtpisa/qtpisa_widget.cpp
	qtpisa/qtpisa_controlpanel.cpp
	qtpisa/qtpisa_intfd3tab.cpp
	qtpisa/qtpisa_stocktab.cpp
	qtpisa/qtpisa_widgetbase.cpp
	qtpisa/qtpisa_core.cpp
	qtpisa/qtpisa_intftab.cpp
	qtpisa/qtpisa_tab.cpp
	qtpisa/qtpisa_datatab.cpp
	qtpisa/qtpisa_mainwindow.cpp
	qtpisa/qtpisa_thread.cpp
	qtpisa/qtpisa_dnldialog.cpp
	qtpisa/qtpisa_mondtab.cpp
	qtpisa/qtpisa_types.cpp)

set(MOC_HEADERS
	qtpisa/qtpisa_asmdtab.h
	qtpisa/qtpisa_helppage.h
	qtpisa/qtpisa_montab.h
	qtpisa/qtpisa_viewer.h
	qtpisa/qtpisa_asmtab.h
	qtpisa/qtpisa_intfd1tab.h
	qtpisa/qtpisa_msgtab.h
	qtpisa/qtpisa_viewsignal.h
	qtpisa/qtpisa_controlpanel.h
	qtpisa/qtpisa_intfd2tab.h
	qtpisa/qtpisa_prefdialog.h
	qtpisa/qtpisa_widget.h
	qtpisa/qtpisa_core.h
	qtpisa/qtpisa_intfd3tab.h
	qtpisa/qtpisa_stocktab.h
	qtpisa/qtpisa_widgetbase.h
	qtpisa/qtpisa_datatab.h
	qtpisa/qtpisa_intftab.h
	qtpisa/qtpisa_tab.h
	qtpisa/qtpisa_mainwindow.h
	qtpisa/qtpisa_thread.h
	qtpisa/qtpisa_dnldialog.h
	qtpisa/qtpisa_mondtab.h)
QT4_WRAP_CPP(MOC_SOURCES ${MOC_HEADERS})

# qtpisa (binary program) in the default location, i.e. top-level build dir,
# would cause name conflict with qtpisa/ directory that contains MOC files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(qtpisa ${QTPISA_SRC} ${MOC_SOURCES} ${RCS} )

target_link_libraries(qtpisa pisalib ${CCP4_LIBRARIES}
                      ${FGL_LIBRARY} ${QTX_LIBRARY} ${QT_LIBRARIES})

install(TARGETS qtpisa DESTINATION bin)

if(APPLE AND INSTALL_MAC_BUNDLE)
  install(DIRECTORY DESTINATION qtpisa.app/Contents/MacOS)
  install(DIRECTORY DESTINATION qtpisa.app/Contents/Resources)
  install(FILES qtpisa/images/qtpisa.icns
          DESTINATION qtpisa.app/Contents/Resources)
  install(FILES qtpisa/Info.plist
          DESTINATION qtpisa.app/Contents)
  install(CODE "
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                      \"../../../bin/qtpisa\"
                      \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/qtpisa.app/Contents/MacOS/qtpisa\")
  ")
endif()

endif(USE_QT)

install(TARGETS pisa DESTINATION bin)
install(TARGETS jspisa DESTINATION bin)
install(FILES pisa.cfg molref/molref.idx molref/molref.rdt
              jspisa.cfg
              pisastore/agents.dat pisastore/asm_params.dat pisastore/intfstats.dat
              pisastore/rcsb_symops.dat pisastore/syminfo_pisa.lib
        DESTINATION share/pisa)
install(DIRECTORY help/ DESTINATION html/pisa)
