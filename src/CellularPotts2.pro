TEMPLATE = app 
GRAPHICS = qt
CONFIG += console 
CONFIG += release
QT += widgets
QT += gui
CONFIG -= debug
CONFIG -= app_bundle
OBJECTS_DIR = model
QMAKE_DISTCLEAN += -r data_film org-data


contains( GRAPHICS, qt ) {
  
}	

TARGET = multisort
MAINFILE = $$join(TARGET, " ", , ".cpp" )
PARAMFILE = parameter-files/parameter_$${TARGET}.cpp

message( $$MAINFILE )
message( $$TARGET )
# Input
HEADERS += model/ca.h \
	   model/hull.h \
           model/cell.h \
           model/conrec.h \
           model/dish.h \
           model/graph.h \
           model/info.h \
           model/misc.h \
           model/output.h \
           model/parameter.h \
           model/parse.h \
           model/pde.h \
           model/random.h \
           model/sqr.h \
           model/sticky.h \
	   model/crash.h \
	   model/warning.h \ 
	   model/storage.h \
	   model/fft.h \
      model/connections.h

        
SOURCES += model/ca.cpp \
	   model/hull.cpp \
           model/cell.cpp \
           model/conrec.cpp \
           model/dish.cpp \
           model/info.cpp \
           model/misc.cpp \
           model/output.cpp \
           model/parse.cpp \
           model/pde.cpp \
           model/random.cpp \
           model/crash.cpp \
           model/warning.cpp \
           model/storage.cpp \
	   model/fft.cpp \	
           model/connections.cpp

SOURCES += $$MAINFILE
SOURCES += $$PARAMFILE
       
#QMAKE_CXXFLAGS_RELEASE += -fexceptions
#QMAKE_CXXFLAGS_DEBUG += -fexceptions
#QMAKE_LFLAGS_RELEASE += -O4
#QMAKE_CXXFLAGS_RELEASE += -O4

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS +=  -fopenmp
LIBS += -fopenmp

contains( GRAPHICS, qt ) {
   message( "Building Qt executable" )
   SOURCES += model/qtgraph.cpp
   HEADERS += model/qtgraph.h
   QMAKE_CXXFLAGS_RELEASE += -DQTGRAPHICS
   QMAKE_CXXFLAGS_DEBUG += -DQTGRAPHICS 
#   QT += qt3support
   unix {
      system(rm $$TARGET.o)
   } 
   win32 {
     QMAKE_LFLAGS += -L "\"C:\Program Files\GnuWin32\lib\"" -lpng -lzdll
     QMAKE_CXXFLAGS += -I "\"C:\Program Files\GnuWin32\include\""
   }
   #LIBS += -lpng -fopenmp
}

contains( GRAPHICS, qt3 ) {
   message( "Building Qt executable" )
   SOURCES += model/qt3graph.cpp
   HEADERS += model/qt3graph.h
   QMAKE_CXXFLAGS_RELEASE += -DQTGRAPHICS
   QMAKE_CXXFLAGS_DEBUG += -DQTGRAPHICS 
   unix {
      system(rm vessel.o)
   } 
   win32 {
     QMAKE_LFLAGS += -L "C:\Program Files\GnuWin32\lib" -lpng -lzdll
     QMAKE_CXXFLAGS += -I "C:\Program Files\GnuWin32\include"
   }
   LIBS += -lpng
}
contains( GRAPHICS, x11 ) {
   !unix {
     error("X11 graphics only available on Unix systems.")
   }
   message("Building X11 executable")
   SOURCES += model/x11graph.cpp
   HEADERS += model/x11graph.h
   QMAKE_CXXFLAGS_RELEASE += -DX11GRAPHICS
   QMAKE_CXXFLAGS_DEBUG += -DX11GRAPHICS 
   unix {
      system(rm vessel.o)
   }
   CONFIG -= qt
   CONFIG += x11
   unix:LIBS += -lpng
}


#The following line was inserted by qt3to4
QT +=  
