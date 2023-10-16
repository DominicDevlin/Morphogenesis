TEMPLATE = app 
GRAPHICS = qt
CONFIG += console 
CONFIG += release
QT += widgets
QT += gui
CONFIG -= debug
CONFIG -= app_bundle


contains( GRAPHICS, qt ) {
  
}	

TARGET = overlap
MAINFILE = $$join(TARGET, " ", , ".cpp" )

message( $$MAINFILE )
message( $$TARGET )
# Input
HEADERS += ca.h \
	   hull.h \
           cell.h \
           conrec.h \
           dish.h \
           graph.h \
           info.h \
           misc.h \
           output.h \
           parameter.h \
           parse.h \
           pde.h \
           random.h \
           sqr.h \
           sticky.h \
       	   crash.h \
	   warning.h \ 
	   storage.h \
	   fft.h \
      connections.h

        
SOURCES += ca.cpp \
	   hull.cpp \
           cell.cpp \
           conrec.cpp \
           dish.cpp \
           info.cpp \
           misc.cpp \
           output.cpp \
           parameter.cpp \
           parse.cpp \
           pde.cpp \
           random.cpp \
           crash.cpp \
           warning.cpp \
           storage.cpp \
	   fft.cpp \	
           connections.cpp

SOURCES += $$MAINFILE
       
#QMAKE_CXXFLAGS_RELEASE += -fexceptions
#QMAKE_CXXFLAGS_DEBUG += -fexceptions
#QMAKE_LFLAGS_RELEASE += -O4
#QMAKE_CXXFLAGS_RELEASE += -O4

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS +=  -fopenmp
LIBS += -fopenmp

contains( GRAPHICS, qt ) {
   message( "Building Qt executable" )
   SOURCES += qtgraph.cpp
   HEADERS += qtgraph.h
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
   SOURCES += qt3graph.cpp
   HEADERS += qt3graph.h
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
   SOURCES += x11graph.cpp
   HEADERS += x11graph.h
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
