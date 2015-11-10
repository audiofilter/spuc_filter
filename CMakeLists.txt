CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "DEBUG")
ENDIF()
PROJECT("plotdes")
Message("Project spuce plotdes")

set(CMAKE_CXX_FLAGS " -std=c++11 ")
set(CMAKE_PREFIX_PATH "/usr/local/opt/qt5")
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

FIND_PACKAGE(Qt5 REQUIRED Gui Core Widgets)


ADD_DEFINITIONS("-Wall -fPIC")
INCLUDE_DIRECTORIES(".." "/usr/local/opt/qt5/include")
LINK_DIRECTORIES("/usr/local/lib")

set(EXECUTABLE_OUTPUT_PATH ..)

set(SOURCES
	make_filter.cpp 
	des_filter.cpp 
)

SET(FILT_SRCS  
	spuce/filters/butterworth_iir.cpp
	spuce/filters/chebyshev_iir.cpp
	spuce/filters/chebyshev2_iir.cpp
	spuce/filters/elliptic_iir.cpp
	spuce/filters/find_roots.cpp
	spuce/filters/fir_coeff.cpp
	spuce/filters/iir_coeff.cpp
	# FIR stuff
	spuce/filters/gaussian_fir.cpp
	spuce/filters/butterworth_fir.cpp
	spuce/filters/raised_cosine_imp.cpp
	spuce/filters/root_raised_cosine_imp.cpp
	spuce/filters/remez_fir.cpp
	spuce/filters/sinc_fir.cpp
	spuce/filters/create_remez_lpfir.cpp
	spuce/filters/fir_inv_dft.cpp
	spuce/filters/transform_fir.cpp
	# Windows
	spuce/filters/window.cpp
	# Other spuce/filters
	spuce/filters/butterworth_allpass.cpp
	spuce/filters/elliptic_allpass.cpp
  spuce/filters/shelf_allpass1.cpp
)

SET(plotdes_SOURCES main.cpp mainwindow.cpp other_freq.cpp qcustomplot.cpp)
SET(plotdes_HEADERS mainwindow.h other_freq.h qcustomplot.h)

INCLUDE_DIRECTORIES(${QT5_INCLUDES})
ADD_DEFINITIONS(${QT5Widgets_DEFINITIONS})

ADD_LIBRARY(filtlib ${SOURCES} ${FILT_SRCS})
ADD_EXECUTABLE(plotdes ${plotdes_SOURCES} ${plotdes_HEADERS_MOC})
TARGET_LINK_LIBRARIES(plotdes filtlib ${QT_LIBRARIES} ${SPUCEFILT_LIBRARY})
QT5_USE_Modules(plotdes Gui)
QT5_USE_Modules(plotdes Core)
QT5_USE_Modules(plotdes Widgets)
QT5_USE_Modules(plotdes PrintSupport)

