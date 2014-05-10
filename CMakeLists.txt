CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "DEBUG")
ENDIF()
PROJECT("plotdes")
Message("Project spuc plotdes")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

FIND_PACKAGE(Qt5Widgets)
FIND_PACKAGE(Qt5Core)
FIND_PACKAGE(Qt5Gui) 

set(SPUCDIR ./spuc/spuc)

ADD_DEFINITIONS("-Wall -fPIC -DSPUCBOOST")
INCLUDE_DIRECTORIES(".." "/usr/local/opt/qt5/include" "${SPUCDIR}/..")
LINK_DIRECTORIES("/usr/local/lib")

set(EXECUTABLE_OUTPUT_PATH ..)

# not sure if needed if libfilt/spuc is available!
set(SOURCES
${SPUCDIR}/root_raised_cosine_imp.cpp 
${SPUCDIR}/shelf_allpass1.cpp 
${SPUCDIR}/iir_coeff.cpp 
${SPUCDIR}/gaussian_fir.cpp 
${SPUCDIR}/fir_inv_dft.cpp 
${SPUCDIR}/elliptic_iir.cpp 
${SPUCDIR}/elliptic_allpass.cpp 
${SPUCDIR}/chebyshev_iir.cpp 
${SPUCDIR}/butterworth_iir.cpp 
${SPUCDIR}/butterworth_fir.cpp 
${SPUCDIR}/butterworth_allpass.cpp 
${SPUCDIR}/root_raised_cosine_quantized.cpp 
${SPUCDIR}/spuc_math.cpp 
${SPUCDIR}/find_roots.cpp 
${SPUCDIR}/create_remez_lpfir.cpp 
${SPUCDIR}/create_remez_fir.cpp 
${SPUCDIR}/fir_coeff.cpp 
${SPUCDIR}/notch_example.cpp 
${SPUCDIR}/remez_fir.cpp 
${SPUCDIR}/raised_cosine_imp.cpp 
make_filter.cpp 
des_filter.cpp 
)

SET(plotdes_SOURCES main.cpp mainwindow.cpp other_freq.cpp qcustomplot.cpp)
SET(plotdes_HEADERS mainwindow.h other_freq.h qcustomplot.h)

INCLUDE_DIRECTORIES(${QT5_INCLUDES})
ADD_DEFINITIONS(${QT5Widgets_DEFINITIONS})

ADD_LIBRARY(filtlib ${SOURCES})
ADD_EXECUTABLE(plotdes ${plotdes_SOURCES} ${plotdes_HEADERS_MOC})
TARGET_LINK_LIBRARIES(plotdes filtlib ${QT_LIBRARIES})
QT5_USE_Modules(plotdes Gui)
QT5_USE_Modules(plotdes Core)
QT5_USE_Modules(plotdes Widgets)
QT5_USE_Modules(plotdes PrintSupport)

