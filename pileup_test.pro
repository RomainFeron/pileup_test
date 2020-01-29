TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -lz -llzma -lbz2 -L../include/htslib/libhts.a

INCLUDEPATH += include/

SOURCES += \
    src/main.cpp

DISTFILES += \

HEADERS += \
