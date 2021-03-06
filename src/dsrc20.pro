TEMPLATE = app
CONFIG += console
CONFIG -= qt

QMAKE_CXXFLAGS += -Wall

# select between boost and c++11
QMAKE_CXXFLAGS += -DUSE_BOOST_THREAD
# QMAKE_CXXFLAGS += -std=c++11

LIBS += -lpthread
LIBS += -lboost_thread
LIBS += -lboost_system

INCLUDEPATH += /usr/include/python2.7

SOURCES += \
    huffman.cpp \
    DsrcFile.cpp \
    FileStream.cpp \
    QualityPositionModeler.cpp \
    QualityRLEModeler.cpp \
    DnaModelerHuffman.cpp \
    RecordsProcessor.cpp \
    TagModeler.cpp \
    BlockCompressor.cpp \
    FastqParser.cpp \
    FastqIo.cpp \
    FastqStream.cpp \
    StdStream.cpp \
    DsrcWorker.cpp \
    DsrcOperator.cpp \
    DsrcIo.cpp \
    RecordsBlockCompressor.cpp

SOURCES += \
    DsrcArchive.cpp \
    FastqFile.cpp \
    DsrcModule.cpp

#SOURCES += ../py/Interface.cpp
#LIBS += -lboost_python -lpython2.7


SOURCES += main.cpp


HEADERS += \
    utils.h \
    huffman.h \
    FileStream.h \
    FastqIo.h \
    DsrcFile.h \
    DataPool.h \
    DsrcIo.h \
    DataQueue.h \
    Buffer.h \
    BitMemory.h \
    FastqParser.h \
    RangeCoder.h \
    QualityModeler.h \
    DnaModeler.h \
    DnaModelerBasicB2.h \
    DnaModelerRCO.h \
    SymbolCoderRC.h \
    QualityPositionModeler.h \
    QualityRLEModeler.h \
    DnaModelerHuffman.h \
    RecordsProcessor.h \
    QualityModelerProxy.h \
    DnaModelerProxy.h \
    Stats.h \
    TagModeler.h \
    BlockCompressor.h \
    Fastq.h \
    QualityOrderModeler.h \
    FastqStream.h \
    QualityEncoder.h \
    DataStream.h \
    StdStream.h \
    DsrcWorker.h \
    DsrcOperator.h \
    Common.h \
    Crc32.h \
    ErrorHandler.h \
    RecordsBlockCompressor.h

HEADERS += \
    ../include/dsrc/DsrcArchive.h \
    ../include/dsrc/FastqFile.h \
    ../include/dsrc/DsrcModule.h \
    ../include/dsrc/Dsrc.h \
    ../include/dsrc/Globals.h \
    ../include/dsrc/FastqRecord.h
