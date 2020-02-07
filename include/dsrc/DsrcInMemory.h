#ifndef H_DSRCINMEMORY
#define H_DSRCINMEMORY

#include "../../src/DsrcFile.h"
#include "../../src/DsrcIo.h"
#include "../../src/Fastq.h"

namespace dsrc{
  namespace comp {
    using namespace core;
    using namespace fq;

    /**
     * Utilizes the implementation of DsrcFileReader to read the binary contents
     * of a .dsrc file into memory and return the uncompressed FASTQ content as
     * a string
     *
     * This allows the contents of a .dsrc file to be uncompressed without
     * needing to write the contents to a file, which can be beneficial when
     * concatenating multiple .dsrc files into a single concatenated .fastq file
     * with limited hard disk space (i.e. no longer needing to uncompress each
     * .dsrc file into temporary .fastq files before concatenation)
     */
    class DsrcInMemory {
      public:
        DsrcInMemory(const std::string& dsrcFilename_);
        ~DsrcInMemory();
        std::string getNextChunk();

      private:
        DsrcFileReader* reader = NULL;
        DsrcDataChunk* dsrcChunk = NULL;
        FastqDataChunk* fastqChunk = NULL;
    };
  }
}

#endif
