#ifndef H_DSRCINMEMORY
#define H_DSRCINMEMORY

#include "../../src/DsrcFile.h"
#include "../../src/DsrcIo.h"
#include "../../src/Fastq.h"

namespace dsrc{
  namespace ext {

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
        bool IsError() const;
        void AddError(const std::string& err_);
        void ClearError();

      private:
        comp::DsrcFileReader* reader = NULL;
        comp::DsrcDataChunk* dsrcChunk = NULL;
        fq::FastqDataChunk* fastqChunk = NULL;
        std::string errorMsg = "";
    };
  }
}

#endif
