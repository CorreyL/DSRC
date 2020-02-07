#include "../include/dsrc/DsrcInMemory.h"

#include "DsrcFile.h"
#include "DsrcIo.h"
#include "Buffer.h"
#include "Fastq.h"

#include <iostream>

namespace dsrc{
  namespace comp {
    using namespace core;
    using namespace fq;

    /**
     * Creates a pointer to the given .dsrc file to read the buffer contents
     * in chunks
     *
     * Allocates memory for dsrcChunk* for the buffer to be loaded into memory,
     * allocates memory for fastqChunk* for the decompressed buffer contents to
     * be loaded into memory
     */
    DsrcInMemory::DsrcInMemory(const std::string& dsrcFilename_) {
      try {
        reader = new DsrcFileReader();
        reader->StartDecompress(dsrcFilename_);
        dsrcChunk = new DsrcDataChunk(DsrcDataChunk::DefaultBufferSize);
        fastqChunk = new FastqDataChunk(FastqDataChunk::DefaultBufferSize);
      } catch (const DsrcException& e_) {
        /**
         * @todo Implement IsError() functionality? Extend from IDsrcOperator?
         */
        std::cerr << e_.what() << std::endl;
        exit(1);
      }
    }

    /**
     * Removes all contents from memory allocated for dsrcChunk and fastqChunk,
     * closes the decompression fileStream, and frees up all memory allocated
     * for each pointer
     */
    DsrcInMemory::~DsrcInMemory() {
      dsrcChunk->Reset();
      fastqChunk->Reset();
      reader->FinishDecompress();
      delete dsrcChunk;
      delete fastqChunk;
      delete reader;
    }
  }
}
