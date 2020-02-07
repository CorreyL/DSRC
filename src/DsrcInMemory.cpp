#include "../include/dsrc/DsrcInMemory.h"

#include "BitMemory.h"
#include "BlockCompressor.h"
#include "Common.h"
#include "DsrcFile.h"
#include "DsrcIo.h"
#include "Buffer.h"
#include "Fastq.h"

#include <iostream>

namespace dsrc{
  namespace comp {
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
        reader = new comp::DsrcFileReader();
        reader->StartDecompress(dsrcFilename_);
        dsrcChunk = new comp::DsrcDataChunk(comp::DsrcDataChunk::DefaultBufferSize);
        fastqChunk = new fq::FastqDataChunk(fq::FastqDataChunk::DefaultBufferSize);
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

    /**
     * Retrieves the next chunk of byte code from the `.dsrc` file and loads it
     * into a buffer, converts the buffer and returns the content as a string
     *
     * This function must be called multiple times in order to obtain whole
     * content of the input file
     */
    std::string DsrcInMemory::getNextChunk(){
      comp::BlockCompressor superblock(
        reader->GetDatasetType(),
        reader->GetCompressionSettings()
      );
      if (reader->ReadNextChunk(dsrcChunk)) {
        core::BitMemoryReader bitMemory(
          dsrcChunk->data.Pointer(),
          dsrcChunk->size
        );
        superblock.Read(bitMemory, *fastqChunk);
        std::string chunkContents = std::string(
          reinterpret_cast<char const*>(fastqChunk->data.Pointer()),
          fastqChunk->size
        );
        dsrcChunk->Reset();
        fastqChunk->Reset();
        return chunkContents;
      }
      return "";
    }
  }
}
