/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
  
#include <boost/python/detail/wrap_python.hpp>		// must be included before <Python>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/overloads.hpp>

#include <Python.h>
#include <string>
#include <sstream>
#include <vector>

#include "../include/dsrc/Globals.h"
#include "../include/dsrc/FastqRecord.h"
#include "../include/dsrc/FastqFile.h"
#include "../include/dsrc/DsrcModule.h"
#include "../include/dsrc/DsrcArchive.h"
#include "../include/dsrc/DsrcInMemory.h"

namespace dsrc
{

namespace py
{

namespace boo = boost::python;
using namespace ext;


// Exception handling
//
class PyException : public DsrcException
{
public:
	PyException(const char* msg_)
		: DsrcException(msg_)
	{}

	PyException(const std::string& msg_)
		: DsrcException(msg_)
	{}

	static void translate(const PyException& e_)		// boost::python requirement
	{
		PyErr_SetString(PyExc_RuntimeError, e_.what());
	}
};

template <class _T>
void TCheckError(_T& obj_)
{
	if (obj_.IsError())
	{
		std::string err = obj_.GetError();
		obj_.ClearError();
		throw PyException(err);
	}
}

template <class _T>
void TCheckPtrError(_T*& obj_)
{
	if (obj_->IsError())
	{
		std::string err = obj_->GetError();
		obj_->ClearError();
		throw PyException(err);
	}
}

// Compression settings
//
struct PyDsrcCompressionSettings : public DsrcCompressionSettings
{
	int GetDnaCompressionLevel() const
	{
		return dnaCompressionLevel;
	}

	void SetDnaCompressionLevel(int v_)
	{
		if (v_ < (int)DsrcCompressionSettings::MinDnaCompressionLevel || v_ > (int)DsrcCompressionSettings::MaxDnaCompressionLevel)
			throw PyException("Invalid DNA compression level specified");
		dnaCompressionLevel = v_;
	}

	int GetQualityCompressionLevel() const
	{
		return qualityCompressionLevel;
	}

	void SetQualityCompressionLevel(int v_)
	{
		if (v_ < (int)DsrcCompressionSettings::MinQualityCompressionLevel || v_ > (int)DsrcCompressionSettings::MaxQualityCompressionLevel)
			throw PyException("Invalid quality compression level specified");
		qualityCompressionLevel = v_;
	}

	bool IsLossyQualityCompression() const
	{
		return lossyQualityCompression;
	}

	void SetLossyQualityCompression(bool v_)
	{
		lossyQualityCompression = v_;
	}

	long int GetTagFieldPreserveMask() const
	{
		return tagPreserveMask;
	}

	void SetTagFieldPreserveMask(long int v_)
	{
		tagPreserveMask = v_;
	}

	int GetFastqBufferSizeMb() const
	{
		return fastqBufferSizeMb;
	}

	void SetFastqBufferSizeMb(int v_)
	{
		if (v_ < (int)DsrcCompressionSettings::MinFastqBufferSizeMB || v_ > (int)DsrcCompressionSettings::MaxFastqBufferSizeMB)
			throw PyException("Invalid FASTQ buffer size specified");
		fastqBufferSizeMb = v_;
	}

	bool IsCrc32Checking() const
	{
		return calculateCrc32;
	}

	void SetCrc32Checking(bool v_)
	{
		calculateCrc32 = v_;
	}
};


// FASTQ interfaces
//
class PyFastqFileRecordsReader
{
public:
	void Open(const std::string &filename_)
	{
		if (!reader.Open(filename_))
			throw PyException("Cannot open file: " + filename_);
	}

	void Close()
	{
		reader.Close();
	}

	bool ReadNextRecord(FastqRecord& rec_)
	{
		return reader.ReadNextRecord(rec_);
	}

private:
	FastqFileRecordsReader reader;
};


class PyFastqFileRecordsWriter
{
public:
	void Open(const std::string &filename_)
	{
		if (!writer.Open(filename_))
			throw PyException("Cannot open file: " + filename_);
	}

	void Close()
	{
		writer.Close();
	}

	void WriteNextRecord(const FastqRecord& rec_)
	{
		writer.WriteNextRecord(rec_);
	}

private:
	FastqFileRecordsWriter writer;
};


// DSRC module interfaces
//
class PyDsrcModule
{
public:
	void Compress(const std::string& inFastqFilename_,
				  const std::string& outDsrcFilename_,
				  const PyDsrcCompressionSettings& compSettings_,
				  uint32 threadsNum_,
				  bool useFastqStdIo_ = false,
				  uint32 qualityOffset_ = 0)
	{
		dsrc.Compress(inFastqFilename_,
					  outDsrcFilename_,
					  compSettings_,
					  threadsNum_,
					  useFastqStdIo_,
					  qualityOffset_);
		TCheckError(dsrc);
	}

	void Decompress(const std::string& inDsrcFilename_,
					const std::string& outFastqFilename_,
					uint32 threadsNum_,
					bool useFastqStdIo_ = false)
	{
		dsrc.Decompress(inDsrcFilename_,
						outFastqFilename_,
						threadsNum_,
						useFastqStdIo_);
		TCheckError(dsrc);
	}

	// we need to overload member functions with default arguments for Python <-> C++ API compatibility
	BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Compress_overload, Compress, 4, 6)
	BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Decompress_overload, Decompress, 3, 4)

private:
	DsrcModule dsrc;
};


// DSRC archive interfaces
//
class PyDsrcArchiveRecordsWriter
{
public:
	void StartCompress(const std::string& dsrcFilename_,
					   const PyDsrcCompressionSettings& compressionSettings_,
					   uint32 qualityOffset_ = 0)
	{
		writer.StartCompress(dsrcFilename_,
							 compressionSettings_,
							 1,						// threads num -- atm only 1
							 qualityOffset_);
		TCheckError(writer);
	}

	void FinishCompress()
	{
		writer.FinishCompress();
		TCheckError(writer);
	}

	void WriteNextRecord(const FastqRecord& rec_)
	{
		writer.WriteNextRecord(rec_);
		TCheckError(writer);
	}

	BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(StartCompress_overload, StartCompress, 2, 3)

private:
	DsrcArchiveRecordsWriter writer;
};


class PyDsrcArchiveRecordsReader
{
public:
	void StartDecompress(const std::string& dsrcFilename_)
	{
		reader.StartDecompress(dsrcFilename_, 1);		// threads num -- atm only 1
		TCheckError(reader);
	}

	void FinishDecompress()
	{
		reader.FinishDecompress();
		TCheckError(reader);
	}

	bool ReadNextRecord(FastqRecord& rec_)
	{
		if (!reader.ReadNextRecord(rec_))
		{
			TCheckError(reader);
			return false;
		}
		return true;
	}

private:
	DsrcArchiveRecordsReader reader;
};

class PyDsrcReadInMemory {
	public:
		void Open(const std::string& inDsrcFilename_) {
			dsrcInMemory = new DsrcInMemory(inDsrcFilename_);
			TCheckPtrError(dsrcInMemory);
		}

		void Close() {
			delete dsrcInMemory;
			TCheckPtrError(dsrcInMemory);
			dsrcInMemory = NULL;
		}

		/**
		 * Returns the next line of the opened .dsrc file
		 *
		 * Retrieves chunks of the .dsrc file by calling ReadNextChunk, and storing
		 * each line of the chunk in a Queue
		 *
		 * The Queue is accessed for the next line until the Queue is empty, at
		 * which point the next chunk is retrieved
		 */
		std::string readline() {
			if (chunk.empty()) {
				std::string lines = ReadNextChunk();
				std::vector<std::string> split_lines = split(lines, '\n');
				if (split_lines.capacity() == 0) {
					// End of file has been reached
					return lines;
				}
				for (
					std::vector<std::string>::iterator it = split_lines.begin();
					it != split_lines.end();
					++it
				) {
					chunk.push(*it + "\n");
				}
			}
			std::string line = chunk.front();
			chunk.pop();
			return line;
		}

		bool closed() {
			return dsrcInMemory == NULL;
		}

		PyDsrcReadInMemory __enter__() {
			return *this;
		}

		/**
		 * Arguments are typed as `void*` to simulate the `NoneType` type in Python,
		 * ensuring argument types match between C++ and Python when the function
		 * is invoked
		 */
		void __exit__no_error(
			void* exceptionType,
			void* exceptionValue,
			void* traceback
		) {
			Close();
		}

		void __exit__error(
			PyObject* exceptionType,
			PyObject* exceptionValue,
			PyObject* traceback
		) {
			Close();
		}


	private:
		DsrcInMemory * dsrcInMemory = NULL;
		std::queue<std::string> chunk;

		/**
		 * C++ does not have a built-in std::string.split() function that tokenizes
		 * strings on a given delimiter
		 *
		 * This implementation is based off of a Stack Overflow implementation:
		 *
		 * https://stackoverflow.com/a/236803/10491481
		 */
		template <typename Out>
		void split(const std::string &s, char delim, Out result) {
			std::istringstream iss(s);
			std::string item;
			while (std::getline(iss, item, delim)) {
				*result++ = item;
			}
		}

		std::vector<std::string> split(const std::string &s, char delim) {
			std::vector<std::string> elems;
			split(s, delim, std::back_inserter(elems));
			return elems;
		}

		std::string ReadNextChunk() {
			std::string chunk = dsrcInMemory->getNextChunk();
			TCheckPtrError(dsrcInMemory);
			return chunk;
		}
};


BOOST_PYTHON_MODULE(pydsrc)
{
	boo::register_exception_translator<PyException>(&PyException::translate);

	boo::class_<FastqRecord>("FastqRecord")
		.def_readwrite("tag", &FastqRecord::tag)
		.def_readwrite("sequence", &FastqRecord::sequence)
		.def_readwrite("plus", &FastqRecord::plus)
		.def_readwrite("quality", &FastqRecord::quality)
	;

	boo::class_<PyFastqFileRecordsReader, boost::noncopyable>("FastqFileRecordsReader")
		.def("Open", &PyFastqFileRecordsReader::Open)
		.def("Close", &PyFastqFileRecordsReader::Close)
		.def("ReadNextRecord", &PyFastqFileRecordsReader::ReadNextRecord)
	;

	boo::class_<PyFastqFileRecordsWriter, boost::noncopyable>("FastqFileRecordsWriter")
		.def("Open", &PyFastqFileRecordsWriter::Open)
		.def("Close", &PyFastqFileRecordsWriter::Close)
		.def("WriteNextRecord", &PyFastqFileRecordsWriter::WriteNextRecord)
	;

	boo::class_<PyDsrcCompressionSettings>("CompressionSettings")
		.add_property("DNACompressionLevel", &PyDsrcCompressionSettings::GetDnaCompressionLevel, &PyDsrcCompressionSettings::SetDnaCompressionLevel)
		.add_property("QualityCompressionLevel", &PyDsrcCompressionSettings::GetQualityCompressionLevel, &PyDsrcCompressionSettings::SetDnaCompressionLevel)
		.add_property("LossyCompression", &PyDsrcCompressionSettings::IsLossyQualityCompression, &PyDsrcCompressionSettings::SetLossyQualityCompression)
		.add_property("TagFieldFilterMask", &PyDsrcCompressionSettings::GetTagFieldPreserveMask, &PyDsrcCompressionSettings::SetTagFieldPreserveMask)
		.add_property("FastqBufferSizeMB", &PyDsrcCompressionSettings::GetFastqBufferSizeMb, &PyDsrcCompressionSettings::SetFastqBufferSizeMb)
		.add_property("Crc32Checking", &PyDsrcCompressionSettings::IsCrc32Checking, &PyDsrcCompressionSettings::SetCrc32Checking)
	;

	boo::class_<FieldMask>("FieldMask")
		.def("AddField", &FieldMask::AddField)
		.def("GetMask", &FieldMask::GetMask)
	;

	boo::class_<PyDsrcModule, boost::noncopyable>("DsrcModule")
		.def("Compress", &PyDsrcModule::Compress, PyDsrcModule::Compress_overload())
		.def("Decompress", &PyDsrcModule::Decompress, PyDsrcModule::Decompress_overload())
	;

	boo::class_<PyDsrcArchiveRecordsWriter, boost::noncopyable>("DsrcArchiveRecordsWriter")
		.def("StartCompress", &PyDsrcArchiveRecordsWriter::StartCompress, PyDsrcArchiveRecordsWriter::StartCompress_overload())
		.def("WriteNextRecord", &PyDsrcArchiveRecordsWriter::WriteNextRecord)
		.def("FinishCompress", &PyDsrcArchiveRecordsWriter::FinishCompress)
	;

	boo::class_<PyDsrcArchiveRecordsReader, boost::noncopyable>("DsrcArchiveRecordsReader")
		.def("StartDecompress", &PyDsrcArchiveRecordsReader::StartDecompress)
		.def("ReadNextRecord", &PyDsrcArchiveRecordsReader::ReadNextRecord)
		.def("FinishDecompress", &PyDsrcArchiveRecordsReader::FinishDecompress)
	;

	boo::class_<PyDsrcReadInMemory>("DsrcReadInMemory")
		.def("open", &PyDsrcReadInMemory::Open)
		.def("readline", &PyDsrcReadInMemory::readline)
		.def("close", &PyDsrcReadInMemory::Close)
		.def("closed", &PyDsrcReadInMemory::closed)
		.def("__enter__", &PyDsrcReadInMemory::__enter__)
		.def("__exit__", &PyDsrcReadInMemory::__exit__no_error)
		.def("__exit__", &PyDsrcReadInMemory::__exit__error)
	;
}

} // namespace py

} // namespace dsrc
