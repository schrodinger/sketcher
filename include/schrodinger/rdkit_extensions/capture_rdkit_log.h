/* -------------------------------------------------------------------------
 * Class to capture RDKit logging
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <ostream>
#include <string>
#include <sstream>

#include <RDGeneral/RDLog.h>
#include <boost/core/noncopyable.hpp>

#include "schrodinger/rdkit_extensions/definitions.h"

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * Captures all messages issued to RDKit error logging
 */
class RDKIT_EXTENSIONS_API CaptureRDErrorLog : private boost::noncopyable
{
  public:
    CaptureRDErrorLog();
    ~CaptureRDErrorLog();

    /**
     * @return captured messages from rdErrorLog
     */
    std::string messages() const;

  private:
    std::stringstream m_messages;
    std::ostream* m_saved_dp_dest;
    boost::logging::RDTeeStream* m_saved_teestream;
};

} // namespace adapter
} // namespace schrodinger
