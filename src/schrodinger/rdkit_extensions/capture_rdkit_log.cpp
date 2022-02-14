#include "schrodinger/rdkit_extensions/capture_rdkit_log.h"

namespace schrodinger
{
namespace rdkit_extensions
{

CaptureRDErrorLog::CaptureRDErrorLog()
    : m_saved_dp_dest(&m_messages), m_saved_teestream(nullptr)
{
    if (rdErrorLog == nullptr) {
        RDLog::InitLogs();
    }
    std::swap(rdErrorLog->dp_dest, m_saved_dp_dest);
    std::swap(rdErrorLog->teestream, m_saved_teestream);
}

CaptureRDErrorLog::~CaptureRDErrorLog()
{
    std::swap(rdErrorLog->dp_dest, m_saved_dp_dest);
    std::swap(rdErrorLog->teestream, m_saved_teestream);
}

std::string CaptureRDErrorLog::messages() const
{
    return m_messages.str();
}

} // namespace rdkit_extensions
} // namespace schrodinger
