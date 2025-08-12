#include "schrodinger/rdkit_extensions/capture_rdkit_log.h"

namespace schrodinger
{
namespace rdkit_extensions
{

CaptureRDErrorLog::CaptureRDErrorLog() :
    m_saved_dp_dest(&m_messages),
    m_saved_teestream(nullptr)
{
    if (rdErrorLog == nullptr) {
        RDLog::InitLogs();
    } else {
        // Make sure at least the error log is active
        // so we can capture something.
        m_error_log_initial_state = rdErrorLog->df_enabled;
        rdErrorLog->df_enabled = true;
    }
    std::swap(rdErrorLog->dp_dest, m_saved_dp_dest);
    std::swap(rdErrorLog->teestream, m_saved_teestream);
}

CaptureRDErrorLog::~CaptureRDErrorLog()
{
    std::swap(rdErrorLog->dp_dest, m_saved_dp_dest);
    std::swap(rdErrorLog->teestream, m_saved_teestream);
    rdErrorLog->df_enabled = m_error_log_initial_state;
}

std::string CaptureRDErrorLog::messages() const
{
    return m_messages.str();
}

} // namespace rdkit_extensions
} // namespace schrodinger
