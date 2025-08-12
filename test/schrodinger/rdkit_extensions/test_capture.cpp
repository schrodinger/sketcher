/* -------------------------------------------------------------------------
 * Tests class schrodinger::rdkit_extensions:: rdkit error log capture
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_MODULE rdkit_extensions_capture

#include <boost/test/unit_test.hpp>

#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/RDGeneral/RDLog.h>

#include "schrodinger/rdkit_extensions/capture_rdkit_log.h"

using namespace schrodinger::rdkit_extensions;

BOOST_AUTO_TEST_CASE(testCaptureErrorLogEnabled)
{
    // SHARED-8613

    std::string bad_smarts{"garbage"};

    // Make sure logging is not redirected before we start
    // (it should be either stderr or not initialized)
    BOOST_REQUIRE(rdErrorLog == nullptr || rdErrorLog->dp_dest == &std::cerr);

    {
        CaptureRDErrorLog outer_capture_log;
        RDLog::LogStateSetter silence_rdkit_logging;

        RDKit::SmartsToMol(bad_smarts);
        BOOST_REQUIRE(outer_capture_log.messages().empty());

        {
            CaptureRDErrorLog capture_log;
            RDKit::SmartsToMol(bad_smarts);

            // Capture enables the error log, so now we expect a message
            std::string err_msg{"SMARTS Parse Error: Failed parsing SMARTS "
                                "'garbage' for input: 'garbage'"};
            BOOST_CHECK(capture_log.messages().find(err_msg) !=
                        std::string::npos);
        }

        RDKit::SmartsToMol(bad_smarts);
        BOOST_REQUIRE(outer_capture_log.messages().empty());
    }

    // We should have initialized && restored the error stream
    BOOST_CHECK(rdErrorLog && rdErrorLog->dp_dest == &std::cerr);
}
