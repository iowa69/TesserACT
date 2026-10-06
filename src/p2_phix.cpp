// Phase 2 W1 F5 (EVAL_PLAN_P2): the registered spike-in set -- Enterobacteria phage phiX174,
// NCBI RefSeq NC_001422.1, complete genome (5,386 bp, circular), the Illumina run-control
// spike-in. Sequence copied from work/backward/finishing/db/phix174.fa (the file the backward
// PhiX census and backward/verify_numbers/phix_scan.py used). md5 of the 5,386 uppercase bases:
// 3332ed720ac7eaa9b3655c06f6b9e196.
#include "p2_emit.h"
#include "p2_phix.h"

namespace ts {
namespace p2 {

const std::string& phix174() {
    // One copy of the sequence in the binary: emit-B's kPhiX174 (p2_phix.h), the same 5,386 bases
    // (md5 3332ed72..., checked by test_p2_emit).
    static const std::string kPhix(kPhiX174);
    return kPhix;
}

}  // namespace p2
}  // namespace ts
