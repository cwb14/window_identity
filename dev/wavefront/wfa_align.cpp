// Batch wavefront aligner for Step 10.
//
// Reads alignment jobs from stdin and writes CIGARs to stdout, one line each:
//     stdin :  id \t target_seq \t query_seq
//     stdout:  id \t status \t score \t cigar        (cigar uses =/X/I/D)
//
// Streaming rather than one process per segment: gene-pair partitioning produces ~16k
// segments, so process startup would dominate. All coordinate arithmetic and PAF emission
// live in Python (bin/paf_emit.py); this binary only aligns.
//
// Two defaults here are deliberate corrections of the previous attempt
// (dev/wavefront/wavefront_global.cpp), which scored ~70% identity on real blocks:
//
//   1. setHeuristicNone() is called explicitly. WFA2-lib's wavefront_aligner_attr_default
//      sets .heuristic.strategy = wf_heuristic_wfadaptive with max_distance_threshold=50,
//      and the C++ ctor copies those defaults, so NOT calling this leaves lossy adaptive
//      wavefront pruning switched on. A 100bp-15kb TE indel needs an excursion of thousands
//      of diagonals; anything >50 off the best diagonal is pruned and the true path dies.
//      Measured cost of leaving it on: -14 identity points at 30kb, ~0 at 5kb -- which is
//      exactly why the old attempt looked fine on small gene pairs and broke on blocks.
//
//   2. MemoryUltralow (BiWFA) is the default. Exact WFA memory is O(s^2) under MemoryHigh;
//      the largest pair here is 8.4Mb x 10.1Mb at ~20% divergence, which would exhaust RAM.
//
// Build:
//   g++ -O3 -std=c++11 wfa_align.cpp ../../../WFA2-lib/bindings/cpp/WFAligner.cpp \
//       -I../../../WFA2-lib -I../../../WFA2-lib/bindings/cpp \
//       ../../../WFA2-lib/build/libwfa2.a -o wfa_align

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include "bindings/cpp/WFAligner.hpp"

using namespace std;
using namespace wfa;

// Stream jobs through an already-configured aligner.
static int run(WFAligner& aligner, const string& heur, const string& span, int free_slack) {
    // Not optional. See the header comment: the library default is a lossy heuristic.
    if (heur == "none") aligner.setHeuristicNone();

    string line;
    while (getline(cin, line)) {
        if (line.empty()) continue;
        size_t t1 = line.find('\t');
        if (t1 == string::npos) continue;
        size_t t2 = line.find('\t', t1 + 1);
        if (t2 == string::npos) continue;

        string id   = line.substr(0, t1);
        string tseq = line.substr(t1 + 1, t2 - t1 - 1);
        string qseq = line.substr(t2 + 1);

        int status;
        if (span == "end2end") {
            status = aligner.alignEnd2End(tseq, qseq);
        } else {
            status = aligner.alignEndsFree(tseq, free_slack, free_slack,
                                           qseq, free_slack, free_slack);
        }
        if (status != 0) {
            cout << id << "\t" << status << "\t0\t\n";
            cout.flush();
            continue;
        }
        cout << id << "\t0\t" << aligner.getAlignmentScore() << "\t"
             << aligner.getCIGAR(true) << "\n";
        cout.flush();
    }
    return 0;
}

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " [options]   (reads id\\tTSEQ\\tQSEQ from stdin)\n"
         << "Penalties:\n"
         << "  -x <int>   mismatch (default 6)\n"
         << "  -O <int>   gap open 1 (default 4)\n"
         << "  -E <int>   gap extend 1 (default 2)\n"
         << "  -o <int>   gap open 2 (default 100, affine2p only)\n"
         << "  -e <int>   gap extend 2 (default 1, affine2p only)\n"
         << "Model:\n"
         << "  -m <str>   affine | affine2p            (default affine2p)\n"
         << "  -s <str>   end2end | endsfree           (default end2end)\n"
         << "  -f <int>   ends-free slack in bp        (default 0)\n"
         << "  -M <str>   high|med|low|ultralow        (default ultralow = BiWFA)\n"
         << "  -H <str>   none|adaptive               (default none)\n"
         << "  -h         help\n";
}

int main(int argc, char** argv) {
    int mismatch = 6, o1 = 4, e1 = 2, o2 = 100, e2 = 1;
    string model = "affine2p", span = "end2end", mem = "ultralow", heur = "none";
    int free_slack = 0;

    int opt;
    while ((opt = getopt(argc, argv, "x:O:E:o:e:m:s:f:M:H:h")) != -1) {
        switch (opt) {
            case 'x': mismatch = atoi(optarg); break;
            case 'O': o1 = atoi(optarg); break;
            case 'E': e1 = atoi(optarg); break;
            case 'o': o2 = atoi(optarg); break;
            case 'e': e2 = atoi(optarg); break;
            case 'm': model = optarg; break;
            case 's': span = optarg; break;
            case 'f': free_slack = atoi(optarg); break;
            case 'M': mem = optarg; break;
            case 'H': heur = optarg; break;
            case 'h': usage(argv[0]); return 0;
            default: usage(argv[0]); return 1;
        }
    }
    if (model != "affine" && model != "affine2p") {
        cerr << "Error: -m must be affine or affine2p\n"; return 1;
    }
    if (span != "end2end" && span != "endsfree") {
        cerr << "Error: -s must be end2end or endsfree\n"; return 1;
    }

    WFAligner::MemoryModel mm;
    if      (mem == "high")     mm = WFAligner::MemoryHigh;
    else if (mem == "med")      mm = WFAligner::MemoryMed;
    else if (mem == "low")      mm = WFAligner::MemoryLow;
    else if (mem == "ultralow") mm = WFAligner::MemoryUltralow;
    else { cerr << "Error: -M must be high|med|low|ultralow\n"; return 1; }

    if (heur != "none" && heur != "adaptive") {
        cerr << "Error: -H must be none or adaptive\n"; return 1;
    }

    // WFAligner's destructor is protected, so the concrete aligner is held by value in
    // main and only a reference is passed around.
    if (model == "affine") {
        WFAlignerGapAffine aligner(mismatch, o1, e1, WFAligner::Alignment, mm);
        return run(aligner, heur, span, free_slack);
    }
    WFAlignerGapAffine2Pieces aligner(mismatch, o1, e1, o2, e2, WFAligner::Alignment, mm);
    return run(aligner, heur, span, free_slack);
}
