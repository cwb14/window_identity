#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <unistd.h>                 // for getopt
#include <iomanip>
#include "bindings/cpp/WFAligner.hpp"

// g++ -O3 -std=c++11 wavefront_global.cpp WFA2-lib/bindings/cpp/WFAligner.cpp     -I./WFA2-lib -I./WFA2-lib/bindings/cpp     ./WFA2-lib/build/libwfa2.a -o wavefront_global

using namespace std;
using namespace wfa;

// Print help/usage message
void print_usage(const char* prog) {
    cerr << "Usage: " << prog << " [options] <fasta1> <fasta2>\n"
         << "Options:\n"
         << "  -x <int>        mismatch penalty (default 6)\n"
         << "  -O <int>        gap opening1 penalty (default 4)\n"
         << "  -E <int>        gap extension1 penalty (default 2)\n"
         << "  -o <int>        gap opening2 penalty (default 100)\n"
         << "  -e <int>        gap extension2 penalty (default 1)\n"
         << "  -c [sam|full]   CIGAR format (default full)\n"
         << "  -h              show this help message\n";
}

// Read the first record from a FASTA file into header and sequence
bool read_fasta(const string& path, string& header, string& seq) {
    ifstream in(path);
    if (!in) return false;
    string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            header = line.substr(1);
            break;
        }
    }
    if (header.empty()) return false;
    seq.clear();
    while (getline(in, line) && !line.empty() && line[0] != '>') {
        seq += line;
    }
    return true;
}

// Compute gap-compressed identity from a full CIGAR (=,X,D,I)
double compute_gc_identity(const string& cigar) {
    long matches = 0, mismatches = 0, gaps = 0;
    int  num = 0;
    for (char c : cigar) {
        if (isdigit(c)) {
            num = num * 10 + (c - '0');
        } else {
            switch (c) {
                case '=': matches    += num; break;
                case 'X': mismatches += num; break;
                case 'D':
                case 'I': gaps       += 1;   break;
                default:               break;
            }
            num = 0;
        }
    }
    return double(matches) / (matches + mismatches + gaps);
}

int main(int argc, char** argv) {
    // default penalties and format
    int mismatch  = 6;
    int go1       = 4;
    int ge1       = 2;
    int go2       = 100;
    int ge2       = 1;
    string format = "full";

    // parse options
    int opt;
    while ((opt = getopt(argc, argv, "x:O:E:o:e:c:h")) != -1) {
        switch (opt) {
            case 'x': mismatch  = stoi(optarg); break;
            case 'O': go1       = stoi(optarg); break;
            case 'E': ge1       = stoi(optarg); break;
            case 'o': go2       = stoi(optarg); break;
            case 'e': ge2       = stoi(optarg); break;
            case 'c': format    = optarg;       break;
            case 'h':
            default:  print_usage(argv[0]); return (opt=='h'?0:1);
        }
    }
    if (format!="full" && format!="sam") {
        cerr << "Error: unknown CIGAR format '" << format << "'\n";
        print_usage(argv[0]);
        return 1;
    }
    if (optind + 2 != argc) {
        print_usage(argv[0]);
        return 1;
    }
    const string file1 = argv[optind];
    const string file2 = argv[optind+1];

    // read FASTA
    string h1, s1, h2, s2;
    if (!read_fasta(file1, h1, s1)) {
        cerr << "Error reading FASTA from " << file1 << "\n";
        return 1;
    }
    if (!read_fasta(file2, h2, s2)) {
        cerr << "Error reading FASTA from " << file2 << "\n";
        return 1;
    }

    // create and run the 2-piece gap-affine aligner
    WFAlignerGapAffine2Pieces aligner(
        mismatch, go1, ge1, go2, ge2,
        WFAligner::Alignment, WFAligner::MemoryHigh
    );
    auto status = aligner.alignEnd2End(s1, s2);
    if (status < 0) {
        cerr << "Alignment failed (status=" << status << ")\n";
        return 1;
    }

    // get requested CIGAR
    string cigar = aligner.getCIGAR(format=="full");

    // output
    cout << file1 << "\t" << file2 << "\t" << cigar;
    if (format=="full") {
        double id = compute_gc_identity(cigar);
        cout << "\t" << fixed << setprecision(6) << id;
    }
    cout << "\n";
    return 0;
}
