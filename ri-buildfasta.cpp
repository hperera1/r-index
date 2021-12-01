#include <cstdio>
#include <algorithm>
#include <chrono>
#include <string>
#include <getopt.h>

#include <iostream>
#include "internal/utils.hpp"
#include "internal/r_index.hpp"

using std::cout;
using std::cerr;

struct ribuild_args {
    std::string input_fname = "";
    std::string outpre = "";
    std::string bwt_alg = "bigbwt";
    uint32_t threads = 0;
    bool acgt = false;
};

void print_help() {
	fprintf(stderr, "ri-buildfasta: builds the r-index from a fasta file. Extension .ri is automatically added to output index file\n");
	fprintf(stderr, "Usage: ri-buildfasta [options] <input_fasta_name>\n");
	fprintf(stderr, "    --output_prefix <basename>\n");
	fprintf(stderr, "    -o <basename>                  use 'basename' as prefix for all index files. Default: basename is the specified input_file_name\n");
	fprintf(stderr, "    --bwt_alg <alg>\n");
	fprintf(stderr, "    -b <alg>                       bwt algorithm (bigbwt, sais), default:bigbwt\n");
    fprintf(stderr, "    --threads <int>\n");
    fprintf(stderr, "    -t <int>                       number of threads for parallelized construction. parallelism not yet supported for gzipped files\n");
    fprintf(stderr, "    --acgt\n");
    fprintf(stderr, "    -a                             strip out all non-ACGT characters from the input (this means Ns will also be removed)\n");
	fprintf(stderr, "    <input_file_name>              input text file.\n");
}

ribuild_args parse_args(int argc, char** argv) {
    int c;
    char* end;
    ribuild_args args;
    static struct option long_options[] {
        {"output_prefix", required_argument, 0, 'o'},
        {"bwt_alg", required_argument, 0, 'b'},
        {"acgt", no_argument, 0, 'a'},
        {"threads", required_argument, 0, 't'}
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "o:b:at:", long_options, &long_index)) != -1) {
        switch (c) {
            case 'o':
                args.outpre = optarg;
                break;
            case 'b':
                args.bwt_alg = optarg;
                break;
            case 'a':
                args.acgt = true;
                break;
            case 't':
                args.threads = std::strtoul(optarg, &end, 10);
                break;
            default:
                print_help();
                exit(1);
                break;
        }
    }

    if (argc - optind < 1) {
        fprintf(stderr, "no argument provided\n");
        exit(1);
    }

    args.input_fname = argv[optind];
    if (args.outpre == "") {
        args.outpre = args.input_fname;
    }
    return args;
}

std::string read_fasta(std::string input_file) {
    kseq_t *seq;
    int l;
    gzFile fp = gzopen(input_file.c_str(), "r");
    seq = kseq_init(fp);
    std::string input("");
    while ((l = kseq_read(seq)) >= 0) {
        std::string line(seq->seq.s);
        std::transform(line.begin(), line.end(), line.begin(), ::toupper);
        input += line;
    }
    kseq_destroy(seq);
    gzclose(fp);
    return input;
}

int main(int argc, char** argv) {
    ribuild_args args(parse_args(argc, argv));
    bool fast = false;
    std::string path(args.outpre + ".ri");
    auto t1 = std::chrono::high_resolution_clock::now();
    std::ofstream out(path);
    out.write((char*) &fast, sizeof(fast));
    ri::r_index<> idx;
    if (args.bwt_alg == "sais") {
        std::string input(read_fasta(args.input_fname));
        std::cout << "building forward index using sais" << std::endl;
        idx.init_sais(input,true);
    } else if (args.bwt_alg == "bigbwt") {
        std::cout << "building forward index using bigbwt" << std::endl;
        idx.init_bigbwt(args.input_fname, args.threads);
        std::string fai_path(args.outpre + ".1.ri");
        std::cout << "generating faidx of sequences, saving to" << fai_path << std::endl;
        ri::build_seqidx(args.input_fname.c_str(), fai_path.c_str());
    } else if (args.bwt_alg == "from_bwt") {
        std::cout << "building forward index from existing bwt in "  << args.input_fname << std::endl;
        idx.from_bwt(args.input_fname);
    } else {
        std::cerr << "\"" << args.bwt_alg << "\" is an invalid bwt algortihm. Please use one of [sais, bigbwt, from_bwt]" << std::endl;
        exit(1);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    ulint total = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1>>>(t2 - t1).count();
    std::cout << "Build time : " << get_time(total) << endl;
    std::cout << "writing forward index to " << path << std::endl;
    idx.serialize(out);
}
