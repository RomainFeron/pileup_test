#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include "htslib/htslib/sam.h"
#include "htslib/faidx.h"


// Simple structure holding all information about an input file
struct inputFile {
    htsFile *sam;  // Main file descriptor (for the alignment file)
    hts_idx_t *idx;  // Index file descriptor
    sam_hdr_t *header;  // Header information read directly from main file
    uint16_t file_n;  // Input file number
};


// Open an alignment file in a format-agnostic way and fill an inputFile object with all the information
int open_input(char *fn_in, inputFile *file, std::string &reference, uint16_t file_n) {

    // Open alignment file and handle opening error
    if ((file->sam = hts_open(fn_in, "r")) == nullptr) {
        std::cerr << "Error opening alignment file <" << fn_in << ">" << std::endl;
        return 1;
    }

    // CRAM files require a reference. Need to add the reference path and reference index path to the file descriptor
    if (file->sam->is_cram) {
        // Add reference file path to file descriptor
        std::string ref_option = "reference=" + reference; // Create the string "reference=<provided/path/to/ref>" to add as option to format in htsFile
        hts_opt_add(reinterpret_cast<hts_opt **>(&file->sam->format.specific), ref_option.c_str());  // Add reference to htsFile
        std::string fai_path = reference + ".fai";  // Create the string "<provided/path/to/ref.fai>"
        if (hts_set_fai_filename(file->sam, fai_path.c_str()) < 0) {  // Set reference index path in file descriptor
            std::cerr << "Warning: index file not found for reference file <" << reference << ">. Indexing reference" << std::endl;
            if (fai_build(reference.c_str()) < 0) {  // Build reference fasta index if missing
                std::cerr << "Error: could not build index for reference file <" << reference << ">" << std::endl;
                return 1;
            }
        }
    }

    // Read file header and handle errors
    if ((file->header = sam_hdr_read(file->sam)) == nullptr) {
        std::cerr << "Error reading header for alignment file <" << fn_in << ">" << std::endl;
        return 1;
    }

    file->idx = sam_index_load(file->sam, fn_in); // Load index for alignment file. Index name is automatically infered from alignment file name

    // Handle error opening index for alignment file
    if (file->idx == nullptr) {
        std::cerr << "Alignment file <" << file->sam->fn << "> is not indexed. Index with 'samtools index " << file->sam->fn << "'" << std::endl;
        return 1;
    }

    file->file_n = file_n;

    return 0;
}


int process_file(inputFile* input, char *region, std::vector<std::vector<uint16_t>>& depths, uint min_qual=0) {

    hts_itr_t *iter = nullptr;
    bam1_t *b = bam_init1();
    int result;

    // sam_itr_querys parses the string given by <contig> to find the region with format `contig:start-end' and returns an iterator
    if ((iter = sam_itr_querys(input->idx, input->header, region)) == nullptr) {
        std::cerr << "Region <" << region << "> not found in index file";
        return 1;
    }

    uint32_t mapping_position = 0;
    uint8_t *sequence = nullptr;
    uint16_t mapping_quality = 0;
    char nucleotide;
    uint32_t *cigar = nullptr;

    // Iterate through all alignments in the specified region
    while ((result = sam_itr_next(input->sam, iter, b)) >= 0) {
        mapping_quality = b->core.qual ;
        if (mapping_quality < min_qual) continue;  // Skip reads with low mapping quality
        mapping_position = static_cast<uint32_t>(b->core.pos);
        sequence = bam_get_seq(b);
        cigar = bam_get_cigar(b);
        for (uint k = 0; k < b->core.n_cigar; ++k) {
          uint op = bam_cigar_op(cigar[k]);
          uint l = bam_cigar_oplen(cigar[k]);
          if (op == BAM_CMATCH) {
            for (uint j = mapping_position; j < mapping_position + l; ++j) {
                nucleotide = seq_nt16_str[bam_seqi(sequence, j - mapping_position)]; // Get nucleotide id from read sequence and convert it to <ATGCN>.
                switch (nucleotide) {
                    case 'A':
                        ++depths[j][input->file_n * 6];
                        break;
                    case 'T':
                        ++depths[j][input->file_n * 6 + 1];
                        break;
                    case 'C':
                        ++depths[j][input->file_n * 6 + 2];
                        break;
                    case 'G':
                        ++depths[j][input->file_n * 6 + 3];
                        break;
                    case 'N':
                        ++depths[j][input->file_n * 6 + 4];
                        break;
                    default:
                        ++depths[j][input->file_n * 6 + 5];
                        break;
                }
            }
            mapping_position += l;
          } else if (op == BAM_CREF_SKIP || op == BAM_CDEL) {
              mapping_position += l;
          }
        }
    }

    // Destroy objects
    hts_itr_destroy(iter);
    bam_destroy1(b);

    if (result < -1) {
        std::cerr << "Error processing region <" << region << "> in file <" << input->sam->fn << "> due to truncated file or corrupt BAM index file";
        return 1;
    }

    return 0;
}


int main(int argc, char *argv[]) {

    if (argc < 3) {
        std::cerr << "Usage: test reference.fa in.<sam|bam|cram> [in2.<sam|bam|cram> ...]";
        return 1;
    }

    int main_return = 0;
    char *contig = nullptr;
    uint32_t contig_len = 0;
    std::vector<std::vector<uint16_t>> depths;
    uint n_files = static_cast<uint>(argc) - 2;  // Number of alignment files to process

    std::string reference = argv[1];

    // Properly open all alignment files with all necessary information (header, indexes, reference ...) and store them in a vector
    std::cout << "#Files";  // Comment line in output with names of all processed alignment files in order
    std::vector<inputFile> input;
    for (uint16_t i=2; i<argc; ++i) {
        inputFile tmp;
        if (open_input(argv[i], &tmp, reference, i - 2) != 0) {
            main_return = 1;
            goto end;
        }
        input.push_back(tmp);
        std::cout << "\t" << argv[i];  // Output alignment file path to comment output string
    }
    std::cout << "\n";

    // Process all alignment files contig by contig to reduce memory usage
    for (int i=0; i<input[0].header->n_targets; ++i) {

        contig = input[0].header->target_name[i];
        contig_len = input[0].header->target_len[i];

        std::cerr << "Processing contig " << contig << " (" << contig_len << " bp)" << std::endl;

        // Depths: {position: [nA, nT, nG, nC, nN, nOther] * number of files}
        depths.resize(0);
        depths.resize(contig_len);
        for (uint k=0; k<contig_len; ++k) {
            depths[k].resize(6 * n_files);
        }

        // Process each alignment file
        for (auto file: input) {
            if (process_file(&file, contig, depths, 0) != 0) {
                main_return = 1;
                goto end;
            }
        }

        // Output depths for this region. Format:
        // - 1 line with format "region=<region>\t<len=<region_length>"
        // - for each position in region (in order), "nA, nT, nG, nC, nN, nOther" for each alignment file, alignment files are tab-separated
        std::cout << "region=" << contig << "\tlen=" << contig_len << "\n";
        for (uint j=0; j<contig_len; ++j) {
            for (uint k=0; k<n_files; ++k) {
                for (uint l=0; l<6; ++l) {
                    std::cout << depths[j][l + 6 * k];
                    if (l < 5) std::cout << ",";
                }
                (k < input.size() - 1) ? std::cout << "\t" : std::cout << "\n";
            }
        }
    }

end:
    for (auto f: input) {  // Destroy all created objects
        if (f.sam) hts_close(f.sam);
        if (f.header) sam_hdr_destroy(f.header);
        if (f.idx) hts_idx_destroy(f.idx);
    }

    return main_return;
}
