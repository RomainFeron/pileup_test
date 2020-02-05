#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>
#include "htslib/htslib/sam.h"
#include "htslib/faidx.h"


struct inputFile {
    htsFile *sam;
    hts_idx_t *idx;
    sam_hdr_t *header;
    uint16_t file_n;
};


int open_input(char *fn_in, inputFile *file, char *reference, uint16_t file_n) {

    // Open alignment file and handle opening error
    if ((file->sam = hts_open(fn_in, "r")) == nullptr) {
        std::cerr << "Error opening alignment file <" << fn_in << ">" << std::endl;
        return 1;
    }

    // CRAM files require a reference
    if (file->sam->is_cram) {
        char *ref = reinterpret_cast<char *>(malloc(10 + strlen(reference) + 1));
        char *fai_path;
        if (!ref) {
            std::cerr << "Error allocating memory for reference file <" << reference << ">" << std::endl;
            return 1;
        }
        sprintf(ref, "reference=%s", reference);  // Creating the string "reference=<provided/path/to/ref>" to add as option to format in htsFile
        hts_opt_add(reinterpret_cast<hts_opt **>(&file->sam->format.specific), ref);  // Add reference to htsFile
        free(ref);
        fai_path = reinterpret_cast<char *>(malloc(5 + strlen(reference)));
        sprintf(fai_path, "%s.fai", reference);
        hts_set_fai_filename(file->sam, fai_path);
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


int process_file(inputFile* input, char *contig, std::vector<std::vector<uint16_t>>& depths, uint min_qual=0) {

    hts_itr_t *iter = nullptr;
    bam1_t *b = bam_init1();
    int result;

    iter = sam_itr_querys(input->idx, input->header, contig); // parse a region in the format like `chr2:100-200'
    if (iter == nullptr) { // region invalid or reference name not found
        std::cerr << "Region <" << contig << "> not found in index file";
        return 1;
    }

    uint32_t pos = 0;
    long read_length = 0;
    uint8_t *q = nullptr;
    uint16_t mapping_quality = 0;
    char qseq;

    while ((result = sam_itr_next(input->sam, iter, b)) >= 0) {
        mapping_quality = b->core.qual ;
        if (mapping_quality < min_qual) continue;
        pos = static_cast<uint32_t>(b->core.pos);
        read_length = b->core.l_qseq;
        q = bam_get_seq(b);
        for(uint16_t p=0; p<read_length ; p++){
            qseq = seq_nt16_str[bam_seqi(q,p)]; //gets nucleotide id and converts them into IUPAC id.
            switch (qseq) {
                case 'A':
                    ++depths[pos + p + input->file_n * 6][0];
                    break;
                case 'T':
                    ++depths[pos + p + input->file_n * 6][1];
                    break;
                case 'G':
                    ++depths[pos + p + input->file_n * 6][2];
                    break;
                case 'C':
                    ++depths[pos + p + input->file_n * 6][3];
                    break;
                case 'N':
                    ++depths[pos + p + input->file_n * 6][4];
                    break;
                default:
                    ++depths[pos + p + input->file_n * 6][5];
                    break;
            }
        }
    }

    // Destroy iterator
    hts_itr_destroy(iter);
    bam_destroy1(b);

    if (result < -1) {
        std::cerr << "Error processing contig <" << contig << "> in file <" << input->sam->fn << "> due to truncated file or corrupt BAM index file";
        return 1;
    }

    return 0;
}

int main(int argc, char *argv[]) {

    int main_return = 0, step_return = 0;
    char *contig = nullptr;
    uint32_t contig_len = 0;
    std::vector<std::vector<uint16_t>> depths;
    uint n_files = static_cast<uint>(argc) - 2;

    if (argc < 3) {
        std::cerr << "Usage: test reference.fa in.<sam|bam|cram> [in2.<sam|bam|cram> ...]";
        return 1;
    }

    std::cerr << "Getting command-line arguments" << std::endl;

    char *reference = argv[1];

    std::vector<inputFile> input;
    for (uint16_t i=2; i<argc; ++i) {
        inputFile tmp;
        step_return = open_input(argv[i], &tmp, reference, i - 2);
        if (step_return != 0) {
            main_return = 1;
            goto end;
        }
        input.push_back(tmp);
    }

    for (int i=0; i<input[0].header->n_targets; ++i) {

        contig = input[0].header->target_name[i];
        contig_len = input[0].header->target_len[i];

        std::cerr << "Processing contig " << contig << " (" << contig_len << " bp)" << std::endl;

        depths.resize(0);
        depths.resize(contig_len);
        for (uint k=0; k<contig_len; ++k) {
            depths[k].resize(6 * n_files);
        }

        for (auto file: input) {
            step_return = process_file(&file, contig, depths);
            if (step_return != 0) {
                main_return = 1;
                goto end;
            }
        }

        std::cout << "contig=" << contig << "\tlen=" << contig_len << "\n";
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
    for (auto f: input) {
        if (f.sam) hts_close(f.sam);
        if (f.header) sam_hdr_destroy(f.header);
        if (f.idx) hts_idx_destroy(f.idx);
    }
    return main_return;
}
