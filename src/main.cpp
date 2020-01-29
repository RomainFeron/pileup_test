#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <iostream>
#include "htslib/htslib/sam.h"

int main(int argc, char *argv[]) {

    int ret = 0;
    samFile *in = nullptr;
    sam_hdr_t *header = nullptr;
    char *fn_in = nullptr;

    std::cout << "Parsing arguments" << std::endl;
    fn_in = argv[1];
    htsFormat format;
    bam1_t *b = bam_init1();

    // First open file to detect format
    if ((in = hts_open(fn_in, "r")) == nullptr) {
        ret = 1;
        goto view_end;
    }

    std::cout << "Detecting format" << std::endl;

    if (not in->is_cram) {
        if (hts_detect_format(in->fp.hfile, &format) < 0) {
            ret = 1;
            goto view_end;
        }
    } else {
        hts_parse_format(&format, "cram");
    }

    // First open file to detect format
    hts_close(in);


    if (argc == 3) {
        std::cout << "Updating format" << std::endl;
        char *ref = reinterpret_cast<char *>(malloc(10 + strlen(argv[2]) + 1));
        hts_opt_add(reinterpret_cast<hts_opt **>(format.specific), ref); // HERE UNDERSTAND WHAT TO DO
    }

    std::cout << "Opening files" << std::endl;

    // open file handlers
    if ((in = hts_open_format(fn_in, "r", &format)) == nullptr) {
        ret = 1;
        goto view_end;
    }

    if ((header = sam_hdr_read(in)) == nullptr) {
        ret = 1;
        goto view_end;
    }

    std::cout << "Reading bam" << std::endl;

    int r;
    while ((r = sam_read1(in, header, b)) >= 0) {
        std::cout << "New" << std::endl;
        long pos = b->core.pos +1; //left most position of alignment in zero based coordianate (+1)
        std::cout << "Pos: " << pos << std::endl;
        std::cout << "chr-tid: " << b->core.tid << std::endl;
        char *chr = header->target_name[b->core.tid] ; //contig name (chromosome)
        std::cout << "chr: " << chr << std::endl;
        long len = b->core.l_qseq; //length of the read.
        std::cout << "len: " << len << std::endl;
        uint8_t *q = bam_get_seq(b); //quality string
        uint32_t q2 = b->core.qual ; //mapping quality
        std::cout << "mapQ: " << q2 << std::endl;
        char *qseq = (char *)malloc(len);
        for(int i=0; i< len ; i++){
            qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
        }
    }
    if (r < -1) {
        fprintf(stderr, "[main_samview] truncated file.\n");
        ret = 1;
    }
    bam_destroy1(b);


view_end:

    // close files, free and return
    if (in) hts_close(in);
//    if (header) sam_hdr_destroy(header);
    return ret;
}
