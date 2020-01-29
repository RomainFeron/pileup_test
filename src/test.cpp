#include <iostream>
#include <string>
#include <htslib/htslib/hts.h>
#include <htslib/htslib/sam.h>

using namespace std;

int maiaan(int argc, char *argv[]) {

    if (argc == 1) {
        std::cout << "Usage: test file.<sam|bam|cram>" << std::endl;
        return 1;
    }

    const std::string sam_file_path = argv[1];
    htsFile* file = hts_open(sam_file_path.c_str(), "r");

    bam_hdr_t *bamHdr = sam_hdr_read(file);
    bam1_t *aln = bam_init1();

//    std::cout << file->is_bin << "\t" << file->is_cram << std::endl;

//    for (auto i=0; i<bamHdr->n_targets; ++i) std::cout << i << "\t" << bamHdr->target_name[i] << std::endl;

    while(sam_read1(file, bamHdr, aln) >= 0){

        std::cout << "New" << std::endl;
        int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
        std::cout << "Pos: " << pos << std::endl;
        std::cout << "chr-tid: " << aln->core.tid << std::endl;
        char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
        std::cout << "chr: " << chr << std::endl;
        uint32_t len = aln->core.l_qseq; //length of the read.
        std::cout << "len: " << len << std::endl;
        uint8_t *q = bam_get_seq(aln); //quality string
        std::cout << "q: " << q << std::endl;
        uint32_t q2 = aln->core.qual ; //mapping quality
        std::cout << "q2: " << q2 << std::endl;
        char *qseq = (char *)malloc(len);
        for(int i=0; i< len ; i++){
            qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
        }
    }

    bam_destroy1(aln);
    hts_close(file);

    return 0;
}
