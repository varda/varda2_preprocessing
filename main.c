#include <stddef.h>     // size_t
#include <stdio.h>      // fprintf
#include <stdlib.h>     // EXIT_*, atoi

#include "htslib/vcf.h"     // hts*, bcf_*


int
main(int argc, char* argv[])
{
    if (2 > argc)
    {
        fprintf(stderr, "Usage: %s THRESHOLD\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    int const threshold = atoi(argv[1]);  // FIXME: unsafe

    htsFile* const restrict fh = bcf_open("-", "r");
    if (NULL == fh)
    {
        fprintf(stderr, "bcf_open() failed\n");
        return EXIT_FAILURE;
    } // if

    bcf_hdr_t* const restrict hdr = bcf_hdr_read(fh);

    if (1 != bcf_hdr_nsamples(hdr))
    {
        fprintf(stderr, "#samples = %d\n", bcf_hdr_nsamples(hdr));
        bcf_hdr_destroy(hdr);
        bcf_close(fh);
        return EXIT_FAILURE;
    } // if

    int nseq = 0;
    char const** seqnames = bcf_hdr_seqnames(hdr, &nseq);
    if (NULL == seqnames)
    {
        fprintf(stderr, "bcf_hdr_seqnames() failed\n");
        bcf_hdr_destroy(hdr);
        bcf_close(fh);
        return EXIT_FAILURE;
    } // if

    bcf1_t* restrict rec = bcf_init();
    if (NULL == rec)
    {
        fprintf(stderr, "bcf_init() failed\n");
        bcf_hdr_destroy(hdr);
        bcf_close(fh);
        return EXIT_FAILURE;
    } // if

    int ndp_arr = 0;
    int* restrict dp = NULL;

    while (0 == bcf_read(fh, hdr, rec))
    {
        int depth = 0;
        if (1 == bcf_get_format_int32(hdr, rec, "DP", &dp, &ndp_arr))
        {
            depth = dp[0];
        } // else

        if (threshold <= depth)
        {
            fprintf(stdout, "%s\t%d\t%d\n", seqnames[rec->rid],
                                            rec->pos,
                                            rec->pos + rec->rlen);
        } // if
    } // while

    free(seqnames);
    free(dp);

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fh);

    return EXIT_SUCCESS;
} // main
