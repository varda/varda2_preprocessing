#include <stddef.h>     // size_t
#include <stdio.h>      // fprintf
#include <stdlib.h>     // EXIT_*, atoi
#include <getopt.h>     // getopt

#include "htslib/vcf.h"     // hts*, bcf_*


int
main(int argc, char* argv[])
{

    extern char *optarg;
    int threshold = 10, merge = 1, distance = 0;
    int nflag = 0, dflag = 0;
    int c, err = 0;

    static char usage[] = "usage: %s [-t threshold] [-d distance] [-n]\n";

    while ((c = getopt(argc, argv, "t:nd:")) != -1) {
        switch (c) {
        case 't':
            threshold = atoi(optarg);
            break;
        case 'n':
            nflag = 1;
            merge = 0;
            break;
        case 'd':
            dflag = 1;
            distance = atoi(optarg);
            break;
        case '?':
            err = 1;
            break;
        }
    }

    if (nflag && dflag) {
        fprintf(stderr, "Both no merge and distance specified!\n");
        err += 1;
    }

    if (err) {
        fprintf(stderr, usage, argv[0]);
        return EXIT_FAILURE;
    }

    htsFile* const fh = bcf_open("-", "r");
    if (NULL == fh)
    {
        fprintf(stderr, "bcf_open() failed\n");
        return EXIT_FAILURE;
    } // if

    bcf_hdr_t* const hdr = bcf_hdr_read(fh);

    if (1 != bcf_hdr_nsamples(hdr))
    {
        fprintf(stderr, "#samples = %d\n", bcf_hdr_nsamples(hdr));
        goto error1;
    } // if


    int nseq = 0;
    char const** seqnames = bcf_hdr_seqnames(hdr, &nseq);
    if (NULL == seqnames)
    {
        fprintf(stderr, "bcf_hdr_seqnames() failed\n");
        goto error1;
    } // if

    bcf1_t* rec = bcf_init();
    if (NULL == rec)
    {
        fprintf(stderr, "bcf_init() failed\n");
        goto error2;
    } // if

    int32_t *dp = NULL;
    int32_t *gt = NULL;

    while (0 == bcf_read(fh, hdr, rec))
    {
        int32_t depth = 0;
        if (1 == bcf_get_format_int32(hdr, rec, "DP", &dp, &(int){0}))
        {
            depth = dp[0];
        } // if

        int ngt_arr = 0;
        int const ploidy = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
        if (0 > ploidy)
        {
            fprintf(stderr, "bcf_get_genotypes failed\n");
            goto error;
        } // if

        if (ploidy != ngt_arr)
        {
            fprintf(stderr, "ploidy != count\n");
            goto error;
        } // if

        if (threshold <= (int)depth)
        {
            fprintf(stdout, "%s\t%d\t%d\t%d\n", seqnames[rec->rid],
                                                rec->pos,
                                                rec->pos + rec->rlen,
                                                ploidy);
        } // if
    } // while

    bcf_destroy(rec);
    free(dp);
    free(gt);
    free(seqnames);
    bcf_hdr_destroy(hdr);
    bcf_close(fh);

    return EXIT_SUCCESS;

error:
    bcf_destroy(rec);
    free(dp);
    free(gt);

error2:
    free(seqnames);

error1:
    bcf_hdr_destroy(hdr);
    bcf_close(fh);

    return EXIT_FAILURE;

} // main
