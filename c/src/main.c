#include <stddef.h>     // size_t
#include <stdio.h>      // fprintf
#include <stdlib.h>     // EXIT_*, atoi
#include <getopt.h>     // getopt
#include <sys/param.h>  // MIN, MAX

#include "htslib/vcf.h"     // hts*, bcf_*


int
main(int argc, char* argv[])
{

    extern char *optarg;
    int threshold = 10, merge = 1, distance = 0;
    int nflag = 0, dflag = 0;
    int c, err = 0;

    static char usage[] = "usage: %s [-t threshold] [-d distance] [-n]\n";

    while ((c = getopt(argc, argv, "t:nd:")) != -1)
    {
        switch (c)
        {
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
        } // switch
    } // while

    if (nflag && dflag)
    {
        fprintf(stderr, "Both no merge and distance specified!\n");
        err += 1;
    } // if

    if (err)
    {
        fprintf(stderr, usage, argv[0]);
        return EXIT_FAILURE;
    } // if

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

    int first = 1;
    int jump = 1;

    int window_start = 0;
    int window_end = 0;
    const char * window_chrom = "";
    int window_ploidy = 0;

    while (0 == bcf_read(fh, hdr, rec)) {

        jump = 0;

        int32_t depth = 0;
        if (1 == bcf_get_format_int32(hdr, rec, "DP", &dp, &(int){0}))
        {
            depth = dp[0];
        } // if

        //
        // If depth is below the threshold, no need to proceed
        //
        if (depth < threshold)
        {
            continue;
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

        //
        // Convenience handles
        //
        const char *chrom = seqnames[rec->rid];
        int start = rec->pos;
        int end = rec->pos + rec->rlen;

        //
        // When we don't merge, just print here and proceed
        //
        if (!merge)
        {
            fprintf(stdout, "%s\t%d\t%d\t%d\n", chrom, start, end, ploidy);
            continue;
        } // if

        //
        // We just started
        //
        if (first)
        {
            window_start = start;
            window_end = end;
            window_chrom = chrom;
            window_ploidy = ploidy;

            // eprint(f"First! c:{window_chrom} s:{start}, w_s={window_start} e:{end} w_e={window_end}")

            first = 0;
            continue;
        } // if

        if (window_chrom != chrom)
        {
            // eprint(f"Chrom changed from {window_chrom} to {chrom}.")
            jump = 1;
        } else if (window_ploidy != ploidy)
        {
            // eprint(f"Ploidy changed from {window_ploidy} to {ploidy}")
            jump = 1;
        } else if (window_end + distance < start)
        {
            // eprint("Gap! (window_end:%d < start:%d)" % (window_end + distance, start))
            jump = 1;
        } // if

        if (jump)
        {
            fprintf(stdout, "%s\t%d\t%d\t%d\n", window_chrom, window_start, window_end, window_ploidy);

            window_start = start;
            window_end = end;
            window_chrom = chrom;
            window_ploidy = ploidy;
        } else
        {
            window_start = MIN(window_start, start);
            window_end = MAX(window_end, end);
            // eprint(f"No jump! s:{start}, w_s={window_start} e:{end} w_e={window_end}")
        } // if

    } // while

    //
    // If the last iteration of the loop was not a jump, we still need to print
    //
    if (merge && !jump)
    {
        fprintf(stdout, "%s\t%d\t%d\t%d\n", window_chrom, window_start, window_end, window_ploidy);
    }

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
