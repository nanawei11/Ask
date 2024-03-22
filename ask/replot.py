################################################################################
# APIs to AmpliconSeeK
################################################################################
import os
import pandas as pd
import covbam
import cndetect
import bpdetect
import bpjoint
import bppair
import ggraph
import grange
import trackplot
import misc
import cbs

#------------------------------------------------------------------------------#
def plot_amplicon(circ_anno, line_anno = None, cn_amp = None,
                  genefile = None, bincnt = None, binsize = 10000,
                  ext = 0.3, fig_dir = 'plot', plot_n = 5,
                  fig_width = 15, fontsize = 12):
    """plot amplicon structure

    Parameters
    ----------
    circ_anno : circular amplicon dataframe
    line_anno : linear amplicon dataframe
    cn_amp : amplified segments dataframe
    genefile : gene annotation file
    bincnt : bin read counts per binsize genomic region
    ext : extent roi by 30 percent
    plot_dir : plot folder
    plot_n : only plot n structures for linear or cn
        Note : already plot all structures for circular
    """

    #--------------------------------------------------------------------------#
    # readin bed12 gene annotation file
    try:
        genebed = trackplot.read_bedf(genefile)
    except:
        print('cannot load')
    genebed = [g for g in genebed if (g.name != '')]

    #--------------------------------------------------------------------------#
    # plot amplified segments
    if (cn_amp is not None and not cn_amp.empty):
        # sort cn_amp by copy number gain
        cn_amp_list = trackplot.get_rois_with_score(cn_amp, binsize * 10)
        cn_amp_list = sorted(cn_amp_list, key=lambda x: x[3], reverse = True)

        # plot
        for tag, roi in enumerate(cn_amp_list[:plot_n]):
            # get sub dataframe
            df_sub = trackplot.get_sub_df(cn_amp, roi)

            # output pdf file name
            fig_out = os.path.join(fig_dir, 'ampseg_' + str(tag) + '.pdf')

            # prepare plot
            rois, segs, links, bgs, beds = trackplot.prepare_plot(
                df_sub, genebed, bincnt, binsize, ext, 'cn')

            # make plot
            trackplot.plot_amp(rois, segs, links, bgs, beds, fig_out,
                               fig_width, fontsize)

    #--------------------------------------------------------------------------#
    # plot circular amplicon
    if (not circ_anno.empty):
        print(circ_anno)
        for tag, circ_sub in circ_anno.groupby(['AmpliconID']):

            # output pdf file name
            try:
                fig_out = os.path.join(fig_dir, 'circular_' + tag + '.pdf')
            except:
                fig_out = os.path.join(fig_dir, 'circular_' + tag[0] + '.pdf')
            # prepare plot
            circ_sub_copy = circ_sub.copy()
            rois, segs, links, bgs, beds = trackplot.prepare_plot(
                circ_sub_copy, genebed, bincnt, binsize, ext, 'circular')

            # make plot
            trackplot.plot_amp(rois, segs, links, bgs, beds, fig_out,
                               fig_width, fontsize)

    #--------------------------------------------------------------------------#
    # plot linear amplicon
    if (line_anno is not None and not line_anno.empty):
        id2plot = list(
            line_anno.groupby('AmpliconID').CN.max().\
                sort_values(ascending = False)[0:plot_n].index)
        for tag in id2plot:
            # output pdf file name
            try:
                fig_out = os.path.join(fig_dir, 'Linear_' + tag + '.pdf')
            except:
                fig_out = os.path.join(fig_dir, 'Linear_' + tag[0] + '.pdf')
            # get the sub dataframe
            line_sub = line_anno[line_anno['AmpliconID'] == tag]

            # prepare plot
            line_sub_copy = line_sub.copy()
            rois, segs, links, bgs, beds = trackplot.prepare_plot(
                line_sub_copy, genebed, bincnt, binsize, ext, 'linear')

            # make plot
            trackplot.plot_amp(rois, segs, links, bgs, beds, fig_out,
                               fig_width, fontsize)

if __name__ == '__main__':
    datapath = '/cluster/home/WeiNa/project/ecDNA/ask_run/data/'
    resultpath = '/cluster/home/WeiNa/project/ecDNA/ask_run/test_result/GBM39_chr7_all'
    resfigpath = '/cluster/home/WeiNa/project/ecDNA/ask_run/test_result/GBM39_chr7_all/GBM39_chr7_all_ask_plot/'
    from pathlib import Path
    import os
    import time
    from pathlib import Path
    import pickle
    #------------------------------------------------------------------------#
    #@ create base dir of output files
    Path(os.path.dirname(resultpath)).mkdir(parents=True, exist_ok=True)

    #@ create base dir of plot files
    Path(resfigpath).mkdir(parents=True, exist_ok=True)

    #@ set output file names
    output_stats     = resultpath + '_ask_stats.tsv'
    output_cnamp     = resultpath + '_ask_amplified_segment.tsv'
    output_cnseg     = resultpath + '_ask_cn_segmentation.tsv'

    output_bpcand    = resultpath + '_ask_breakpoint.tsv'
    output_bpduo     = resultpath + '_ask_breakpoint_pair_raw.tsv'
    output_bppair    = resultpath + '_ask_breakpoint_pair.tsv'
    output_circ      = resultpath + '_ask_amplicon_circular.tsv'
    output_line      = resultpath + '_ask_amplicon_linear.tsv'
    output_clip      = resultpath + '_ask_clip_count.bedgraph'
    output_bincnt    = resultpath + '_ask_bin_count.tsv'
    output_binnorm   = resultpath + '_ask_bin_count_norm.tsv'
    output_align     = resultpath + '_ask_breakpoint_pair_alignment.tsv'

    output_pdat_1    = resultpath + '_ask_step1.pdat'
    output_pdat_2    = resultpath + '_ask_step2.pdat'
    output_pdat_3    = resultpath + '_ask_step3.pdat'
    output_pdat_4    = resultpath + '_ask_step4.pdat'


    #-------------------------------------------------------------------------#
    #-------------------------------------------------------------------------#

    # blacklist file
    blfile = os.path.join(datapath,  'hg38_blacklist.bed')

    # gene annotation file
    genefile = os.path.join(
        datapath, 'hg38_refgene.bed12.gz')

    # genome size file (currently not in use)
    gsfile = os.path.join(
        datapath, 'hg38.genome')

    # cancer gene file
    cgfile = os.path.join(
        datapath, 'Census_all_20200624_14_22_39.tsv')

    bamfile = os.path.join(
        datapath, 'GBM39_chr7.bam')

    # bias file (GC + mappability)
    biasfile = os.path.join(
            datapath,
        'hg38_bias.bed.gz')

    binsize = 10000
    mapq = 20
    nmmax = 1
    segmode = 'standard'
    # segmode = 'bias'
    sub_binsize = 1000
    segquant = 0.5
    mincn = 5
    bpcount = 2
    juncread = 5
    max_n_breakpoint = 500
    cleanbp_cutoff = 0.2
    sub_binsize = 1000
    only_keep_cleanbp = True

    startTime = time.time()
    filenames_ = ['GBM39_chr7_all_ask_step1.pdat',
                    'GBM39_chr7_all_ask_step2.pdat',
                    'GBM39_chr7_all_ask_step3.pdat',
                    'GBM39_chr7_all_ask_step4.pdat']

    filepath_ ='/cluster/home/WeiNa/project/ecDNA/ask_run/test_result/GBM39_chr7_all/'
    data = {}
    for fn in filenames_:
        with open(filepath_+fn, 'rb') as f:
            data[fn] = pickle.load(f)
    circ_anno = data['GBM39_chr7_all_ask_step4.pdat'][0]
    line_anno = data['GBM39_chr7_all_ask_step4.pdat'][1]
    cn_amp = data['GBM39_chr7_all_ask_step2.pdat'][0]
    bin_norm = data['GBM39_chr7_all_ask_step2.pdat'][2]
    #--------------------------------------------------------------------------#
    # step 5 : plot amplicons
    #--------------------------------------------------------------------------#
    # run module∂
    outfig = resfigpath

    # try:
    plot_amplicon(circ_anno, line_anno, cn_amp, genefile,
                        bin_norm, binsize = binsize,
                        fig_dir = outfig, plot_n = 5,
                        fig_width = 15, fontsize = 12, ext = 0.3)
    # output log
    print('plot amplicons - done \
        - {0:0.1f} seconds\n'.format(time.time() - startTime))
    # except:
        # print('plot amplicons - error \
        #     - {0:0.1f} seconds\n'.format(time.time() - startTime))


    #--------------------------------------------------------------------------#
    # output stats
    #--------------------------------------------------------------------------#
    # output statistics table
    # with open(output_stats, 'w') as f:
    #     f.write('# of amplified segments: ' + str(len(cn_amp)) + '\n')
    #     f.write('# of candidate breakpoints: ' + str(len(bp_cand_stats)) + '\n')
    #     f.write('# of breakpoint pairs: ' + str(len(bp_duo)) + '\n')
    #     f.write('# of circular amplicons: '
    #         + str(circ_anno['AmpliconID'].nunique()) + '\n')
    #     f.write('# of linear amplicons: '
    #         + str(line_anno['AmpliconID'].nunique()) + '\n')
    #
    # print('All - done - {} seconds\n'.format(time.time() - startTime))
