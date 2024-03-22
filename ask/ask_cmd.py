#!/usr/bin/env python
import argparse
import os
import time
from pathlib import Path
import pickle
import pandas as pd
import ask


#------------------------------------------------------------------------------#
def get_args():
    """ input arguments
    """
    parser = argparse.ArgumentParser(
        description='Seek amplicons from ChIP-seq data')

    parser.add_argument('-i', '--bamfile', required=True,
        help="input bam file, also need the index file in the same folder, \
              which is generated by: samtools index <in.bam>. \
              Recommended to be deduplicated by: samtools markdup")
    parser.add_argument('-g', '--genome', default='hg38',
        help="genome work on (hg19 or hg38, default: hg38)")
    parser.add_argument('-o', '--outprefix',
        help="output file prefix, could include folder")
    parser.add_argument('-b', '--bwfile',
        help="input bigwig file, not currently in use")  # ?????????
    parser.add_argument('-w', '--binsize', type=int, default=10000,
        help="bin size for copy number segmentation")
    parser.add_argument('-c', '--mincn', type=float, default=5,
        help="minimal absolute copy number of bins to be considered \
              as amplicons; if None, cutoff set as mean +- 2*sd")
    parser.add_argument('-k', '--mapq', type=int, default=20,
        help="minimal mapping quality (include)")
    parser.add_argument('-l', '--nmmax', type=int, default=1,
        help="maximal number of mismatches (include)")

    parser.add_argument('-n', '--bpcount', type=int, default=5,
        help="# clip reads for a breakpoint to include in analysis")
    parser.add_argument('-m', '--juncread', type=int, default=10,
        help="# clip reads for a breakpoint pair to include in analysis")

    parser.add_argument('--max_n_breakpoint', type=int, default=500,
        help="max # of breakpoints to included in the analysis, if the \
            sequence depth is too high, the clip points will be too many \
            in this case, only the top N breakpoints will be analyzed \
            if tie @ N all breakpoints <= counts @ N will be removed")

    parser.add_argument('-d', '--segmode', default='standard',
        help='''
            segmentation mode --
            'standard' : calculate raw read counts in target bins,
                         used for data with low coverage bias,
                         such as ChIP-seq input, WGS.
            'bias' : calculate smoothed read counts in target bins,
                     used for data with coverage bias,
                     such as ChIP-seq test, ATAC-seq, Exome-seq.
            ''')
    parser.add_argument('-s', '--sub-binsize', type=int, default=1000,
        help="sub bin size, for bias segmentation mode only")
    parser.add_argument('-q', '--segquant', default=0.5,
        help="calculate robust read counts in target genomic bins \
              by using quantiles of read counts in smaller bins, \
              for bias segmentation mode only")

    parser.add_argument('--cleanbp-cutoff', type=float, default=0.2,
        help="percentage cutoff of the CleanBP (breakpoint), \
              defined as the ClipDepth/Depth >= perc")  

    parser.add_argument('--only-keep-cleanbp', default=True,
        help="remove non CleanBP for segmentation")

    parser.add_argument('--SA_with_nm', action='store_true', required=False,
        help="only consider the reads with SA nm greater than 1")
    
    parser.add_argument('--subseg', action='store_true', required=False,
        help="infer the sub segment on amp")
    
    parser.add_argument('--num_chrom', default = None,
        help="infer the sub segment on amp")
    
    parser.add_argument('--method', default = 'both',
        help="infer the sub segment on amp")
    
    parser.add_argument('--rm_amp_df', action='store_true', required = False,
        help="infer the sub segment on amp")
    
    # parser.add_argument('--seg_len', type=int, default=9000000,
    # help="infer the sub segment on amp")

    parser.add_argument('--run_from_pdat', type = int, default = 0,
        help="run ask from saved pdat, resume run from k + 1 step")   #?????????????
    
    parser.add_argument('--knn', type = int, default = 3,
        help="run ask from saved pdat, resume run from k + 1 step")   #?????????????

    parser.add_argument('--first_filter', action='store_true', required=False,
        help="infer the sub segment on amp")
    
    parser.add_argument('--second_filter', action='store_true', required=False,
        help="infer the sub segment on amp")
    
    parser.add_argument('--third_filter', action='store_true', required=False,
        help="infer the sub segment on amp")
    
    return parser.parse_args()


#------------------------------------------------------------------------------#
def process_args(args):
    """ post process of arguments
    """
    #@ output file prefix
    if args.outprefix is None:
        args.outprefix = args.bamfile.replace('.bam', '')

    #@ output plot folder
    args.outfig = args.outprefix + '_ask_plot'

    #@ set genome related files
    if args.genome == 'hg19':
        # blacklist file
        args.blfile = os.path.join(
            os.path.dirname(__file__), '..', 'data', 'hg19_blacklist.bed')  #

        # gene annotation file
        args.genefile = os.path.join(
            os.path.dirname(__file__), '..', 'data', 'hg19_refgene_process.bed12')
        
        # super encanher annotation file
        args.sefile = os.path.join(
            os.path.dirname(__file__), '..', 'data', 'se_hg19_sort.bed')

        # genome size file (currently not in use)
        args.gsfile = os.path.join(
            os.path.dirname(__file__), '..', 'data', 'hg19.genome')

        # cancer gene file
        args.cgfile = os.path.join(
            os.path.dirname(__file__), '..', 'data',
                'Census_all_20200624_14_22_39.tsv')

        # bias file (GC + mappability)
        args.biasfile = os.path.join(
            os.path.dirname(__file__), '..', 'data',
                'hg19_bias.bed.gz')  # ????????????

    elif args.genome == 'hg38':
        # blacklist file
        args.blfile = os.path.join(
            os.path.dirname(__file__), '..', 'data', 'hg38_blacklist.bed')

        # gene annotation file
        args.genefile = os.path.join(
            os.path.dirname(__file__), '..', 'data', 'hg38_refgene_process.bed12')

        # super encanher annotation file
        args.sefile = os.path.join(
            os.path.dirname(__file__), '..', 'data', 'se_hg38_sort.bed')

        # genome size file (currently not in use)
        args.gsfile = os.path.join(
            os.path.dirname(__file__), '..', 'data', 'hg38.genome')

        # cancer gene file
        args.cgfile = os.path.join(
            os.path.dirname(__file__), '..', 'data',
                'Census_all_20200624_14_22_39.tsv')

        # bias file (GC + mappability)
        args.biasfile = os.path.join(
            os.path.dirname(__file__), '..', 'data',
                'hg38_bias.bed.gz')

    else:
        args.blfile = None
        args.genefile = None
        args.gsfile = None
        args.cgfile = None
        print('Invalid genome build provided -- \
            No blacklist filtering applied, no gene annotation applied.')

    #@ create base dir of output files
    Path(os.path.dirname(args.outprefix)).mkdir(parents=True, exist_ok=True)

    #@ create base dir of plot files
    Path(args.outfig).mkdir(parents=True, exist_ok=True)

    #@ set output file names
    args.output_stats     = args.outprefix + '_ask_stats.tsv'
    args.output_cnamp     = args.outprefix + '_ask_amplified_segment.tsv'
    args.output_cnseg     = args.outprefix + '_ask_cn_segmentation.tsv'

    args.output_bpcand    = args.outprefix + '_ask_breakpoint.tsv'
    args.output_bpduo     = args.outprefix + '_ask_breakpoint_pair_raw.tsv'
    args.output_bppair    = args.outprefix + '_ask_breakpoint_pair.tsv'
    args.output_seg       = args.outprefix + '_ask_breakpoint_seg.tsv'
    args.output_circ      = args.outprefix + '_ask_amplicon_circular.tsv'
    args.output_line      = args.outprefix + '_ask_amplicon_linear.tsv'
    args.output_clip      = args.outprefix + '_ask_clip_count.bedgraph'
    args.output_bincnt    = args.outprefix + '_ask_bin_count.tsv'
    args.output_binnorm   = args.outprefix + '_ask_bin_count_norm.tsv'
    args.output_align_dir     = args.outprefix + '_ask_junctionseq'
    args.output_circ_stat = args.outprefix + '_ask_amplicon_circular_stat.tsv'

    args.output_pdat_1    = args.outprefix + '_ask_step1.pdat'
    args.output_pdat_2    = args.outprefix + '_ask_step2.pdat'
    args.output_pdat_3    = args.outprefix + '_ask_step3.pdat'
    args.output_pdat_4    = args.outprefix + '_ask_step4.pdat'

    return args


#------------------------------------------------------------------------------#
def main():
    #--------------------------------------------------------------------------#
    # set start time
    #--------------------------------------------------------------------------#
    startTime = time.time()


    #--------------------------------------------------------------------------#
    # get arguments
    #--------------------------------------------------------------------------#
    args = get_args()
    args = process_args(args)


    #--------------------------------------------------------------------------#
    # load pdat
    #--------------------------------------------------------------------------#
    if int(args.run_from_pdat) >= 1 and os.path.isfile(args.output_pdat_1):
        try:
            bp_all, clip_bg, bin_count, bp_duo_bam = pd.read_pickle(args.output_pdat_1)
        except:
            bp_all, clip_bg, bin_count = pd.read_pickle(args.output_pdat_1)

    if int(args.run_from_pdat) >= 2 and os.path.isfile(args.output_pdat_2):
        # with open(args.output_pdat_2, "rb") as f:
        cn_amp, cn_seg, bin_norm = pd.read_pickle(args.output_pdat_2)
        # with open('/cluster/home/WeiNa/project/ecDNA/result/ChIP-seq-ask2/ask3_new/ask3_both_nm1_seg1_GBM39_1013.txt','a') as f:
        #     f.write(f'step2 data already read.\n')        
        # print(f'step2 data already read.')
    if int(args.run_from_pdat) >= 3 and os.path.isfile(args.output_pdat_3):
        # with open(args.output_pdat_3, "rb") as f:
        bp_duo, bp_cand_stats = pd.read_pickle(args.output_pdat_3)
        # with open('/cluster/home/WeiNa/project/ecDNA/result/ChIP-seq-ask2/ask3_new/ask3_both_nm1_seg1_GBM39_1013.txt','a') as f:
        #     f.write(f'step3 data already read.\n')     
    if int(args.run_from_pdat) >= 4 and os.path.isfile(args.output_pdat_4):
        # with open(args.output_pdat_4, "rb") as f:
        circ_anno, line_anno, bp_pair, seg = pd.read_pickle(args.output_pdat_4)


    #--------------------------------------------------------------------------#
    # step 1 : process bam file to generate bin counts and clip reads
    #--------------------------------------------------------------------------#
    if 'bp_all' not in locals(): # and not os.path.exists(args.output_pdat_1)
        # run module
        bp_all, clip_bg, bin_count, bp_duo_bam = ask.process_alignment(
            args.bamfile, gsfile = args.gsfile, binsize = args.binsize,
            mapq = args.mapq, nmmax = args.nmmax,
            mode = args.segmode, sub_binsize = args.sub_binsize,
            seg_robost_quant = args.segquant, SA_with_nm = args.SA_with_nm)

        # impossible = bp_all['left'].empty and bp_all['right'].empty
        # assert impossible==True "This data do not have split reads, therefore not be deteced ecDNA!"

        # output results
        bin_count.to_csv(args.output_bincnt, sep='\t', index=False)
        clip_bg.to_csv(args.output_clip, sep='\t', index=False, header=False)

        # save pdat
        pdat = [bp_all, clip_bg, bin_count, bp_duo_bam]
        with open(args.output_pdat_1, "wb") as f:
            pickle.dump(pdat, f)

        # output log
        print('alignment process - done \
            - {0:0.1f} seconds\n'.format(time.time() - startTime))
    # else:
    #     with open(args.output_pdat_1, 'rb') as f:
    #         data = pickle.load(f)
    # bp_all = data[0]
    # clip_bg = data[1]
    # bin_count = data[2]


    #--------------------------------------------------------------------------#
    # step 2 : identify amplified segments
    #--------------------------------------------------------------------------#
    if 'cn_amp' not in locals(): #and not os.path.exists(args.output_pdat_2)
        # run module
        cn_amp, cn_seg, bin_norm = \
            ask.detect_amplified_segment(
                bin_count, bp_all,
                args.blfile, args.genefile, args.gsfile, args.cgfile,
                biasfile = args.biasfile,
                binsize = args.binsize, min_cn = args.mincn)

        # output results
        cn_amp.to_csv(args.output_cnamp, sep='\t', index=False)
        cn_seg.to_csv(args.output_cnseg, sep='\t', index=False)
        bin_norm.to_csv(args.output_binnorm, sep='\t', index=False)

        # save pdat
        pdat = [cn_amp, cn_seg, bin_norm]
        with open(args.output_pdat_2, "wb") as f:
            pickle.dump(pdat, f)

        # output log
        print('detect amplified segments - done \
            - {0:0.1f} seconds\n'.format(time.time() - startTime))
    # else:
    #     with open(args.output_pdat_2, 'rb') as f:
    #         data = pickle.load(f)
    #     cn_amp = data[0]
    #     cn_seg = data[1]
        # bin_count = data[2]
    
    #--------------------------------------------------------------------------#
    # step 3 : identify breakpoint pairs
    #--------------------------------------------------------------------------#
    if 'bp_duo' not in locals():
        # run module
        bp_duo, bp_cand_stats = ask.detect_bp_pair(
            args.bamfile, bp_duo_bam,  bp_all, cn_amp, args.blfile,
            binsize = args.binsize, bp_min_clip = args.bpcount,
            mapq = args.mapq, nmmax = args.nmmax,
            max_n_bp = args.max_n_breakpoint,
            clean_bp_perc = args.cleanbp_cutoff,
            only_keep_clean_bp = args.only_keep_cleanbp,
            bp_pair_info_file = args.output_align, num_chrom = args.num_chrom,
            way = args.method)  # bp_pair_info_file can remove

        # output results
        bp_cand_stats.to_csv(args.output_bpcand, sep='\t', index=False)
        bp_duo.to_csv(args.output_bpduo, sep='\t', index=False)

        # save pdat
        pdat = [bp_duo, bp_cand_stats]
        with open(args.output_pdat_3, "wb") as f:
            pickle.dump(pdat, f)

        # output log
        print('detect breakpoint pairs - done \
            - {0:0.1f} seconds\n'.format(time.time() - startTime))


    #--------------------------------------------------------------------------#
    # step 4 : construct amplicons
    #--------------------------------------------------------------------------#
    if 'circ_anno' not in locals():
        # run module
        circ_anno, line_anno, bp_pair, seg, circ_score = ask.construct_amplicon(
            bp_duo, bp_cand_stats, cn_amp, bin_norm, args.genefile, args.sefile, args.cgfile, 
            segment_restrict = False,
            min_junc_cnt = args.juncread, subseg = args.subseg, rm_amp_df = args.rm_amp_df,
            first_filter = args.first_filter, second_filter = args.second_filter, third_filter = args.third_filter,
            knn = args.knn, saveseg = args.output_seg, savebp_bp_pair = args.output_bppair)

        # output results
        circ_anno.to_csv(args.output_circ, sep='\t', index=False)
        line_anno.to_csv(args.output_line, sep='\t', index=False)
        bp_pair.to_csv(args.output_bppair, sep='\t', index=False)
        seg.to_csv(args.output_seg, sep='\t', index=False)

        # # output breakpoint pair alignment for all bp in bp_duo
        # ask.output_bppair_alignment(bp_duo, args.bamfile, args.output_align)

        # output breakpoint pair alignment
        # output_bppair_alignment(bp_pair, bamfile, output_align_dir, circ_anno)
        ask.output_bppair_alignment(bp_pair, args.bamfile, args.output_align_dir, circ_anno)

        # save pdat
        pdat = [circ_anno, line_anno, bp_pair, seg]
        with open(args.output_pdat_4, "wb") as f:
            pickle.dump(pdat, f)

        # output log
        print('construct amplicons - done \
            - {0:0.1f} seconds\n'.format(time.time() - startTime))

    if not circ_anno.empty:
        try:
            # circ_stat = ask.score_circ(circ_anno, bin_norm, binsize = 10000)
            circ_score.to_csv(args.output_circ_stat, sep='\t', index=False)
        except:
            print('ask.score_circ - error \
            - {0:0.1f} seconds\n'.format(time.time() - startTime))
    #--------------------------------------------------------------------------#
    # step 5 : plot amplicons
    #--------------------------------------------------------------------------#
    # run module
    try:
        ask.plot_amplicon(circ_anno, line_anno, cn_amp, args.genefile, args.sefile, 
                          bin_norm, binsize = args.binsize,
                          fig_dir = args.outfig, plot_n = 5,
                          fig_width = 15, fontsize = 12, ext = 0.3)
        # output log
        print('plot amplicons - done \
            - {0:0.1f} seconds\n'.format(time.time() - startTime))
    except:
        print('plot amplicons - error \
            - {0:0.1f} seconds\n'.format(time.time() - startTime))


    #--------------------------------------------------------------------------#
    # output stats
    #--------------------------------------------------------------------------#
    # output statistics table
    with open(args.output_stats, 'w') as f:
        f.write('# of amplified segments: ' + str(len(cn_amp)) + '\n')
        f.write('# of candidate breakpoints: ' + str(len(bp_cand_stats)) + '\n')
        f.write('# of breakpoint pairs: ' + str(len(bp_duo)) + '\n')
        f.write('# of circular amplicons: '
            + str(circ_anno['AmpliconID'].nunique()) + '\n')
        f.write('# of linear amplicons: '
            + str(line_anno['AmpliconID'].nunique()) + '\n')

    print('All - done - {} seconds\n'.format(time.time() - startTime))
    if args.subseg:
        print(f'seg is {args.subseg}')
    else:
        print(f'seg is {args.subseg}')
        print(0)


#------------------------------------------------------------------------------#
if __name__ == '__main__':
    main()
