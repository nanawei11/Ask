################################################################################
# APIs to AmpliconSeeK
################################################################################
import os
import sys
import pandas as pd
import numpy as np
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
import getbppair
import time
import logger

logg = logger.Logger()

#------------------------------------------------------------------------------#
def process_alignment(bamfile, gsfile = None, binsize = 10000,
                      mapq = 1, nmmax = 3,
                      mode = 'standard', sub_binsize = 1000,
                      seg_robost_quant = 0.5, SA_with_nm = True):  # ******change
    """from bam file to amplified segments and breakpoint list

    Parameters
    ----------
    bamfile : input bam file with index in the same folder
        better bwa mem aligned, coordinate sorted and duplicates marked
    gsfile : genome size file (not yet in use)
    binsize : genomic bin size to count reads in
    mapq : minimal mapping quality (include)
    nmmax : maximal number of mismatches (include)
    mode : bin read count mode
        'standard' : calculate total read counts in target bins
                     used for data with low bias, such as Input
        'bias' : calculate robust statistics in target bins
                 used for data with bias coverage, such as ChIP
    sub_binsize : smaller bin size within target bin size to
        calculate robust statistics
    seg_robost_quant : quantile as robust statistics
        if None, if None auto determine for different library size
    """
    filename = os.path.splitext(os.path.basename(bamfile))
    logg.info(f'process_alignment: {filename}')
    # extract all clip reads
    bp_all, reads_with_SA = covbam.clip_from_bam(bamfile, mapq = mapq, nmmax = nmmax, \
                                                 SA_with_nm = SA_with_nm) # ******change

    # output bedgraph
    clip_bg = covbam.clip2bedgraph(bp_all)

    # get read counts in genomic bins
    if (mode == 'standard'):
        bin_count = covbam.region_count(bamfile, gsfile,
            binsize = binsize, mapq = mapq, nmmax = nmmax)
    else:
        bin_count = covbam.region_count_2pass(bamfile, gsfile,
            binsize = binsize, mapq = mapq, nmmax = nmmax,
            sub_binsize = sub_binsize, q = seg_robost_quant)
    print(f'possible bp_pair more than {len(reads_with_SA)}')
    if len(reads_with_SA) > 0:
        bp_duo_bam =  getbppair.get_bp_pair(reads_with_SA)
    else:
        bp_duo_bam = pd.DataFrame()
    
    logg.info(f'process_alignment finished')
    ## return
    return bp_all, clip_bg, bin_count, bp_duo_bam


#------------------------------------------------------------------------------#
def detect_amplified_segment(bin_count, bp_all, blfile = None, genefile = None,
                             gsfile = None, cgfile = None, biasfile = None,
                             binsize = 10000,
                             min_cn = None, std_scale = 8,
                             min_segsize_scale = 5,
                             min_nt = 5):
    """detect amplified segments and breakpoint list

    Parameters
    ----------
    bin_count : bin counts
    blfile : blacklist file in bed format
        include regions having universal high read counts in any samples
    genefile : gene annotation file in bed12 format
    gsfile : genome size file (not yet in use)
    cgfile : cancer gene annotation file, gene symbol in 1st column
    binsize : genomic bin size to count reads in
    min_cn : minimal copy number to treat bin as amplification
        score based on assuming no whole genome duplication
        (None - use mean + std_scale * std as cutoff)
    std_scale : determin min_cn by mean + std_scale * std
    min_segsize_scale: minimal segment size factor
        min_segsize_scale * binsize
    min_nt : minimal nt required to be matched to bridge two breakend
    """
    logg.info(f'detect_amplified_segment ....')
    # get the bias in binsize and remove bias from count data
    if (biasfile is not None):
        dfm = cndetect.get_bin_bias(bin_count, biasfile, binsize)
        bin_norm = cndetect.norm_glm(dfm)
    else:
        bin_norm = bin_count

    # set blacklist bins to diploid CN
    if (blfile is not None):
        bin_norm = cndetect.bin_filter_by_blacklist(
            bin_norm, blfile, binsize = binsize)

    # smoothing on CN column
    bin_norm = cndetect.smooth_cn(bin_norm)

    # segmentation
    cn_seg = cbs.cbs(bin_norm, binsize = binsize)

    # get amplicon
    cn_amp_raw = cbs.cbs_amplicon(cn_seg, min_cn = min_cn)

    # filter out small spikes (fold of binsize)
    cn_amp_raw = cndetect.segment_filter_by_size(cn_amp_raw, \
        binsize = binsize, fold = min_segsize_scale)

    # cn_seg boundary refine
    cn_amp = cndetect.finemap_cn_segment_boundary(cn_amp_raw, bp_all, \
        binsize = binsize)

    # annotation and output
    if (genefile is not None):
        cn_amp = cndetect.segment_annotation(cn_amp, genefile)
        if (cgfile is not None):
            cn_amp = misc.map_cancergene(cn_amp, cgfile)

    logg.info(f'detect_amplified_segment finished')
    ## return
    return cn_amp, cn_seg, bin_norm


#-----------------------------------------------------------------------------#
def compute_depth(bin_norm, cn_amp, cn_seg):
    table_df1 = bin_norm.agg({'Count': ['sum', 'mean'],
                        'CN': ['sum', 'mean']}).T
    table_df1.columns = ['bin_norm_sum','bin_norm_mean']
    table_df2 = cn_amp.agg({'Count': ['sum', 'mean'],
                            'CN': ['sum', 'mean']}).T
    table_df2.columns = ['cn_amp_sum','cn_amp_mean']
    table_df3 = cn_seg.agg({'Count': ['sum', 'mean'],
                            'CN': ['sum', 'mean']}).T
    table_df3.columns = ['cn_seg_sum','cn_seg_mean']
    table_df = pd.concat([table_df1, table_df2, table_df3], axis=1)
    return table_df

#-----------------------------------------------------------------------------#
def detect_bp_pair(bamfile, bp_duo_bam, bp_all, cn_amp,
                   blfile = None, bp_min_clip = 10,
                   binsize = 10000, max_n_bp = 500, clean_bp_perc = 0.2,
                   only_keep_clean_bp = True, min_nt = 5,
                   mapq = 1, nmmax = 3, bp_pair_info_file = None, 
                   num_chrom = None, num_stratify = None, way = 'both',
                   multi_process = False, num_cores = 8):
    """ from breakpoint list to breakpoint pairs

    Parameters
    ----------
    bamfile : input bam file with index in the same folder
        better bwa mem aligned, coordinate sorted and duplicates marked
    bp_all : all breakpoints in a dict
        dict{'left' = pd_df, 'right' = pd_df}
    cn_amp : amplified segments in pd_df
    blfile : blacklist file in bed format
        include regions having universal high read counts in any samples
    bp_min_clip : minimal clipped reads of breakpoint to include
    binsize : genomic bin size
    max_n_bp : max # of breakpoints to included in the analysis
        if the sequence depth is too high, the clip points will be too many
        in this case, only the top N breakpoints will be analyzed
        if tie @ N all breakpoints <= counts @ N will be removed
    clean_bp_perc : percentage cutoff of the "CleanBP"
        CleanBP (breakpoint) is defined as the ClipDepth/InDepth >= perc
        - calculate the sequence depth on the breakpoint (InDepth)
        - # of clip reads (ClipDepth)
    only_keep_clean_bp: remove non CleanBP
    min_nt : minimal nt required to be matched to bridge two breakend
    mapq : minimal mapping quality (include)
    nmmax : maximal number of mismatches (include)

    """
    logg.info(f'detect_bp_pair ....')
    print(cn_amp.iloc[:,:6])
    cn_amp = cn_amp[['Chrom', 'Start', 'End', 'CN', 'ClipLeft', 'ClipRight', 'Gene', 'CancerGene']]
    
    if way == 'ask2' or way == 'both':
        # get bp pairs from bp_cand_stats
        if not bp_duo_bam.empty:                
            print(f'get the number of bppair from bamfile  {bp_duo_bam.shape}')
            # filter breakpoints in blacklist regions
            if (blfile is not None):
                bp_duo_bam = bpdetect.bp_filter_by_blacklist(bp_duo_bam, blfile)
            print(f'get the number of bppair after filter_by_blacklist {bp_duo_bam.shape}')

            # oxnly include the breakpoints in amplified segments
            bp_duo_bam = bpdetect.bp_filter_by_amplicon(bp_duo_bam, cn_amp, binsize)
            print(f'get the number of bppair after filter_by_amplicon {bp_duo_bam.shape}')

            # filter when the bp pair is very large
            if bp_duo_bam.shape[0] >= 20000:
                bp_duo_bam = bp_duo_bam[bp_duo_bam["Count"] > 1]
            bp_cand_bam =  bpdetect.refine_bp_cand_bam(bp_duo_bam, bp_all)  # bp_duo_bam[bp_duo_bam["Count"] > 1]
        else:
            bp_cand_bam = pd.DataFrame()
    else:
        bp_duo_bam, bp_cand_bam = pd.DataFrame(), pd.DataFrame()   

    if way == 'both' or way == 'ask':
        # get the breakpoint candidates >= "n" read counts
        bp_cand = bpdetect.bp_filter_by_n(bp_all, n = bp_min_clip)
        print(f'bp_cand after filter breakpoints  >= "n" read counts: left {bp_cand["left"].shape},  right {bp_cand["right"].shape}')

        # filter breakpoints in blacklist regions
        if (blfile is not None):
            bp_cand = bpdetect.bp_filter_by_blacklist(bp_cand, blfile, isdf = False)
        print(f'bp_cand after filter breakpoints in blacklist regions: left {bp_cand["left"].shape},  right {bp_cand["right"].shape}')

        # only include the breakpoints in amplified segments
        bp_cand = bpdetect.bp_filter_by_amplicon(bp_cand, cn_amp, binsize, isdf = False)
        print(f'bp_cand after filter breakpoints in amplified segments: left {bp_cand["left"].shape},  right {bp_cand["right"].shape}')
    
        # in this case, only the top N breakpoints will be analyzed
        bp_cand_df = bpdetect.bp_top_n(bp_cand, topn = max_n_bp, num_chrom = num_chrom, num_stratify = num_stratify)
        print(f'bp_cand_df.shape: {bp_cand_df.shape}')

        # put bp candidates into one table  ----------------------(consuming time )
        bp_cand_all = bpdetect.bp_seq_depth(bp_cand_df, \
                                            bamfile, perc = clean_bp_perc, mapq = mapq, nmmax = nmmax, isdf = False)
        # remove noise
        if (only_keep_clean_bp):
            bp_cand_stats = bp_cand_all[bp_cand_all['CleanBP'] == True]
        else:
            bp_cand_stats = bp_cand_all

        # print(f'bp_duo_bam.shape: {bp_duo_bam.shape}')
        # bp_duo_bam.to_csv('bp_duo_bam_' + way +'.csv')
    else:
        bp_cand_stats = pd.DataFrame()

    print(f'get bp pairs')
    # get bp pairs
    bp_duo = bpjoint.get_bppair(bamfile, bp_cand_stats, bp_duo_bam, min_nt = min_nt,
                                multi_process = multi_process, num_cores = num_cores)  # consuming time

    print(f'bp_duo.shape: {bp_duo.shape}')

    # merge bp cand
    try:
        bp_cand_stats.drop(columns = 'Group', inplace=True)
    except:
        print(f'bp_cand_stats no Group column or empty')  

    if not bp_cand_bam.empty:
        bp_cand_stats = pd.concat([bp_cand_stats, bp_cand_bam]).drop_duplicates(['Chrom', 'Coord', 'Clip'], keep='first')
        
    # get bp pairs from paired reads
    bp_cn = bpjoint.bp_cn_boundary(cn_amp)
    bp_duo_pe = bpjoint.get_bppair_peread(bamfile, bp_cn)

    # merge bp pairs
    bp_duo = pd.concat([bp_duo_pe, bp_duo])
    logg.info(f'detect_bp_pair finished')
    return bp_duo, bp_cand_stats


#------------------------------------------------------------------------------#
def output_bppair_alignment(bp_pair, bamfile, output_align_dir, circ_anno):
    """output breakpoint pair alignments
    """
    if (circ_anno.empty == False):
        bpjoint.ouput_alignment(bp_pair, bamfile, output_align_dir, circ_anno)

def find_circ(bp_pair, bp_fine, seg, genefile, cgfile, sefile):
    '''
    search circ based DFS.
    '''
    # create a new ggraph
    gg = ggraph.GGraph()
    
    # build ggraph from breakpoints
    gg.build_ggraph_from_bp(bp_pair, bp_fine, seg) 

    # keep thr first top 4
    second_find = False
    if (seg.shape[0] > 400) & (bp_pair.shape[0] > 400):
        print(f'seg.shape: {seg.shape}, bp_pair.shape{bp_pair.shape}')
        for i, v in gg.graph.items():
            if len(v) >= 4:
                count_neigh = [-gg.vertex[i][4] for i in v]
                gg.graph[i] = [v[i] for i in np.argsort(count_neigh)[:4]]
        second_find = True

    # get all circles in the ggraph
    all_circ_path = gg.get_all_path_contain_circle() 
    #--------------------------------------------------------------------------#
    # remove redundant
    circ_uniq = gg.path_unique(all_circ_path)

    print(f'get representitive circles in the ggraph')
    # choose a representitive circle from each circle cluster
    circ = gg.get_representitive_path(circ_uniq)

    # make interpretable circular amplicon
    circ_df = gg.make_amplicon_df(circ)

    # annotate circular amplicon
    if (genefile is not None):
        circ_anno = circ_df.assign(
            Gene = grange.GRange.map_gene(circ_df, genefile)['Gene'])
        if (cgfile is not None):
            circ_anno = misc.map_cancergene(circ_anno, cgfile)
    else:
        circ_anno = circ_df
    
    if (sefile is not None):
        circ_anno = grange.GRange.map_senchancer(circ_anno, sefile)
    return gg, all_circ_path, circ_anno, second_find

#------------------------------------------------------------------------------#
def construct_amplicon(bp_duo, bp_cand_stats, cn_amp, bin_norm, 
                       genefile = None, sefile = None, cgfile = None, 
                       min_junc_cnt = 10, bpp_min_dist = 200,
                       segment_restrict = True, subseg = True, rm_amp_df = False,
                       first_filter = False, second_filter = False, third_filter = False, knn = 3,
                       saveseg = None, savebp_bp_pair = None):  # ******change
    """construct amplicon structure

    Parameters
    ----------
    bp_duo : final breakpoint pairs used to construct amplicon
    bp_cand_stats: breakpoint candidate stats table
    cn_amp : copy number amplified segments
    genefile : gene annotation file
    cgfile : cancer gene file with first column gene symbol
    min_junc_cnt : minimal supporting read counts of breakpoint pairs
    bpp_min_dist : remove breakpoint pairs with distance <= N
    segment_restrict : True - remove non CleanBP
        False - use non CleanBP as alternative segment stop
    """
    logg.info(f'construct_amplicon ...')
    #--------------------------------------------------------------------------#
    # return empty construct if no bp pairs detected
    if (bp_duo.empty):
        colnames = ['Chrom', 'Start', 'End', 'Strand',
                    'SplitCount', 'CN', 'AmpliconID', 'Gene']
        colnames1 = ['AmpliconID', 'Chrom','Start', 'End', 'Seg_num', 'Length', 'SplitCount_sum', 'SplitCount_mean', 
                     'SplitCount_std', 'CN_sum', 'CN_mean', 'CN_std', 
                     'Left_CN', 'Right_CN', 'Gene_num', 'Cancergene_num', 'SE_num',
                     'F1','F2','F3','F4','Score']
        circ_anno = pd.DataFrame(columns = colnames)
        line_anno = pd.DataFrame(columns = colnames)
        circ_score = pd.DataFrame(columns = colnames1)
        return circ_anno, line_anno, bp_duo, cn_amp, circ_score

    #--------------------------------------------------------------------------#
    cn_amp = cn_amp[['Chrom', 'Start', 'End', 'CN', 'ClipLeft', 'ClipRight', 'Gene', 'CancerGene']]

    print(f'filter bp duo: bp_duo.shape {bp_duo.shape}')
    
    # estimate the thred for bp_pair
    if bp_duo.shape[0] >= 1000:
        _, min_junc_cnt, _, _ = bppair.estimaed_thred(bp_duo)
    else:
        min_junc_cnt = 2
    # remove breakpoint pairs in a short distance or have few split reads
    bp_pair = bppair.bp_pair_filter(bp_duo, bpp_min_dist, min_junc_cnt)
    print(f'filter bp_pair by thred: {min_junc_cnt}: bp_pair.shape {bp_duo.shape}')
    # if subseg and bpp_max_dist:
    #     print(f'only subseg is {subseg} == True')
    #     ind = (bp_pair['Chrom1']==bp_pair['Chrom2']) & (np.abs(bp_pair['Coord2'] - bp_pair['Coord1']) > bpp_max_dist)
    #     bp_pair = bp_pair[~ind]

    # breakpoint stats table
    print(f'bp_refine: bp_pair.shape {bp_pair.shape}')
    bp_fine, cn_amp_merged = bppair.bp_refine(bp_pair, bp_cand_stats, cn_amp)

    seg = bppair.get_segment(bp_fine, bp_pair, cn_amp_merged, 
                             restrict = segment_restrict, subseg = subseg, rm_amp_df = rm_amp_df, knn = knn) # ******change 
    # add cn to segments
    seg = bppair.add_cn(seg, cn_amp, bin_norm)
    # if saveseg:
    #     seg.to_csv(saveseg, sep='\t', index=False)
    # if savebp_bp_pair:
    #     bp_pair.to_csv(savebp_bp_pair, sep='\t', index=False)

    #--------------------------------------------------------------------------#
    print(f'get all circles in the ggraph')
    gg, all_circ_path, circ_anno, second_find = find_circ(bp_pair, bp_fine, seg, genefile, cgfile, sefile)

    # compute fretures of  circs and  take the better circs based on scores
    if not circ_anno.empty:
        circ_score = ggraph.add_stats_circ(circ_anno, bin_norm, binsize = 10000)
    
        result = (circ_score.groupby(['Chrom','Start','End'])
                    .apply(lambda x: x.nlargest(1, 'Score'))
                    .reset_index(drop=True))
        select_circ = result['AmpliconID'][(result['Score'] >= result['Score'].quantile(0.85)) & (result['F1'] >= 1)]
        print(f'select_circ: {select_circ}')

        # search the circ again on specfic chroms 
        circ_new = pd.DataFrame()
        if second_find:
            cand_serch_circ = circ_score[circ_score['AmpliconID'].isin(select_circ)][['AmpliconID', 'Chrom', 'Start', 'End']]
            select_chrom = cand_serch_circ['Chrom'].tolist()
            
            bp_pair1 = bp_pair[(bp_pair['Chrom1'].isin(select_chrom)) & (bp_pair['Chrom2'].isin(select_chrom))]
            seg1 = seg[seg['Chrom'].isin(select_chrom)]
            bp_fine1 = bp_fine[bp_fine['Chrom'].isin(select_chrom)]
            _, _, circ_anno1, _ = find_circ(bp_pair1, bp_fine1, seg1, genefile, cgfile, sefile)
            
            # merge the circ of two search
            query_circ_circ = np.unique(circ_anno['AmpliconID'])
            query_circ_list = []
            for i in query_circ_circ:
                query_circ_list_ = []
                for i in circ_anno[ circ_anno['AmpliconID'] == i ].itertuples():
                    query_circ_list_.append(i[1:7])
                query_circ_list.append(query_circ_list_)

            circ_cand = []
            circ_id = np.unique(circ_anno1['AmpliconID'])
            circ_num = len(query_circ_circ)
            for i in circ_id:
                circ_cand_ = []
                circ_df_sub = circ_anno1[ circ_anno1['AmpliconID'] == i ]
                for j in circ_df_sub.itertuples():
                    circ_cand_.append(j[1:7])
                if circ_cand_ not in query_circ_list:
                    circ_reindex = 'circ_' + str(circ_num)
                    circ_df_sub['AmpliconID'] = circ_reindex
                    circ_cand.append( circ_df_sub )
                    circ_num += 1
            if circ_cand:
                circ_new = pd.concat(circ_cand)
                # compute fretures of  circs and  take the better circs based on scores
                circ_new_score = ggraph.add_stats_circ(circ_new, bin_norm, binsize = 10000)
            else:
                circ_new = pd.DataFrame()

        if not circ_new.empty:
            circ_anno = pd.concat([circ_anno, circ_new])
            circ_score = pd.concat([circ_score, circ_new_score])
    else:
        circ_score = pd.DataFrame()

    #--------------------------------------------------------------------------#
    # remove redundant
    line_uniq = gg.path_unique(all_circ_path, type = 'linear')
    print(f'get representitive line in the ggraph')
    # choose a representitive one
    line = gg.get_representitive_path(line_uniq, type = 'linear')

    # make interpretable amplicon
    line_df = gg.make_amplicon_df(line, type = 'linear')

    # annotate circular amplicon
    if (genefile is not None):
        line_anno = line_df.assign(
            Gene = grange.GRange.map_gene(line_df, genefile)['Gene'])
        if (cgfile is not None):
            line_anno = misc.map_cancergene(line_anno, cgfile)
    else:
        line_anno = line_df
    
    if (sefile is not None):
        line_anno = grange.GRange.map_senchancer(line_anno, sefile)

    logg.info(f'construct_amplicon finished')
    return circ_anno, line_anno, bp_pair, seg, circ_score

#------------------------------------------------------------------------------#
def plot_amplicon(circ_anno, line_anno = None, cn_amp = None,
                  genefile = None, sefile = None, bincnt = None, binsize = 10000,
                  ext = 0.3, fig_dir = 'plot', plot_n = 5,
                  fig_width = 15, fontsize = 12):
    """plot amplicon structure

    Parameters
    ----------
    circ_anno : circular amplicon dataframe
    line_anno : linear amplicon dataframe
    cn_amp : amplified segments dataframe
    genefile : gene annotation file
    sefile: super enhancer annotation file *********************  no set
    bincnt : bin read counts per binsize genomic region
    ext : extent roi by 30 percent
    plot_dir : plot folder
    plot_n : only plot n structures for linear or cn
        Note : already plot all structures for circular
    """

    #--------------------------------------------------------------------------#
    # readin bed12 gene annotation file
    if (genefile is not None):
        genebed = trackplot.read_bedf(genefile)
        genebed = [g for g in genebed if (g.name != '')]
    
    # readin bed4 super encancer annotation file
    if (sefile is not None):
        sebed = trackplot.read_bedf(sefile)
        sebed = [g for g in sebed if (g.name != '')]

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
