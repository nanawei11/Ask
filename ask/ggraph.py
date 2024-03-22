################################################################################
# GGraph class: search for circular amplicons
################################################################################
#------------------------------------------------------------------------------#
from collections import Counter, defaultdict
import numpy as np
import pandas as pd


#------------------------------------------------------------------------------#
import misc



#------------------------------------------------------------------------------#
# GGraph class: search for circular amplicons
#------------------------------------------------------------------------------#
class GGraph:

    # all nodes must be 0, 1, 2, 3, 4,...
    #--------------------------------------------------------------------------#
    def __init__(self):
        self.graph = defaultdict(list)
        self.vertex = defaultdict(list)
        self.edges = defaultdict(list)
        self.node_dict = {}
        self.single_segment_loop = []

    #--------------------------------------------------------------------------#
    def add_edge(self, a, b, weight):
        if b not in self.graph[a]:
            self.graph[a].append(b)
        if a not in self.graph[b]:
            self.graph[b].append(a)
        self.edges[(a, b)] = weight
        self.edges[(b, a)] = weight

    #--------------------------------------------------------------------------#
    def add_vertex(self, vertex_key, vertex_value):
        self.vertex[vertex_key] = vertex_value

    #--------------------------------------------------------------------------#
    def get_nodes(self):
        return list(self.graph.keys())
    
    #--------------------------------------------------------------------------#
    def dfs_search_circle(self, v, visited, parent, path, allpath, nodebefore):
        """
        depth-first search for circles in a graph

        Note: make visited and path have the same nodes
        """

        # mark current node as visited
        visited[v] = True

        # store current node to list
        path.append(v)

        # loop for adjacent nodes of the current node
        for v_adj in self.graph[v]:
            if v_adj not in nodebefore:
                # edge type is different from the last one
                # one bp_pair edge and one segment edge
                if (parent == -1 \
                    or self.edges[(parent, v)][0] != self.edges[(v, v_adj)][0]):
                    if (visited[v_adj] == False): # recur if the node is not visited
                        self.dfs_search_circle(v_adj, visited, v, path, allpath, nodebefore)
                    # cycle detected if the adjacent node is visited
                    # and is not the parent of current node
                    elif (parent != -1 and parent != v_adj): #and v_adj in path
                        path.append(v_adj) # add v_adj, it's the loop start
                        # store when loop detected
                        allpath.append(tuple(path)) # must use tuple! or can't run
                        path.pop() # remove v_adj from path
                        # file3 = '/cluster/home/WeiNa/project/ecDNA/result/ChIP-seq-ask2/ask3_1020/ask3_both_knn3_mr10_GBM39_1020_nofilter_circ/all_path.txt'
                        # with open(file3,'a') as f:
                        #     f.write(f'circ\t{tuple(path)}\n')
                else:
                    if (len(self.graph[v]) == 1): # end of a path, linear amplicon
                        # store when to the branch end
                        allpath.append(tuple(path)) # must use tuple! or can't run
                        # file3 = '/cluster/home/WeiNa/project/ecDNA/result/ChIP-seq-ask2/ask3_1020/ask3_both_knn3_mr10_GBM39_1020_nofilter_circ/all_path.txt'
                        # with open(file3,'a') as f:
                        #     f.write(f'linear\t{tuple(path)}\n')

        # remove current edge from path and mark it as unvisited
        path.pop()
        visited[v] = False

        return allpath


    #--------------------------------------------------------------------------#
    def get_all_path_contain_circle(self):
        """
        get all paths containing circle
        also include paths not containing any circle
        """

        # init onjects
        path = []
        allpath = []
        node = sorted(self.get_nodes())

        # mark all nodes as not vsisited
        visited = [False]*(max(node) + 1)

        # record nodes has already detected circle
        nodebefore = []
        # detect circle in different subgraph
        # say until all nodes are visisted
        for i in node:
            if visited[i] == False:
                self.dfs_search_circle(i, visited, -1, path, allpath, nodebefore)
            # nodebefore.append(i)
        return allpath


    #--------------------------------------------------------------------------#
    # amplicon processing functions
    #--------------------------------------------------------------------------#
    def path_unique(self, circ, type = 'circular'):
        """
        get unique paths of given list of paths

        type: circular or linear
        """
        # get all circle and unique
        if (type == 'circular'):
            circ = misc.unique([row[row.index(row[-1]):-1] for row in circ])
            circ = [x for x in circ if x] # remove empty path
        elif (type == 'linear'):
            circ = [row for row in circ if (row[-1] not in row[:-1])]

        # keep ones both start and end with inner
        circ = [row for row in circ \
            if (self.edges[row[:2]][0] == 'inner' and \
                self.edges[row[-2:]][0] == 'inner')]

        if (type == 'circular'):
            # rotate path to smallest in the first
            circ = [self.rotate_path(i) for i in circ]

            # invert path to make second one smaller than last one
            circ = misc.unique([
                self.invert_path(i) if (i[1] > i[-1]) else i for i in circ])
        elif (type == 'linear'):
            circ = misc.unique([
                i[::-1] if (i[0] > i[-1]) else i for i in circ])

        return circ

    #--------------------------------------------------------------------------#
    def path_longest(self, circ, type = 'circular', check_chrom = False):
        """
        get the longest path

        if multiple longest path,
        choose the one with max edge weight

        if multiple max edge weight,
        choose the first one
        """
        # determine whether include the last to first loop
        if (type == 'circular'):
            last_idx = 1
        elif (type == 'linear'):
            last_idx = 0

        if type == 'linear':
            # get all longest
            cc = [len(row) for row in circ]
            idx = misc.all_max_index(cc)
            longest = [circ[i] for i in idx]

            # get the first one with max edge weight
            oop = []
            for row in longest:
                op = []
                for v in zip(row, row[1:] + row[:last_idx]):
                    if (self.edges[v][0] == 'outer'):
                        op.append(self.edges[v][1])
                oop.append(np.sum(op))
                circ_ref = longest[np.argmax(oop)]
        else:
            # get all longest
            cc = [len(row) for row in circ]
            circ_ref = []
            for i in set(cc):
                paths = [circ[ind] for ind, v in enumerate(cc) if v == i]
                # get the first one with max edge weight
                oop = []
                cand_ind = []
                for ind, row in enumerate(paths):
                    if check_chrom:
                        # find the all path is from same chrom 
                        chr = set([self.vertex[v][0] for v in row])
                        if len(chr) == 1:
                            cand_ind.append(ind)
                    op = []
                    for v in zip(row, row[1:] + row[:last_idx]):
                        if (self.edges[v][0] == 'outer'):
                            op.append(self.edges[v][1]) 
                    oop.append(np.sum(op))
                cand_ind.append(np.argmax(oop))
                circ_ref += [paths[t] for t in set(cand_ind)]
            # print(circ_ref)
        
        return circ_ref

    #--------------------------------------------------------------------------#
    def get_representitive_path(self, circ, type = 'circular', first_filter = False, second_filter = False, third_filter = False):
        """
        get representitive non-overlapping paths
        """
        paths = circ.copy()
        circ_repr = []
        if type == 'circular':
            t = 0
            while(paths):
                if t == 0:
                    check_chrom = True
                else:
                    check_chrom = False
                # search for longest path
                longest = self.path_longest(paths, type = type, check_chrom = check_chrom)
                circ_repr += longest
                # search non-intersected path
                # tf = [not bool(set(longest).intersection(row)) for row in paths]
                tf = []
                for row in paths:
                    tf_ = True
                    for c in longest:
                        if len(set(c).intersection(row))/len(row) > 0.7:
                            tf_ = False
                            break
                    tf.append(tf_)
                print(len(tf), sum(tf))
                # tf = [len(set(longest).intersection(row))/len(row) >= 0.4 for row in paths]
                paths = [i for (i, v) in zip(paths, tf) if v]
                t += 1
                if sum(tf) == len(paths):
                    break

            # # filter circ_repr
            # last_idx = 1
            # circ_cls = []
            # for row in circ_repr:
            #     op = []
            #     for v in zip(row, row[1:] + row[:last_idx]):
            #         if (self.edges[v][0] == 'outer'):
            #             op.append(self.edges[v][1])
            #     circ_cls.append([row[0], row[1], row[-2], row[-1], row, np.mean(op)])
                # circ_cls[row[:2] + row[-2:]].append(np.mean(op)) 
            
            # circ_cls_df = pd.DataFrame(circ_cls, columns = ['P1','P2','P3','P4', 'Path', 'Weight'])
            # if first_filter:
            #     circ_cls_df = circ_cls_df.groupby(['P1','P2','P3','P4']).apply(lambda x: x.loc[x['Weight'].idxmax()]).reset_index(drop=True)
            #     print(f'first_filter after: circ_cls_df.shape ({circ_cls_df.shape})')
            # if second_filter:
            #     circ_cls_df = circ_cls_df.groupby(['P1','P2','P4']).apply(lambda x: x.loc[x['Weight'].idxmax()]).reset_index(drop=True)
            #     print(f'second_filter after: circ_cls_df.shape ({circ_cls_df.shape})')
            # if third_filter:
            #     circ_cls_df = circ_cls_df.groupby(['P1','P3','P4']).apply(lambda x: x.loc[x['Weight'].idxmax()]).reset_index(drop=True)
            #     print(f'third_filter after: circ_cls_df.shape ({circ_cls_df.shape})')
            # circ_repr = circ_cls_df['Path'].tolist()
            
            # circ_repr = circ_cls_df.groupby(['P1','P3','P4']).apply(lambda x: x.loc[x['Weight'].idxmax()])['Path'].tolist()
        else:
            while(paths):
                # search for longest path
                longest = self.path_longest(paths, type = type)
                circ_repr.append(longest)
                # search non-intersected path
                tf = [not bool(set(longest).intersection(row)) for row in paths]
                paths = [i for (i, v) in zip(paths, tf) if v]
    
        return circ_repr

    #--------------------------------------------------------------------------#
    def make_amplicon_df(self, circ, type = 'circular'):
        """
        convert ggraph to interpretable amplicon dataframe
        """
        # determine whether include the last to first loop
        if (type == 'circular'):
            last_idx = 1
            tag = 'circ_'
        elif (type == 'linear'):
            last_idx = 0
            tag = 'line_'

        op = []
        for idx in range(len(circ)):
            row = circ[idx]
            op_seg = []
            op_bpp = []
            for v in zip(row, row[1:] + row[:last_idx]):
                if (self.edges[v][0] == 'inner'):
                    op_seg.append(self.vertex[v[0]][0:3] \
                        + self.vertex[v[1]][0:3] + self.edges[v])
                else:
                    op_bpp.append(self.edges[v])

            # put zero count for the last one of the linear amplicon
            if (type == 'linear'):
                op_bpp.append((('outer', 0)))

            for i in range(len(op_seg)):
                row_seg = op_seg[i]
                row_bpp = op_bpp[i]

                if (row_seg[2] == 'L'):
                    op.append([row_seg[0], row_seg[1], row_seg[4], \
                        '+', row_bpp[1], row_seg[7], tag + str(idx)])
                else:
                    op.append([row_seg[0], row_seg[4], row_seg[1], \
                        '-', row_bpp[1], row_seg[7], tag + str(idx)])

        colnames = ['Chrom', 'Start', 'End', 'Strand', \
            'SplitCount', 'CN', 'AmpliconID']
        df = pd.DataFrame(op, columns = colnames)

        # also add single segment loops
        if (type == 'circular'):
            node_in_path = [v for path in circ for v in path]
            ssl = [row for row in self.single_segment_loop \
                if (self.node_dict[(row[0], row[1], 'L')] not in node_in_path)]

            # add index
            ssl_df = []
            idx = len(circ)
            for row in ssl:
                ssl_df.append(row + ['circ_' + str(idx)])
                idx += 1

            # make dataframe
            ssl_df = pd.DataFrame(ssl_df, columns = colnames)
            df = pd.concat([df, ssl_df])

        return df
    
    @staticmethod
    def add_stats_circ(circ_anno, bin_norm, binsize):
        '''
        compute the stats for each circ
        '''
        # the length of each segs  
        circ_anno['Length'] = circ_anno['End'] - circ_anno['Start']

        # group by AmpliconID
        grouped = circ_anno.groupby('AmpliconID')

        # add stats
        result = grouped.agg({
            'Chrom': 'count', 
            'Length': 'sum',
            'SplitCount': ['sum','mean','var'],
            'CN': ['sum','mean','var']
        }).reset_index()

        if 'Gene' in circ_anno.columns:
            result['Gene_num'] = grouped['Gene'].apply(lambda x: x.str.split(';').str.len().sum())[0]
        else:
            result['Gene_count'] = 0

        if 'CancerGene' in circ_anno.columns:  
            result['CancerGene'] = grouped['CancerGene'].apply(lambda x: x.str.split(';').str.len().sum())[0]
        else:
            result['cancergenes_count'] = 0
            
        if 'SE' in circ_anno.columns:
            result['SE_count'] = grouped['SE'].apply(lambda x: x.str.split(';').str.len().sum())[0]
        else:
            result['SE_count'] = 0
            
        result.columns = ['AmpliconID', 'Seg_num', 'Length', 'SplitCount_sum', 'SplitCount_mean', 'SplitCount_Var', 
                        'CN_sum', 'CN_mean', 'CN_var', 'Gene_num', 'Cancergene_num', 'SE_num']

        # 使用transform获取首尾行
        start = grouped.min()[['Chrom', 'Start']]
        end = grouped.max()['End']
        result[['Chrom','Start', 'End']] = pd.concat([start, end], axis = 1).values
        CN_ext = {}
        # add the mean CN of +/- 100000 of start/end
        for chrom, start, end in zip(start['Chrom'], start['Start'], end):
            if (chrom, start, end) not in CN_ext:
                left = bin_norm[(bin_norm['Chrom'] == chrom) & 
                                (bin_norm['Coord'] < start - binsize) &  
                                (bin_norm['Coord'] > (start - binsize * 10))]['CN'].mean()
                right = bin_norm[(bin_norm['Chrom'] == chrom) & 
                                (bin_norm['Coord'] > end + binsize) &  
                                (bin_norm['Coord'] < (end + binsize * 10))]['CN'].mean()
                CN_ext[(chrom, start, end)] = [left, right]
            else:
                left, right = CN_ext[(chrom, start, end)]
            if 'Circ_ext_CN' not in CN_ext:
                CN_ext['Circ_ext_CN'] = [(left, right)]
            else:
                CN_ext['Circ_ext_CN'].append((left, right))
        result[['Left_CN','Right_CN']]   = pd.DataFrame(CN_ext['Circ_ext_CN'])
        rearrange = ['AmpliconID', 'Chrom','Start', 'End', 'Seg_num', 'Length', 'SplitCount_sum', 'SplitCount_mean', 'SplitCount_Var', 'CN_sum', 'CN_mean', 'CN_var', 'Left_CN', 'Right_CN', 'Gene_num', 'Cancergene_num', 'SE_num']
        return result[rearrange]

    #--------------------------------------------------------------------------#
    def build_ggraph_from_bp(self, bp_pair, bp_fine, seg):
        """
        build ggraph from breakpoint data
        """

        for row in bp_fine.itertuples():
            self.add_vertex(row[0], row[1:])

        for row in bp_fine.itertuples():
            self.node_dict[row[1:4]] = row[0]

        for row in bp_pair.itertuples():
            if row[1:4] in self.node_dict and row[4:7] in self.node_dict:
                a = self.node_dict[row[1:4]]
                b = self.node_dict[row[4:7]]
                w = ('outer', row[7])
                self.add_edge(a, b, w)

        for row in seg.itertuples():
            if (row[1], row[2], 'L') in self.node_dict and \
                (row[1], row[3], 'R') in self.node_dict:
                a = self.node_dict[(row[1], row[2], 'L')]
                b = self.node_dict[(row[1], row[3], 'R')]
                w = ('inner', row[4])
                # if already exist, it's a single segment loop
                if (a, b) in self.edges:
                    self.single_segment_loop.append([row[1], row[2], row[3], \
                        '+', self.edges[(a, b)][1], row[4]])
                    # don't add inner in this case
                else:
                    self.add_edge(a, b, w)


    #--------------------------------------------------------------------------#
    @staticmethod
    def invert_path(path):
        path = path[::-1]
        return path[-1:] + path[:-1]

    @staticmethod
    def rotate_path(path):
        i = path.index(min(path))
        return path[i:]+path[:i]

    @staticmethod
    def is_new_path(path, paths):
        return not path in paths


    #--------------------------------------------------------------------------#
    # currently not in use
    #--------------------------------------------------------------------------#
    def dfs_connected(self, v, visited, visited_list):
        """
        depth-first search for connected nodes (not in use currently)
        """

        # mark current node as visited
        visited[v] = True

        # store current node to list
        visited_list.append(v)

        # recur for adjacent nodes of the current node
        for v_adj in self.graph[v]:
            if (visited[v_adj] == False): # add to list if not visited
                visited_list = self.dfs_connected(v_adj, visited, visited_list)

        return visited_list

    #--------------------------------------------------------------------------#
    def get_connected_subgraph(self):
        """
        get connected subgraph nodes (not in use currently)
        """

        # init objects
        node = sorted(self.get_nodes())
        visited = []
        op = []

        # mark all nodes as not vsisited
        visited = [False]*(max(node) + 1)

        # loop to search for all subgraphs
        for v in node:
            if (visited[v] == False):
                visited_list = []
                op.append(self.dfs_connected(v, visited, visited_list))
        return op


def add_stats_circ(circ_anno, bin_norm, binsize):
    '''
    Compute the features for each circ and score for each circ based on these features. 
    '''
    # the length of each segs 
    circ_anno['Length'] = circ_anno['End'] - circ_anno['Start']

    # group by AmpliconID
    grouped = circ_anno.groupby('AmpliconID')

    # add stats for each circ
    circ_stat = grouped.agg({
        'Chrom': 'count', 
        'Length': 'sum',
        'SplitCount': ['sum','mean','std'],
        'CN': ['sum','mean', 'std']
    }).reset_index()

    if 'Gene' in circ_anno.columns:
        circ_stat['Gene_num'] = grouped['Gene'].apply(lambda x: x.str.split(';').str.len().sum()).reset_index(drop=True)
    else:
        circ_stat['Gene_count'] = 0

    if 'CancerGene' in circ_anno.columns:
        try:
            circ_stat['CancerGene_count'] = grouped['CancerGene'].apply(lambda x: x.str.split(';').str.len().sum()).reset_index(drop=True)
        except:
            print(f'CancerGene_count is error; set is 0 {circ_anno["CancerGene"]}.')
            circ_stat['CancerGene_count'] = 0
    else:
        circ_stat['CancerGene_count'] = 0
        
    if 'SE' in circ_anno.columns:
        circ_stat['SE_count'] = grouped['SE'].apply(lambda x: x.str.split(';').str.len().sum()).reset_index(drop=True)
    else:
        circ_stat['SE_count'] = 0
        
    circ_stat.columns = ['AmpliconID', 'Seg_num', 'Length', 'SplitCount_sum', 'SplitCount_mean', 'SplitCount_std', 
                    'CN_sum', 'CN_mean', 'CN_std', 'Gene_num', 'Cancergene_num', 'SE_num']

    start = grouped.head(1)[['AmpliconID', 'Chrom', 'Start']].reset_index(drop=True)
    end = grouped.tail(1)[['Chrom', 'End']].reset_index(drop=True)
    circ_interval = pd.concat([start, end], axis = 1)
    CN_ext = {}
    # add the mean CN of +/- 100000 of start/end
    for row in circ_interval.itertuples():
        if row[2:] not in CN_ext:
            left = bin_norm[(bin_norm['Chrom'] == row[2]) & 
                            (bin_norm['Coord'] < (row[3] - 2 * binsize)) &  
                            (bin_norm['Coord'] > (row[3] - binsize * 10))]['CN'].mean()
            right = bin_norm[(bin_norm['Chrom'] == row[4]) & 
                            (bin_norm['Coord'] > (row[5] + 2 * binsize)) &  
                            (bin_norm['Coord'] < (row[5] + binsize * 10))]['CN'].mean()
            CN_ext[row[2:]] = [left, right]
        else:
            left, right = CN_ext[row[2:]]
            
        if 'Circ_ext_CN' not in CN_ext:
            CN_ext['Circ_ext_CN'] = [(left, right)]
        else:
            CN_ext['Circ_ext_CN'].append((left, right))

    circ_interval[['Left_CN','Right_CN']]   = pd.DataFrame(CN_ext['Circ_ext_CN'])
    circ_stat = pd.merge(circ_stat, circ_interval.iloc[:,[0,1,2,4,5,6]], on = 'AmpliconID')
    rearrange = ['AmpliconID', 'Chrom','Start', 'End', 'Seg_num', 'Length', 'SplitCount_sum', 'SplitCount_mean', 'SplitCount_std', 'CN_sum', 'CN_mean', 'CN_std', 'Left_CN', 'Right_CN', 'Gene_num', 'Cancergene_num', 'SE_num']
    circ_stat = circ_stat[rearrange]

    circ_stat.iloc[:,2:] = np.log(circ_stat.iloc[:,2:].astype(float) + 1)

    ## integrate feature
    circ_stat['F1'] = circ_stat['SplitCount_mean']/(circ_stat['SplitCount_std'] + 0.5) + circ_stat['SplitCount_mean']/(circ_stat['SplitCount_mean']).sum()
    circ_stat['F2'] = circ_stat['CN_mean']/(circ_stat['CN_std'] + 0.5)  + circ_stat['CN_mean']/(circ_stat['CN_mean']).sum()
    circ_stat['F3'] = circ_stat['CN_mean']/(circ_stat['Right_CN'])
    circ_stat['F4'] = circ_stat['CN_mean']/(circ_stat['Left_CN'])
    circ_stat = circ_stat.fillna(0)   
    coef_ = np.array([0.01, 1.5,  1.5,  0.5, 0.5])
    # Score the circ
    score = circ_stat.loc[circ_stat['Gene_num'] > 0, ['Seg_num', 'F1', 'F2', 'F3','F4']].values @ coef_ 
    circ_stat['Score'] = 0
    circ_stat['Score'][circ_stat['Gene_num'] > 0]  = score
   
    return circ_stat