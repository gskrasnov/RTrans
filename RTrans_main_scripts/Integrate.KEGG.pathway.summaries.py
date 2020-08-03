import sys,os,optparse
import xlsxwriter
from xlsxwriter.utility import xl_range

use_glob = True
try:
    import glob
except ImportError:
    use_glob = False


import math

def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    if not N:
        return None

    N = sorted(N)
    k = (len(N)-1) * (float(percent)/100)
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1

def color(RGB):
    R,G,B = RGB
    Rhex=hex(R)[2:]
    if len(Rhex)==1:
        Rhex='0'+Rhex
    Ghex=hex(G)[2:]
    if len(Ghex)==1:
        Ghex='0'+Ghex
    Bhex=hex(B)[2:]
    if len(Bhex)==1:
        Bhex='0'+Bhex
    return '#%s%s%s'%(Rhex,Ghex,Bhex)

def Gradient(Points,Position):
    Position = min(1,max(0,Position))
    if Position == 1:
        return Points[-1]
    FragmentSize = (1/(len(Points)-1))
    FragmentNumber = int(Position / FragmentSize)

    InFragmentPosition = (Position%FragmentSize)/FragmentSize

    RColor=int(Points[FragmentNumber][0] + (Points[FragmentNumber+1][0] - Points[FragmentNumber][0])*InFragmentPosition)
    GColor=int(Points[FragmentNumber][1] + (Points[FragmentNumber+1][1] - Points[FragmentNumber][1])*InFragmentPosition)
    BColor=int(Points[FragmentNumber][2] + (Points[FragmentNumber+1][2] - Points[FragmentNumber][2])*InFragmentPosition)

    return RColor,GColor,BColor

OverexpressionMaxColor = (226,101,0)
DownregulationMaxColor = (29,136,234)

OverexpressionMaxColor = (237,128,38)
DownregulationMaxColor = (29,136,234)

Bar_and_Background_colors_difference = 0.1

OverexpressionMaxColor_bg = Gradient([OverexpressionMaxColor,(255,255,255)],Bar_and_Background_colors_difference)
DownregulationMaxColor_bg = Gradient([DownregulationMaxColor,(255,255,255)],Bar_and_Background_colors_difference)

parser = optparse.OptionParser()
parser.add_option('--input-file-mask', dest='Mask',default=None, help='')
parser.add_option('--input-file-names', dest='FileNames',default=None, help='')
parser.add_option('--out-tsv-prefix', dest='Out_tsv_prefix',default=None, help='')
parser.add_option('--out-excel', dest='Out_Excel_FileName',default=None, help='')
parser.add_option('--heatmax-max-percentile', dest='Heatmap_Max_Percentile',default='95', help='')
parser.add_option('--bar-max-percentile', dest='Bar_Max_Percentile',default='100', help='')
parser.add_option('--bar-max-value-override', dest='Bar_Max_Value_override',default=None, help='')
parser.add_option('--combine-node-count-and-perc', dest='Combine_Node_count_and_perc',default='no', help='')
opts, args  = parser.parse_args(sys.argv[1:])


if opts.Mask != None and not use_glob:
    print('Cannot use mask-based file search since "glob" python module cannot be imported')
    exit(127)
    
if (opts.FileNames == None) == (opts.Mask == None):
    print('either File names or file mask must be specified')
    exit()

if opts.FileNames == None:
    FileNames = glob.glob(opts.Mask)
else:
    FileNames = opts.FileNames.replace('"','').replace("'",'').split(';')
    if 'win' in sys.platform: opts.FileNames = opts.FileNames.replace('/','\\')

Models = [os.path.split(x)[1].split(', KEGG pathways DE')[0] for x in FileNames]

opts.Heatmap_Max_Percentile = float(opts.Heatmap_Max_Percentile)
opts.Bar_Max_Percentile = float(opts.Bar_Max_Percentile)
opts.Bar_Max_Value_override = float(opts.Bar_Max_Value_override)
opts.Combine_Node_count_and_perc = opts.Combine_Node_count_and_perc.casefold() in ['yes','on','y']

class TPathwayInfo:
    def __init__(self,id,name,threshold,total_nodes):
        self.id = id
        self.name = name
        self.nodes_count_by_model = dict()
        self.nodes_percentage_by_model = dict()
        self.threshold = threshold
        self.total_nodes = total_nodes

PathwayInfos_by_threshold_by_pathway_id = dict()
Pathway_names_by_pathway_id = dict()
Total_nodes_by_pathway_id = {}

all_thresholds = []

for n in range(len(FileNames)):
    FN = FileNames[n]
    M = Models[n]
    F = open(FN)
    header = F.readline().rstrip().split('\t')
    col_indici = dict([(header[x],x+1) for x in range(2,len(header))])

    gene_count_col_index_by_threshold = {}
    gene_percent_col_index_by_threshold = {}

    for x in list(col_indici.keys()):
        thr = float(x.replace(', %',"").replace('nodes with LogFC > ','').replace('nodes with LogFC < ',''))
        if ', %' in x:
            gene_percent_col_index_by_threshold[thr] = col_indici[x]
        else:  gene_count_col_index_by_threshold[thr] = col_indici[x]

        if not thr in PathwayInfos_by_threshold_by_pathway_id: PathwayInfos_by_threshold_by_pathway_id[thr] = dict()

    current_thresholds = sorted(PathwayInfos_by_threshold_by_pathway_id)

    for L in F.readlines():
        C = L.rstrip().split('\t')
        pathway_id = C[0]
        pathway_name = C[1]
        total_nodes = int(C[2])
        Total_nodes_by_pathway_id[pathway_id] = total_nodes

        Pathway_names_by_pathway_id[pathway_id] = pathway_name
        for thr in current_thresholds:
            if not pathway_id in PathwayInfos_by_threshold_by_pathway_id[thr]:
                PathwayInfos_by_threshold_by_pathway_id[thr][pathway_id] = TPathwayInfo(pathway_id,pathway_name,thr,total_nodes)
            PathwayInfos_by_threshold_by_pathway_id[thr][pathway_id].nodes_count_by_model[M] = C[gene_count_col_index_by_threshold[thr]]
            PathwayInfos_by_threshold_by_pathway_id[thr][pathway_id].nodes_percentage_by_model[M] = C[gene_percent_col_index_by_threshold[thr]]

    all_thresholds += current_thresholds



all_thresholds = list(sorted(set(all_thresholds)))

if opts.Out_tsv_prefix != None:
    for thr in all_thresholds:
        if thr < 0: ins = 'lt_-'
        else: ins = 'gt_+'
        with open(opts.Out_tsv_prefix + "__LogFC_%s%g__node_count.tsv"%(ins,thr),'w') as out:
            out.write('Pathway ID\tPathway name\tTotal nodes\t' + '\t'.join(Models) + '\n')
            for id in sorted(Pathway_names_by_pathway_id):
                out.write("%s\t%s\t%d\t"%(id,Pathway_names_by_pathway_id[id],Total_nodes_by_pathway_id[id]))
                for M in Models:
                    if not id in PathwayInfos_by_threshold_by_pathway_id[thr]: text = 'na'
                    elif not M in PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_count_by_model: text = 'na'
                    else:  text = PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_count_by_model[M]
                    out.write(text + '\t')
                out.write('\n')

        with open(opts.Out_tsv_prefix + "__LogFC_%s%g__node_percentage.tsv"%(ins,thr),'w') as out:
            out.write('Pathway ID\tPathway name\tTotal nodes\t' + '\t'.join(Models) + '\n')
            for id in sorted(Pathway_names_by_pathway_id):
                out.write("%s\t%s\t%d\t"%(id,Pathway_names_by_pathway_id[id],Total_nodes_by_pathway_id[id]))
                for M in Models:
                    if not id in PathwayInfos_by_threshold_by_pathway_id[thr]: text = 'na'
                    elif not M in PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_percentage_by_model: text = 'na'
                    else:  text = PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_percentage_by_model[M]
                    out.write(text + '\t')
                out.write('\n')

if opts.Out_Excel_FileName != None:
    Workbook = xlsxwriter.Workbook(opts.Out_Excel_FileName)

    bold = Workbook.add_format({'bold': True, 'italic': False})
    italic = Workbook.add_format({'bold': False, 'italic': True})
    bold_italic = Workbook.add_format({'bold': True, 'italic': True})

    # zero_format = Workbook.add_format()
    zero_format = Workbook.add_format({'font_color': '#888888','bg_color':'#ffffff'})

    simple_wrap = Workbook.add_format({'bold': False, 'italic': False, 'align': 'center', 'valign': 'vcenter'})
    simple_wrap.set_text_wrap()
    bold_wrap = Workbook.add_format({'bold': True, 'italic': False, 'align': 'center', 'valign': 'vcenter'})
    bold_wrap.set_text_wrap()
    italic_wrap = Workbook.add_format({'bold': False, 'italic': True, 'align': 'center', 'valign': 'vcenter'})
    italic_wrap.set_text_wrap()
    bold_italic_wrap = Workbook.add_format({'bold': True, 'italic': True, 'align': 'center', 'valign': 'vcenter'})
    bold_italic_wrap.set_text_wrap()

    format0 = Workbook.add_format()
    format0.set_num_format('0')
    format1 = Workbook.add_format()
    format1.set_num_format('0.0')
    format2 = Workbook.add_format()
    format2.set_num_format('0.00')

    border_format = Workbook.add_format({'bold': True, 'italic': True, 'border': True, 'align': 'center', 'valign': 'vcenter'})
    border_format.set_text_wrap()

    all_thresholds_abs_values = list(sorted(set([abs(thr) for thr in all_thresholds])))

    for thr in all_thresholds_abs_values:
        sheet = Workbook.add_worksheet('|LogFC| > %g'%thr)
        sheet.write(1, 0, 'KEGG ID', italic_wrap)
        sheet.set_column(0, 0, 9)
        sheet.write(1, 1, 'pathway name', italic_wrap)
        sheet.set_column(1, 1, 42)
        sheet.write(1, 2, 'nodes count in pathway', italic_wrap)
        sheet.set_column(3, 3, 5)
        start_ColN = 4

        positive_LogFC_node_count_column_coords = []  # format: (start_col, end_col)
        negative_LogFC_node_count_column_coords = []
        positive_LogFC_node_perc_column_coords = []
        negative_LogFC_node_perc_column_coords = []

        if thr in all_thresholds:
            if len(Models) > 1:
                sheet.merge_range(0, start_ColN, 0, start_ColN + len(Models) - 1, 'nodes with LogFC > %g'%thr, bold_wrap)
            else:
                sheet.write(0, start_ColN, 'nodes with LogFC > %g'%thr, bold_wrap)
            positive_LogFC_node_count_column_coords.append((start_ColN, start_ColN + len(Models) - 1))
            for x in range(len(Models)):  sheet.write(1, start_ColN + x, Models[x],simple_wrap)
            start_ColN += len(Models)
            sheet.set_column(start_ColN, start_ColN, 5)
            start_ColN += 1

        if (-1)*thr in all_thresholds:
            if len(Models) > 1:
                sheet.merge_range(0, start_ColN, 0, start_ColN + len(Models) - 1, 'nodes with LogFC < %g'%thr, bold_wrap)
            else:
                sheet.write(0, start_ColN, 'nodes with LogFC < %g'%thr, bold_wrap)
            negative_LogFC_node_count_column_coords.append((start_ColN, start_ColN + len(Models) - 1))
            for x in range(len(Models)):  sheet.write(1, start_ColN + x, Models[x],simple_wrap)
            start_ColN += len(Models)
            sheet.set_column(start_ColN, start_ColN, 5)
            start_ColN += 1

        if thr in all_thresholds and not opts.Combine_Node_count_and_perc:
            if len(Models) > 1:
                sheet.merge_range(0, start_ColN, 0, start_ColN + len(Models) - 1, 'perc. nodes with LogFC > %g' % thr, bold_wrap)
            else:
                sheet.write(0, start_ColN, 'perc. nodes with LogFC > %g' % thr, bold_wrap)
            positive_LogFC_node_perc_column_coords.append((start_ColN, start_ColN + len(Models) - 1))
            for x in range(len(Models)):  sheet.write(1, start_ColN + x, Models[x],simple_wrap)
            start_ColN += len(Models)
            sheet.set_column(start_ColN, start_ColN, 5)
            start_ColN += 1

        if (-1)*thr in all_thresholds and not opts.Combine_Node_count_and_perc:
            if len(Models) > 1:
                sheet.merge_range(0, start_ColN, 0, start_ColN + len(Models) - 1, 'perc. nodes with LogFC < %g' % thr, bold_wrap)
            else:
                sheet.write(0, start_ColN, 'perc. nodes with LogFC < %g' % thr, bold_wrap)
            negative_LogFC_node_perc_column_coords.append((start_ColN, start_ColN + len(Models) - 1))
            for x in range(len(Models)):  sheet.write(1, start_ColN + x, Models[x],simple_wrap)
            start_ColN += len(Models)
            sheet.set_column(start_ColN, start_ColN, 5)
            start_ColN += 1

        all_node_counts_values = []
        all_node_perc_values = []
        for pathway_n in range(len(Pathway_names_by_pathway_id)):
            id = sorted(Pathway_names_by_pathway_id)[pathway_n]
            sheet.write(2 + pathway_n, 0,id)
            sheet.write(2 + pathway_n, 1,Pathway_names_by_pathway_id[id])
            sheet.write_number(2 + pathway_n, 2,Total_nodes_by_pathway_id[id])

            start_ColN = 4
            if thr in all_thresholds:
                for x in range(len(Models)):
                    M = Models[x]
                    if not id in PathwayInfos_by_threshold_by_pathway_id[thr]:  sheet.write(2 + pathway_n, start_ColN + x, 'na')
                    elif not M in PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_count_by_model:  sheet.write(2 + pathway_n, start_ColN + x, 'na')
                    else:
                        try:
                            node_count = int(PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_count_by_model[M])
                            if node_count == 0:  sheet.write_number(2 + pathway_n, start_ColN + x, node_count,zero_format)
                            else:  sheet.write_number(2 + pathway_n, start_ColN + x, node_count)
                            all_node_counts_values.append(node_count)
                        except ValueError:
                            sheet.write(2 + pathway_n, start_ColN + x, PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_count_by_model[M])

                start_ColN += len(Models) + 1

            if (-1)*thr in all_thresholds:
                for x in range(len(Models)):
                    M = Models[x]
                    if not id in PathwayInfos_by_threshold_by_pathway_id[-thr]:  sheet.write(2 + pathway_n, start_ColN + x, 'na')
                    elif not M in PathwayInfos_by_threshold_by_pathway_id[-thr][id].nodes_count_by_model:  sheet.write(2 + pathway_n, start_ColN + x, 'na')
                    else:
                        try:
                            node_count = int(PathwayInfos_by_threshold_by_pathway_id[-thr][id].nodes_count_by_model[M])
                            if node_count == 0:  sheet.write_number(2 + pathway_n, start_ColN + x, node_count,zero_format)
                            else:  sheet.write_number(2 + pathway_n, start_ColN + x, node_count)
                            all_node_counts_values.append(node_count)
                        except ValueError:
                            sheet.write(2 + pathway_n, start_ColN + x, PathwayInfos_by_threshold_by_pathway_id[-thr][id].nodes_count_by_model[M])

                start_ColN += len(Models) + 1

            if thr in all_thresholds:
                for x in range(len(Models)):
                    M = Models[x]
                    if not id in PathwayInfos_by_threshold_by_pathway_id[thr]:
                        if not opts.Combine_Node_count_and_perc:  sheet.write(2 + pathway_n, start_ColN + x, 'na')
                    elif not M in PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_percentage_by_model:
                        if not opts.Combine_Node_count_and_perc:  sheet.write(2 + pathway_n, start_ColN + x, 'na')
                    else:
                        try:
                            node_perc = float(PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_percentage_by_model[M])
                            if not opts.Combine_Node_count_and_perc:
                                if node_perc == 0:  sheet.write_number(2 + pathway_n, start_ColN + x, node_perc,zero_format)
                                else:  sheet.write_number(2 + pathway_n, start_ColN + x, node_perc,format1)
                            all_node_perc_values.append(node_perc)
                        except ValueError:
                            if not opts.Combine_Node_count_and_perc:  sheet.write(2 + pathway_n, start_ColN + x, PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_percentage_by_model[M])

                start_ColN += len(Models) + 1

            if (-1)*thr in all_thresholds:
                for x in range(len(Models)):
                    M = Models[x]
                    if not id in PathwayInfos_by_threshold_by_pathway_id[-thr]:
                        if not opts.Combine_Node_count_and_perc:  sheet.write(2 + pathway_n, start_ColN + x, 'na')
                    elif not M in PathwayInfos_by_threshold_by_pathway_id[-thr][id].nodes_percentage_by_model:
                        if not opts.Combine_Node_count_and_perc:  sheet.write(2 + pathway_n, start_ColN + x, 'na')
                    else:
                        try:
                            node_perc = float(PathwayInfos_by_threshold_by_pathway_id[-thr][id].nodes_percentage_by_model[M])
                            if not opts.Combine_Node_count_and_perc:
                                if node_perc == 0:  sheet.write_number(2 + pathway_n, start_ColN + x, node_perc,zero_format)
                                else:  sheet.write_number(2 + pathway_n, start_ColN + x, node_perc,format1)
                            all_node_perc_values.append(node_perc)
                        except ValueError:
                            if not opts.Combine_Node_count_and_perc:  sheet.write(2 + pathway_n, start_ColN + x, PathwayInfos_by_threshold_by_pathway_id[-thr][id].nodes_percentage_by_model[M])

                start_ColN += len(Models) + 1

        pathway_count = len(Pathway_names_by_pathway_id)
        sheet.conditional_format(2, 2, 1 + pathway_count,2, {'type': 'data_bar', 'bar_color': '#ffb910',
                                                                         'min_type': 'num', 'max_type': 'percentile',
                                                                         'min_value': 0, 'max_value': 95})

        node_count__graph_max_value = max(1, percentile(all_node_counts_values,opts.Heatmap_Max_Percentile))
        if opts.Bar_Max_Value_override != None and opts.Bar_Max_Value_override != 0:
            node_perc__graph_max_value = percentile(all_node_perc_values,opts.Bar_Max_Percentile)
        else:
            node_perc__graph_max_value = opts.Bar_Max_Value_override

        if not opts.Combine_Node_count_and_perc:
            for coords in positive_LogFC_node_count_column_coords:
                col_start,col_end = coords
                sheet.conditional_format(2,col_start,2 + pathway_count,col_end, {'type': '2_color_scale',
                                                          'min_color': "#ffffff",'max_color': "#de6f16",
                                                          'min_type': 'num','max_type': 'num',
                                                          'min_value': 0, 'max_value': node_count__graph_max_value+0.01})

            for coords in negative_LogFC_node_count_column_coords:
                col_start,col_end = coords
                sheet.conditional_format(2,col_start,2 + pathway_count,col_end, {'type': '2_color_scale',
                                                          'min_color': "#ffffff",'max_color': "#1985e8",
                                                          'min_type': 'num','max_type': 'num',
                                                          'min_value': 0, 'max_value': node_count__graph_max_value+0.01})

            for coords in positive_LogFC_node_perc_column_coords:
                col_start,col_end = coords
                sheet.conditional_format(2,col_start,2 + pathway_count,col_end, {'type': 'data_bar', 'bar_color': '#f48226',
                                                                     'min_type': 'num', 'max_type': 'num',
                                                                     'min_value': 0, 'max_value': 100})

            for coords in negative_LogFC_node_perc_column_coords:
                col_start,col_end = coords
                sheet.conditional_format(2,col_start,2 + pathway_count,col_end, {'type': 'data_bar', 'bar_color': '#3794e9',
                                                                     'min_type': 'num', 'max_type': 'num',
                                                                     'min_value': 0, 'max_value': 100})

        if opts.Combine_Node_count_and_perc:
            for pathway_n in range(len(Pathway_names_by_pathway_id)):
                id = sorted(Pathway_names_by_pathway_id)[pathway_n]

                start_ColN = 4
                if thr in all_thresholds:
                    for x in range(len(Models)):
                        M = Models[x]
                        if not id in PathwayInfos_by_threshold_by_pathway_id[thr]: continue
                        elif not M in PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_count_by_model: continue
                        else:
                            try:
                                node_count = int(PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_count_by_model[M])
                                node_perc = float(PathwayInfos_by_threshold_by_pathway_id[thr][id].nodes_percentage_by_model[M])

                                if node_count == 0 or node_perc == 0:
                                    sheet.write_number(2 + pathway_n, start_ColN + x, node_count,zero_format)
                                    continue

                                current_bg_color = color(Gradient([(255,255,255),OverexpressionMaxColor_bg],float(node_count)/float(node_count__graph_max_value)))
                                current_bar_color = color(Gradient([(255,255,255),OverexpressionMaxColor],float(node_count)/float(node_count__graph_max_value)))
                                current_format = Workbook.add_format({'bg_color': current_bg_color})

                                sheet.write_number(2 + pathway_n, start_ColN + x, node_count,current_format)
                                sheet.conditional_format(2 + pathway_n, start_ColN + x, 2 + pathway_n, start_ColN + x,
                                                     {'type': 'data_bar', 'bar_color': current_bar_color,
                                                      'min_type': 'num', 'max_type': 'num',
                                                      'min_value': 0, 'max_value': node_count/node_perc*100})

                            except ValueError:
                                continue

                    start_ColN += len(Models) + 1

                if -thr in all_thresholds:
                    for x in range(len(Models)):
                        M = Models[x]
                        if not id in PathwayInfos_by_threshold_by_pathway_id[-thr]: continue
                        elif not M in PathwayInfos_by_threshold_by_pathway_id[-thr][id].nodes_count_by_model: continue
                        else:
                            try:
                                node_count = int(PathwayInfos_by_threshold_by_pathway_id[-thr][id].nodes_count_by_model[M])
                                node_perc = float(PathwayInfos_by_threshold_by_pathway_id[-thr][id].nodes_percentage_by_model[M])

                                if node_count == 0 or node_perc == 0:
                                    sheet.write_number(2 + pathway_n, start_ColN + x, node_count,zero_format)
                                    continue

                                current_bg_color = color(Gradient([(255,255,255),DownregulationMaxColor_bg],float(node_count)/float(node_count__graph_max_value)))
                                current_bar_color = color(Gradient([(255,255,255),DownregulationMaxColor],float(node_count)/float(node_count__graph_max_value)))
                                current_format = Workbook.add_format({'bg_color': current_bg_color})

                                sheet.write_number(2 + pathway_n, start_ColN + x, node_count,current_format)
                                sheet.conditional_format(2 + pathway_n, start_ColN + x, 2 + pathway_n, start_ColN + x,
                                                     {'type': 'data_bar', 'bar_color': current_bar_color,
                                                      'min_type': 'num', 'max_type': 'num',
                                                      'min_value': 0, 'max_value': node_count/node_perc*100})

                            except ValueError:
                                continue

                    start_ColN += len(Models) + 1

        sheet.write(1 + pathway_count + 2,0,"* KEGG nodes represent subfamilies of genes/proteins (e.g. node 'Ldh' is a composition of LdhA and LdhB enzymes)",italic)
        sheet.write(1 + pathway_count + 3,0,"   or single gene/protein (e.g. KRAS, BRAF), and vise versa, one protein can participate metabolic pathway is different places.",italic)
        sheet.write(1 + pathway_count + 4,0,"   In such cases, one gene/protein corresponds to several KEGG nodes",italic)
        sheet.write(1 + pathway_count + 5,0,"   (for an example, see pathway hsa00140; genes HSD3B1 (16 nodes), UGT2B11 (8 nodes))",italic)
        sheet.write(1 + pathway_count + 7,0,"* cell background indicates the number of affected nodes in pathway",italic)
        sheet.write(1 + pathway_count + 8,0,"* bar indicates the percentage of affected nodes in pathway",italic)



    Workbook.close()

    # sheet.write(1,2,'Gene ID',bold_wrap)
