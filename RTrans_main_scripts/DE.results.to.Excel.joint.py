__author__ = 'George'
import sys,math,os,math
import xlsxwriter
import argparse
from xlsxwriter.utility import xl_range

###
# usage: python3 Convert.TRans.ouput.to.xlsx out.xlsx file1.txt file2.txt

class TGO_Term_LogFCs():
    def __init__(self,ID,Name,min_P,min_FDR,gene_max_P,gene_min_CPM,Genes_count):
        self.ID =ID
        self.Name = Name
        self.min_P = min_P
        self.min_FDR = min_FDR
        self.gene_max_P = gene_max_P
        self.gene_min_CPM = gene_min_CPM
        self.Genes_count = Genes_count
        self.LogFCs = []

class TGO_Term_Info():
    def __init__(self,ID,Name):
        self.ID =ID
        self.Name = Name
        self.LogFCs_by_combination_code = {}



def Convert_ColorPoints_to_DecTuple(Points):
    Points_adj = []
    for p in Points:
        if type(p) is str:
            p = p.lstrip('#')
            p = tuple(int(p[i:i + 2], 16) for i in (0, 2, 4))
        Points_adj.append(p)

    return(Points_adj)

def Gradient(Points,Position):
    Position = min(1,max(0,Position))
    Points = Convert_ColorPoints_to_DecTuple(Points)

    if Position == 1:
        return Points[-1]
    FragmentSize = (1/(len(Points)-1))
    FragmentNumber = int(Position / FragmentSize)

    InFragmentPosition = (Position%FragmentSize)/FragmentSize

    RColor=int(Points[FragmentNumber][0] + (Points[FragmentNumber+1][0] - Points[FragmentNumber][0])*InFragmentPosition)
    GColor=int(Points[FragmentNumber][1] + (Points[FragmentNumber+1][1] - Points[FragmentNumber][1])*InFragmentPosition)
    BColor=int(Points[FragmentNumber][2] + (Points[FragmentNumber+1][2] - Points[FragmentNumber][2])*InFragmentPosition)

    return RColor,GColor,BColor

def TwoD_Gradient(Points_sets, PositionX, PositionY):
    PositionX = min(1,max(0,PositionX))
    PositionY = min(1,max(0,PositionY))

    ## an example:     Points_sets = [[(128,128,128), (255,0,0)], [(128,128,128), (0,255,0)], [(128,128,128), (0,0,255)]]


    if PositionX < 1:
        FragmentSizeX = (1/(len(Points_sets)-1))
        FragmentNumberX = int(PositionX / FragmentSizeX)

        InFragmentPositionX = (PositionX % FragmentSizeX)/FragmentSizeX

        Points_high = Points_sets[FragmentNumberX + 1]
        Points_low = Points_sets[FragmentNumberX]

        Color_high = Gradient(Points_high, PositionY)
        Color_low = Gradient(Points_low, PositionY)

        RColor=int(Color_low[0] + (Color_high[0] - Color_low[0])*InFragmentPositionX)
        GColor=int(Color_low[1] + (Color_high[1] - Color_low[1])*InFragmentPositionX)
        BColor=int(Color_low[2] + (Color_high[2] - Color_low[2])*InFragmentPositionX)
    else:
        return Gradient(Points_sets[-1], PositionY)

    return RColor,GColor,BColor

TwoD_ColorPointsSet1 = [
    ['434ff6','ffffff', 'ddd730'],
    ['00c8ca','ffffff', 'd6a51e'],
    ['32d026','ffffff', 'd6532b']
    ]

TwoD_ColorPointsSet2 = [
    ['a0a6f7','ffffff', 'e8e598'],
    ['6beff0','ffffff', 'e0bd5f'],
    ['32d026','ffffff', 'd6532b']
    ]

TwoD_ColorPointsSet3 = [
    ['d3e1fc','ffffff', 'fae9d8'],
    ['2568f1','a0a0a0', 'df7208']
    ]

Color_gradient_schemas = [TwoD_ColorPointsSet1, TwoD_ColorPointsSet2, TwoD_ColorPointsSet3]


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


# x = TwoD_Gradient(Color_gradient_schemas[0], 0.5, 0.5)
# print(color(x))
# exit()

def write_number__mod(sheet, rowN, colN, text, format=None, bypass_if_NaN = True,
    pvalue_mode = False, pvalue_formats = None):
    if (pvalue_formats is None) == pvalue_mode:
        raise('write_number__mod: pvalue_formats and pvalue_mode conflict')
    if pvalue_mode and format != None:
        raise('write_number__mod: format should not be defined in pvalue_mode')

    try:
        d = float(text)
        if pvalue_mode:
            if not math.isnan(d) and not math.isinf(d):
                digits_count = int((-1)*math.log10(d)) + 1
                if digits_count < len(pvalue_formats):
                    format = pvalue_formats[digits_count]

        if format != None:
            if (math.isnan(d) or math.isinf(d)):
                if bypass_if_NaN: return
                sheet.write(rowN,colN,str(text),format)
            else: sheet.write_number(rowN,colN,d,format)
        else:
            if (math.isnan(d) or math.isinf(d)):
                if bypass_if_NaN:  return
                sheet.write(rowN,colN,str(text))
            else: sheet.write_number(rowN,colN,d)
    except ValueError:
        if bypass_if_NaN:
            return
        if format != None:
            sheet.write(rowN,colN,str(text),format)
        else:
            sheet.write(rowN,colN,str(text))


def FormatLogFC_Cell(sheet,FC,RowN,ColN,coeff1 = 1.0, coeff2 = 1.0):
    if math.isnan(FC): return
    OverexpressionMaxColor = (226,101,0)
    DownregulationMaxColor = (29,136,234)

    LogFC = math.log2(FC)
    if LogFC >= 0:
        C = color(Gradient(((255,255,255),OverexpressionMaxColor), min(1,(coeff1*LogFC + 0.1)/2.5)))
        sheet.conditional_format(RowN,ColN,RowN,ColN, {'type': 'data_bar','bar_color': C,
                                         'min_type':'num','max_type':'num',
                                         'min_value':0,'max_value':7/coeff2})
    else:
        C = color(Gradient(((255,255,255),DownregulationMaxColor), min(1,(coeff1*LogFC*(-1) + 0.1)/2.5)))
        sheet.conditional_format(RowN,ColN,RowN,ColN, {'type': 'data_bar','bar_color': C,
                                         'min_type':'num','max_type':'num',
                                         'min_value':LogFC*2,'max_value':LogFC*2+7/coeff2})


def FormatLogCPM_Cell(sheet,logCPM,RowN,ColN, max_value = 8):
    if math.isnan(logCPM): return
    maxColor = (220,181,0)

    if logCPM <= 0: return
    C = color(Gradient(((255,255,255),maxColor), min(1,(logCPM + 0.1)/max_value)))
    sheet.conditional_format(RowN,ColN,RowN,ColN, {'type': 'data_bar','bar_color': C,
                                     'min_type':'num','max_type':'num',
                                     'min_value':0,'max_value':max_value})

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

class TDE_Info():
    def __init__(self):
        self.LogFC = float('NaN')
        self.LogCPM = float('NaN')
        self.LR = float('NaN')
        self.F = float('NaN')
        self.GLM_P = float('NaN')
        self.GLM_FDR = float('NaN')
        self.Score = float('NaN')
        self.raw_mean = None
        self.ReadCounts_Passed = None
        self.ttest_p = float('NaN')
        self.ttest_FDR = float('NaN')
        self.Wilcoxon_p = float('NaN')
        self.Wilcoxon_FDR = float('NaN')
        self.Mann_Wh_p = float('NaN')
        self.Mann_Wh_FDR = float('NaN')
        self.Spearman_r = float('NaN')
        self.Spearman_p = float('NaN')
        self.Pearson_r = float('NaN')
        self.Pearson_p = float('NaN')

class TGene_Info():
    def __init__(self,GeneID,mode='edger'):
        self.GeneID = GeneID
        self.Symbol = ''
        self.mode = mode
        self.Biotype = ''
        self.RefSeq_summary = ''
        self.Name = ''
        self.CPM_by_SampleName = dict()
        self.DE_Info_by_GLM = dict()
        self.Places_in_rating = []

    def Calculate_Summary_Rating(self, method = 'rating', ttest_missing_Score_multiplier = 0.7):
        if method == 'rating':
            if len(self.Places_in_rating) == 0:
                self.avg_Place_in_rating = float('NaN')
                self.Overall_score = 0
                return
            self.avg_Place_in_rating = sum(self.Places_in_rating)/len(self.Places_in_rating)
            self.Overall_score = (1000 / (self.avg_Place_in_rating + 50)) * ((len(self.Places_in_rating))**2)/1000
        elif method == 'individual_scores':
            all_DE_Scores = []
            for DE in self.DE_Info_by_GLM.values():
                if not math.isnan(DE.Score):
                    if math.isnan(DE.ttest_p):   all_DE_Scores.append(DE.Score**0.6 * ttest_missing_Score_multiplier)
                    else:   all_DE_Scores.append(DE.Score**0.6)
            self.Overall_score = (sum(all_DE_Scores)/len(self.DE_Info_by_GLM))**(1/0.6)

        else:
            raise('Calculate_Summary_Rating: incorrect method ' + method)

        return self.Overall_score


def ToBool(word, exit_if_error = True):
    word = word.casefold()
    if word in ['yes','y','on', 't', 'true']: return  True
    elif word in ['no','n','off', 'f', 'false']: return  False
    if exit_if_error:
        print('Incorrect input "%s"'%word)
        exit(127)
    else:
        return None

def to_float(value):
  try:    return float(value)
  except ValueError:   return float('NaN')

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Creating Excel workbooks from RTrans GO-centric expression prfiles')
    parser.add_argument('-m','--mode', dest='Mode', nargs='?',action='store', required=False,default='RTrans', help='')
    # parser.add_argument('--meta-r-model-names', dest='Mode', nargs='*',action='store', required=True,default='Meta-R', help='')
    parser.add_argument('-i','--in', dest='Input_FileNames', nargs='*',action='store', required=True,default=None, help='')
    parser.add_argument('-n','--in-names', dest='Sheet_Names', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('-cpm','--cpm-table', dest='Input_CPM_Table_FileNames', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('-o','--out-excel', dest='Workbook_FileName', nargs='?',action='store', required=True,default=None, help='')
    parser.add_argument('--logfc', dest='insert_LogFC', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--logfc-colname', dest='LogFC_colname', nargs='?',action='store', required=False,default='auto', help='')
    parser.add_argument('--pvalue-colname', dest='PValue_colname', nargs='?',action='store', required=False,default='auto', help='')
    parser.add_argument('--fdr-colname', dest='FDR_colname', nargs='?',action='store', required=False,default='auto', help='')
    parser.add_argument('--logcpm-colname', dest='LogCPM_colname', nargs='?',action='store', required=False,default='auto', help='')
    parser.add_argument('--logcpm', dest='insert_LogCPM', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--logcpm-common', dest='insert_LogCPM_common', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--cpms', dest='insert_CPMs', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--lr', dest='insert_LR', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--p', dest='insert_pvalue', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--fdr', dest='insert_FDR', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--score', dest='insert_Score', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--spearmanr', dest='insert_spearman_r', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--spearmanp', dest='insert_spearman_p', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--pearsonr', dest='insert_pearson_r', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--pearsonp', dest='insert_pearson_p', nargs='?',action='store', required=False,default='yes', help='')

    parser.add_argument('--mannp', dest='insert_mann_wh_p', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--wilcoxonp', dest='insert_wilcoxon_p', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--ttestp', dest='insert_ttest_p', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--mannfdr', dest='insert_mann_wh_FDR', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--wilcoxonfdr', dest='insert_wilcoxon_FDR', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--ttestfdr', dest='insert_ttest_FDR', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--readcountsok', dest='insert_read_counts_ok', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--rawmean', dest='insert_raw_mean', nargs='?',action='store', required=False,default='yes', help='')

    parser.add_argument('-s', '--sort-genes-by-overall-score', dest='Sort_Genes_by_Overall_Score', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--max-formatted-cells', dest='max_formatted_cells', nargs='?',action='store', required=False,default='35000', help='')
    parser.add_argument('--simple-logfc-format-mode', dest='simple_logfc_format_mode', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--exclude-genes-without-de-info', dest='exclude_genes_without_de_info', nargs='?',action='store', required=False,default='no', help='')

    parser.add_argument('--ttest-missing-score-multiplier', dest='ttest_missing_Score_multiplier', nargs='?',action='store', required=False,default='0.7', help='')
    args = parser.parse_args()

    for x in ['insert_LogFC', 'insert_LR','insert_pvalue', 'insert_FDR', 'insert_Score', 'insert_spearman_r',
              'insert_spearman_p', 'insert_pearson_r', 'insert_pearson_p','insert_CPMs','insert_LogCPM','insert_LogCPM_common',
              'simple_logfc_format_mode','exclude_genes_without_de_info', 'Sort_Genes_by_Overall_Score', 'insert_wilcoxon_p',
              'insert_ttest_p', 'insert_wilcoxon_FDR' , 'insert_ttest_FDR' , 'insert_read_counts_ok', 'insert_raw_mean',
              'insert_mann_wh_p', 'insert_mann_wh_FDR']:
        setattr(args,x, ToBool(getattr(args,x)))

    args.max_formatted_cells = int(args.max_formatted_cells)
    args.ttest_missing_Score_multiplier = float(args.ttest_missing_Score_multiplier)

    args.Mode = args.Mode.casefold()
    if not args.Mode in ['meta-r', 'rtrans']:
        print('args.Mode should be either Meta-R or RTrans')
        exit(127)

    
    Input_FileNames_Pack = []
    while True:
        if '+' in args.Input_FileNames:
            plus_pos = args.Input_FileNames.index('+')
            Input_FileNames_Pack.append(args.Input_FileNames[:plus_pos])
            args.Input_FileNames = args.Input_FileNames[plus_pos + 1 : ]
        else:
            Input_FileNames_Pack.append(args.Input_FileNames)
            break

    if args.Input_CPM_Table_FileNames is None:
        args.Input_CPM_Table_FileNames = [None] * len(Input_FileNames_Pack)

    if args.Sheet_Names is None:
        args.Sheet_Names = [None] * len(Input_FileNames_Pack)

    # print(Input_FileNames_Pack)
    # exit()

    if len(Input_FileNames_Pack) != len(args.Input_CPM_Table_FileNames) or len(args.Input_CPM_Table_FileNames) != len(args.Sheet_Names):
        # print(Input_FileNames_Pack)
        # print(args.Input_CPM_Table_FileNames)
        # print(args.Sheet_Names)
        print('Lengths of Input_FileNames_Pack, args.Input_CPM_Table_FileNames, args.Sheet_Names should be equal (args.Input_FileNames can be separated with "+" symbol)')
        exit(127)


    Workbook = xlsxwriter.Workbook(args.Workbook_FileName)

    # GO_Sheet = Workbook.add_worksheet('GO-centric DE profiles')
    bold = Workbook.add_format({'bold': True, 'italic': False})
    italic = Workbook.add_format({'bold': False, 'italic': True})
    bold_italic = Workbook.add_format({'bold': True, 'italic': True})
    # bold_italic.set_text_wrap()
    pvalue_formats = []
    for digits in range(5):
        pvalue_formats.append(Workbook.add_format())
        pvalue_formats[-1].set_num_format('0.%s'%('0'*digits))
    

    format0 = Workbook.add_format()
    format0.set_num_format('0')
    format1 = Workbook.add_format()
    format1.set_num_format('0.0')
    format2 = Workbook.add_format()
    format2.set_num_format('0.00')
    format_center = Workbook.add_format({'align': 'center'})
    format1_center = Workbook.add_format({'align': 'center'})
    format1_center.set_num_format('0.0')

    headers_format = Workbook.add_format({
        'bold': True,
        'italic': False,
        'border': True,
        'align': 'center',
        'valign': 'vcenter'})
    headers_format.set_text_wrap()

    secondary_headers_format = Workbook.add_format({
        'bold': False,
        'italic': False,
        'border': True,
        'align': 'center',
        'valign': 'vcenter'})
    secondary_headers_format.set_text_wrap()


    merge_format = Workbook.add_format({
        'bold': True,
        'italic': True,
        'border': True,
        'align': 'center',
        'valign': 'vcenter'})
    merge_format.set_text_wrap()


    sheet_N = 0
    for (Input_FileNames, Input_CPM_Table_FileName, Sheet_Name) in zip(Input_FileNames_Pack, args.Input_CPM_Table_FileNames, args.Sheet_Names):
        #print((Input_FileNames, Input_CPM_Table_FileName, Sheet_Name))
        sheet_N += 1

        Genes_by_ID = {}
        all_sample_names = []

        if not Input_CPM_Table_FileName is None and args.insert_CPMs:
            if os.path.exists(Input_CPM_Table_FileName):
                F = open(Input_CPM_Table_FileName,'r')
                sample_names = F.readline().rstrip().split('\t')
                all_sample_names = list(sample_names)
                for L in F.readlines():
                    C = L.rstrip().split('\t')
                    GeneID = C[0]
                    Genes_by_ID[GeneID] = TGene_Info(GeneID)
                    CPMs = [to_float(x) for x in C[1:]]
                    if len(CPMs) != len(sample_names):
                        print('%s: strings with sample names and CPMs do not match (cell counts %d vs %d)'%(Input_CPM_Table_FileName,
                            len(sample_names),len(CPMs)))
                        exit(127)

                    for x in range(len(CPMs)):
                        Genes_by_ID[GeneID].CPM_by_SampleName[sample_names[x]] = CPMs[x]
            else:
                print('Warning! CPM table file %s is not found'%(Input_CPM_Table_FileName))
                args.insert_CPMs = False


        GLM_predictor_values_by_predictor_by_SampleName = dict()

        # P_CPM_combinations = []
        # P_CPM_combination_codes = []
        # GO_terms_list = []
        # min_abs_LogFC = 0.4
        # trunc_sheet_number = 1


        Any_RefSeq_Summary_is_present = False
        LR_is_present_by_model = dict()
        F_is_present_by_model = dict()
        raw_mean_is_present_by_model = dict()

        all_Models = []
        for FileName in Input_FileNames:
            #'Age_edgeR_combined.tsv'
            if args.Mode == 'meta-r':
                # 'DTA results, Phenol Fri. 100x vs Fri. nt [complete]'
                base = os.path.split(FileName)[1]
                model = base[:base.rfind(' [')].replace('DTA results, ', '')

            elif 'edger' in FileName:
                try:  model = os.path.split(FileName.split('_edger_combined')[-2])[1]
                except ValueError:
                    print('Incorrect file name. Cannot detect type/number of top genes. Correct format: Age, GSEA for top-40 downreg genes.txt')
                    exit()
            elif 'deseq' in FileName:
                try:  model = os.path.split(FileName.split('_deseq_combined')[-2])[1]
                except ValueError:
                    print('Incorrect file name. Cannot detect type/number of top genes. Correct format: Age, GSEA for top-40 downreg genes.txt')
                    exit()
            elif 'spearman' in FileName:
                try:  model = os.path.split(FileName.split('_spearman_combined')[-2])[1]
                except ValueError:
                    print('Incorrect file name. Cannot detect type/number of top genes. Correct format: Age, GSEA for top-40 downreg genes.txt')
                    exit()
            elif 'pearson' in FileName:
                try:  model = os.path.split(FileName.split('_pearson_combined')[-2])[1]
                except ValueError:
                    print('Incorrect file name. Cannot detect type/number of top genes. Correct format: Age, GSEA for top-40 downreg genes.txt')
                    exit()
            elif 'combined_cor' in FileName:
                try:  model = os.path.split(FileName.split('_combined_cor_combined')[-2])[1]
                except ValueError:
                    print('Incorrect file name. Cannot detect type/number of top genes. Correct format: Age, GSEA for top-40 downreg genes.txt')
                    exit()
            elif 'external_data' in FileName:
                try:  model = os.path.split(FileName.split('_external_data')[-2])[1]
                except ValueError:
                    print('Incorrect file name. Cannot detect type/number of top genes. Correct format: Age, GSEA for top-40 downreg genes.txt')
                    exit()
            else:
                print('Incorrect file name. Cannot detect type/number of top genes. Correct format: Age, GSEA for top-40 downreg genes.txt')
                exit()

            all_Models.append(model)
            File = open(FileName,'r')
            L = File.readline()
            L = L.replace('\r','').replace('\n','').replace('"','')
            C = L.split('\t')
            first_Line_cells = C
            if args.Mode == 'rtrans':
                Stats_cells_count = 2 + max([C.index(x) for x in ['p (Pearson)', 'p (Spearman)', 'PValue', 'FDR', 'LogCPM', 'LogFC', 'logFC', 'TBA LogFC',
                    'pearson p-value', 'spearman p-value', 'Pearson r', 'Spearman r', 'Wilcoxon p-value', 'Wilcoxon FDR', 'Mann-Wh. FDR', 't-test FDR', 'Mann-Wh. p-value', 't-test p-value'] if x in C])
                CPM_cells_count = len(C) - Stats_cells_count + 1
            else:
                Stats_cells_count = 2 + max([C.index(x) for x in ['pearson p-value', 'spearman p-value', 'pearson rs', 'pearson rs', 'Mann-Wh. FDR', 't-test FDR', 'Mann-Wh. p-value', 't-test p-value'] if x in C])
                CPM_cells_count = len(C) - Stats_cells_count + 1

            cell_index_by_name = {}
            for n in range(Stats_cells_count-1):
                cell_index_by_name[C[n]] = n + 1


            # if args.Mode == 'rtrans':
            if args.LogFC_colname in ['auto', '{auto}']:
                current_LogFC_colname = None
                possible_LogFC_colnames = ['TBA LogFC', 'trimmed LogFC', 'Trimmed LogFC', 'tLogFC', 'LogFC', 'logFC']
                for cn in possible_LogFC_colnames:
                    if cn in cell_index_by_name:
                        current_LogFC_colname = cn
                if current_LogFC_colname is None:
                    print('Cannot find LogFC column name. Available only:')
                    print(cell_index_by_name)
                    continue
            else:
                current_LogFC_colname = args.LogFC_colname

            if args.PValue_colname in ['auto', '{auto}']:
                current_PValue_colname = None
                possible_PValue_colnames = ["p (QLF test)", "p (LR test)", "p (ex. test)", 'PValue', 'p-value', 'P-Value', 'Pv', "P", "p"]
                for cn in possible_PValue_colnames:
                    if cn in cell_index_by_name:
                        current_PValue_colname = cn
                if current_PValue_colname is None and args.Mode != 'meta-r':
                    print('Cannot find PValue column name. Available only:')
                    print(cell_index_by_name)
                    continue
            else:
                current_PValue_colname = args.PValue_colname

            if args.FDR_colname in ['auto', '{auto}']:
                current_FDR_colname = None
                possible_FDR_colnames = ["FDR (QLF test)", "FDR (LR test)", "FDR (ex. test)", 'FDR' ,'fdr',
                   'FDR (QLF test paired)', 'FDR (LR test paired)', 'FDR (Wilcoxon paired)']
                for cn in possible_FDR_colnames:
                    if cn in cell_index_by_name:
                        current_FDR_colname = cn
                if current_FDR_colname is None and args.Mode != 'meta-r':
                    print('Cannot find FDR column name. Available only:')
                    print(cell_index_by_name)
                    continue
            else:
                current_FDR_colname = args.FDR_colname

            if args.LogCPM_colname in ['auto', '{auto}']:
                current_LogCPM_colname = None
                possible_LogCPM_colnames = ["LogCPM", "logCPM"]
                for cn in possible_LogCPM_colnames:
                    if cn in cell_index_by_name:
                        current_LogCPM_colname = cn
                if current_LogCPM_colname is None:
                    print('Cannot find LogCPM column name. Available only:')
                    print(cell_index_by_name)
                    continue
            else:
                current_LogCPM_colname = args.LogCPM_colname

                
                

            if args.Mode == 'rtrans' and any([(x not in cell_index_by_name) for x in ("Symbol", "Biotype", "Name", "Score", current_LogFC_colname, current_LogCPM_colname, current_PValue_colname, current_FDR_colname, "Spearman r", "p (Spearman)", "Pearson r", "p (Pearson)")]):
                print('Incorrect file %s'%FileName)
                print('The following columns should present: "Symbol", "Biotype", "Name", "RefSeq_Summary", "Score", "%s", "%s", "%s", "%s", "Spearman r", "p (Spearman)", "Pearson r", "p (Pearson)"'%(
                    current_LogFC_colname, current_LogCPM_colname, current_PValue_colname, current_FDR_colname))
                print('Available only:')
                print(cell_index_by_name)
                continue
            elif args.Mode == 'meta-r' and any([(x not in cell_index_by_name) for x in ("Score", current_LogCPM_colname, current_LogFC_colname)]):
                print('Incorrect file %s'%FileName)
                print('The following columns should present: "Score", "%s", "%s"'%(current_LogCPM_colname, current_LogFC_colname))
                print('Available only:')
                print(cell_index_by_name)
                continue

            MetaR_to_RTrans_colnames = {'LogCPM':'logCPM', 'spearman rs':'Spearman r', 'spearman p-value':'p (Spearman)',
                                        'pearson rs':'Pearson r', 'pearson p-value':'p (Pearson)', 'Mann-Wh. p-value':'Mann-Wh. P',
                                        'Mann-Wh. FDR':'Mann-Wh. FDR', 't-test p-value':'t-test P' }

            for x in MetaR_to_RTrans_colnames:
                if x in cell_index_by_name:
                    cell_index_by_name[MetaR_to_RTrans_colnames[x]] = cell_index_by_name[x]

            src_raw_mean_cols = [cell_index_by_name[x] for x in cell_index_by_name if ('raw mean' in x)]

            LR_is_present_by_model[model] = 'LR' in cell_index_by_name
            F_is_present_by_model[model] = 'F' in cell_index_by_name
            raw_mean_is_present_by_model[model] = len(src_raw_mean_cols) > 0

            # sheet_name = model.replace('[','(').replace(']',')').replace('?','_').replace(':','(d)').replace('*','(m)').replace('\\','_').replace('/','_')
            # if len(sheet_name) > 30:
            #     sheet_name = sheet_name[:25] + '...%d'%trunc_sheet_number
            #     trunc_sheet_number += 1


            # sheet = Workbook.add_worksheet(sheet_name)

            RefSeq_Summary_is_present = ('RefSeq_Summary' in cell_index_by_name)

            #
            # sheet.write(0,0,'Gene ID',bold)
            # sheet.set_column(0,0,19)
            # sheet.write(0,1,'Symbol',bold)
            # sheet.write(0,2,'Biotype',bold)
            # sheet.set_column(2,2,18)
            # sheet.write(0,3,'Name',bold)
            # sheet.set_column(3,3,35)
            # if RefSeq_Summary_is_present:  sheet.write(0,4,'Summary',bold)
            # sheet.write(0,4 + RefSeq_Summary_is_present,'Score',bold)
            # sheet.write(0,5 + RefSeq_Summary_is_present,'LogFC',bold)
            # sheet.write(0,6 + RefSeq_Summary_is_present,'LogCPM',bold)
            # sheet.write(0,7 + RefSeq_Summary_is_present,'LR',bold)
            # sheet.write(0,8 + RefSeq_Summary_is_present,'GLM P',bold)
            # sheet.write(0,9 + RefSeq_Summary_is_present,'GLM FDR',bold)
            # sheet.write(0,10 + RefSeq_Summary_is_present,'Spear. r',bold)
            # sheet.write(0,11 + RefSeq_Summary_is_present,'Spear. P',bold)
            # sheet.write(0,12 + RefSeq_Summary_is_present,'Pearson r',bold)
            # sheet.write(0,13 + RefSeq_Summary_is_present,'Pearson P',bold)
            # sheet.write(0,14 + RefSeq_Summary_is_present,'CPM:',italic)

            # for n in range(CPM_cells_count):
            #     sheet.write(0,15 + RefSeq_Summary_is_present + n,C[13 + RefSeq_Summary_is_present + n])

            stringN = 0
            for L in open(FileName,'r').readlines()[1:]:
                stringN += 1
                L = L.replace('\r','').replace('\n','').replace('"','')
                C = L.split('\t')
                if C[cell_index_by_name[current_LogFC_colname]] == '' and args.Mode == 'rtrans': ## predictors line
                    pred_name = C[0]
                    for n in range(Stats_cells_count,Stats_cells_count+CPM_cells_count):
                        if not pred_name in GLM_predictor_values_by_predictor_by_SampleName:
                            GLM_predictor_values_by_predictor_by_SampleName[pred_name] = dict()
                        sample_name = first_Line_cells[n-1]
                        try:
                            GLM_predictor_values_by_predictor_by_SampleName[pred_name][sample_name] = float(C[n])
                        except ValueError:
                            if not sample_name in GLM_predictor_values_by_predictor_by_SampleName[pred_name]:
                                GLM_predictor_values_by_predictor_by_SampleName[pred_name][sample_name] = float('NaN')

                        if not sample_name in all_sample_names:
                            all_sample_names.append(sample_name)

                    continue

                ## casual line
                GeneID = C[0]
                if not GeneID in Genes_by_ID:   Genes_by_ID[GeneID] = TGene_Info(GeneID)
                Gene = Genes_by_ID[GeneID]

                Gene.Places_in_rating.append(stringN)

                if args.Mode == 'rtrans':
                    Gene.Symbol = C[cell_index_by_name['Symbol']]
                    Gene.Biotype = C[cell_index_by_name['Biotype']]
                    Gene.Name = C[cell_index_by_name['Name']]
                    if Gene.RefSeq_summary == '' and RefSeq_Summary_is_present:
                        Gene.RefSeq_summary = C[cell_index_by_name['RefSeq_Summary']]

                DE = TDE_Info()
                Gene.DE_Info_by_GLM[model] = DE
                # DE.Score = C[cell_index_by_name['Score']]

                if args.insert_Score:
                    try:   DE.Score = float(C[cell_index_by_name['Score']])
                    except ValueError: pass

                if args.insert_LogFC:
                    try:   DE.LogFC = float(C[cell_index_by_name[current_LogFC_colname]])
                    except ValueError: pass

                # if args.insert_LogCPM
                try:   DE.LogCPM = float(C[cell_index_by_name[current_LogCPM_colname]])
                except ValueError: pass

                if args.Mode == 'meta-r' and args.insert_raw_mean:
                    raws = []
                    for x in src_raw_mean_cols:
                        try:  r = float(C[x])
                        except ValueError:
                            raws = None
                            break
                        if math.isnan(r):
                            raws = None
                            break
                        raws.append(r)

                    if raws is not None:
                        if len(raws) > 0:
                            DE.raw_mean = sum(raws)/len(raws)


                if args.Mode == 'rtrans':
                    try:
                        if 'LR' in cell_index_by_name:  DE.LR = float(C[cell_index_by_name['LR']])
                    except ValueError: pass

                    try:
                        if 'F' in cell_index_by_name:  DE.F = float(C[cell_index_by_name['F']])
                    except ValueError: pass

                    try:   DE.GLM_P = float(C[cell_index_by_name[current_PValue_colname]])
                    except ValueError: pass

                    try:   DE.GLM_FDR = float(C[cell_index_by_name[current_FDR_colname]])
                    except ValueError: pass

                if args.insert_spearman_r:
                    try:   DE.Spearman_r = float(C[cell_index_by_name['Spearman r']])
                    except ValueError: pass

                if args.insert_spearman_p:
                    try:   DE.Spearman_p = float(C[cell_index_by_name['p (Spearman)']])
                    except ValueError: pass

                if args.insert_pearson_r:
                    try:   DE.Pearson_r = float(C[cell_index_by_name['Pearson r']])
                    except ValueError: pass

                if args.insert_pearson_p:
                    try:   DE.Pearson_p = float(C[cell_index_by_name['p (Pearson)']])
                    except ValueError: pass

                if args.Mode == 'meta-r':
                    if args.insert_ttest_p:
                        try:   DE.ttest_p = float(C[cell_index_by_name['t-test P']])
                        except: pass

                    if args.insert_ttest_FDR:
                        try:   DE.ttest_FDR = float(C[cell_index_by_name['t-test FDR']])
                        except: pass

                    if args.insert_wilcoxon_p:
                        try:   DE.Wilcoxon_p = float(C[cell_index_by_name['Wilcoxon P']])
                        except: pass

                    if args.insert_wilcoxon_FDR:
                        try:   DE.Wilcoxon_FDR = float(C[cell_index_by_name['Wilcoxon FDR']])
                        except: pass

                    if args.insert_mann_wh_p:
                        try:   DE.Mann_Wh_p = float(C[cell_index_by_name['Mann-Wh. P']])
                        except: pass

                    if args.insert_mann_wh_FDR:
                        try:   DE.Mann_Wh_FDR = float(C[cell_index_by_name['Mann-Wh. FDR']])
                        except: pass

                    if args.insert_read_counts_ok:
                        try:   DE.ReadCounts_Passed = ToBool(C[cell_index_by_name['read counts OK']], exit_if_error = False)
                        except: pass


                # if abs(DE.LogFC) > 0.3:  FormatLogFC_Cell(sheet,2**logFC,stringN,5 + RefSeq_Summary_is_present,coeff = 1.4)

                # for n in range(CPM_cells_count):
                #     try:  Gene.CPM_by_SampleName[] = float(C[14 + RefSeq_Summary_is_present + n])
                #     except ValueError:  Gene.CPM_by_SampleName[] = float('NaN')
                #
                # for n in range(CPM_cells_count):
                #     sheet.write_number(stringN,15 + RefSeq_Summary_is_present + n,float(C[14 + RefSeq_Summary_is_present + n]))

            Any_RefSeq_Summary_is_present = Any_RefSeq_Summary_is_present or RefSeq_Summary_is_present



        Gene_ID_sequence = [Genes_by_ID[x] for x in sorted(Genes_by_ID.keys())]
        for Gene in Gene_ID_sequence:
            if args.Mode == 'meta-r':
                Gene.Calculate_Summary_Rating(method = 'individual_scores', ttest_missing_Score_multiplier = args.ttest_missing_Score_multiplier)
            else:
                Gene.Calculate_Summary_Rating()

        if args.Sort_Genes_by_Overall_Score:
            Gene_ID_sequence.sort(key = (lambda x: x.Overall_score),reverse=True)


        # header_nonmerge_format = Workbook.add_format({
        #     'bold': True,
        #     'italic': True,
        #     'border': True,
        #     'align': 'center',
        #     'valign': 'vcenter'})
        # header_nonmerge_format.set_text_wrap()


        LogFC_cols = []
        LogCPM_cols = []
        LR_cols = []
        P_cols = []
        FDR_cols = []
        raw_mean_cols = []
        Score_cols = []
        Rs_cols = []
        Filter_cols = []

        if Sheet_Name is None:
            sheet = Workbook.add_worksheet('Joint data (%d)'%sheet_N)
        else:
            sheet = Workbook.add_worksheet(Sheet_Name)

        # sheet.set_row(1, 1, 30)

        if args.Mode == 'rtrans':
            sheet.write(1,0,'Gene ID',bold)
            sheet.set_column(0,0,18)
            sheet.write(1,1,'Symbol',bold)
            sheet.write(1,2,'Biotype',bold)
            sheet.set_column(2,2,18)
            sheet.write(1,3,'Name',bold)
            sheet.set_column(3,3,30)
            if args.insert_LogCPM_common:
                sheet.write(1,4,'LogCPM (mean)',bold)
                LogCPM_cols.append(4)
                # if args.Input_CPM_Table_FileName != None:
                #     Gene
            if Any_RefSeq_Summary_is_present:
                sheet.write(1, 4 + args.insert_LogCPM_common, 'RefSeq Summary', bold)
                sheet.set_column(4 + args.insert_LogCPM_common, 4 + args.insert_LogCPM_common, 12)

            ColN = 5 + Any_RefSeq_Summary_is_present + args.insert_LogCPM_common
        else:
            sheet.write(1,0,'Taxon',bold)
            sheet.set_column(0,0,18)
            ColN = 1


        sheet.set_column(ColN,ColN,2)
        if args.Mode == 'rtrans':
            GLM_info_columns_count = sum([getattr(args,x) for x in ['insert_LogFC', 'insert_LR', 'insert_LogCPM','insert_pvalue',
                                                                    'insert_FDR', 'insert_Score','insert_spearman_r','insert_spearman_p',
                                                                    'insert_pearson_r', 'insert_pearson_p']])
        else:
            GLM_info_columns_count = sum([getattr(args,x) for x in ['insert_LogFC', 'insert_LogCPM','insert_wilcoxon_p','insert_mann_wh_p',
                                                                    'insert_ttest_p', 'insert_wilcoxon_FDR', 'insert_mann_wh_FDR', 'insert_ttest_FDR',
                                                                    'insert_read_counts_ok',
                                                                    'insert_Score','insert_spearman_r','insert_spearman_p',
                                                                    'insert_pearson_r', 'insert_pearson_p', 'insert_raw_mean']])



        for model in all_Models:
            ColN += 1
            if GLM_info_columns_count > 1:
                sheet.merge_range(0,ColN,0,ColN + GLM_info_columns_count-1,model,merge_format)
            else:
                sheet.write(0,ColN,model,merge_format)

            if args.insert_LogFC:   sheet.write(1,ColN,current_LogFC_colname,bold); LogFC_cols.append(ColN); ColN += 1
            if args.insert_LogCPM:   sheet.write(1,ColN,'LogCPM',bold);  LogCPM_cols.append(ColN); ColN += 1
            if args.Mode == 'rtrans':
                if args.insert_LR and LR_is_present_by_model[model]:    sheet.write(1,ColN,'LR test',bold); LR_cols.append(ColN); ColN += 1
                elif args.insert_LR and F_is_present_by_model[model]:    sheet.write(1,ColN,'F test',bold); LR_cols.append(ColN); ColN += 1
                if args.insert_pvalue:   sheet.write(1,ColN,'P-value',bold);  P_cols.append(ColN); ColN += 1
                if args.insert_FDR:   sheet.write(1,ColN,'FDR',bold);  FDR_cols.append(ColN); ColN += 1
            if args.Mode == 'meta-r':
                if args.insert_raw_mean and raw_mean_is_present_by_model[model]:   sheet.write(1,ColN,'raw mean',bold);  raw_mean_cols.append(ColN); ColN += 1
            if args.insert_Score:   sheet.write(1,ColN,'Score',bold);  Score_cols.append(ColN); ColN += 1
            if args.Mode == 'meta-r':
                if args.insert_read_counts_ok:   sheet.write(1,ColN,'reads count OK',bold);  Filter_cols.append(ColN); ColN += 1
                if args.insert_ttest_p:   sheet.write(1,ColN,'t-test P',bold);  P_cols.append(ColN); ColN += 1
                if args.insert_ttest_FDR:   sheet.write(1,ColN,'t-test FDR',bold);  FDR_cols.append(ColN); ColN += 1
                if args.insert_mann_wh_p:   sheet.write(1,ColN,'Mann-Wh. P',bold);  P_cols.append(ColN); ColN += 1
                if args.insert_mann_wh_FDR:   sheet.write(1,ColN,'Mann-Wh. FDR',bold);  FDR_cols.append(ColN); ColN += 1
                if args.insert_wilcoxon_p:   sheet.write(1,ColN,'Wilcoxon P',bold);  P_cols.append(ColN); ColN += 1
                if args.insert_wilcoxon_FDR:   sheet.write(1,ColN,'Wilcoxon FDR',bold);  FDR_cols.append(ColN); ColN += 1
            if args.insert_spearman_r:   sheet.write(1,ColN,'Spearman r',bold);  Rs_cols.append(ColN); ColN += 1
            if args.insert_spearman_p:   sheet.write(1,ColN,'p (Spearman)',bold);  P_cols.append(ColN); ColN += 1
            if args.insert_pearson_r:   sheet.write(1,ColN,'Pearson r',bold);  Rs_cols.append(ColN); ColN += 1
            if args.insert_pearson_p:   sheet.write(1,ColN,'p (Pearson)',bold);  P_cols.append(ColN); ColN += 1
            sheet.set_column(ColN,ColN,3)

        all_predictors = []
        if args.Mode == 'rtrans':
            sheet.set_column(ColN+1,ColN+1,3)
            sheet.write(1,ColN,'CPMs',italic)
            ColN += 2
            CPMs_start_col = ColN
            all_sample_names.sort()

            # Max_LogFC_formatted_cells_per_model = int(float(args.max_formatted_cells)/float(len(all_Models)))

            if args.insert_CPMs:
                all_predictors = sorted(GLM_predictor_values_by_predictor_by_SampleName)
                if len(all_predictors) > 0:
                    for x in range(len(all_predictors)):
                        sheet.write(2+x, 0, all_predictors[x],italic)
                for sample_n in range(len(all_sample_names)):
                    sample_name = all_sample_names[sample_n]
                    sheet.write(1,CPMs_start_col + sample_n,sample_name,italic)
                    for pred_n in range(len(all_predictors)):
                        pred = all_predictors[pred_n]
                        if sample_name not in GLM_predictor_values_by_predictor_by_SampleName[pred]: continue
                        val = GLM_predictor_values_by_predictor_by_SampleName[pred][sample_name]
                        value = GLM_predictor_values_by_predictor_by_SampleName[pred][sample_name]
                        if not math.isnan(val):
                            sheet.write(2 + pred_n,CPMs_start_col + sample_n, value)

            RowN = 2 + len(all_predictors)
        else:
            RowN = 2

        formatted_cells_count = 0
        for Gene in Gene_ID_sequence:
            if args.exclude_genes_without_de_info and len(Gene.DE_Info_by_GLM) == 0:
                continue

            if args.Mode == 'rtrans':
                sheet.write(RowN,0,Gene.GeneID)
                sheet.write(RowN,1,Gene.Symbol)
                sheet.write(RowN,2,Gene.Biotype)
                sheet.write(RowN,3,Gene.Name)
                if args.insert_LogCPM_common:
                    meanCPM = mean(list(Gene.CPM_by_SampleName.values()))
                    if meanCPM > 0:
                        sheet.write_number(RowN,4,math.log2(meanCPM),format2)
                    else:
                        sheet.write(RowN, 4, 'nan')

                if Any_RefSeq_Summary_is_present:
                    sheet.write(RowN,4+args.insert_LogCPM_common,Gene.RefSeq_summary)
                    sheet.write(RowN,6+args.insert_LogCPM_common,' ')
                ColN = 5 + Any_RefSeq_Summary_is_present + args.insert_LogCPM_common
            
            else:
                sheet.write(RowN, 0, Gene.GeneID)
                ColN = 1


            for model in all_Models:
                ColN += 1
                if model not in Gene.DE_Info_by_GLM:
                    ColN += GLM_info_columns_count
                    continue
                DE = Gene.DE_Info_by_GLM[model]
                if args.insert_LogFC:
                    sheet.write_number(RowN,ColN,DE.LogFC,format2)
                    if formatted_cells_count < args.max_formatted_cells and not args.simple_logfc_format_mode:
                        if abs(DE.LogFC) > 0.3:
                            FormatLogFC_Cell(sheet, 2**DE.LogFC, RowN, ColN, coeff1 = 1.4, coeff2 = 1.4)
                            formatted_cells_count += 1
                    ColN += 1

                if args.insert_LogCPM:
                    write_number__mod(sheet, RowN, ColN, DE.LogCPM, format1)
                    ColN += 1

                if args.Mode == 'rtrans':
                    if args.insert_LR and LR_is_present_by_model[model]:
                        write_number__mod(sheet, RowN, ColN, DE.LR, format1)
                        ColN += 1
                    elif args.insert_LR and F_is_present_by_model[model]:
                        write_number__mod(sheet, RowN, ColN, DE.F, format1)
                        ColN += 1

                    if args.insert_pvalue:
                        write_number__mod(sheet,RowN,ColN,DE.GLM_P, pvalue_mode = True, pvalue_formats = pvalue_formats)
                        ColN += 1
                    if args.insert_FDR:
                        write_number__mod(sheet,RowN,ColN,DE.GLM_FDR, pvalue_mode = True, pvalue_formats = pvalue_formats)
                        ColN += 1
                if args.Mode == 'meta-r':
                    if args.insert_raw_mean:
                        if DE.raw_mean is not None:
                            write_number__mod(sheet, RowN, ColN, DE.raw_mean, format1)
                            ColN += 1

                if args.insert_Score:
                    write_number__mod(sheet, RowN, ColN, DE.Score, format1)
                    ColN += 1

                if args.Mode == 'meta-r':
                    if args.insert_read_counts_ok:
                        if ReadCounts_Passed is not None:
                            if ReadCounts_Passed:
                                sheet.write(RowN, ColN, 'yes')
                        else:
                            sheet.write(RowN, ColN, 'no data')
                    if args.insert_ttest_p:
                        write_number__mod(sheet,RowN,ColN,DE.ttest_p, pvalue_mode = True, pvalue_formats = pvalue_formats)
                        ColN += 1
                    if args.insert_ttest_FDR:
                        write_number__mod(sheet,RowN,ColN,DE.ttest_FDR, pvalue_mode = True, pvalue_formats = pvalue_formats)
                        ColN += 1
                    if args.insert_mann_wh_p:
                        write_number__mod(sheet,RowN,ColN,DE.Mann_Wh_p, pvalue_mode = True, pvalue_formats = pvalue_formats)
                        ColN += 1
                    if args.insert_mann_wh_FDR:
                        write_number__mod(sheet,RowN,ColN,DE.Mann_Wh_FDR, pvalue_mode = True, pvalue_formats = pvalue_formats)
                        ColN += 1
                    if args.insert_wilcoxon_p:
                        write_number__mod(sheet,RowN,ColN,DE.Wilcoxon_p, pvalue_mode = True, pvalue_formats = pvalue_formats)
                        ColN += 1
                    if args.insert_wilcoxon_FDR:
                        write_number__mod(sheet,RowN,ColN,DE.Wilcoxon_FDR, pvalue_mode = True, pvalue_formats = pvalue_formats)
                        ColN += 1


                if args.insert_spearman_r:
                    write_number__mod(sheet,RowN,ColN,DE.Spearman_r,format2)
                    ColN += 1
                if args.insert_spearman_p:
                    write_number__mod(sheet,RowN,ColN,DE.Spearman_p, pvalue_mode = True, pvalue_formats = pvalue_formats)
                    ColN += 1
                if args.insert_pearson_r:
                    write_number__mod(sheet, RowN, ColN, DE.Pearson_r, format2)
                    ColN += 1
                if args.insert_pearson_p:
                    write_number__mod(sheet, RowN, ColN, DE.Pearson_p, pvalue_mode = True, pvalue_formats = pvalue_formats)
                    ColN += 1

            ColN += 2
            CPMs_start_col = ColN

            if args.insert_CPMs:
                for sample_n in range(len(all_sample_names)):
                    sample_name = all_sample_names[sample_n]
                    if sample_name not in Gene.CPM_by_SampleName: continue
                    if not math.isnan(Gene.CPM_by_SampleName[sample_name]):
                        sheet.write(RowN,CPMs_start_col + sample_n,Gene.CPM_by_SampleName[sample_name],format2)

            RowN += 1

        Last_Row = RowN

        if args.simple_logfc_format_mode:
            if args.Mode == 'meta-r':   max_value = 5.1
            else:  max_value = 5.1
            for ColN in LogFC_cols:
                sheet.conditional_format(2 + len(all_predictors), ColN, Last_Row, ColN,
                                         {'type': '3_color_scale',
                                                          'min_color': "#1170f1",'mid_color': "#ffffff",'max_color': "#ec4a18",
                                                          'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                          'min_value': -max_value, 'mid_value': 0.0, 'max_value': max_value})

        if args.Mode == 'meta-r':
            for ColN in Score_cols:
                sheet.conditional_format(1, ColN, Last_Row, ColN, {'type': '3_color_scale',
                                                                   'min_color': "#ffffff", 'mid_color': "#ffe08d",
                                                                   'max_color': "#9ece49",
                                                                   'min_type': 'num', 'mid_type': 'num',
                                                                   'max_type': 'num',
                                                                   'min_value': 2, 'mid_value': 12,
                                                                   'max_value': 25})
        
        for ColN in Rs_cols:
            sheet.conditional_format(1,ColN,Last_Row,ColN, {'type': '3_color_scale',
                                                      'min_color': "#1170f1",'mid_color': "#ffffff",'max_color': "#ec4a18",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': -0.98, 'mid_value': 0.0, 'max_value': 0.98})
                

        for ColN in LogCPM_cols:
            if args.Mode == 'meta-r':   max_value = 20
            else:  max_value = 9
            sheet.conditional_format(2 + len(all_predictors),ColN,Last_Row,ColN, {'type': 'data_bar','bar_color': '#ffb910',
                                                                                  'min_type':'num','max_type':'num',
                                                                                  'min_value':0,'max_value':max_value})

        for ColN in P_cols + FDR_cols:
            sheet.conditional_format(2 + len(all_predictors),ColN,Last_Row,ColN, {'type': '3_color_scale',
                                                      'min_color': "#9ece49",'mid_color': "#ffe08d",'max_color': "#ffffff",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': 0, 'mid_value': 0.0002, 'max_value': 0.07})

        # Light red fill with dark red text.
        format_lincRNA = Workbook.add_format({'bg_color':   '#c3d594'})

        # Light yellow fill with dark yellow text.
        format_antisense = Workbook.add_format({'bg_color':   '#d9c188'})

        # Green fill with dark green text.
        format_pseudogene = Workbook.add_format({'bg_color':   '#a0b4c7'})

        if args.Mode == 'rtrans':
            sheet.conditional_format(2 + len(all_predictors),2,Last_Row,2, {'type':'text',
                                               'criteria': 'containing', 'value':    'lincRNA', 'format':   format_lincRNA})
            sheet.conditional_format(2 + len(all_predictors),2,Last_Row,2, {'type':'text',
                                               'criteria': 'containing', 'value':    'antisense', 'format':   format_antisense})
            sheet.conditional_format(2 + len(all_predictors),2,Last_Row,2, {'type':'text',
                                               'criteria': 'containing', 'value':    'pseudogene', 'format':   format_pseudogene})

    Workbook.close()


