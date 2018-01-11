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

ScoreColorPoints = [(225,225,225),(0,0,0)]
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

def write_number__mod(sheet,rowN,colN,text,format=None):
    try:
        d = float(text)
        if format != None:
            if math.isnan(d) or math.isinf(d): sheet.write(rowN,colN,text,format)
            else: sheet.write_number(rowN,colN,d,format)
        else:
            if math.isnan(d) or math.isinf(d): sheet.write(rowN,colN,text)
            else: sheet.write_number(rowN,colN,d)
    except:
        if format != None:
            sheet.write(rowN,colN,text,format)
        else:
            sheet.write(rowN,colN,text)


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

def FormatLogFC_Cell(sheet,FC,RowN,ColN,coeff = 1.0):
    if math.isnan(FC): return
    OverexpressionMaxColor = (226,101,0)
    DownregulationMaxColor = (29,136,234)

    LogFC = math.log2(FC)
    if LogFC >= 0:
        C = color(Gradient(((255,255,255),OverexpressionMaxColor), min(1,(coeff*LogFC + 0.1)/2.5)))
        sheet.conditional_format(RowN,ColN,RowN,ColN, {'type': 'data_bar','bar_color': C,
                                         'min_type':'num','max_type':'num',
                                         'min_value':0,'max_value':7/coeff})
    else:
        C = color(Gradient(((255,255,255),DownregulationMaxColor), min(1,(coeff*LogFC*(-1) + 0.1)/2.5)))
        sheet.conditional_format(RowN,ColN,RowN,ColN, {'type': 'data_bar','bar_color': C,
                                         'min_type':'num','max_type':'num',
                                         'min_value':LogFC*2,'max_value':LogFC*2+7/coeff})


def FormatLogCPM_Cell(sheet,logCPM,RowN,ColN):
    if math.isnan(logCPM): return
    maxColor = (220,181,0)

    if logCPM <= 0: return
    C = color(Gradient(((255,255,255),maxColor), min(1,(logCPM + 0.1)/8)))
    sheet.conditional_format(RowN,ColN,RowN,ColN, {'type': 'data_bar','bar_color': C,
                                     'min_type':'num','max_type':'num',
                                     'min_value':0,'max_value':9})

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

class TDE_Info():
    def __init__(self):
        self.LogFC = float('NaN')
        self.LogCPM = float('NaN')
        self.LR = float('NaN')
        self.GLM_P = float('NaN')
        self.GLM_FDR = float('NaN')
        self.Score = float('NaN')
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

    def Calculate_Summary_Rating(self):
        if len(self.Places_in_rating) == 0:
            self.avg_Place_in_rating = float('NaN')
            self.Overall_score = 0
            return
        self.avg_Place_in_rating = sum(self.Places_in_rating)/len(self.Places_in_rating)
        self.Overall_score = (1000 / (self.avg_Place_in_rating + 50)) * ((len(self.Places_in_rating))**2)/1000
        return self.Overall_score


def ToBool(word):
    word = word.casefold()
    if word in ['yes','y','on']: return  True
    elif word in ['no','n','off']: return  False
    print('Incorrect input "%s"'%word)
    exit(127)

def to_float(value):
  try:    return float(value)
  except ValueError:   return float('NaN')

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Creating Excel workbooks from RTrans GO-centric expression prfiles')
    parser.add_argument('-i','--in', dest='Input_FileNames', nargs='*',action='store', required=True,default=None, help='')
    parser.add_argument('-cpm','--cpm-table', dest='Input_CPM_Table_FileName', nargs='?',action='store', required=False,default=None, help='')
    parser.add_argument('-o','--out-excel', dest='Workbook_FileName', nargs='?',action='store', required=True,default=None, help='')
    parser.add_argument('--logfc', dest='insert_LogFC', nargs='?',action='store', required=False,default='yes', help='')
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
    parser.add_argument('--max-formatted-cells', dest='max_formatted_cells', nargs='?',action='store', required=False,default='35000', help='')
    parser.add_argument('--simple-logfc-format-mode', dest='simple_logfc_format_mode', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--exclude-genes-without-de-info', dest='exclude_genes_without_de_info', nargs='?',action='store', required=False,default='no', help='')
    args = parser.parse_args()

    for x in ['insert_LogFC', 'insert_LR','insert_pvalue', 'insert_FDR', 'insert_Score', 'insert_spearman_r',
              'insert_spearman_p', 'insert_pearson_r', 'insert_pearson_p','insert_CPMs','insert_LogCPM','insert_LogCPM_common',
              'simple_logfc_format_mode','exclude_genes_without_de_info']:
        setattr(args,x, ToBool(getattr(args,x)))

    args.max_formatted_cells = int(args.max_formatted_cells)

    Genes_by_ID = {}
    all_sample_names = []

    if not args.Input_CPM_Table_FileName is None and args.insert_CPMs:
        if os.path.exists(args.Input_CPM_Table_FileName):
            F = open(args.Input_CPM_Table_FileName,'r')
            sample_names = F.readline().rstrip().split('\t')
            all_sample_names = list(sample_names)
            for L in F.readlines():
                C = L.rstrip().split('\t')
                GeneID = C[0]
                Genes_by_ID[GeneID] = TGene_Info(GeneID)
                CPMs = [to_float(x) for x in C[1:]]
                if len(CPMs) != len(sample_names):
                    print('%s: strings with sample names and CPMs do not match (cell counts %d vs %d)'%(args.Input_CPM_Table_FileName,
                        len(sample_names),len(CPMs)))
                    exit(127)

                for x in range(len(CPMs)):
                    Genes_by_ID[GeneID].CPM_by_SampleName[sample_names[x]] = CPMs[x]
        else:
            print('Warning! CPM table file %s is not found'%(args.Input_CPM_Table_FileName))
            args.insert_CPMs = False


    GLM_predictor_values_by_predictor_by_SampleName = dict()

    # P_CPM_combinations = []
    # P_CPM_combination_codes = []
    # GO_terms_list = []
    # min_abs_LogFC = 0.4
    # trunc_sheet_number = 1


    Any_RefSeq_Summary_is_present = False

    all_Models = []
    for FileName in args.Input_FileNames:
        #'Age_edgeR_combined.tsv'
        if 'edger' in FileName:
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
        Stats_cells_count = C.index('Pearson P') + 2
        CPM_cells_count = len(C) - Stats_cells_count + 1
        cell_index_by_name = {}
        for n in range(Stats_cells_count-1):
            cell_index_by_name[C[n]] = n + 1

        if any([(x not in cell_index_by_name) for x in ("Gene name", "Biotype", "description", "Score", "logFC", "logCPM", "LR", "PValue", "FDR", "Spearman r", "Spearman P", "Pearson r", "Pearson P")]):
            print('Incorrect file %s'%FileName)
            print('The following columns should present: "Gene name", "Biotype", "description", "RefSeq_Summary", "Score", "logFC", "logCPM", "LR", "PValue", "FDR", "Spearman r", "Spearman P", "Pearson r", "Pearson P"')
            print('Available only:')
            print(cell_index_by_name)
            continue

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
            if C[cell_index_by_name['logFC']] == '': ## predictors line
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

            Gene.Symbol = C[cell_index_by_name['Gene name']]
            Gene.Biotype = C[cell_index_by_name['Biotype']]
            Gene.Name = C[cell_index_by_name['description']]
            if Gene.RefSeq_summary == '' and RefSeq_Summary_is_present:
                Gene.RefSeq_summary = C[cell_index_by_name['RefSeq_Summary']]

            DE = TDE_Info()
            Gene.DE_Info_by_GLM[model] = DE
            DE.Score = C[cell_index_by_name['Score']]

            try:   DE.Score = float(C[cell_index_by_name['Score']])
            except ValueError: pass

            try:   DE.LogFC = float(C[cell_index_by_name['logFC']])
            except ValueError: pass

            try:   DE.LogCPM = float(C[cell_index_by_name['logCPM']])
            except ValueError: pass

            try:   DE.LR = float(C[cell_index_by_name['LR']])
            except ValueError: pass

            try:   DE.GLM_P = float(C[cell_index_by_name['PValue']])
            except ValueError: pass

            try:   DE.GLM_FDR = float(C[cell_index_by_name['FDR']])
            except ValueError: pass

            try:   DE.Spearman_r = float(C[cell_index_by_name['Spearman r']])
            except ValueError: pass

            try:   DE.Spearman_p = float(C[cell_index_by_name['Spearman P']])
            except ValueError: pass

            try:   DE.Pearson_r = float(C[cell_index_by_name['Pearson r']])
            except ValueError: pass

            try:   DE.Pearson_p = float(C[cell_index_by_name['Pearson P']])
            except ValueError: pass


            # if abs(DE.LogFC) > 0.3:  FormatLogFC_Cell(sheet,2**logFC,stringN,5 + RefSeq_Summary_is_present,coeff = 1.4)

            # for n in range(CPM_cells_count):
            #     try:  Gene.CPM_by_SampleName[] = float(C[14 + RefSeq_Summary_is_present + n])
            #     except ValueError:  Gene.CPM_by_SampleName[] = float('NaN')
            #
            # for n in range(CPM_cells_count):
            #     sheet.write_number(stringN,15 + RefSeq_Summary_is_present + n,float(C[14 + RefSeq_Summary_is_present + n]))

        Any_RefSeq_Summary_is_present = Any_RefSeq_Summary_is_present or RefSeq_Summary_is_present



    Gene_ID_sequence = [Genes_by_ID[x] for x in sorted(Genes_by_ID.keys())]
    for Gene in Gene_ID_sequence:  Gene.Calculate_Summary_Rating()
    Gene_ID_sequence.sort(key = (lambda x: x.Overall_score),reverse=True)


    Workbook = xlsxwriter.Workbook(args.Workbook_FileName)

    # GO_Sheet = Workbook.add_worksheet('GO-centric DE profiles')
    bold = Workbook.add_format({'bold': True, 'italic': False})
    italic = Workbook.add_format({'bold': False, 'italic': True})
    bold_italic = Workbook.add_format({'bold': True, 'italic': True})
    # bold_italic.set_text_wrap()
    format0 = Workbook.add_format()
    format0.set_num_format('0')
    format1 = Workbook.add_format()
    format1.set_num_format('0.0')
    format2 = Workbook.add_format()
    format2.set_num_format('0.00')
    format_center = Workbook.add_format({'align': 'center'})

    merge_format = Workbook.add_format({
        'bold': True,
        'italic': True,
        'border': True,
        'align': 'center',
        'valign': 'vcenter'})
    merge_format.set_text_wrap()

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
    Score_cols = []
    Rs_cols = []

    sheet = Workbook.add_worksheet('Joint data')

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
        sheet.write(1,4+args.insert_LogCPM_common,'RefSeq Summary',bold)
        sheet.set_column(4+args.insert_LogCPM_common,4+args.insert_LogCPM_common,12)

    ColN = 5+Any_RefSeq_Summary_is_present+args.insert_LogCPM_common

    sheet.set_column(ColN,ColN,2)
    GLM_info_columns_count = sum([getattr(args,x) for x in ['insert_LogFC', 'insert_LR', 'insert_LogCPM','insert_pvalue',
                                                            'insert_FDR', 'insert_Score','insert_spearman_r','insert_spearman_p',
                                                            'insert_pearson_r', 'insert_pearson_p']])



    for model in all_Models:
        ColN += 1
        if GLM_info_columns_count > 1:
            sheet.merge_range(0,ColN,0,ColN + GLM_info_columns_count-1,model,merge_format)
        else:
            sheet.write(0,ColN,model,merge_format)

        if args.insert_LogFC:   sheet.write(1,ColN,'LogFC',bold); LogFC_cols.append(ColN); ColN += 1
        if args.insert_LogCPM:   sheet.write(1,ColN,'LogCPM',bold);  LogCPM_cols.append(ColN); ColN += 1
        if args.insert_LR:   sheet.write(1,ColN,'LR',bold); LR_cols.append(ColN); ColN += 1
        if args.insert_pvalue:   sheet.write(1,ColN,'GLM p',bold);  P_cols.append(ColN); ColN += 1
        if args.insert_FDR:   sheet.write(1,ColN,'GLM FDR',bold);  FDR_cols.append(ColN); ColN += 1
        if args.insert_Score:   sheet.write(1,ColN,'Score',bold);  Score_cols.append(ColN); ColN += 1
        if args.insert_spearman_r:   sheet.write(1,ColN,'Spearman r',bold);  Rs_cols.append(ColN); ColN += 1
        if args.insert_spearman_p:   sheet.write(1,ColN,'Spearman p',bold);  P_cols.append(ColN); ColN += 1
        if args.insert_pearson_r:   sheet.write(1,ColN,'Pearson r',bold);  Rs_cols.append(ColN); ColN += 1
        if args.insert_pearson_p:   sheet.write(1,ColN,'Pearson p',bold);  P_cols.append(ColN); ColN += 1
        sheet.set_column(ColN,ColN,3)

    sheet.set_column(ColN+1,ColN+1,3)
    sheet.write(1,ColN,'CPMs',italic)
    ColN += 2
    CPMs_start_col = ColN
    all_sample_names.sort()

    # Max_LogFC_formatted_cells_per_model = int(float(args.max_formatted_cells)/float(len(all_Models)))

    all_predictors = []
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
    formatted_cells_count = 0
    for Gene in Gene_ID_sequence:
        if args.exclude_genes_without_de_info and len(Gene.DE_Info_by_GLM) == 0:
            continue

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
                        FormatLogFC_Cell(sheet,2**DE.LogFC,RowN,ColN,coeff = 1.4)
                        formatted_cells_count += 1
                ColN += 1

            if args.insert_LogCPM:
                if not math.isnan(args.insert_LogCPM):   sheet.write_number(RowN,ColN,DE.LogCPM,format1)
                else:   sheet.write(RowN,ColN,'nan')
                ColN += 1
            if args.insert_LR:
                if not math.isnan(DE.LR):   sheet.write_number(RowN,ColN,DE.LR,format1)
                else:   sheet.write(RowN,ColN,'nan')
                ColN += 1
            if args.insert_pvalue:
                if not math.isnan(DE.GLM_P):   sheet.write_number(RowN,ColN,DE.GLM_P)
                else:   sheet.write(RowN,ColN,'nan')
                ColN += 1
            if args.insert_FDR:
                if not math.isnan(DE.GLM_FDR):  sheet.write_number(RowN,ColN,DE.GLM_FDR)
                else:  sheet.write(RowN,ColN,'nan')
                ColN += 1
            if args.insert_Score:
                if not math.isnan(DE.Score):   sheet.write_number(RowN,ColN,DE.Score,format2)
                else:   sheet.write(RowN,ColN,'nan')
                ColN += 1
            if args.insert_spearman_r:
                if not math.isnan(DE.Spearman_r):   sheet.write_number(RowN,ColN,DE.Spearman_r,format2)
                else:   sheet.write(RowN,ColN,'nan')
                ColN += 1
            if args.insert_spearman_p:
                if not math.isnan(DE.Spearman_p):  sheet.write_number(RowN,ColN,DE.Spearman_p)
                else:  sheet.write(RowN,ColN,'nan')
                ColN += 1
            if args.insert_pearson_r:
                if not math.isnan(DE.Pearson_r):  sheet.write_number(RowN,ColN,DE.Pearson_r,format2)
                else:  sheet.write(RowN, ColN, 'nan')
                ColN += 1
            if args.insert_pearson_p:
                if not math.isnan(DE.Pearson_p):  sheet.write_number(RowN,ColN,DE.Pearson_p)
                else:  sheet.write(RowN,ColN,'nan')
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
        for ColN in LogFC_cols:
            sheet.conditional_format(2 + len(all_predictors), ColN, Last_Row, ColN,
                                     {'type': '3_color_scale',
                                                      'min_color': "#1170f1",'mid_color': "#ffffff",'max_color': "#ec4a18",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': -5.1, 'mid_value': 0.0, 'max_value': 5.1})

    for ColN in LogCPM_cols:
        sheet.conditional_format(2 + len(all_predictors),ColN,Last_Row,ColN, {'type': 'data_bar','bar_color': '#ffb910',
                                                                              'min_type':'num','max_type':'num',
                                                                              'min_value':0,'max_value':9})

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

    sheet.conditional_format(2 + len(all_predictors),2,Last_Row,2, {'type':'text',
                                       'criteria': 'containing', 'value':    'lincRNA', 'format':   format_lincRNA})
    sheet.conditional_format(2 + len(all_predictors),2,Last_Row,2, {'type':'text',
                                       'criteria': 'containing', 'value':    'antisense', 'format':   format_antisense})
    sheet.conditional_format(2 + len(all_predictors),2,Last_Row,2, {'type':'text',
                                       'criteria': 'containing', 'value':    'pseudogene', 'format':   format_pseudogene})

    Workbook.close()


