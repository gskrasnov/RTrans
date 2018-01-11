__author__ = 'George'
import sys,math,os,math
import xlsxwriter
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


if __name__ == '__main__':

    if len(sys.argv) < 3:
        print('Too few args. syntax: out.xlsx file1.tsv file2.tsv')
        print(sys.argv)
        exit()


    Workbook = xlsxwriter.Workbook(sys.argv[1])
    # GO_Sheet = Workbook.add_worksheet('GO-centric DE profiles')
    bold = Workbook.add_format({'bold': True, 'italic': False})
    italic = Workbook.add_format({'bold': False, 'italic': True})
    format0 = Workbook.add_format()
    format0.set_num_format('0')
    format1 = Workbook.add_format()
    format1.set_num_format('0.0')
    format2 = Workbook.add_format()
    format2.set_num_format('0.00')
    format_center = Workbook.add_format({'align': 'center'})

    text_wrap = Workbook.add_format()
    text_wrap.set_text_wrap()
    text_wrap.set_align('center')
    text_wrap.set_align('vcenter')
    merge_format = Workbook.add_format({
        'bold': 1,
        'border': 1,
        'align': 'center',
        'valign': 'vcenter'})


    # P_CPM_combinations = []
    # P_CPM_combination_codes = []
    # GO_terms_list = []
    # min_abs_LogFC = 0.4
    trunc_sheet_number = 1
    for FileName in sys.argv[2:]:
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
            try:
                model = os.path.split(FileName.split('_external_data')[-2])[1]
            except ValueError:
                print('Incorrect file name. Cannot detect type/number of top genes. Correct format: Age, GSEA for top-40 downreg genes.txt')
                exit()
        else:
            print('Incorrect file name. Cannot detect type/number of top genes. Correct format: Age, GSEA for top-40 downreg genes.txt')
            exit()

        File = open(FileName,'r')
        L = File.readline()
        L = L.replace('\r','').replace('\n','').replace('"','')
        C = L.split('\t')
        if 'Pearson P' in C:
            Stats_cells_count = C.index('Pearson P') + 2
            CPM_cells_count = len(C) - Stats_cells_count + 1
        else:
            CPM_cells_count = 0
            Stats_cells_count = len(C) + 1

        cell_index_by_name = {}
        for n in range(Stats_cells_count-1):
            cell_index_by_name[C[n]] = n + 1

        mandatory_cols = ("Gene name", "Biotype", "description", "Score", "logFC", "PValue", "FDR")
        if any([(x not in cell_index_by_name) for x in mandatory_cols]):
            print('Incorrect file %s'%FileName)
            print('The following columns should present: %s'%(str(mandatory_cols)))
            print('Available only:')
            print(cell_index_by_name)
            continue

        sheet_name = model.replace('[','(').replace(']',')').replace('?','_').replace(':','(d)').replace('*','(m)').replace('\\','_').replace('/','_')
        if len(sheet_name) > 30:
            sheet_name = sheet_name[:25] + '...%d'%trunc_sheet_number
            trunc_sheet_number += 1
        sheet = Workbook.add_worksheet(sheet_name)
        Summary_is_present = 'RefSeq_Summary' in cell_index_by_name
        LogCPM_is_present = "logCPM" in cell_index_by_name
        LR_is_present = "LR" in cell_index_by_name
        Spearman_r_is_present = "Spearman r" in cell_index_by_name
        Spearman_P_is_present = "Spearman P" in cell_index_by_name
        Pearson_r_is_present = "Pearson r" in cell_index_by_name
        Pearson_P_is_present = "Pearson P" in cell_index_by_name

        excel_Col_index_by_Col_name = dict()

        ColN = 0
        sheet.write(0,ColN,'Gene ID',bold)
        excel_Col_index_by_Col_name['Gene ID'] = ColN
        sheet.set_column(ColN,ColN,19)
        ColN += 1

        sheet.write(0,ColN,'Symbol',bold)
        excel_Col_index_by_Col_name['Symbol'] = ColN
        ColN += 1

        sheet.write(0,ColN,'Biotype',bold)
        excel_Col_index_by_Col_name['Biotype'] = ColN
        sheet.set_column(ColN,ColN,18)
        ColN += 1

        sheet.write(0,ColN,'Name',bold)
        excel_Col_index_by_Col_name['Name'] = ColN
        sheet.set_column(ColN,ColN,35)
        ColN += 1

        if Summary_is_present:
            sheet.write(0,ColN,'Summary',bold)
            excel_Col_index_by_Col_name['Summary'] = ColN
            ColN += 1

        sheet.write(0, ColN, 'Score', bold)
        excel_Col_index_by_Col_name['Score'] = ColN
        ColN += 1

        sheet.write(0, ColN, 'LogFC', bold)
        excel_Col_index_by_Col_name['LogFC'] = ColN
        ColN += 1

        if LogCPM_is_present:
            sheet.write(0, ColN, 'LogCPM', bold)
            excel_Col_index_by_Col_name['LogCPM'] = ColN
            ColN += 1

        if LR_is_present:
            sheet.write(0, ColN, 'LR', bold)
            excel_Col_index_by_Col_name['LR'] = ColN
            ColN += 1

        sheet.write(0, ColN, 'GLM P', bold)
        excel_Col_index_by_Col_name['GLM P'] = ColN
        ColN += 1

        sheet.write(0, ColN, 'GLM FDR', bold)
        excel_Col_index_by_Col_name['GLM FDR'] = ColN
        ColN += 1

        if Spearman_r_is_present:
            sheet.write(0, ColN, 'Spear. r', bold)
            excel_Col_index_by_Col_name['Spear. r'] = ColN
            ColN += 1

        if Spearman_P_is_present:
            sheet.write(0, ColN, 'Spear. P', bold)
            excel_Col_index_by_Col_name['Spear. P'] = ColN
            ColN += 1

        if Pearson_r_is_present:
            sheet.write(0, ColN, 'Pearson r', bold)
            excel_Col_index_by_Col_name['Pearson r'] = ColN
            ColN += 1

        if Pearson_P_is_present:
            sheet.write(0, ColN, 'Pearson P', bold)
            excel_Col_index_by_Col_name['Pearson P'] = ColN
            ColN += 1

        if CPM_cells_count > 0:
            sheet.write(0, ColN, 'CPM:', bold)
            excel_Col_index_by_Col_name['CPM:'] = ColN
            ColN += 1
        CPM_start_ColN = ColN

        for n in range(CPM_cells_count):
            sheet.write(0,CPM_start_ColN + n,C[13 + Summary_is_present + n])

        stringN = 0
        for L in open(FileName,'r').readlines()[1:]:
            stringN += 1
            L = L.replace('\r','').replace('\n','').replace('"','')
            C = L.split('\t')
            if C[cell_index_by_name['logFC']] == '': ## predictors line
                sheet.write(stringN,0,C[0],italic)
                if CPM_cells_count > 0:
                    for n in range(Stats_cells_count,Stats_cells_count+CPM_cells_count):
                        try:  sheet.write_number(stringN,n+1,float(C[n]),italic)
                        except:  sheet.write(stringN,n+1,C[n],italic)
                continue

            ## casual line
            sheet.write(stringN,0,C[0],italic)
            sheet.write(stringN,excel_Col_index_by_Col_name['Symbol'],C[cell_index_by_name['Gene name']])
            sheet.write(stringN,excel_Col_index_by_Col_name['Biotype'],C[cell_index_by_name['Biotype']])
            sheet.write(stringN,excel_Col_index_by_Col_name['Name'],C[cell_index_by_name['description']])
            if Summary_is_present: sheet.write(stringN,excel_Col_index_by_Col_name['Summary'],C[cell_index_by_name['RefSeq_Summary']])

            # print(C[cell_index_by_name['Score']])
            write_number__mod(sheet,stringN,excel_Col_index_by_Col_name['Score'],C[cell_index_by_name['Score']])

            try:
                logFC = float(C[cell_index_by_name['logFC']])
                sheet.write_number(stringN,excel_Col_index_by_Col_name['LogFC'],logFC,format2)
                if abs(logFC) > 0.3:  FormatLogFC_Cell(sheet,2**logFC,stringN,excel_Col_index_by_Col_name['LogFC'],coeff = 1.4)
            except ValueError:  sheet.write(stringN,excel_Col_index_by_Col_name['LogFC'],C[cell_index_by_name['logFC']])

            if LogCPM_is_present:
                try:
                    logCPM = float(C[cell_index_by_name['logCPM']])
                    sheet.write_number(stringN,excel_Col_index_by_Col_name['LogCPM'],logCPM,format2)
                    # FormatLogCPM_Cell(sheet,logCPM,stringN,6 + Summary_is_present)
                except ValueError:  sheet.write(stringN,excel_Col_index_by_Col_name['LogCPM'],C[cell_index_by_name['logCPM']])

            if LR_is_present:
                write_number__mod(sheet,stringN,excel_Col_index_by_Col_name['LR'],C[cell_index_by_name['LR']],format2)

            write_number__mod(sheet,stringN,excel_Col_index_by_Col_name['GLM P'],C[cell_index_by_name['PValue']])

            write_number__mod(sheet,stringN,excel_Col_index_by_Col_name['GLM FDR'],C[cell_index_by_name['FDR']])

            if Spearman_r_is_present:
                try:  sheet.write_number(stringN,excel_Col_index_by_Col_name['Spear. r'],float(C[cell_index_by_name['Spearman r']]),format2)
                except ValueError:  sheet.write(stringN,excel_Col_index_by_Col_name['Spear. r'],C[cell_index_by_name['Spearman r']])
            if Spearman_P_is_present:
                try:  sheet.write_number(stringN,excel_Col_index_by_Col_name['Spear. P'],float(C[cell_index_by_name['Spearman P']]))
                except ValueError:  sheet.write(stringN,excel_Col_index_by_Col_name['Spear. P'],C[cell_index_by_name['Spearman P']])
            if Pearson_r_is_present:
                try:  sheet.write_number(stringN,excel_Col_index_by_Col_name['Pearson r'],float(C[cell_index_by_name['Pearson r']]),format2)
                except ValueError:  sheet.write(stringN,excel_Col_index_by_Col_name['Pearson r'],C[cell_index_by_name['Pearson r']])
            if Pearson_P_is_present:
                try:  sheet.write_number(stringN,excel_Col_index_by_Col_name['Pearson P'],float(C[cell_index_by_name['Pearson P']]))
                except ValueError:  sheet.write(stringN,excel_Col_index_by_Col_name['Pearson P'],C[cell_index_by_name['Pearson P']])

            for n in range(CPM_cells_count):
                sheet.write_number(stringN,CPM_start_ColN + n,float(C[14 + Summary_is_present + n]))

        for col_name in ('GLM P','GLM FDR','Spear. P','Pearson P'):
            if not col_name in excel_Col_index_by_Col_name: continue

            sheet.conditional_format(1,excel_Col_index_by_Col_name[col_name],1+stringN,excel_Col_index_by_Col_name[col_name], {'type': '3_color_scale',
                                                      'min_color': "#9ece49",'mid_color': "#ffe08d",'max_color': "#ffffff",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': 0, 'mid_value': 0.0002, 'max_value': 0.07})

        if LogCPM_is_present:
            sheet.conditional_format(1,excel_Col_index_by_Col_name['LogCPM'],1+stringN,excel_Col_index_by_Col_name['LogCPM'], {'type': 'data_bar','bar_color': '#ffb910',
                                                                                                 'min_type':'num','max_type':'num',
                                                                                                 'min_value':0,'max_value':9})
    Workbook.close()
