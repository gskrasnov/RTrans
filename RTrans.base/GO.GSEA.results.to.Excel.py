__author__ = 'George'
import sys,math
import xlsxwriter,os
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
    min_abs_LogFC = 0.4
    for FileName in sys.argv[2:]:
        'Age, GO GSEA for top-40 downreg genes'
        signature = FileName.replace('GO GSEA for top-','\t').replace(' gene','\t').split('\t')[-2]
        maxGenes,geneTypes = signature.split(' ')
        try:  maxGenes = int(maxGenes)
        except ValueError:
            print('Incorrect file name. Cannot detect type/number of top genes. Correct format: Age, GSEA for top-40 downreg genes.txt')
            continue

        if not os.path.exists(FileName):
            print('File %s does not exist'%FileName)
            continue

        sheet = Workbook.add_worksheet('top %d %s'%(maxGenes,geneTypes))
        sheet.write(0,0,'GO id',text_wrap)
        sheet.set_column(0,0,14)
        sheet.write(0,1,'Term',text_wrap)
        sheet.set_column(1,1,35)
        sheet.write(0,2,'Genes in genome (background)',text_wrap)
        sheet.write(0,3,'Genes in top-%d %s'%(maxGenes,geneTypes),text_wrap)
        sheet.write(0,4,'Expected genes in top-%d %s'%(maxGenes,geneTypes),text_wrap)
        sheet.write(0,5,'Rank in P.Fisher.elim',text_wrap)
        sheet.write(0,6,'P.Fisher',text_wrap)
        sheet.write(0,7,'P.Fisher.elim',text_wrap)
        sheet.write(0,8,'FDR.Fisher',text_wrap)
        sheet.write(0,9,'FDR.Fisher.elim',text_wrap)
        sheet.write(0,10,'Description',text_wrap)
        sheet.set_column(10,10,20)
        sheet.write(0,11,'Avg.LogFC for top-10 upreg',text_wrap)
        sheet.write(0,12,'Avg.LogFC for top-10 downreg',text_wrap)
        sheet.write(0,13,'Upreg genes (LogFC > 0.4)',text_wrap)
        sheet.write(0,14,'Downreg genes (LogFC < -0.4)',text_wrap)
        stringN = 1
        for L in open(FileName,'r').readlines()[1:]:
            L = L.replace('\r','').replace('\n','').replace('"','')
            C = L.split('\t')
            if len(C) >= 15:
                GO_id,Term,Genes_in_genome,Genes_in_top,Expected_genes_in_top,Rank_in_Fisher_elim,P_Fisher,P_Fisher_elim,FDR_Fisher,FDR_Fisher_elim,genes_simple,Description,genes_all,LogFCs_all = C[1:17]
            elif len(C) >= 13:
                GO_id,Term,Genes_in_genome,Genes_in_top,Expected_genes_in_top,Rank_in_Fisher_elim,P_Fisher,P_Fisher_elim,FDR_Fisher,FDR_Fisher_elim,genes_simple,Description,genes_all,LogFCs_all = C[1:15]
            else:
                print('File %s, string %d has incorrect format'%(FileName,stringN+1))
                continue

            Genes_in_genome = int(Genes_in_genome)
            Genes_in_top =int(Genes_in_top)
            Expected_genes_in_top = float(Expected_genes_in_top)
            Rank_in_Fisher_elim = int(Rank_in_Fisher_elim)
            P_Fisher = float(P_Fisher)
            P_Fisher_elim = float(P_Fisher_elim)
            FDR_Fisher = float(FDR_Fisher)
            FDR_Fisher_elim = float(FDR_Fisher_elim)

            genes_all = genes_all.split(', ')
            LogFCs_all = [float(x) for x in LogFCs_all.split(', ')]
            if len(LogFCs_all) != len(genes_all):
                print('File %s, string %d has incorrect format'%(FileName,stringN+1))
                continue
            pairs = [(genes_all[x],LogFCs_all[x]) for x in range(len(genes_all))]

            top_up_pairs = [x for x in pairs if x[1] > min_abs_LogFC]
            top_down_pairs = list(reversed([x for x in pairs if x[1] < (-1)*min_abs_LogFC]))

            sheet.write(stringN,0,GO_id,italic)
            sheet.write(stringN,1,Term)
            sheet.write_number(stringN,2,Genes_in_genome)
            sheet.write_number(stringN,3,Genes_in_top)
            sheet.write_number(stringN,4,Expected_genes_in_top)
            sheet.write_number(stringN,5,Rank_in_Fisher_elim)
            sheet.write_number(stringN,6,P_Fisher)
            sheet.write_number(stringN,7,P_Fisher_elim)
            sheet.write_number(stringN,8,FDR_Fisher)
            sheet.write_number(stringN,9,FDR_Fisher_elim)
            sheet.write(stringN,10,Description)
            avg = sum(LogFCs_all[:10])/(min(10,len(LogFCs_all)))
            sheet.write_number(stringN,11,avg)
            FormatLogFC_Cell(sheet,2**avg,stringN,11,coeff=2.5)
            avg = sum(LogFCs_all[-10:])/(min(10,len(LogFCs_all)))
            sheet.write_number(stringN,12,avg)
            FormatLogFC_Cell(sheet,2**avg,stringN,12,coeff=2.5)
            sheet.write(stringN,13,', '.join([x[0] for x in top_up_pairs]))
            sheet.write(stringN,14,', '.join([x[0] for x in top_down_pairs]))


            stringN += 1

        sheet.conditional_format(1,6,1+stringN,9, {'type': '3_color_scale',
                                                  'min_color': "#8cc031",'mid_color': "#ffe08d",'max_color': "#ffffff",
                                                  'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                  'min_value': 0, 'mid_value': 0.0002, 'max_value': 0.07})

    Workbook.close()
