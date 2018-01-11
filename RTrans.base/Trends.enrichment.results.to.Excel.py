__author__ = 'George'
import sys,math
import xlsxwriter,os
from xlsxwriter.utility import xl_range

###

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

    if len(sys.argv) < 4:
        print('Too few args. syntax: db out.xlsx file1.tsv file2.tsv')
        print(sys.argv)
        exit()


    DB_name = sys.argv[1]
    Workbook = xlsxwriter.Workbook(sys.argv[2])
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
    for FileName in sys.argv[3:]:
        # 'Age, KEGG GSEA for top-40 downreg genes'
        # signature = FileName.replace('KEGG trends enrichment for top-','\t').replace(' gene','\t').split('\t')[-2]
        # maxGenes,geneTypes = signature.split(' ')
        # try:  maxGenes = int(maxGenes)
        # except ValueError:
        #     print('Incorrect file name. Cannot detect type/number of top genes. Correct format: Age, GSEA for top-40 downreg genes.txt')
        #     continue
        #
        # if not os.path.exists(FileName):
        #     print('File %s does not exist'%FileName)
        #     continue

        sheet_up = Workbook.add_worksheet('%s trends enrich., UP'%DB_name)
        sheet_down = Workbook.add_worksheet('%s trends enrich., DOWN'%DB_name)
        for sheet in (sheet_up, sheet_down):
            sheet.write(0,0,'%s id'%DB_name,text_wrap)
            sheet.set_column(0,0,14)
            sheet.write(0,1,'Pathway Name',text_wrap)
            sheet.set_column(1,1,35)
            sheet.write(0, 2,'genes in pathway',text_wrap)
            sheet.write(0, 3,'Norm. Enrich. Score',text_wrap)
            sheet.write(0, 4,'p-value',text_wrap)
            sheet.write(0, 5,'FDR',text_wrap)
            sheet.write(0, 6,'q-value',text_wrap)
            sheet.write(0, 7,'LEA: tags, %',text_wrap)
            sheet.write(0, 8,'LEA: list, %',text_wrap)
            sheet.write(0, 9,'LEA: signal, %',text_wrap)
            sheet.write(0, 10,'the most differentially expressed genes')

        for (sheet, direction) in [(sheet_up, 'up'), (sheet_down, 'down')]:
            stringN = 1
            for L in open(FileName,'r').readlines()[1:]:
                L = L.replace('\r','').replace('\n','').replace('"','')
                C = L.split('\t')
                if len(C) == 12:
                    KEGG_id, Pathway_Name, setSize, enrichmentScore, NES, pvalue, p_adjusted, qvalue, rank, leading_edge, core_enrichment = C[1:]
                else:
                    print('[Trends.enrichment.results.to.Excel.py] File %s, string %d has incorrect format (1)'%(FileName,stringN+1));   continue

                try:
                    setSize = int(setSize)
                    enrichmentScore = float(enrichmentScore)
                    NES = float(NES)
                    pvalue = float(pvalue)
                    p_adjusted = float(p_adjusted)
                    qvalue = float(qvalue)
                except ValueError:
                    print('[Trends.enrichment.results.to.Excel.py] File %s, string %d has incorrect format (2)' % (FileName, stringN + 1))
                    continue

                if direction == 'up' and NES < 0: continue
                if direction == 'down' and NES > 0: continue

                try:
                    'tags=27%, list=20%, signal=22%'
                    leading_edge = leading_edge.replace(' ','').replace('tags=','').replace('list=','').replace('signal=','').replace('%','')
                    tags, list, signal = [int(x) for x in leading_edge.split(',')]
                except:
                    print('[Trends.enrichment.results.to.Excel.py] File %s, string %d has incorrect format (3)' % (FileName, stringN + 1))
                    continue

                core_enrichment = core_enrichment.split('/')

                sheet.write(stringN,0,KEGG_id,italic)
                sheet.write(stringN,1,Pathway_Name)
                sheet.write_number(stringN,2,setSize)
                sheet.write_number(stringN,3,NES,format2)
                # FormatLogFC_Cell(sheet, 2**NES, stringN, 3, coeff=1.3)
                sheet.write_number(stringN,4,pvalue)
                sheet.write_number(stringN,5,p_adjusted)
                sheet.write_number(stringN,6,qvalue)
                sheet.write_number(stringN,7,tags)
                sheet.write_number(stringN,8,list)
                sheet.write_number(stringN,9,signal)
                sheet.write(stringN,10,', '.join(core_enrichment))

                stringN += 1

            # enrichment score
            if direction == 'up':
                sheet.conditional_format(1, 3, 1 + stringN, 3, {'type': 'data_bar', 'bar_color': '#e26500',
                                                                 'min_type': 'num', 'max_type': 'num',
                                                                 'min_value': 0, 'max_value': 4.0})
            else:
                sheet.conditional_format(1, 3, 1 + stringN, 3, {'type': 'data_bar', 'bar_color': '#1a9fec',
                                                                'min_type': 'num', 'max_type': 'num',
                                                                'min_value': 0, 'max_value': -4.0})
            # p,q, FDR
            sheet.conditional_format(1, 4, 1+stringN, 6, {'type': '3_color_scale',
                                                      'min_color': "#8cc031",'mid_color': "#ffe08d",'max_color': "#ffffff",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': 0, 'mid_value': 0.001, 'max_value': 0.07})

            # gene count
            sheet.conditional_format(1, 2, 1 + stringN, 2, {'type': 'data_bar', 'bar_color': '#ffb910',
                                                             'min_type': 'num', 'max_type': 'percentile',
                                                             'min_value': 0, 'max_value': 97})

            # tags
            sheet.conditional_format(1, 7, 1 + stringN, 7, {'type': 'data_bar', 'bar_color': '#e6da84',
                                                             'min_type': 'num', 'max_type': 'num',
                                                             'min_value': 0, 'max_value': 100})
            # list
            sheet.conditional_format(1, 8, 1 + stringN, 8, {'type': 'data_bar', 'bar_color': '#b7e691',
                                                             'min_type': 'num', 'max_type': 'num',
                                                             'min_value': 0, 'max_value': 100})
            # signal
            sheet.conditional_format(1, 9, 1 + stringN, 9, {'type': 'data_bar', 'bar_color': '#aadeed',
                                                             'min_type': 'num', 'max_type': 'num',
                                                             'min_value': 0, 'max_value': 100})


    Workbook.close()
