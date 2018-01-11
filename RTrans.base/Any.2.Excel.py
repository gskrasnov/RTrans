__author__ = 'George'
import sys,math,os,math
import xlsxwriter
import argparse,glob
from xlsxwriter.utility import xl_range


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
    parser = argparse.ArgumentParser(description='Creating Excel workbooks from txt|tsv LogFC, logCMP tables')
    parser.add_argument('-i','--in', dest='Input_FileNames', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('-d','--dir', dest='Input_Dir', nargs='?',action='store', required=False,default=None, help='')
    parser.add_argument('-ext','--extensions', dest='file_extensions', nargs='?',action='store', required=False,default=['txt','tsv'], help='')
    parser.add_argument('-m','--max-depth', dest='max_depth', nargs='?',action='store', required=False,default='10000', help='')
    parser.add_argument('-r','--recursive', dest='Recursive', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('-o','--out-excel', dest='Workbook_FileName', nargs='?',action='store', required=False,default=None, help='')
    parser.add_argument('--one-book', dest='Generate_SingleBook', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--max-formatted-cells', dest='max_formatted_cells', nargs='?',action='store', required=False,default='35000', help='')
    args = parser.parse_args()

    for x in ['Generate_SingleBook','Recursive']:
        setattr(args,x, ToBool(getattr(args,x)))

    args.max_formatted_cells = int(args.max_formatted_cells)
    args.max_depth = int(args.max_depth)

    if (args.Input_Dir == None) and (args.Input_FileNames == None):
        print('Please specify input directory or input files')
        exit(127)


    if args.Input_FileNames != None:  all_src_FileNames = args.Input_FileNames
    else:  all_src_FileNames = []

    if args.Input_Dir != None:
        for ext in args.file_extensions:
            files = glob.glob('%s/**/*.%s'%(args.Input_Dir,ext), recursive=True)
            for f in files:
                depth = f.replace('\\','/').count('/') - args.Input_Dir.replace('\\','/').count('/')
                if depth <= args.max_depth:
                    all_src_FileNames += [f]

    print('Total %d src-files were found'%(len(all_src_FileNames)))

    file_N = 0
    print('\n')
    for FN in all_src_FileNames:
        file_N += 1
        sys.stdout.write('\rProcessing file %d of %d...'%(file_N,len(all_src_FileNames)))
        f = open(FN,'r')
        header_cells_count = f.readline().count('\t') + 1
        real_cells_count = f.readline().count('\t') + 1
        if header_cells_count == real_cells_count - 1:  row_names_is_present = True
        elif header_cells_count == real_cells_count:  row_names_is_present = False
        else:  print('Incorrect header and body cells count (%d and %d)'%(header_cells_count,real_cells_count));  exit(127)
        f.close()

        f = open(FN,'r')
        header_cells = f.readline().replace('\r','').replace('\n','').split('\t')
        if row_names_is_present:  header_cells = ['entry'] + header_cells

        header_cells_cf = [x.casefold().replace(',',' ').replace(';',' ').replace('.',' ') for x in header_cells]
        logFC_cols = [x for x in range(len(header_cells)) if any([y == 'logfc' for y in header_cells_cf[x].split(' ')])]
        # print(logFC_cols)
        # print(header_cells_cf)
        # for x in range(len(header_cells)):
        #     print(header_cells_cf[x].split(' '))
        logCPM_cols = [x for x in range(len(header_cells)) if any([y == 'logcpm' for y in header_cells_cf[x].split(' ')])]
        CPM_cols = [x for x in range(len(header_cells)) if any([y == 'cpm' for y in header_cells_cf[x].split(' ')]) and not any([y == 'logcpm' for y in header_cells_cf[x].split(' ')])]
        P_value_cols = [x for x in range(len(header_cells)) if any([y == 'p' or y == 'p.value' or y == 'p.adjust' or y == 'pvalue' or y == 'p_value' or y == 'qvalue' or y == 'q.value' or y == 'q_value' for y in header_cells_cf[x].split(' ')])]
        FDR_cols = [x for x in range(len(header_cells)) if any([y == 'fdr' for y in header_cells_cf[x].split(' ')])]
        Biotype_cols = [x for x in range(len(header_cells)) if 'biotype' in header_cells_cf[x]]
        Correlation_r_cols = [x for x in range(len(header_cells)) if
                              (any([y == 'spearman' for y in header_cells_cf[x].split(' ')]) or
                               any([y == 'pearson' for y in header_cells_cf[x].split(' ')]) or
                               any([('corr' in y) for y in header_cells_cf[x].split(' ')])) and
                              (any([y == 'r' for y in header_cells_cf[x].split(' ')]) or any([y == 'rs' for y in header_cells_cf[x].split(' ')]))]
        Score_cols = [x for x in range(len(header_cells)) if 'score' in header_cells_cf[x]]



        Workbook = xlsxwriter.Workbook(FN[:FN.rfind('.')]+ '.xlsx')

        # GO_Sheet = Workbook.add_worksheet('GO-centric DE profiles')
        bold = Workbook.add_format({'bold': True, 'italic': False})
        bold_ww = Workbook.add_format({'bold': True, 'italic': False})
        bold_ww.set_text_wrap()
        italic = Workbook.add_format({'bold': False, 'italic': True})
        bold_italic = Workbook.add_format({'bold': True, 'italic': True})
        bold_italic_ww = Workbook.add_format({'bold': True, 'italic': True})
        bold_italic_ww.set_text_wrap()
        format0 = Workbook.add_format()
        format0.set_num_format('0')
        format1 = Workbook.add_format()
        format1.set_num_format('0.0')
        format2 = Workbook.add_format()
        format2.set_num_format('0.00')
        format_center = Workbook.add_format({'align': 'center'})

        headers_format = Workbook.add_format({
            'bold': True,
            'italic': False,
            'border': True,
            'align': 'center',
            'valign': 'vcenter'})
        headers_format.set_text_wrap()

        merge_format = Workbook.add_format({
            'bold': True,
            'italic': True,
            'border': True,
            'align': 'center',
            'valign': 'vcenter'})
        merge_format.set_text_wrap()

        trunc_sheet_number = 1
        base = os.path.split(FN)[-1]
        sheet_name = base.replace('[','(').replace(']',')').replace('?','_').replace(':','(d)').replace('*','(m)').replace('\\','_').replace('/','_')
        if len(sheet_name) > 30:
            sheet_name = sheet_name[:25] + '...%d'%trunc_sheet_number
            trunc_sheet_number += 1
        sheet = Workbook.add_worksheet(sheet_name)

        for x in range(len(header_cells)):
            sheet.write(0,x,header_cells[x],headers_format)

        RowN = 1
        for L in f.readlines():
            C = L.rstrip().split('\t')
            for x in range(len(C)):
                if (x in logFC_cols) or (x in logCPM_cols) or (x in CPM_cols) or (x in Correlation_r_cols):
                    write_number__mod(sheet,RowN,x,C[x],format2)
                else:  write_number__mod(sheet,RowN,x,C[x])
            RowN += 1

        Last_Row = RowN

        for ColN in logFC_cols:
            sheet.conditional_format(1,ColN,Last_Row,ColN, {'type': '3_color_scale',
                                                      'min_color': "#1170f1",'mid_color': "#ffffff",'max_color': "#ec4a18",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': -5.1, 'mid_value': 0.0, 'max_value': 5.1})
        for ColN in Correlation_r_cols:
            sheet.conditional_format(1,ColN,Last_Row,ColN, {'type': '3_color_scale',
                                                      'min_color': "#1170f1",'mid_color': "#ffffff",'max_color': "#ec4a18",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': -0.98, 'mid_value': 0.0, 'max_value': 0.98})

        for ColN in logCPM_cols:
            sheet.conditional_format(1,ColN,Last_Row,ColN, {'type': 'data_bar','bar_color': '#ffb910',
                                                                                  'min_type':'num','max_type':'num',
                                                                                  'min_value':0,'max_value':9})

        for ColN in P_value_cols + FDR_cols:
            sheet.conditional_format(1,ColN,Last_Row,ColN, {'type': '3_color_scale',
                                                      'min_color': "#9ece49",'mid_color': "#ffe08d",'max_color': "#ffffff",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': 0, 'mid_value': 0.0002, 'max_value': 0.07})

        # Light red fill with dark red text.
        format_lincRNA = Workbook.add_format({'bg_color':   '#c3d594'})

        # Light yellow fill with dark yellow text.
        format_antisense = Workbook.add_format({'bg_color':   '#d9c188'})

        # Green fill with dark green text.
        format_pseudogene = Workbook.add_format({'bg_color':   '#a0b4c7'})

        for ColN in Biotype_cols:
            sheet.conditional_format(1,ColN,Last_Row,ColN, {'type':'text',
                                               'criteria': 'containing', 'value':    'lincRNA', 'format':   format_lincRNA})
            sheet.conditional_format(1,ColN,Last_Row,ColN, {'type':'text',
                                               'criteria': 'containing', 'value':    'antisense', 'format':   format_antisense})
            sheet.conditional_format(1,ColN,Last_Row,ColN, {'type':'text',
                                               'criteria': 'containing', 'value':    'pseudogene', 'format':   format_pseudogene})

        Workbook.close()
    print('\rCompleted                               ')


