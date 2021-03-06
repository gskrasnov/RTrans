__author__ = 'George'
import sys,math,os
import xlsxwriter
from xlsxwriter.utility import xl_range
import argparse

def write_number__mod(sheet, rowN, colN, text, format=None, bypass_if_NaN = False,
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


def ResizeArray(src_array,desired_length):
    initial_length = len(src_array)
    if initial_length == desired_length:
        return src_array

    res_array = [0.0]*desired_length  #numpy.zeros(desired_length,dtype=float)
    part_size =  initial_length / desired_length
    for part_number in range(desired_length):
        start_coord = part_number*part_size
        end_coord = start_coord + part_size
        values = []
        weights = []

        start_fragment_weight = int(min(start_coord + 1,end_coord)) - start_coord
        if start_fragment_weight > 0:
            values.append(src_array[int(start_coord)])
            weights.append(start_fragment_weight)

        end_fragment_weight = end_coord - int(max(end_coord,start_coord+1))
        if end_fragment_weight > 0:
            values.append(src_array[min(initial_length-1, int(end_coord))])
            weights.append(end_fragment_weight)


        if start_fragment_weight <= 0 and end_fragment_weight <= 0 :
            if int(start_coord) != int(end_coord):
                print('smth strange...')
            values.append(src_array[int(start_coord)])
            weights.append(1)

        for x in range(int(start_coord)+1,int(end_coord)):
            values.append(src_array[x])
            weights.append(1)

        res_array[part_number] = sum([values[x]*weights[x] for x in range(len(values))])/sum(weights)
    return res_array

class TGO_Term_LogFCs():
    def __init__(self,ID,Name,min_P_up,min_FDR_up,min_P_down,min_FDR_down,model,gene_max_P,gene_min_CPM,Genes_count):
        self.ID =ID
        self.ID_plus = ID
        self.Name = Name
        self.min_P_up = min_P_up
        self.min_FDR_up = min_FDR_up
        self.min_P_down = min_P_down
        self.min_FDR_down = min_FDR_down
        self.model = model
        self.gene_max_P = gene_max_P
        self.gene_min_CPM = gene_min_CPM
        self.Genes_count = Genes_count
        self.LogFCs = []

class TGO_Term_Info():
    def __init__(self,ID,Name):
        self.ID =ID
        self.ID_plus = ID
        self.Name = Name
        self.LogFCs_by_combination_code = {}

ScoreColorPoints = [(225,225,225),(0,0,0)]
ScoreColorPoints_up = [(225,225,225),(255,25,0)]
ScoreColorPoints_down = [(225,225,225),(0,89,255)]

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



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Creating Excel workbooks from RTrans GO-centric expression profiles')
    parser.add_argument('-i','--in', dest='Input_FileNames', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--models', dest='Input_ModelNames', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('-o','--out-excel', dest='Workbook_FileName', nargs='?',action='store', required=True,default=None, help='')
    parser.add_argument('-s','--insert-spaces', dest='Insert_spaces', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('-l','--sparklines-axis-limit', dest='SparklinesAxisLimits', nargs='?',action='store', required=False,default='2', help='')
    parser.add_argument('-max','--maximal-genes-count', dest='Max_Genes_Count', nargs='?',action='store', required=False,default='120', help='')
    parser.add_argument('-tol','--maximal-genes-count-tolerance-percents', dest='Max_Genes_Count_tolerance_percents', nargs='?',action='store', required=False,default='30', help='')
    parser.add_argument('-db','--database', dest='Database', nargs='?',action='store', required=False,default='GO', help='')
    args = parser.parse_args()


    # if len(sys.argv) < 4:
    #     print('Too few args. syntax: <insert.spaces?> out.xlsx file1.tsv file2.tsv')
    #     print(sys.argv)
    #     exit()

    GO_Term_Info_by_ID = {}

    args.Insert_spaces = args.Insert_spaces.casefold() in ['yes','y','on']
    if not args.Max_Genes_Count is None:  args.Max_Genes_Count = int(args.Max_Genes_Count)
    args.Max_Genes_Count_tolerance_percents = int(args.Max_Genes_Count_tolerance_percents)
    args.SparklinesAxisLimits = float(args.SparklinesAxisLimits)
    Workbook = xlsxwriter.Workbook(args.Workbook_FileName)

    GO_Sheet = Workbook.add_worksheet('%s-centric DE profiles'%args.Database)
    bold = Workbook.add_format({'bold': True, 'italic': False})
    bold.set_bg_color('#ffffff')
    italic = Workbook.add_format({'bold': False, 'italic': True})
    italic.set_bg_color('#ffffff')
    # format0 = Workbook.add_format()
    # format0.set_num_format('0')
    # format1 = Workbook.add_format()
    # format1.set_num_format('0.0')
    # format2 = Workbook.add_format()
    # format2.set_num_format('0.00')

    text_wrap = Workbook.add_format({'bold': True, 'italic': False})
    text_wrap.set_text_wrap()
    text_wrap.set_align('center')
    text_wrap.set_align('vcenter')
    text_wrap.set_bg_color('#ffffff')

    Casual_Format = Workbook.add_format()
    Casual_Format.set_bg_color('#ffffff')

    pvalue_formats = []
    for digits in range(5):
        pvalue_formats.append(Workbook.add_format())
        pvalue_formats[-1].set_num_format('0.%s' % ('0' * digits))


    if(len(args.Input_ModelNames) != len(args.Input_FileNames)):
        print('args.Input_FileNames and args.Input_ModelNames lengths should be equal')
        exit(127)


    P_CPM_combinations = []
    P_CPM_combination_codes = []
    GO_terms_list = []

    for f_N in range(len(args.Input_FileNames)):
        FileName = args.Input_FileNames[f_N]
        model = args.Input_ModelNames[f_N]
        #'Age, custom GO terms gene-centric info, ' '0' '0.01'
        tmp = os.path.split(FileName)[1].replace('.tsv','').replace(', P.lt.','\t').replace('logCPM.gt.','\t').replace(', GO-DE info','\t').replace(', KEGG-DE info','\t').replace(', Reactome-DE info','\t').split('\t')
        # print(tmp)
        ### model = tmp[0]
        gene_min_CPM = float(tmp[1])
        gene_max_P = float(tmp[2])
        P_CPM_combinations.append((model,gene_max_P,gene_min_CPM))
        P_CPM_combination_codes.append('%s/%g/%g'%(model,gene_max_P,gene_min_CPM))

        for L in open(FileName,'r').readlines()[1:]:
            L = L.replace('\r','').replace('\n','')
            C = L.split('\t')
            GO_ID = C[1]
            GO_name = C[2]

            try: min_P_up = float(C[3])
            except ValueError: min_P_up = 1

            try: min_FDR_up = float(C[4])
            except ValueError: min_FDR_up = 1

            try: min_P_down = float(C[5])
            except ValueError: min_P_down = 1

            try: min_FDR_down = float(C[6])
            except ValueError: min_FDR_down = 1

            Genes_count = int(C[7])
            if Genes_count == 0: continue

            if not GO_ID in GO_Term_Info_by_ID:
                GO_Term_Info_by_ID[GO_ID] = TGO_Term_Info(GO_ID,GO_name)
                GO_terms_list.append(GO_ID)
            gt = GO_Term_Info_by_ID[GO_ID]
            lgfcs = TGO_Term_LogFCs(GO_ID,GO_name,min_P_up,min_FDR_up,min_P_down,min_FDR_down,model,gene_max_P,gene_min_CPM,Genes_count)
            gt.LogFCs_by_combination_code[P_CPM_combination_codes[-1]] = lgfcs
            try:  lgfcs.LogFCs = [float(x) for x in C[8:]]
            except ValueError:
                print('Incorrect file %s'%FileName)
                print('String: ')
                print(C)
                exit()
            # if len(lgfcs.LogFCs) != Genes_count:
            #     print('(Error 1) Incorrect file %s'%FileName)
            #     exit()

            if (not args.Max_Genes_Count is None) and (len(lgfcs.LogFCs) > args.Max_Genes_Count*(100 + args.Max_Genes_Count_tolerance_percents)/100):
                lgfcs.LogFCs = ResizeArray(lgfcs.LogFCs,args.Max_Genes_Count)
                lgfcs.ID_plus = lgfcs.ID + ' [resized to %d]'%(args.Max_Genes_Count)

    order = list(range(len(P_CPM_combinations)))
    order.sort(key=(lambda x: P_CPM_combinations[x][1]), reverse=True)

    P_CPM_combinations = [P_CPM_combinations[x] for x in order]
    P_CPM_combination_codes = [P_CPM_combination_codes[x] for x in order]

    GO_Sheet.write(0,0,'%s ID'%args.Database,bold)
    GO_Sheet.set_column(0,0,14)
    GO_Sheet.write(0,1,'%s name'%args.Database,bold)
    GO_Sheet.set_column(1,1,25)
    GO_Sheet.write(0,2,'gene count',bold)
    GO_Sheet.set_column(2,2,10)
    GO_Sheet.write(0,3,'min p across tests (up)',text_wrap)
    GO_Sheet.set_column(3,3,10)
    GO_Sheet.write(0,4,'min FDR across tests (up)',text_wrap)
    GO_Sheet.set_column(4,4,10)
    GO_Sheet.write(0,5,'min p across tests (down)',text_wrap)
    GO_Sheet.set_column(5,5,10)
    GO_Sheet.write(0,6,'min FDR across tests (down)',text_wrap)
    GO_Sheet.set_column(6,6,10)
    GO_Sheet.write(0,7,'gen.sel.:',italic)
    for n_Comb in range(len(P_CPM_combinations)):
        if P_CPM_combinations[n_Comb][1] == 1:
            GO_Sheet.write(0, 8 + 2 * n_Comb, '~ %s; logCPM > %g' % (P_CPM_combinations[n_Comb][0], P_CPM_combinations[n_Comb][2]), text_wrap)
        else:
            GO_Sheet.write(0,8 + 2*n_Comb,'~ %s; p < %g, logCPM > %g'%P_CPM_combinations[n_Comb],text_wrap)

        GO_Sheet.set_column(8 + 2*n_Comb,8 + 2*n_Comb,11)
        GO_Sheet.set_column(8 + 2*n_Comb + 1,8 + 2*n_Comb + 1,2)

        # start_ColN = 8 + 2*len(P_CPM_combinations)
        # GO_Sheet.write(0,start_ColN + n_Comb,'p-val. for ~ %s; p < %g, CPM > %g'%P_CPM_combinations[n_Comb],text_wrap)
        # GO_Sheet.set_column(start_ColN + n_Comb,start_ColN + n_Comb,12)

    LogFC_Sheet_by_Combination_codes = {}
    for n_GO in range(len(GO_terms_list)):
        ID = GO_terms_list[n_GO]
        gt = GO_Term_Info_by_ID[ID]
        GO_Sheet.set_row(1 + n_GO*(1+args.Insert_spaces),18)
        if args.Insert_spaces:  GO_Sheet.set_row(2 + n_GO*(1+args.Insert_spaces),4)
        GO_Sheet.write(1 + n_GO*(1+args.Insert_spaces),0,gt.ID,italic)
        GO_Sheet.write(1 + n_GO*(1+args.Insert_spaces),1,gt.Name,Casual_Format)
        # GO_Sheet.write(1 + n_GO,2,'-'.min(gt.Genes_count)
        Genes_count_list = sorted(set([x.Genes_count for x in gt.LogFCs_by_combination_code.values()]))

        min_P_up_list = sorted(set([x.min_P_up for x in gt.LogFCs_by_combination_code.values()]))
        min_FDR_up_list = sorted(set([x.min_FDR_up for x in gt.LogFCs_by_combination_code.values()]))
        min_P_down_list = sorted(set([x.min_P_down for x in gt.LogFCs_by_combination_code.values()]))
        min_FDR_down_list = sorted(set([x.min_FDR_down for x in gt.LogFCs_by_combination_code.values()]))

        # if len(Genes_count_list) == 1:  GO_Sheet.write_number(1 + n_GO*(1+args.Insert_spaces),2,Genes_count_list[0],Casual_Format)
        # else:  GO_Sheet.write(1 + n_GO*(1+args.Insert_spaces),2,str(Genes_count_list[0]) + ' ... ' + str(Genes_count_list[-1]),Casual_Format)

        Genes_count_list_avg = int(sum(Genes_count_list)/len(Genes_count_list) + 0.5)
        GO_Sheet.write_number(1 + n_GO*(1+args.Insert_spaces),2,Genes_count_list_avg,Casual_Format)

        # if len(min_P_list) == 1:  GO_Sheet.write_number(1 + n_GO,3,min_P_list[0],Casual_Format)
        # else:  GO_Sheet.write(1 + n_GO*(1+args.Insert_spaces),3,str(min_P_list[0]) + ' ... ' + str(min_P_list[-1]),Casual_Format)
        # GO_Sheet.write_number(1 + n_GO*(1+args.Insert_spaces),3,min_P_up_list[0],Casual_Format)
        # GO_Sheet.write_number(1 + n_GO*(1+args.Insert_spaces),5,min_P_down_list[0],Casual_Format)
        write_number__mod(GO_Sheet, 1 + n_GO*(1+args.Insert_spaces), 3, min_P_up_list[0], pvalue_mode=True, pvalue_formats=pvalue_formats)
        write_number__mod(GO_Sheet, 1 + n_GO*(1+args.Insert_spaces), 5, min_P_down_list[0], pvalue_mode=True, pvalue_formats=pvalue_formats)


        # if len(min_FDR_list) == 1:  GO_Sheet.write_number(1 + n_GO*(1+args.Insert_spaces),4,min_FDR_list[0],Casual_Format)
        # else:  GO_Sheet.write(1 + n_GO*(1+args.Insert_spaces),4,str(min_FDR_list[0]) + ' ... ' + str(min_FDR_list[-1]),Casual_Format)
        # GO_Sheet.write_number(1 + n_GO*(1+args.Insert_spaces),4,min_FDR_up_list[0],Casual_Format)
        # GO_Sheet.write_number(1 + n_GO*(1+args.Insert_spaces),6,min_FDR_down_list[0],Casual_Format)
        write_number__mod(GO_Sheet, 1 + n_GO*(1+args.Insert_spaces), 4, min_FDR_up_list[0], pvalue_mode=True, pvalue_formats=pvalue_formats)
        write_number__mod(GO_Sheet, 1 + n_GO*(1+args.Insert_spaces), 6, min_FDR_down_list[0], pvalue_mode=True, pvalue_formats=pvalue_formats)

        for n_Comb in range(len(P_CPM_combinations)):
            Comb_code = P_CPM_combination_codes[n_Comb]
            if not Comb_code in gt.LogFCs_by_combination_code: continue
            lgfcs = gt.LogFCs_by_combination_code[Comb_code]

            spark_src_sheet_Name = 'src_sheet_%d'%(n_Comb+1)
            if not Comb_code in LogFC_Sheet_by_Combination_codes:
                spark_src_sheet = Workbook.add_worksheet(spark_src_sheet_Name)
                LogFC_Sheet_by_Combination_codes[Comb_code] = spark_src_sheet
                spark_src_sheet.write(0,0,'%s ID'%args.Database)
                spark_src_sheet.write(0,2,'Log FCs')
            else:  spark_src_sheet = LogFC_Sheet_by_Combination_codes[Comb_code]
            spark_src_sheet.write(1 + n_GO,0,lgfcs.ID_plus,italic)
            for n_LgFC in range(len(lgfcs.LogFCs)):
                spark_src_sheet.write_number(1 + n_GO,1 + n_LgFC,lgfcs.LogFCs[n_LgFC])

            GO_Sheet.add_sparkline(1 + n_GO*(1+args.Insert_spaces),8 + 2*n_Comb, {'range': '%s!%s'%(spark_src_sheet_Name,xl_range(1 + n_GO, 1, 1 + n_GO, len(lgfcs.LogFCs))),
                                            'type': 'column','max': args.SparklinesAxisLimits + 0.1, 'min': (-1)*args.SparklinesAxisLimits - 0.01,
                                            'negative_points':True,
                                            'series_color':'#ff6e2e',
                                            'negative_color':'#327fff'})


            GO_Sheet.write(1 + n_GO*(1+args.Insert_spaces),8 + 2*n_Comb + 1,' ',Casual_Format)
            if lgfcs.min_P_up > 0.05 and lgfcs.min_P_down > 0.05:
                GO_Sheet.write(1 + n_GO*(1+args.Insert_spaces),8 + 2*n_Comb,' ',Casual_Format)
                continue
            CurrentFormat = Workbook.add_format()
            CurrentFormat.set_align('center')
            CurrentFormat.set_border(5)   # 5 -bold, 6- double, 8 - dashed,

            ### scale: p = 0.0005....0.07 - from black to white
            if lgfcs.min_P_down > 0.05: ### case with GSEA for upreg genes
                C = color(Gradient(ScoreColorPoints_up, ((-1)*math.log10(lgfcs.min_P_up) - 1.3)/4))
                CurrentFormat.set_border_color(C)
                GO_Sheet.write(1 + n_GO*(1+args.Insert_spaces),8 + 2*n_Comb,' ',CurrentFormat)
            elif lgfcs.min_P_up > 0.05: ### case with GSEA for downreg genes
                C = color(Gradient(ScoreColorPoints_down, ((-1)*math.log10(lgfcs.min_P_down) - 1.3)/4))
                CurrentFormat.set_border_color(C)
                GO_Sheet.write(1 + n_GO*(1+args.Insert_spaces),8 + 2*n_Comb,' ',CurrentFormat)
            else:

                score_down = ((-1)*math.log10(lgfcs.min_P_down) - 1.3)
                score_up = ((-1)*math.log10(lgfcs.min_P_up) - 1.3)
                if score_up > score_down: pos = 0.5 + (score_up/(score_down) - 1)/8
                else: pos = 0.5 - (score_down/(score_up) - 1)/8
                up_down_ColorPoints = (ScoreColorPoints_down[1],(0,0,0),ScoreColorPoints_up[1])
                maxValue = Gradient(up_down_ColorPoints,pos)

                C = color(Gradient(((255,255,255),maxValue), max(score_down,score_up)/4))
                CurrentFormat.set_border_color(C)
                GO_Sheet.write(1 + n_GO*(1+args.Insert_spaces),8 + 2*n_Comb,' ',CurrentFormat)

    GO_Sheet.conditional_format(1, 3, 1 + len(GO_terms_list) * (1 + args.Insert_spaces), 6, {'type': '3_color_scale',
                                                                                             'min_color': "#8cc031",
                                                                                             'mid_color': "#ffe08d",
                                                                                             'max_color': "#ffffff",
                                                                                             'min_type': 'num',
                                                                                             'mid_type': 'num',
                                                                                             'max_type': 'num',
                                                                                             'min_value': 0,
                                                                                             'mid_value': 0.0002,
                                                                                             'max_value': 0.07})
    for n in range(len(P_CPM_combinations)*2 + 8):
        GO_Sheet.write(1 + len(GO_terms_list) * (1 + args.Insert_spaces), n, ' ', Casual_Format)


    Workbook.close()