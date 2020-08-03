__author__ = 'George'
import sys,math,os,math,shlex
import xlsxwriter
import argparse,glob
import statistics
from xlsxwriter.utility import xl_range

from Excel_formats import Create_workbook_Formats, Create_workbook_Formats_whitespace

ScoreColorPoints = [(225,225,225),(0,0,0)]


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

TwoD_ColorPointsSet4 = [
    ['437cf0','ffffff', 'ee8926'],
    ['437cf0','ffffff', 'ee8926']
    ]

Color_gradient_schemas = [TwoD_ColorPointsSet1, TwoD_ColorPointsSet2, TwoD_ColorPointsSet3, TwoD_ColorPointsSet4]

Correlation_color_points = ['5f95e8', 'ffffff', 'fc6467']
LogFC_color_points = ['5ebeff', 'ffffff', 'ffa95e']

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

def PValue_to_str(pvalue):
    try:  pvalue = float(pvalue)
    except:  pvalue = float('NaN')

    if math.isnan(pvalue) or math.isinf(pvalue):
        return('NaN')

    if pvalue <= 0:  return '0'

    digits_count = int((-1)*math.log10(pvalue) - 0.001) + 1
    if digits_count <= 4:
        return eval('"%%.%df"%%pvalue'%digits_count)
    else:
        return '%.2g'%pvalue

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


def WriteFormatted_Rank_Cell(workbook, sheet, RowN, ColN, rank, upper_part = 17, lower_part = 17, extent = 1.25, default_format = None):

    try:  rank = float(rank)
    except:   return

    if math.isnan(rank):  return

    # upper_ColorPoints = [(255,255,255), (244,99,101)]
    # lower_ColorPoints = [(255,255,255), (83,140,210)]
    upper_ColorPoints = [(255,255,255), (244,108,68)]
    lower_ColorPoints = [(255,255,255), (86,126,234)]

    rank = max(0, min(100, rank))
    if (rank >= lower_part) and (rank <= 100 - upper_part):
        if default_format is None:
            CurrentFormat = workbook.add_format()
            CurrentFormat.set_align('center')
            CurrentFormat.set_num_format('0.0')
            CurrentFormat.set_bg_color('#ffffff')
            CurrentFormat.set_font_size(8)
        else:
            CurrentFormat = default_format

    elif rank < lower_part:
        BgC = color(Gradient(lower_ColorPoints,   ((lower_part - rank) / lower_part)**extent   ))
        CurrentFormat = workbook.add_format()
        CurrentFormat.set_align('center')
        CurrentFormat.set_num_format('0.0')
        CurrentFormat.set_bg_color(BgC)
        CurrentFormat.set_font_size(8)

    else:
        BgC = color(Gradient(upper_ColorPoints,   ((rank - 100 + upper_part) / upper_part)**extent   ))
        CurrentFormat = workbook.add_format()
        CurrentFormat.set_align('center')
        CurrentFormat.set_num_format('0.0')
        CurrentFormat.set_bg_color(BgC)
        CurrentFormat.set_font_size(8)

    sheet.write_number(RowN, ColN, rank, CurrentFormat)



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
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-ip','--in-pvalue-tables', dest='Pv_Input_FileNames', nargs='*',action='store', required=True,default=None, help='')
    parser.add_argument('-ir','--in-rs-tables', dest='Rs_Input_FileNames', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('-if','--in-fdr-tables', dest='FDR_Input_FileNames', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('-o','--out-excel', dest='Workbook_FileName', nargs='?',action='store', required=False,default=None, help='')
    parser.add_argument('--one-book', dest='Generate_SingleBook', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--insert-rs-sheets', dest='Insert_Rs_sheets', nargs='?',action='store', required=False,default='yes', help='')    
    parser.add_argument('--insert-pvalue-sheets', dest='Insert_PValue_sheets', nargs='?',action='store', required=False,default='yes', help='')    
    parser.add_argument('--insert-fdr-sheets', dest='Insert_FDR_sheets', nargs='?',action='store', required=False,default='yes', help='')    
    parser.add_argument('--rs-type', dest='rs_type', nargs='?',action='store', required=False, default='rs', help='')    

    parser.add_argument('-canno', '--cols-annotation-table', dest='cols_annotation_tables', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('-ranno', '--rows-annotation-table', dest='rows_annotation_tables', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--forced-cols', dest='Forced_ColTypes_by_ColNumber', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--forced-col-types-by-names', dest='Forced_ColTypes_by_ColName', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--forced-col-types-by-names-quote-symbol', dest='Forced_ColTypes_by_ColName__quote_symbol', nargs='?',action='store', required=False,default='&', help='')
    parser.add_argument('--forced-col-widths-by-names', dest='Forced_ColWidths_by_ColName', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--forced-col-widths-by-names-quote-symbol', dest='Forced_ColWidths_by_ColName__quote_symbol', nargs='?',action='store', required=False,default='&', help='')
    parser.add_argument('--forced-col-widths', dest='Forced_ColWidths', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--col-full-names', dest='ColFullNames', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--disable-coltype-autoassign', dest='Disable_ColType_autoassign', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--max-chars-in-cell', dest='max_char_in_cell', nargs='?',action='store', required=False, default='unlimited', help='')
    parser.add_argument('--simple-logfc-format', dest='Simple_LogFC_format', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--bidirectional-bars-in-logfc-cells', dest='Bidirectional_bars_in_LogFC_cells', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--add-cpm-cells-formatting', dest='Add_cpm_cells_formatting', nargs='?',action='store', required=False,default='yes', help='')

    parser.add_argument('--sheet-names', dest='Sheet_names', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--rs-col-width', dest='rs_col_width', nargs='?',action='store', required=False,default='5.5', help='')
    parser.add_argument('--rs-row-width', dest='rs_row_width', nargs='?',action='store', required=False,default='14.3', help='')
    parser.add_argument('--first-col-italic', dest='first_col_italic', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--first-col-width', dest='first_col_width', nargs='?',action='store', required=False,default='20', help='')
    parser.add_argument('--first-row-width', dest='first_row_width', nargs='?',action='store', required=False,default='100', help='')

    parser.add_argument('--min-font-size', dest='min_font_size', nargs='?',action='store', required=False,default='6.5', help='')
    parser.add_argument('--max-font-size', dest='max_font_size', nargs='?',action='store', required=False,default='10', help='')
    parser.add_argument('--min-pvalue', dest='min_pvalue', nargs='?',action='store', required=False,default='1e-5', help='')
    parser.add_argument('--max-pvalue', dest='max_pvalue', nargs='?',action='store', required=False,default='0.05', help='')
    parser.add_argument('--pvalue-threshold-for-font-color', dest='pvalue_threshold_for_font_color', nargs='?',action='store', required=False,default='0.05', help='')
    parser.add_argument('--font-color-lower-than-pvalue-threshold', dest='font_color_lower_than_pvalue_threshold', nargs='?',action='store', required=False,default='#aaaaaa', help='')
    parser.add_argument('--pvalue-scale-transform', dest='pvalue_scale_transform', nargs='?',action='store', required=False,default='(-math.log(x))**0.8', help='')

    parser.add_argument('--insert-spacer-between-cells', dest='Insert_spasers_between_cells', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--add-pvalues-to-corr-cells', dest='Add_PValues_to_corr_cells', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--add-fdr-to-corr-cells', dest='Add_FDRvalues_to_corr_cells', nargs='?',action='store', required=False,default='auto', help='')
    parser.add_argument('--whitespace', dest='create_White_background', nargs='?',action='store', required=False, default='no', help='')
    parser.add_argument('-v', '--verbose', dest='verbose', nargs='?',action='store', required=False,default='no', help='')
    # parser.add_argument('--exclude-cols', dest='exclude_Cols', nargs='*',action='store', required=False,default=None, help='')
    args = parser.parse_args()


    for x in ['Generate_SingleBook','first_col_italic', 'Insert_spasers_between_cells',
              'create_White_background', 'verbose', 'Insert_Rs_sheets', 'Insert_PValue_sheets', 'Insert_FDR_sheets', 'Add_PValues_to_corr_cells',
              'Simple_LogFC_format', 'Bidirectional_bars_in_LogFC_cells', 'Disable_ColType_autoassign', 'Add_cpm_cells_formatting']:
        #print(x)
        setattr(args,x, ToBool(getattr(args,x)))

    if args.Add_FDRvalues_to_corr_cells == 'auto':
        args.Add_FDRvalues_to_corr_cells = args.Add_PValues_to_corr_cells
    else:
        args.Add_FDRvalues_to_corr_cells = ToBool(args.Add_FDRvalues_to_corr_cells)

    if args.first_col_width != None:   args.first_col_width = float(args.first_col_width)
    if args.first_row_width != None:   args.first_row_width = float(args.first_row_width)
    if args.rs_col_width != None:   args.rs_col_width = float(args.rs_col_width)
    if args.rs_row_width != None:   args.rs_row_width = float(args.rs_row_width)
    
    if args.rs_type not in ('rs', 'logfc'):
        print('Unknown rs_type. SHould be either rs or logfc')
        exit(127)

    args.min_font_size = float(args.min_font_size)
    args.max_font_size = float(args.max_font_size)

    args.min_pvalue = float(args.min_pvalue)
    args.max_pvalue = float(args.max_pvalue)
    min_pvalue = min(args.min_pvalue, args.max_pvalue)
    max_pvalue = max(args.min_pvalue, args.max_pvalue)

    args.pvalue_threshold_for_font_color = float(args.pvalue_threshold_for_font_color)

    if args.create_White_background:
        Create_workbook_Formats_uni = Create_workbook_Formats_whitespace
    else:
        Create_workbook_Formats_uni = Create_workbook_Formats


    def Get_Trans_pv(x):
        return eval(args.pvalue_scale_transform)

    if args.Rs_Input_FileNames is not None:
        if len(args.Rs_Input_FileNames) != len(args.Pv_Input_FileNames):
            print('The number of input rs and Pvalue tables must be equal')
            exit(127)

    if args.FDR_Input_FileNames is not None:
        if len(args.Pv_Input_FileNames) != len(args.FDR_Input_FileNames):
            print('The number of input Pvalue and FDR tables must be equal')
            exit(127)
    
    if args.Sheet_names is not None:
        if len(args.Sheet_names) != len(args.Pv_Input_FileNames):
            print('Count of the forced sheet names (%d) does not correspond to the number of input files (%d)'%(len(args.Sheet_names), len(args.Pv_Input_FileNames)))
            exit(127)



    if args.Forced_ColTypes_by_ColName != None:
        args.Forced_ColTypes_by_ColName = [x.replace(args.Forced_ColTypes_by_ColName__quote_symbol, '"') for x in args.Forced_ColTypes_by_ColName]
    
    if args.Forced_ColWidths_by_ColName != None:
        args.Forced_ColWidths_by_ColName = [x.replace(args.Forced_ColWidths_by_ColName__quote_symbol, '"') for x in args.Forced_ColWidths_by_ColName]


    allowed_forced_col_types = ['cpm', 'cpm2', 'cpm_array', 'logfc_array', 'logfc_array2', 'logfc', 'logfc_narrow', 'logcpm', 'p', 'fdr', 'filter', 'score', 'spacer', 'corr', 'means', 'raw_means',
     'biotype', 'index', 'index1', 'index2', 'index1_err', 'index2_err', 'counts', 'gene_symbol', 'gene_name', 'rank_1', 'rank_2', 'taxon', 'color',
     'gradual_30', 'ffpm', 'percentage', 'delta_logfc', 'percents_updown', 'warm_gradient', 'cold_gradient', 'gray_gradient', 'hidden']

    Forced_ColNumbers_by_ColType_dict = dict()
    all_forced_cols = []
    if args.Forced_ColTypes_by_ColNumber != None:
        # Forced_ColNumbers_by_ColType_dict = dict()
        for expr in args.Forced_ColTypes_by_ColNumber:
            if not ':' in expr:
                print('[1] Incorrect args.Forced_ColTypes_by_ColNumber argument (%s). Should looks like CPM:1,2,3,5  LogFC:4,8,9   ...'%args.Forced_ColTypes_by_ColNumber)
                exit(127)
            expr = expr.casefold().split(':')
            if len(expr) != 2:
                print('[2] Incorrect args.Forced_ColTypes_by_ColNumber argument (%s). Should looks like CPM:1,2,3,5  LogFC:4,8,9   ...'%args.Forced_ColTypes_by_ColNumber)
                exit(127)

            try:
                cols = [int(x) for x in expr[1].split(',')]
            except ValueError:
                print('[3] Incorrect args.Forced_ColTypes_by_ColNumber argument (%s). Should looks like CPM:1,2,3,5  LogFC:4,8,9   ...'%args.Forced_ColTypes_by_ColNumber)
                exit(127)

            expr[0] = expr[0].casefold()
            if expr[0] in allowed_forced_col_types:
                Forced_ColNumbers_by_ColType_dict[expr[0]] = cols
                all_forced_cols += cols
            else:
                print('[4] Incorrect columns type %s in args.Forced_ColTypes_by_ColNumber argument'%expr[0])
                exit(127)

    Forced_ColTypes_by_ColNames = dict()
    Forced_ColNames_by_ColTypes_dict = dict()
    all_forced_col_names = []
    if args.Forced_ColTypes_by_ColName != None:
        # Forced_ColNumbers_by_ColType_dict = dict()
        for expr in args.Forced_ColTypes_by_ColName:
            if not ':' in expr:
                print('[5] Incorrect args.Forced_ColTypes_by_ColName argument (%s). Should looks like CPM:"sample1 CPM","sample2 CPM",.... '%args.Forced_ColTypes_by_ColName)
                exit(127)
            expr = expr.split(':')
            if len(expr) != 2:
                print('[6] Incorrect args.Forced_ColTypes_by_ColName argument (%s). Should looks like CPM:"sample1 CPM","sample2 CPM",.... '%args.Forced_ColTypes_by_ColName)
                exit(127)

            my_splitter = shlex.shlex(expr[1], posix=True)
            my_splitter.whitespace = ','
            my_splitter.whitespace_split = True
            try:
                col_names = list(my_splitter)
            except:
                print('[7] Incorrect args.Forced_ColTypes_by_ColName argument (%s). Should looks like CPM:"sample1 CPM","sample2 CPM",.... '%args.Forced_ColTypes_by_ColName)
                print('Incorrect field %s'%expr[1])
                exit(127)

            expr[0] = expr[0].casefold()
            if expr[0] in allowed_forced_col_types:
                Forced_ColNames_by_ColTypes_dict[expr[0]] = col_names
                for cn in col_names:
                    Forced_ColTypes_by_ColNames[cn] = expr[0]
                all_forced_col_names += col_names
            else:
                print('[8] Incorrect columns type %s in args.Forced_ColTypes_by_ColName argument'%expr[0])
                exit(127)

    Forced_ColWidths_by_ColNames = dict()

    if args.Forced_ColWidths_by_ColName != None:
        # Forced_ColNumbers_by_ColType_dict = dict()
        for expr in args.Forced_ColWidths_by_ColName:
            if not ':' in expr:
                print('[9] Incorrect args.Forced_ColWidths_by_ColName argument (%s). Should looks like 20:"LogFC","trimmed LogFC",.... '%args.Forced_ColWidths_by_ColName)
                exit(127)
            expr = expr.split(':')
            if len(expr) != 2:
                print('[10] Incorrect args.Forced_ColWidths_by_ColName argument (%s). Should looks like 20:"LogFC","trimmed LogFC",.... '%args.Forced_ColWidths_by_ColName)
                exit(127)

            my_splitter = shlex.shlex(expr[1], posix=True)
            my_splitter.whitespace = ','
            my_splitter.whitespace_split = True
            try:
                col_names = list(my_splitter)
            except:
                print('[11] Incorrect args.Forced_ColWidths_by_ColName argument (%s). Should looks like 20:"LogFC","trimmed LogFC",.... '%args.Forced_ColWidths_by_ColName)
                print('Incorrect field %s'%expr[1])
                exit(127)

            expr[0] = expr[0].casefold()
            try:  expr[0] = float(expr[0])
            except ValueError:
                print('[12] Incorrect args.Forced_ColWidths_by_ColName argument (%s). Should looks like 20:"LogFC","trimmed LogFC",.... '%args.Forced_ColWidths_by_ColName)
                exit(127)

            for cn in col_names:
                Forced_ColWidths_by_ColNames[cn] = expr[0]

    if args.Forced_ColWidths != None:
        for expr in args.Forced_ColWidths:
            if not ':' in expr:
                print('[24] Incorrect args.Forced_ColWidths argument (%s). Should looks like "TargetCPM:20" "ControlCPM:20"  ...'%args.Forced_ColWidths)
                exit(127)
            expr = expr.split(':')
            if len(expr) != 2:
                print('[25] Incorrect args.Forced_ColWidths argument (%s). Should looks like "TargetCPM:20" "ControlCPM:20"   ...'%args.Forced_ColWidths)
                exit(127)

            try:  expr[1] = float(expr[1])
            except ValueError:
                print('[26] Incorrect args.Forced_ColWidths argument (%s). Should looks like "TargetCPM:20" "ControlCPM:20"   .... '%args.Forced_ColWidths)
                exit(127)

            Forced_ColWidths_by_ColNames[expr[0]] = expr[1]


    ColFullNames_dict = dict()
    if args.ColFullNames != None:
        for expr in args.ColFullNames:
            if not ':' in expr:
                print('[22] Incorrect args.ColFullNames argument (%s). Should looks like "TargetCPM:Target CPM" "ControlCPM:Control CPM"  ...'%args.ColFullNames)
                exit(127)
            expr = expr.split(':')
            if len(expr) != 2:
                print('[23] Incorrect args.ColFullNames argument (%s). Should looks like "TargetCPM:Target CPM" "ControlCPM:Control CPM"   ...'%args.ColFullNames)
                exit(127)
            ColFullNames_dict[expr[0]] = expr[1]


    if args.Generate_SingleBook:
        if args.Workbook_FileName is None:
            print('Excel workbook filename should be specified when Generate_SingleBook==True')
            exit(127)

        Workbook = xlsxwriter.Workbook(args.Workbook_FileName)

        (bold, bold_light, bold_ww, italic, bold_italic, bold_italic_ww, format0, format1, format2, format1_narrow, format2_narrow, format1_center, format2_center,
            format2_center_bold, format_center, format0_center, format_perc, default_Rank_format, headers_format, small_headers_format,
            secondary_headers_format, merge_format, rotated_format, rotated_format_underline, rotated_format_center, whitespace_format,
            whitespace_format__right_border, whitespace_format__bottom_border,
            whitespace_format1, percentage_format0, whitespace_format__header, empty_logfc_format, sparklines_cell_format, sparklines_cell_format_no_bg,
            sparklines_cell_format_white_bg, pvalue_formats, format_quasi_invisible,
            format_quasi_invisible_left_border, format_quasi_invisible_right_border, format_quasi_invisible_left_right_borders,
            heatmap_cell_format_left_border, heatmap_cell_format_right_border, heatmap_cell_format_left_right_borders) = Create_workbook_Formats_uni(Workbook)

        # (bold, bold_light, bold_ww, italic, bold_italic, bold_italic_ww, format0, format1, format2, format1_narrow, format2_narrow, format1_center, format2_center,
        #     format2_center_bold, format_center, format0_center, format_perc, default_Rank_format, headers_format, small_headers_format,
        #     secondary_headers_format, merge_format, rotated_format, rotated_format_underline, rotated_format_center, whitespace_format,
        #     whitespace_format__right_border, whitespace_format__bottom_border,
        #     whitespace_format1, percentage_format0, whitespace_format__header, empty_logfc_format, sparklines_cell_format, sparklines_cell_format_no_bg,
        #     sparklines_cell_format_white_bg, pvalue_formats) = Create_workbook_Formats(Workbook)
        
        c_formats__by_font_size = dict()
        c_formats__light__by_font_size = dict()
        c_formats__bold__by_font_size = dict()

    if args.FDR_Input_FileNames is None:
        args.Insert_FDR_sheets = False

    if args.cols_annotation_tables is not None:
        if len(args.cols_annotation_tables) != len(args.Pv_Input_FileNames):
            print('The length of args.cols_annotation_tables should be equal to args.Pv_Input_FileNames')
            exit(127)

    if args.rows_annotation_tables is not None:
        if len(args.rows_annotation_tables) != len(args.Pv_Input_FileNames):
            print('The length of args.rows_annotation_tables should be equal to args.Pv_Input_FileNames')
            exit(127)

    if args.max_char_in_cell == 'unlimited':  args.max_char_in_cell = None
    if args.max_char_in_cell is not None:
        args.max_char_in_cell = int(args.max_char_in_cell)


    file_N = 0
    trunc_sheet_number = 1

    for file_N in range(len(args.Pv_Input_FileNames)):
        if args.verbose:  sys.stdout.write('\rProcessing file %d of %d...'%(file_N + 1, len(all_src_FileNames)))

        pv_lines = [x.rstrip('\r\n').split('\t') for x in open(args.Pv_Input_FileNames[file_N], 'r').readlines()]
        fdr_lines, rs_lines = None, None
        if args.Rs_Input_FileNames is not None:
            rs_lines = [x.rstrip('\r\n').split('\t') for x in open(args.Rs_Input_FileNames[file_N], 'r').readlines()]
        if args.FDR_Input_FileNames is not None:
            fdr_lines = [x.rstrip('\r\n').split('\t') for x in open(args.FDR_Input_FileNames[file_N], 'r').readlines()]

        if rs_lines is not None:
            if len(rs_lines) != len(pv_lines):
                print('The number of rs and pvalue lines (%d and %d) does not match for files %s and %s'%(
                    len(rs_lines), len(pv_lines), args.Rs_Input_FileNames[file_N], args.Pv_Input_FileNames[file_N]))
                exit(127)

        if fdr_lines is not None:
            if len(pv_lines) != len(fdr_lines):
                print('The number of pvalue and FDR lines (%d and %d) does not match for files %s and %s'%(
                    len(pv_lines), len(fdr_lines), args.Pv_Input_FileNames[file_N], args.FDR_Input_FileNames[file_N]))
                exit(127)

        if len(pv_lines) < 2:
            print('File %s: the number of lines is lower than 2. Proceeding to the next pair'%args.Pv_Input_FileNames[file_N])
            continue

        if len(pv_lines[0]) == len(pv_lines[1]) - 1:
            pv_lines[0] = [''] + pv_lines[0]
        elif len(pv_lines[0]) != len(pv_lines[1]):
            print('The number of header cells does not correspond to the number of body cells for file %s'%args.Pv_Input_FileNames[file_N])
            exit(127)

        if len(set([len(x) for x in pv_lines])) != 1:
            print('The numbers of cells among strings are not equal for file %s:'%(args.Pv_Input_FileNames[file_N], sorted(set([len(x) for x in pv_lines]))))
            exit(127)

        if rs_lines is not None:
            if len(rs_lines[0]) == len(rs_lines[1]) - 1:
                rs_lines[0] = [''] + rs_lines[0]
            elif len(rs_lines[0]) != len(rs_lines[1]):
                print('The number of header cells does not correspond to the number of body cells for file %s'%args.Rs_Input_FileNames[file_N])
                exit(127)

            if len(set([len(x) for x in rs_lines])) != 1:
                print('The numbers of cells among strings are not equal for file %s:'%(args.Rs_Input_FileNames[file_N], sorted(set([len(x) for x in rs_lines]))))
                exit(127)

        if fdr_lines is not None:
            if len(fdr_lines[0]) == len(fdr_lines[1]) - 1:
                fdr_lines[0] = [''] + fdr_lines[0]
            elif len(fdr_lines[0]) != len(fdr_lines[1]):
                print('The number of header cells does not correspond to the number of body cells for file %s'%args.FDR_Input_FileNames[file_N])
                exit(127)

            if len(set([len(x) for x in fdr_lines])) != 1:
                print('The numbers of cells among strings are not equal for file %s:'%(args.FDR_Input_FileNames[file_N], sorted(set([len(x) for x in fdr_lines]))))
                exit(127)

        header_cells = pv_lines[0]
        header_cells = [C_el[1:-1] if ((C_el.startswith('"') and C_el.endswith('"')) or (C_el.startswith("'") and C_el.endswith("'"))) else C_el for C_el in header_cells]

        
        if args.cols_annotation_tables is not None:
            cols_annotation_lines = [x.rstrip('\r\n').split('\t') for x in open(args.cols_annotation_tables[file_N], 'r').readlines()]
        else:
            cols_annotation_lines = []

        if args.rows_annotation_tables is not None:
            rows_annotation_lines = [x.rstrip('\r\n').split('\t') for x in open(args.rows_annotation_tables[file_N], 'r').readlines()]
        else:  rows_annotation_lines = []

        if len(cols_annotation_lines) > 1:
            cols_annos_count = len(cols_annotation_lines[1])
        else:  cols_annos_count = 0

        if len(rows_annotation_lines) > 1:
            rows_annos_count = len(rows_annotation_lines[1])
        else:  rows_annos_count = 0
        # print(cols_annotation_lines[1])
        # print(rows_annotation_lines[1])
        
        if not args.Generate_SingleBook:
            Workbook = xlsxwriter.Workbook(args.Workbook_FileName)
            c_formats__by_font_size = dict()
            c_formats__light__by_font_size = dict()
            c_formats__bold__by_font_size = dict()
            very_good_pvalue_format = Workbook.add_format(
                                    {'align': 'center', 'bold': False, 'valign': 'vcenter',
                                    'font_size': font_size, 'num_format': '0.00', 'font_color': '#000000'})


            (bold, bold_light, bold_ww, italic, bold_italic, bold_italic_ww, format0, format1, format2, format1_narrow, format2_narrow, format1_center, format2_center,
                format2_center_bold, format_center, format0_center, format_perc, default_Rank_format, headers_format, small_headers_format,
                secondary_headers_format, merge_format, rotated_format, rotated_format_underline, rotated_format_center, whitespace_format,
                whitespace_format__right_border, whitespace_format__bottom_border,
                whitespace_format1, percentage_format0, whitespace_format__header, empty_logfc_format, sparklines_cell_format, sparklines_cell_format_no_bg,
                sparklines_cell_format_white_bg, pvalue_formats, format_quasi_invisible,
                format_quasi_invisible_left_border, format_quasi_invisible_right_border, format_quasi_invisible_left_right_borders,
                heatmap_cell_format_left_border, heatmap_cell_format_right_border, heatmap_cell_format_left_right_borders) = Create_workbook_Formats_uni(Workbook)

            # (bold, bold_light, bold_ww, italic, bold_italic, bold_italic_ww, format0, format1, format2, format1_narrow, format2_narrow, format1_center, format2_center,
            #     format2_center_bold, format_center, format0_center, format_perc, default_Rank_format, headers_format, small_headers_format,
            #     secondary_headers_format, merge_format, rotated_format, rotated_format_underline, rotated_format_center, whitespace_format,
            #     whitespace_format__right_border, whitespace_format__bottom_border,
            #     whitespace_format1, percentage_format0, whitespace_format__header, empty_logfc_format, sparklines_cell_format, sparklines_cell_format_no_bg,
            #     sparklines_cell_format_white_bg, pvalue_formats) = Create_workbook_Formats(Workbook)

            
        format_empty = Workbook.add_format(
                        {'align': 'center', 'bold': False, 'valign': 'vcenter', 'text_wrap': True,
                        'font_size': args.min_font_size, 'num_format': '0.00',
                        'font_color': '#000000'})

        if args.Sheet_names is None:
            sheet_name_src = os.path.split(args.Pv_Input_FileNames[file_N])[-1]
        else:
            sheet_name_src = args.Sheet_names[file_N]

        sheet_name = sheet_name_src.replace('[', '(').replace(']', ')').replace('?', '_').replace(':', '(d)').replace('*','(m)').replace('\\', '_').replace('/', '_')


        if len(sheet_name) > 30:
            sheet_name = sheet_name[:23] + '...%d'%trunc_sheet_number
            if args.Generate_SingleBook:
                trunc_sheet_number += 1

        sheet_name_rs__FDR = None
        if args.rs_type == 'rs':
            if args.Add_PValues_to_corr_cells:
                sheet_name_rs = sheet_name + ', r+p'
            else:
                sheet_name_rs = sheet_name + ', r'
            if args.Add_FDRvalues_to_corr_cells:
                sheet_name_rs__FDR = sheet_name + ', r+FDR'

        elif args.rs_type == 'logfc':
            if args.Add_PValues_to_corr_cells:
                sheet_name_rs = sheet_name + ', LogFC+p'
            else:
                sheet_name_rs = sheet_name + ', LogFC'
            if args.Add_FDRvalues_to_corr_cells:
                sheet_name_rs__FDR = sheet_name + ', LogFC+FDR'

        else:
            raise('Incorrect args.rs_type')
        sheet_name_pv = sheet_name + ', p'
        sheet_name_fdr = sheet_name + ', FDR'

        sheet_name_rs = sheet_name_rs[:31]
        sheet_name_pv = sheet_name_pv[:31]
        sheet_name_fdr = sheet_name_fdr[:31]
        if sheet_name_rs__FDR is not None:
            sheet_name_rs__FDR = sheet_name_rs__FDR[:31]

        # if args.rs_type == 'rs' or args.Add_PValues_to_corr_cells:
        if args.Insert_Rs_sheets:  sheet_rs = Workbook.add_worksheet(sheet_name_rs)
        else:  sheet_rs = None
        
        if args.Insert_PValue_sheets:  sheet_pv = Workbook.add_worksheet(sheet_name_pv)
        else:  sheet_pv = None

        if sheet_name_rs__FDR is not None:   sheet_rs__FDR = Workbook.add_worksheet(sheet_name_rs__FDR)
        else:  sheet_rs__FDR = None

        # else:
        #     if args.Insert_PValue_sheets:  sheet_pv = Workbook.add_worksheet(sheet_name_pv)
        #     else:  sheet_pv = None
        #     if args.Insert_Rs_sheets:  sheet_rs = Workbook.add_worksheet(sheet_name_rs)
        #     else:  sheet_rs = None

        if args.Insert_FDR_sheets:  sheet_fdr = Workbook.add_worksheet(sheet_name_fdr)
        else:  sheet_fdr = None

        if args.rs_row_width is not None and sheet_rs is not None:
            for RowN in range(len(pv_lines) + 1):
                sheet_rs.set_row(cols_annos_count + RowN, args.rs_row_width, whitespace_format)

        if args.rs_row_width is not None and sheet_rs__FDR is not None:
            for RowN in range(len(pv_lines) + 1):
                sheet_rs__FDR.set_row(cols_annos_count + RowN, args.rs_row_width, whitespace_format)

        for sheet in [sheet_rs, sheet_rs__FDR, sheet_pv, sheet_fdr]:
            if sheet is None:  continue

            if args.rs_col_width is not None:
                sheet.set_column(rows_annos_count + 1, rows_annos_count + len(header_cells), args.rs_col_width)
            if args.first_col_width != None:  sheet.set_column(0, 0, args.first_col_width)
            if args.first_row_width != None:  sheet.set_row(0, args.first_row_width, whitespace_format__header)

            for x in range(len(header_cells)):
                sheet.write(0, rows_annos_count + x, header_cells[x], rotated_format_center)
            
        
        min_transformed_pvalue = min(Get_Trans_pv(min_pvalue), Get_Trans_pv(max_pvalue))
        max_transformed_pvalue = max(Get_Trans_pv(min_pvalue), Get_Trans_pv(max_pvalue))

        
        ###### writing main cor matrix

        for RowN in range(1, len(pv_lines)):
            rs_values, pv_values, fdr_values = None, None, None
            if rs_lines is not None:  rs_values = rs_lines[RowN]
            if pv_lines is not None:  pv_values = pv_lines[RowN]
            if fdr_lines is not None:  fdr_values = fdr_lines[RowN]

            for (sheet, values) in zip([sheet_rs, sheet_rs__FDR, sheet_pv, sheet_fdr], [rs_values, rs_values, pv_values, fdr_values]):
                if sheet is None:  continue
                if values is None:  continue
                if args.first_col_italic:  sheet.write(cols_annos_count + RowN, 0, values[0], italic)
                else:  sheet.write(cols_annos_count + RowN, 0, rs_values[0])
            
            for ColN in range(1, len(header_cells)):
                try:  rs = float(rs_values[ColN])
                except ValueError:  rs = float('NaN')

                try:  pv = float(pv_values[ColN])
                except ValueError:
                    if pv_values[ColN].startswith('<'):  pv = 0
                    else:  pv = float('NaN')

                try:  fdr = float(fdr_values[ColN])
                except ValueError:
                    if fdr_values[ColN].startswith('<'):  fdr = 0
                    else:  fdr = float('NaN')


                current_formats = []
                for c_pv in [pv, fdr]:
                    if not math.isnan(c_pv):
                        trans_pv = Get_Trans_pv(max(c_pv, min_pvalue))
                        ratio = min(1, max(0, (trans_pv - min_transformed_pvalue) / (max_transformed_pvalue - min_transformed_pvalue)))
                        font_size = args.min_font_size + ratio*(args.max_font_size - args.min_font_size)
                        font_size = int(font_size * 10 + 0.5)/10
                        
                        make_bold_font = (ratio == 1)
                        if math.isnan(rs):
                            bg_color = '#ffffff'
                        else:
                            if args.rs_type == 'rs':
                                bg_color = color(Gradient(Correlation_color_points, (rs + 1)/2))
                            elif args.rs_type == 'logfc':
                                bg_color = color(Gradient(LogFC_color_points, (rs + 4)/8))

                        if c_pv > args.pvalue_threshold_for_font_color:
                            if args.Add_PValues_to_corr_cells:
                                current_format = Workbook.add_format(
                                    {'align': 'center', 'bold': make_bold_font, 'valign': 'vcenter', 'text_wrap': True,
                                    'font_size': font_size, 'num_format': '0.00', 'bg_color': bg_color,
                                    'font_color': args.font_color_lower_than_pvalue_threshold})
                            else:
                                if font_size not in c_formats__light__by_font_size:
                                    c_formats__light__by_font_size[font_size] = Workbook.add_format(
                                        {'align': 'center', 'bold': False, 'valign': 'vcenter',
                                        'font_size': font_size, 'num_format': '0.00', 'font_color':args.font_color_lower_than_pvalue_threshold})
                                current_format = c_formats__light__by_font_size[font_size]
                        else:
                            if args.Add_PValues_to_corr_cells:
                                current_format = Workbook.add_format(
                                    {'align': 'center', 'bold': make_bold_font, 'valign': 'vcenter', 'text_wrap': True,
                                    'font_size': font_size, 'num_format': '0.00', 'bg_color': bg_color,
                                    'font_color': '#000000'})
                            else:
                                if make_bold_font:
                                    if font_size not in c_formats__bold__by_font_size:
                                        c_formats__bold__by_font_size[font_size] = Workbook.add_format(
                                            {'align': 'center', 'bold': True, 'valign': 'vcenter',
                                            'font_size': font_size, 'num_format': '0.00', 'font_color': '#000000'})
                                    current_format = c_formats__bold__by_font_size[font_size]                                    
                                else:
                                    if font_size not in c_formats__by_font_size:
                                        c_formats__by_font_size[font_size] = Workbook.add_format(
                                            {'align': 'center', 'bold': False, 'valign': 'vcenter',
                                            'font_size': font_size, 'num_format': '0.00', 'font_color': '#000000'})
                                    current_format = c_formats__by_font_size[font_size]
                    else:
                        current_format = format2
                    current_formats.append(current_format)


                if args.Add_PValues_to_corr_cells:
                    if sheet_rs is not None:
                        if not math.isnan(rs) and not math.isnan(pv):
                            sheet_rs.write(cols_annos_count + RowN, rows_annos_count + ColN, '%.2f\n%s'%(rs, PValue_to_str(pv)), current_formats[0])
                        elif math.isnan(rs) and not math.isnan(pv):
                            sheet_rs.write(cols_annos_count + RowN, rows_annos_count + ColN, '-\n%s'%(PValue_to_str(pv)), current_formats[0])
                        elif math.isnan(rs) and math.isnan(pv):
                            sheet_rs.write(cols_annos_count + RowN, rows_annos_count + ColN, '-', format_empty)

                    if sheet_rs__FDR is not None:
                        if not math.isnan(rs) and not math.isnan(fdr):
                            sheet_rs__FDR.write(cols_annos_count + RowN, rows_annos_count + ColN, '%.2f\n%s'%(rs, PValue_to_str(fdr)), current_formats[1])
                        elif math.isnan(rs) and not math.isnan(fdr):
                            sheet_rs__FDR.write(cols_annos_count + RowN, rows_annos_count + ColN, '-\n%s'%(PValue_to_str(fdr)), current_formats[1])
                        elif math.isnan(rs) and math.isnan(fdr):
                            sheet_rs__FDR.write(cols_annos_count + RowN, rows_annos_count + ColN, '-', format_empty)
                            
                else:
                    if sheet_rs is not None:
                        if math.isnan(rs) or math.isinf(rs):
                            sheet_rs.write(cols_annos_count + RowN, rows_annos_count + ColN, '-', format_empty)
                        else:
                            sheet_rs.write_number(cols_annos_count + RowN, rows_annos_count + ColN, rs, current_formats[0])

                if sheet_pv is not None:
                    if math.isnan(pv):
                        sheet_pv.write(cols_annos_count + RowN, rows_annos_count + ColN, '-', format_empty)
                    else:
                        write_number__mod(sheet_pv, cols_annos_count + RowN, rows_annos_count + ColN, pv, pvalue_mode=True, pvalue_formats = pvalue_formats)

                if sheet_fdr is not None:
                    if math.isnan(fdr):
                        sheet_fdr.write(cols_annos_count + RowN, rows_annos_count + ColN, '-', format_empty)
                    else:
                        write_number__mod(sheet_fdr, cols_annos_count + RowN, rows_annos_count + ColN, fdr, pvalue_mode=True, pvalue_formats = pvalue_formats)


        if not args.Add_PValues_to_corr_cells:
            for sheet_x in (sheet_rs, sheet_rs__FDR):
                if sheet_x is None:  continue
                if args.rs_type == 'rs':
                    sheet_x.conditional_format(1 + cols_annos_count, 1 + rows_annos_count, cols_annos_count + len(pv_lines), rows_annos_count + len(header_cells),
                       {'type': '3_color_scale',
                           'min_color': '#' + Correlation_color_points[0], 'mid_color': '#' + Correlation_color_points[1],
                           'max_color': '#' + Correlation_color_points[2],
                           'min_type': 'num', 'mid_type': 'num',
                           'max_type': 'num',
                           'min_value': -0.8, 'mid_value': 0,
                           'max_value': 0.8})

                elif args.rs_type == 'logfc':
                    sheet_x.conditional_format(1 + cols_annos_count, 1 + rows_annos_count, cols_annos_count + len(pv_lines), rows_annos_count + len(header_cells),
                       {'type': '3_color_scale',
                           'min_color': '#' + LogFC_color_points[0], 'mid_color': '#' + LogFC_color_points[1],
                           'max_color': '#' + LogFC_color_points[2],
                           'min_type': 'num', 'mid_type': 'num',
                           'max_type': 'num',
                           'min_value': -4.1, 'mid_value': 0,
                           'max_value': 4.1})
            
        if sheet_pv is not None:
            sheet_pv.conditional_format(1 + cols_annos_count, 1 + rows_annos_count, cols_annos_count + len(pv_lines), rows_annos_count + len(header_cells),
               {'type': '3_color_scale',
                'min_color': "#9ece49",'mid_color': "#ffe08d",'max_color': "#ffffff",
                'min_type': 'num','mid_type': 'num','max_type': 'num',
                'min_value': 0, 'mid_value': 0.0002, 'max_value': 0.07})
                
        if sheet_fdr is not None:
            sheet_fdr.conditional_format(1 + cols_annos_count, 1 + rows_annos_count, len(pv_lines) + cols_annos_count, len(header_cells) + rows_annos_count,
               {'type': '3_color_scale',
                'min_color': "#9ece49",'mid_color': "#ffe08d",'max_color': "#ffffff",
                'min_type': 'num','mid_type': 'num','max_type': 'num',
                'min_value': 0, 'mid_value': 0.0002, 'max_value': 0.07})

        for sheet_x in (sheet_rs, sheet_rs__FDR):
            if sheet_x is None:  continue

            sheet_x.write(cols_annos_count + len(pv_lines) + 1, rows_annos_count + 1, 'Font size and color depending on p-value (or FDR):', italic)
            c_format = Workbook.add_format(
                {'align': 'center', 'bold': False, 'valign': 'vcenter',
                'font_size': args.min_font_size, 'num_format': '0.00', 'font_color':args.font_color_lower_than_pvalue_threshold})
            sheet_x.write(cols_annos_count + len(pv_lines) + 2, rows_annos_count + 1, '>' + PValue_to_str(args.pvalue_threshold_for_font_color), c_format)

            parts_N = 7
            part_size = (math.log2(max_pvalue) - math.log2(min_pvalue)) / parts_N
            for pn in range(parts_N + 1):
                current_pvalue = 2**(math.log2(max_pvalue) - part_size * pn)
                trans_pv = Get_Trans_pv(max(current_pvalue, min_pvalue))
                # print(trans_pv)
                # print(current_pvalue)
                ratio = min(1, max(0, (trans_pv - min_transformed_pvalue) / (max_transformed_pvalue - min_transformed_pvalue)))
                font_size = args.min_font_size + ratio*(args.max_font_size - args.min_font_size)
                font_size = int(font_size * 10 + 0.5)/10
                
                c_format = Workbook.add_format(
                    {'align': 'center', 'bold': False, 'valign': 'vcenter',
                    'font_size': font_size, 'num_format': '0.00', 'font_color':'#000000'})
                sheet_x.write(cols_annos_count + len(pv_lines) + 2, rows_annos_count + 2 + pn, PValue_to_str(current_pvalue), c_format)

            c_format = Workbook.add_format(
                {'align': 'center', 'bold': True, 'valign': 'vcenter',
                'font_size': args.max_font_size, 'num_format': '0.00', 'font_color':'#000000'})
            sheet_x.write(cols_annos_count + len(pv_lines) + 2, rows_annos_count + 3 + parts_N, '<' + PValue_to_str(min_pvalue), c_format)

        
        ###### adding columns and rows annotations

        for mode, annotation_lines, annotation_tables__file_names in (['cols', cols_annotation_lines, args.cols_annotation_tables], ['rows', rows_annotation_lines, args.rows_annotation_tables]):
            if len(annotation_lines) == 0:  continue
            
            header_cells_count = len(annotation_lines[0])
            real_cells_count = len(annotation_lines[1])
            if header_cells_count == real_cells_count - 1:  row_names_is_present = True
            elif header_cells_count == real_cells_count:  row_names_is_present = False
            else:
                if real_cells_count == 1 or real_cells_count == 0:  continue
                print('file %s: Incorrect header and body cells count (%d and %d)'%(annotation_tables__file_names[file_N], header_cells_count, real_cells_count))
                exit(127)
            
            header_cells = annotation_lines[0]
            header_cells = [C_el[1:-1] if ((C_el.startswith('"') and C_el.endswith('"')) or (C_el.startswith("'") and C_el.endswith("'"))) else C_el for C_el in header_cells]

            Cols_count = len(header_cells)
            header_cells_cf = [x.casefold().replace(',',' ').replace(';',' ').replace('.',' ') for x in header_cells]

            if row_names_is_present:
                annotation_lines = [annotation_lines[0]] + [x[1:] for x in annotation_lines[1:]]

            ##### assigning column/row types

            if 'logfc' in Forced_ColNumbers_by_ColType_dict:  logFC_cols = Forced_ColNumbers_by_ColType_dict['logfc']
            elif 'logfc' in Forced_ColNames_by_ColTypes_dict:  logFC_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['logfc'])]
            elif not args.Disable_ColType_autoassign:  logFC_cols = [x for x in range(len(header_cells)) if any([y == 'logfc' for y in header_cells_cf[x].split(' ')]) and x not in all_forced_cols]
            else:  logFC_cols = []

            if 'logfc_narrow' in Forced_ColNumbers_by_ColType_dict:  logFC_narrow_cols = Forced_ColNumbers_by_ColType_dict['logfc_narrow']
            elif 'logfc_narrow' in Forced_ColNames_by_ColTypes_dict:  logFC_narrow_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['logfc_narrow'])]
            else:  logFC_narrow_cols = []        

            if 'logcpm' in Forced_ColNumbers_by_ColType_dict:   logCPM_cols = Forced_ColNumbers_by_ColType_dict['logcpm']
            elif 'logcpm' in Forced_ColNames_by_ColTypes_dict: logCPM_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['logcpm'])]
            elif not args.Disable_ColType_autoassign:    logCPM_cols = [x for x in range(len(header_cells)) if any([y == 'logcpm' for y in header_cells_cf[x].split(' ')]) and x not in all_forced_cols]
            else:  logCPM_cols = []

            if 'p' in Forced_ColNumbers_by_ColType_dict:   P_value_cols = Forced_ColNumbers_by_ColType_dict['p']
            elif 'p' in Forced_ColNames_by_ColTypes_dict: P_value_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['p'])]
            elif not args.Disable_ColType_autoassign:   P_value_cols = [x for x in range(len(header_cells)) if x not in all_forced_cols and any([y == 'p' or y == 'p.value' or y == 'p.adjust' or y == 'pvalue' or y == 'p-value' or y == 'p_value' or y == 'qvalue' or y == 'q.value' or y == 'q_value' for y in header_cells_cf[x].split(' ')])]
            else:  P_value_cols = []

            if 'fdr' in Forced_ColNumbers_by_ColType_dict:   FDR_cols = Forced_ColNumbers_by_ColType_dict['fdr']
            elif 'fdr' in Forced_ColNames_by_ColTypes_dict: FDR_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['fdr'])]
            elif not args.Disable_ColType_autoassign:    FDR_cols = [x for x in range(len(header_cells)) if x not in all_forced_cols and any([y == 'fdr' for y in header_cells_cf[x].split(' ')])]
            else:  FDR_cols = []

            if 'biotype' in Forced_ColNumbers_by_ColType_dict:   Biotype_cols = Forced_ColNumbers_by_ColType_dict['biotype']
            elif 'biotype' in Forced_ColNames_by_ColTypes_dict: Biotype_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['biotype'])]
            elif not args.Disable_ColType_autoassign:    Biotype_cols = [x for x in range(len(header_cells)) if x not in all_forced_cols and 'biotype' in header_cells_cf[x]]
            else:  Biotype_cols = []

            if 'index' in Forced_ColNumbers_by_ColType_dict:   Index_cols = Forced_ColNumbers_by_ColType_dict['index']
            elif 'index' in Forced_ColNames_by_ColTypes_dict: Index_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['index'])]
            else:  Index_cols = []

            if 'index1' in Forced_ColNumbers_by_ColType_dict:   Index1_cols = Forced_ColNumbers_by_ColType_dict['index1']
            elif 'index1' in Forced_ColNames_by_ColTypes_dict: Index1_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['index1'])]
            else:  Index1_cols = []

            if 'index2' in Forced_ColNumbers_by_ColType_dict:   Index2_cols = Forced_ColNumbers_by_ColType_dict['index2']
            elif 'index2' in Forced_ColNames_by_ColTypes_dict: Index2_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['index2'])]
            else:  Index2_cols = []

            if 'index1_err' in Forced_ColNumbers_by_ColType_dict:   Index1Err_cols = Forced_ColNumbers_by_ColType_dict['index1_err']
            elif 'index1_err' in Forced_ColNames_by_ColTypes_dict: Index1Err_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['index1_err'])]
            else:  Index1Err_cols = []

            if 'index2_err' in Forced_ColNumbers_by_ColType_dict:   Index2Err_cols = Forced_ColNumbers_by_ColType_dict['index2_err']
            elif 'index2_err' in Forced_ColNames_by_ColTypes_dict: Index2Err_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['index2_err'])]
            else:  Index2Err_cols = []

            if 'gradual_30' in Forced_ColNumbers_by_ColType_dict:   gradual_30_cols = Forced_ColNumbers_by_ColType_dict['gradual_30']
            elif 'gradual_30' in Forced_ColNames_by_ColTypes_dict: gradual_30_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['gradual_30'])]
            else:  gradual_30_cols = []

            if 'delta_logfc' in Forced_ColNumbers_by_ColType_dict:   delta_logfc_cols = Forced_ColNumbers_by_ColType_dict['delta_logfc']
            elif 'delta_logfc' in Forced_ColNames_by_ColTypes_dict: delta_logfc_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['delta_logfc'])]
            else:  delta_logfc_cols = []

            if 'percents_updown' in Forced_ColNumbers_by_ColType_dict:   percents_updown_cols = Forced_ColNumbers_by_ColType_dict['percents_updown']
            elif 'percents_updown' in Forced_ColNames_by_ColTypes_dict: percents_updown_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['percents_updown'])]
            else:  percents_updown_cols = []        

            if 'ffpm' in Forced_ColNumbers_by_ColType_dict:   FFPM_cols = Forced_ColNumbers_by_ColType_dict['ffpm']
            elif 'ffpm' in Forced_ColNames_by_ColTypes_dict: FFPM_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['ffpm'])]
            else:  FFPM_cols = []
            
            if 'counts' in Forced_ColNumbers_by_ColType_dict:   Counts_cols = Forced_ColNumbers_by_ColType_dict['counts']
            elif 'counts' in Forced_ColNames_by_ColTypes_dict: Counts_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['counts'])]
            else:  Counts_cols = []

            if 'gene_symbol' in Forced_ColNumbers_by_ColType_dict:   Gene_Symbol_cols = Forced_ColNumbers_by_ColType_dict['gene_symbol']
            elif 'gene_symbol' in Forced_ColNames_by_ColTypes_dict: Gene_Symbol_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['gene_symbol'])]
            else:  Gene_Symbol_cols = []

            if 'gene_name' in Forced_ColNumbers_by_ColType_dict:   Gene_Name_cols = Forced_ColNumbers_by_ColType_dict['gene_name']
            elif 'gene_name' in Forced_ColNames_by_ColTypes_dict: Gene_Name_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['gene_name'])]
            else:  Gene_Name_cols = []

            if 'rank_1' in Forced_ColNumbers_by_ColType_dict:   Rank_1_cols = Forced_ColNumbers_by_ColType_dict['rank_1']
            elif 'rank_1' in Forced_ColNames_by_ColTypes_dict: Rank_1_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['rank_1'])]
            else:  Rank_1_cols = []
            
            if 'rank_2' in Forced_ColNumbers_by_ColType_dict:   Rank_2_cols = Forced_ColNumbers_by_ColType_dict['rank_2']
            elif 'rank_2' in Forced_ColNames_by_ColTypes_dict: Rank_2_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['rank_2'])]
            else:  Rank_2_cols = []

            if 'percentage' in Forced_ColNumbers_by_ColType_dict:   percentage_cols = Forced_ColNumbers_by_ColType_dict['percentage']
            elif 'percentage' in Forced_ColNames_by_ColTypes_dict: percentage_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['percentage'])]
            else:  percentage_cols = []

            if 'taxon' in Forced_ColNumbers_by_ColType_dict:   Taxon_cols = Forced_ColNumbers_by_ColType_dict['taxon']
            elif 'taxon' in Forced_ColNames_by_ColTypes_dict: Taxon_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['taxon'])]
            else:  Taxon_cols = []

            if 'color' in Forced_ColNumbers_by_ColType_dict:   Color_cols = Forced_ColNumbers_by_ColType_dict['color']
            elif 'color' in Forced_ColNames_by_ColTypes_dict: Color_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['color'])]
            else:  Color_cols = []

            if 'warm_gradient' in Forced_ColNumbers_by_ColType_dict:   warm_gradient_cols = Forced_ColNumbers_by_ColType_dict['warm_gradient']
            elif 'warm_gradient' in Forced_ColNames_by_ColTypes_dict: warm_gradient_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['warm_gradient'])]
            else:  warm_gradient_cols = []

            if 'cold_gradient' in Forced_ColNumbers_by_ColType_dict:   cold_gradient_cols = Forced_ColNumbers_by_ColType_dict['cold_gradient']
            elif 'cold_gradient' in Forced_ColNames_by_ColTypes_dict: cold_gradient_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['cold_gradient'])]
            else:  cold_gradient_cols = []

            if 'gray_gradient' in Forced_ColNumbers_by_ColType_dict:   gray_gradient_cols = Forced_ColNumbers_by_ColType_dict['gray_gradient']
            elif 'gray_gradient' in Forced_ColNames_by_ColTypes_dict: gray_gradient_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['gray_gradient'])]
            else:  gray_gradient_cols = []


            if 'hidden' in Forced_ColNumbers_by_ColType_dict:   Hidden_cols = Forced_ColNumbers_by_ColType_dict['hidden']
            elif 'hidden' in Forced_ColNames_by_ColTypes_dict: Hidden_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['hidden'])]
            else:  Hidden_cols = []


            if 'corr' in Forced_ColNumbers_by_ColType_dict:   Correlation_r_cols = Forced_ColNumbers_by_ColType_dict['corr']
            elif 'corr' in Forced_ColNames_by_ColTypes_dict: Correlation_r_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['corr'])]
            elif not args.Disable_ColType_autoassign:     Correlation_r_cols = [x for x in range(len(header_cells)) if
                                  x not in all_forced_cols and
                                  (any([y == 'spearman' for y in header_cells_cf[x].split(' ')]) or
                                   any([y == 'pearson' for y in header_cells_cf[x].split(' ')]) or
                                   any([('corr' in y) for y in header_cells_cf[x].split(' ')])) and
                                  (any([y == 'r' for y in header_cells_cf[x].split(' ')]) or any([y == 'rs' for y in header_cells_cf[x].split(' ')]))]
            else:  Correlation_r_cols = []

            if 'score' in Forced_ColNumbers_by_ColType_dict:  Score_cols = Forced_ColNumbers_by_ColType_dict['score']
            elif 'score' in Forced_ColNames_by_ColTypes_dict: Score_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['score'])]
            elif not args.Disable_ColType_autoassign:   Score_cols = [x for x in range(len(header_cells)) if x not in all_forced_cols and 'score' in header_cells_cf[x]]
            else:  Score_cols = []

            if 'cpm2' in Forced_ColNumbers_by_ColType_dict:    CPM_cols2 = Forced_ColNumbers_by_ColType_dict['cpm2']
            elif 'cpm2' in Forced_ColNames_by_ColTypes_dict:  CPM_cols2 = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['cpm2'])]
            else:  CPM_cols2 = []

            if 'logfc_array' in Forced_ColNumbers_by_ColType_dict:    LogFC_array_cols = Forced_ColNumbers_by_ColType_dict['logfc_array']
            elif 'logfc_array' in Forced_ColNames_by_ColTypes_dict:  LogFC_array_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['logfc_array'])]
            else:  LogFC_array_cols = []

            if 'logfc_array2' in Forced_ColNumbers_by_ColType_dict:    LogFC_array2_cols = Forced_ColNumbers_by_ColType_dict['logfc_array2']
            elif 'logfc_array2' in Forced_ColNames_by_ColTypes_dict:  LogFC_array2_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['logfc_array2'])]
            else:  LogFC_array2_cols = []


            # if args.Metagenome_mode:
            #     if 'means' in Forced_ColNumbers_by_ColType_dict:   means_CPM_cols = Forced_ColNumbers_by_ColType_dict['means']
            #     elif 'means' in Forced_ColNames_by_ColTypes_dict: means_CPM_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['means'])]
            #     elif not args.Disable_ColType_autoassign:   means_CPM_cols = [x for x in range(len(header_cells)) if x not in all_forced_cols and
            #                 ('group'in header_cells_cf[x].split(' ') and 'mean'in header_cells_cf[x].split(' ') and 'raw mean' not in header_cells_cf[x] ) and not any(
            #                     [y == 'logcpm' for y in header_cells_cf[x].split(' ')])]
            #     else:  means_CPM_cols = []

            #     if 'raw_means' in Forced_ColNumbers_by_ColType_dict:    raw_means_CPM_cols = Forced_ColNumbers_by_ColType_dict['raw_means']
            #     elif 'raw_means' in Forced_ColNames_by_ColTypes_dict: raw_means_CPM_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['raw_means'])]
            #     elif not args.Disable_ColType_autoassign:   raw_means_CPM_cols = [x for x in range(len(header_cells)) if x not in all_forced_cols and
            #                 ('group'in header_cells_cf[x].split(' ') and 'raw mean'in header_cells_cf[x]) and not any(
            #                     [y == 'logcpm' for y in header_cells_cf[x].split(' ')])]
            #     else:  raw_means_CPM_cols = []

            #     if len(means_CPM_cols) != len(raw_means_CPM_cols) and len(raw_means_CPM_cols) > 0:
            #         print('File %s: the number of "group means" and "group raw means" columns is not equal: %d and %d'%(FN, len(means_CPM_cols), len(raw_means_CPM_cols)))
            #         exit(131)

            #     last_stat_col = max([0] + Correlation_r_cols + FDR_cols + P_value_cols + logCPM_cols + Score_cols)

            #     if 'cpm_array' in Forced_ColNumbers_by_ColType_dict:  CPM_array_cols = Forced_ColNumbers_by_ColType_dict['cpm_array']
            #     elif 'cpm_array' in Forced_ColNames_by_ColTypes_dict: CPM_array_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['cpm_array'])]
            #     elif not args.Disable_ColType_autoassign:  CPM_array_cols = list(range(last_stat_col + 1, len(header_cells)))
            #     else:  CPM_array_cols = []

            #     if 'cpm' in Forced_ColNumbers_by_ColType_dict:    CPM_cols = Forced_ColNumbers_by_ColType_dict['cpm']
            #     elif 'cpm' in Forced_ColNames_by_ColTypes_dict: CPM_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['cpm'])]
            #     else:   CPM_cols = means_CPM_cols + CPM_array_cols

            # else:
            if 'cpm' in Forced_ColNumbers_by_ColType_dict:    CPM_cols = Forced_ColNumbers_by_ColType_dict['cpm']
            elif 'cpm' in Forced_ColNames_by_ColTypes_dict: CPM_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['cpm'])]
            elif not args.Disable_ColType_autoassign:   CPM_cols = [x for x in range(len(header_cells)) if x not in all_forced_cols and any([y == 'cpm' for y in header_cells_cf[x].split(' ')]) and not any([y == 'logcpm' for y in header_cells_cf[x].split(' ')])]
            else:  CPM_cols = []

            if 'means' in Forced_ColNumbers_by_ColType_dict:   means_CPM_cols = Forced_ColNumbers_by_ColType_dict['means']
            elif 'means' in Forced_ColNames_by_ColTypes_dict: means_CPM_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['means'])]
            else: means_CPM_cols = []

            if 'raw_means' in Forced_ColNumbers_by_ColType_dict:    raw_means_CPM_cols = Forced_ColNumbers_by_ColType_dict['raw_means']
            elif 'raw_means' in Forced_ColNames_by_ColTypes_dict: raw_means_CPM_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['raw_means'])]
            else:   raw_means_CPM_cols = []

            if 'cpm_array' in Forced_ColNumbers_by_ColType_dict:  CPM_array_cols = Forced_ColNumbers_by_ColType_dict['cpm_array']
            elif 'cpm_array' in Forced_ColNames_by_ColTypes_dict: CPM_array_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['cpm_array'])]
            else:   CPM_array_cols = CPM_cols

            #######################

            if 'cpm' not in Forced_ColNumbers_by_ColType_dict and 'cpm' not in Forced_ColNames_by_ColTypes_dict and not args.Disable_ColType_autoassign:
                CPM_cols = [x for x in CPM_cols if x not in all_forced_cols and not header_cells_cf[x].startswith('(spacer ')]
                
            if 'cpm2' not in Forced_ColNumbers_by_ColType_dict and 'cpm2' not in Forced_ColNames_by_ColTypes_dict and not args.Disable_ColType_autoassign:
                CPM_cols2 = [x for x in CPM_cols2 if x not in all_forced_cols and not header_cells_cf[x].startswith('(spacer ')]
                
            if 'cpm_array' not in Forced_ColNumbers_by_ColType_dict and 'cpm_array' not in Forced_ColNames_by_ColTypes_dict and not args.Disable_ColType_autoassign:
                CPM_array_cols = [x for x in CPM_array_cols if x not in all_forced_cols and not header_cells_cf[x].startswith('(spacer ')]

            if 'filter' in Forced_ColNumbers_by_ColType_dict:  Filter_cols = Forced_ColNumbers_by_ColType_dict['filter']
            elif 'filter' in Forced_ColNames_by_ColTypes_dict: Filter_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['filter'])]
            elif not args.Disable_ColType_autoassign:   Filter_cols = [x for x in range(len(header_cells)) if x not in all_forced_cols and 'read counts ok' == header_cells_cf[x]]
            else:  Filter_cols = []

            if 'spacer' in Forced_ColNumbers_by_ColType_dict:   Spacer_cols = Forced_ColNumbers_by_ColType_dict['spacer']
            elif 'spacer' in Forced_ColNames_by_ColTypes_dict: Spacer_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['spacer'])]
            elif not args.Disable_ColType_autoassign:    Spacer_cols = [x for x in range(len(header_cells)) if x not in all_forced_cols and header_cells_cf[x].startswith('(spacer ')]
            else:   Spacer_cols = []



            #### writing header (annotations)

            if mode == 'rows':
                def RowsPos(RowN, ColN):  return RowN + cols_annos_count
                def ColsPos(RowN, ColN):  return 1 + ColN
                def RowsPos__header(RowN, ColN):  return RowN
                def ColsPos__header(RowN, ColN):  return 1 + ColN
            
            if mode == 'cols':
                def RowsPos(RowN, ColN):  return 1 + ColN
                def ColsPos(RowN, ColN):  return RowN + rows_annos_count
                def RowsPos__header(RowN, ColN):  return 1 + ColN
                def ColsPos__header(RowN, ColN):  return RowN

            for sheet in [sheet_rs, sheet_rs__FDR, sheet_pv, sheet_fdr]:
                if sheet is None:  continue
                if mode == 'cols':
                    if rows_annos_count > 0:
                        sheet.set_column(rows_annos_count, rows_annos_count, 1.6)
                        for x in range(len(rows_annotation_lines) - 1):
                            sheet.write(cols_annos_count + 1 + x, rows_annos_count, '', whitespace_format__right_border)
                if mode == 'rows':
                    if cols_annos_count > 0:
                        sheet.set_row(cols_annos_count, 8)
                        for x in range(len(cols_annotation_lines) - 1):
                            sheet.write(cols_annos_count, rows_annos_count + 1 + x, '', whitespace_format__bottom_border)

                for x in range(len(header_cells)):
                    if header_cells[x] in ColFullNames_dict:
                        cell_full_name = ColFullNames_dict[header_cells[x]]
                    else:
                        cell_full_name = header_cells[x]

                    if x in Spacer_cols: continue
                    elif x in logFC_narrow_cols:
                        sheet.write(RowsPos(0, x), ColsPos(0, x), cell_full_name, rotated_format)
                        if mode == 'cols':   sheet.set_column(1 + x, 1 + x, 3.5)
                        if mode == 'rows':   sheet.set_row(1 + x, 3.5)
                    elif x in raw_means_CPM_cols:
                        sheet.write(RowsPos(0, x), ColsPos(0, x),'raw mean',secondary_headers_format)
                    elif x in Rank_1_cols or x in Rank_2_cols:
                        sheet.write(RowsPos__header(0, x), ColsPos__header(0, x), cell_full_name, small_headers_format)
                    else:
                        sheet.write(RowsPos__header(0, x), ColsPos__header(0, x), cell_full_name, headers_format)
                
                RowN = 1

                Max_ColN = 0
                for L in annotation_lines[1:]:
                    C = [C_el[1:-1] if ((C_el.startswith('"') and C_el.endswith('"')) or (C_el.startswith("'") and C_el.endswith("'"))) else C_el for C_el in L]
                    CPM_values = []
                    CPM_values2 = []
                    Max_ColN = max(Max_ColN, len(C))
                    for x in range(len(C)):
                        if x in Spacer_cols: continue
                        if C[x].casefold() in ['na','nan']:
                            if x in logFC_narrow_cols:
                                sheet.write(RowsPos(RowN, x), ColsPos(RowN, x), '', empty_logfc_format)
                            continue

                        if args.max_char_in_cell is not None:
                            if len(C[x]) > args.max_char_in_cell:
                                C[x] = '%s .... [truncated]'%C[x][ : max(1, args.max_char_in_cell - 17)]

                        if (x in P_value_cols) or (x in FDR_cols):
                            write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x], pvalue_mode=True, pvalue_formats=pvalue_formats)
                        elif (x in logFC_cols) or (x in CPM_cols) or (x in Correlation_r_cols) or (x in CPM_cols2):
                            write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x], format2)
                        elif x in logFC_narrow_cols:
                            if ConvertableToNumeric(C[x]):
                                write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x], format2_narrow)
                            else:
                                sheet.write(RowsPos(RowN, x), ColsPos(RowN, x), C[x], empty_logfc_format)

                        elif (x in Score_cols) or (x in raw_means_CPM_cols) or (x in Index1_cols) or (x in Index1Err_cols) or (x in LogFC_array_cols) or (x in LogFC_array2_cols):
                            write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x], format1_center)
                        elif (x in Index2_cols) or (x in Index2Err_cols):
                            write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x], format2_center)
                        elif x in Index_cols or x in gradual_30_cols or x in FFPM_cols:
                            write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x], format_center)
                        elif x in logCPM_cols:
                            write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x], format1)
                        elif x in percentage_cols:
                            write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x], whitespace_format1)
                        elif x in gray_gradient_cols or x in warm_gradient_cols or x in cold_gradient_cols:
                            write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x], format0_center)
                        elif x in delta_logfc_cols:
                            write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x], whitespace_format1)
                        elif x in percents_updown_cols:
                            write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x], percentage_format0)
                        elif x in Taxon_cols:
                            sheet.write(RowsPos(RowN, x), ColsPos(RowN, x), C[x], italic)
                        elif x in Color_cols:
                            if not C[x].startswith('#') or len(C[x]) != 7:
                                print('Incorrect color format %s'%C[x])
                                exit(127)
                            CurrentFormat = Workbook.add_format()
                            CurrentFormat.set_bg_color(C[x])
                            sheet.write(RowsPos(RowN, x), ColsPos(RowN, x), ' ', CurrentFormat)
                        elif x in Filter_cols:
                            cf = C[x].casefold()
                            if 'true' == cf or 'ok' == cf or 'passed' == cf:
                                sheet.write(RowsPos(RowN, x), ColsPos(RowN, x),'yes', format_center)
                        elif x in Counts_cols:
                            write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x], format0)
                        elif x in Rank_1_cols or x in Rank_2_cols:
                            WriteFormatted_Rank_Cell(Workbook, sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x],
                                upper_part = 17, lower_part = 17, extent = 1.7, default_format = default_Rank_format)
                        else:  write_number__mod(sheet, RowsPos(RowN, x), ColsPos(RowN, x), C[x])

                        if not args.Simple_LogFC_format and not args.Bidirectional_bars_in_LogFC_cells and x in logFC_cols:
                            FormatLogFC_Cell(sheet, 2**float(C[x]), RowsPos(RowN, x), ColsPos(RowN, x), coeff1 = 7/args.LogFC_cell_abs_max, coeff2 = 7/args.LogFC_cell_abs_max)

                        if x in CPM_cols:
                            try:  CPM_values.append(float(C[x]))
                            except ValueError:  pass
                        
                        if x in CPM_cols2:
                            try:  CPM_values2.append(float(C[x]))
                            except ValueError:  pass

                    CPM_values_pack = [CPM_values, CPM_values2]
                    CPM_cols_pack = [CPM_cols, CPM_cols2]

                    for current_CPM_values, current_CPM_cols in zip(CPM_values_pack, CPM_cols_pack):
                        if args.Add_cpm_cells_formatting and len(current_CPM_values) > 0:
                            avg_current_CPM = sum(current_CPM_values)/len(current_CPM_values)
                            # median_current_CPM = statistics.median(current_CPM_values)
                            current_CPM_values_mod = [x + 1 + args.cpm_cells_formatting__const_add for x in current_CPM_values]
                            gmean_current_CPM = 1
                            for x in current_CPM_values_mod:
                                gmean_current_CPM = gmean_current_CPM*(x/(avg_current_CPM + 1 + args.cpm_cells_formatting__const_add))
                            gmean_current_CPM = (gmean_current_CPM**(1/len(current_CPM_values_mod)))*(avg_current_CPM  + 1 + args.cpm_cells_formatting__const_add)

                            if gmean_current_CPM == 0:
                                continue

                            current_LogFCs = [math.log2((x + args.cpm_cells_formatting__const_add + 1) / (gmean_current_CPM)) for x in current_CPM_values]
                            current_LogFCs = [min(max(x, (-1)*args.cpm_cells_formatting__LogFC_max), args.cpm_cells_formatting__LogFC_max) for x in current_LogFCs]
                            current_LogFCs_color_coords = [0.5 + x / args.cpm_cells_formatting__LogFC_max / 2 for x in current_LogFCs]

                            cpm_color_coord = (avg_current_CPM/args.cpm_cells_formatting__cpm_max)**args.cpm_cells_formatting__cpm_degree

                            # print(CPM_cols)
                            for n in range(len(current_CPM_cols)):
                                col_n = current_CPM_cols[n]
                                try:  val = float(C[col_n])
                                except:  continue

                                CurrentFormat = Workbook.add_format()
                                CurrentFormat.set_align('center')

                                if math.isnan(val):   BgC = '#ffffff'
                                else:
                                    try:
                                        BgC = color(  TwoD_Gradient(Color_gradient_schemas[args.cpm_cells_formatting__gradient_color_schema - 1],
                                                                    cpm_color_coord, current_LogFCs_color_coords[n]) )
                                    except:
                                        print('Incorrect CPM columns. Please check')
                                        print('n = %d, len(current_LogFCs_color_coords) = %d'%(n, len(current_LogFCs_color_coords)))
                                        try:  print('val = %g, current_LogFCs_color_coords[n] = %g, cpm_color_coord = %g'%(val, current_LogFCs_color_coords[n], cpm_color_coord))
                                        except:  pass
                                        exit(127)

                                # if C[0] == 'Aquificae':
                                #     print('Aquificae')
                                #     print(current_LogFCs)
                                #     print(gmean_current_CPM)
                                #     print(args.cpm_cells_formatting__const_add)
                                # BgC = '#aaaaaa'

                                CurrentFormat.set_bg_color(BgC)
                                CurrentFormat.set_border(0)
                                CurrentFormat.set_num_format('0.' + '0'*args.CPM_digits)
                                write_number__mod(sheet, RowsPos(RowN, col_n), ColsPos(RowN, col_n), val, format=CurrentFormat)

                            # args.cpm_cells_formatting__const_add
                            # args.cpm_cells_formatting__LogFC_max
                            # args.cpm_cells_formatting__cpm_max
                            # args.cpm_cells_formatting__cpm_degree
                            # args.cpm_cells_formatting__gradient_color_schema


                    RowN += 1

                Last_Row = RowN

                if args.Simple_LogFC_format:
                    for ColN in logFC_cols + logFC_narrow_cols:
                        sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                                                                  'min_color': "#1170f1",'mid_color': "#ffffff",'max_color': "#ec4a18",
                                                                  'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                                  'min_value': -4.1, 'mid_value': 0.0, 'max_value': 4.1})
                if args.Bidirectional_bars_in_LogFC_cells:
                    for ColN in logFC_cols + logFC_narrow_cols:
                        sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': 'data_bar',
                                                                  'bar_color': '#f79254',
                                                                  'bar_negative_color':'#59c2ff',
                                                                  'min_type':'num','max_type':'num',
                                                                  'min_value':-3.1,'max_value': 3.1,
                                                                  'min_length':0, 'max_length':100})

                for ColN in Correlation_r_cols:
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                                                              'min_color': "#5c86f1",'mid_color': "#ffffff",'max_color': "#ee6c44",
                                                              'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                              'min_value': -0.98, 'mid_value': 0.0, 'max_value': 0.98})

                for ColN in Counts_cols + Index_cols + Index1_cols + Index2_cols:
                    # sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                    #                                           'min_color': "#f2dc49",'mid_color': "#9af762",'max_color': "#7fa5f2",
                    #                                           'min_type': 'percent','mid_type': 'percent','max_type': 'percent',
                    #                                           'min_value': 0, 'mid_value': 50, 'max_value': 100})

                    # sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                    #                                           'min_color': "#f96f45",'mid_color': "#ffffff",'max_color': "#71e27c",
                    #                                           'min_type': 'percentile','mid_type': 'percentile','max_type': 'percentile',
                    #                                           'min_value': 0, 'mid_value': 50, 'max_value': 98})

                    # sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                    #                                           'min_color': "#5c86f1",'mid_color': "#ffffff",'max_color': "#ee6c44",
                    #                                           'min_type': 'num','mid_type': 'num','max_type': 'num',
                    #                                           'min_value': 0.0, 'mid_value': 50.0, 'max_value': 100.0})

                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                                                              'min_color': "#ffffff",'mid_color': "#afff91",'max_color': "#a1cdff",
                                                              'min_type': 'percentile','mid_type': 'percentile','max_type': 'percentile',
                                                              'min_value': 0, 'mid_value': 50, 'max_value': 98})
                
                for ColN in gradual_30_cols:
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                                                              'min_color': "#ffffff",'mid_color': "#ffe08d",'max_color': "#9ece49",
                                                              'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                              'min_value': 0, 'mid_value': 14, 'max_value': 30})

                for ColN in FFPM_cols:
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                                                              'min_color': "#ffffff",'mid_color': "#ffe08d",'max_color': "#9ece49",
                                                              'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                              'min_value': 0, 'mid_value': 2, 'max_value': 5})

                for ColN in delta_logfc_cols:
                    #'min_color': "#06ade0",'mid_color': "#ffffff",'max_color': "#efb00e",
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': 'data_bar',
                                                              'bar_color': '#efb00e',
                                                              'bar_negative_color':'#06ade0',
                                                              'bar_border_color': '#bc8803',
                                                              'bar_negative_border_color': '#0387af',
                                                              'min_type':'num','max_type':'num',
                                                              'min_value':-3.1,'max_value': 3.1,
                                                              'min_length':0, 'max_length':100})
                    if ColN in Hidden_cols:  continue
                    if mode == 'cols':  sheet.set_column(ColN, ColN, 15)

                for ColN in percents_updown_cols:
                    #'min_color': "#06ade0",'mid_color': "#ffffff",'max_color': "#efb00e",
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                                                              'min_color': "#317bf5",'mid_color': "#ffffff",'max_color': "#ec692f",
                                                              'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                              'min_value': -1, 'mid_value': 0, 'max_value': 4})
                    
                    # sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': 'data_bar',
                    #                                           'bar_color': '#efb00e',
                    #                                           'bar_negative_color':'#06ade0',
                    #                                           'bar_border_color': '#bc8803',
                    #                                           'bar_negative_border_color': '#0387af',
                    #                                           'min_type':'num','max_type':'num',
                    #                                           'min_value': -1.0,'max_value': 3.0,
                    #                                           'min_length':0, 'max_length':100})

                    # if ColN in Hidden_cols:  continue
                    # sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 15)

                if mode == 'cols':
                    # if args.Insert_spasers_between_rows:
                    #     for ColN in range(0, Cols_count*2):
                    #         if ColN in Hidden_cols:  continue
                    #         sheet.set_column(ColN, ColN, 9, whitespace_format)

                    for ColN in Rank_1_cols:
                        sheet.set_column(ColN, ColN, 2.5)
                    
                    for ColN in Rank_2_cols:
                        sheet.set_column(ColN, ColN, 3.5)

                    for ColN in Taxon_cols:
                        sheet.set_column(ColN, ColN, 35)
                        

                for ColN in logCPM_cols:
                    if args.Metagenome_mode:  max_value = 20
                    else:  max_value = 10
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': 'data_bar','bar_color': '#d89c0d',
                                                                                          'min_type':'num','max_type':'num',
                                                                                          'min_value':0,'max_value': max_value})
                for ColN in percentage_cols:
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': 'data_bar','bar_color': '#1cbaff',
                                                                                          'min_type':'num','max_type':'num',
                                                                                          'min_value':0,'max_value': 100})

                for ColN in warm_gradient_cols:
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                                                              'min_color': "#ffffff",'mid_color': "#ffeb84",'max_color': "#f8696b",
                                                              'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                              'min_value': 0, 'mid_value': 31, 'max_value': 82})

                for ColN in cold_gradient_cols:
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                                                              'min_color': "#ffffff",'mid_color': "#85ebe6",'max_color': "#3c86e8",
                                                              'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                              'min_value': 0, 'mid_value': 32, 'max_value': 90})

                for ColN in gray_gradient_cols:
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                                                              'min_color': "#ffffff",'mid_color': "#d9d9d9",'max_color': "#808080",
                                                              'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                              'min_value': 0, 'mid_value': 32, 'max_value': 90})

                for ColN in P_value_cols + FDR_cols:
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                                                              'min_color': "#9ece49",'mid_color': "#ffe08d",'max_color': "#ffffff",
                                                              'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                              'min_value': 0, 'mid_value': 0.0002, 'max_value': 0.07})

                for ColN in LogFC_array_cols + LogFC_array2_cols:
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                                                              'min_color': "#06ade0",'mid_color': "#ffffff",'max_color': "#efb00e",
                                                              'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                              'min_value': -args.LogFC_array_abs_max, 'mid_value': 0.0, 'max_value': args.LogFC_array_abs_max})

                # if args.Metagenome_mode:
                #     for ColN in Score_cols:
                #         sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type': '3_color_scale',
                #                                                            'min_color': "#ffffff", 'mid_color': "#ffe08d",
                #                                                            'max_color': "#9ece49",
                #                                                            'min_type': 'num', 'mid_type': 'num',
                #                                                            'max_type': 'num',
                #                                                            'min_value': 2, 'mid_value': 12,
                #                                                            'max_value': 25})

                # Light red fill with dark red text.
                format_lincRNA = Workbook.add_format({'bg_color':   '#c3d594'})

                # Light yellow fill with dark yellow text.
                format_antisense = Workbook.add_format({'bg_color':   '#d9c188'})

                # Green fill with dark green text.
                format_pseudogene = Workbook.add_format({'bg_color':   '#a0b4c7'})

                for ColN in Biotype_cols:
                    if ColN in Hidden_cols:  continue
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type':'text',
                                                       'criteria': 'containing', 'value':    'lincRNA', 'format':   format_lincRNA})
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type':'text',
                                                       'criteria': 'containing', 'value':    'lncRNA', 'format':   format_lincRNA})
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type':'text',
                                                       'criteria': 'containing', 'value':    'antisense', 'format':   format_antisense})
                    sheet.conditional_format(RowsPos(1, ColN), ColsPos(1, ColN), RowsPos(Last_Row, ColN), ColsPos(Last_Row, ColN), {'type':'text',
                                                       'criteria': 'containing', 'value':    'pseudogene', 'format':   format_pseudogene})
                    if mode == 'cols':  sheet.set_column(ColN, ColN, 18)


                if mode == 'cols':  
                    for ColN in Gene_Name_cols:
                        if ColN in Hidden_cols:  continue
                        sheet.set_column(ColN, ColN, 35)

                    # if args.Metagenome_mode:  sheet.set_column(WSparkLine_Offset(args, 0), WSparkLine_Offset(args, 0), 19)

                    for ColN in warm_gradient_cols + cold_gradient_cols + gray_gradient_cols:
                        if ColN in Hidden_cols:  continue
                        sheet.set_column(ColN, ColN, 5.5)

                    for ColN in Spacer_cols:
                        if ColN in Hidden_cols:  continue
                        sheet.set_column(ColN, ColN, 2)

                    for ColN in Hidden_cols:
                        sheet.set_column(ColN, ColN, None, None, {'hidden': True})

                    for x in range(len(header_cells)):
                        if header_cells[x] in Forced_ColWidths_by_ColNames:
                            sheet.set_column(x, x, Forced_ColWidths_by_ColNames[header_cells[x]])


        
        if not args.Generate_SingleBook:
            Workbook.close()

    if args.Generate_SingleBook:
        Workbook.close()

    if args.verbose:
        print('\rExcel workbook generation completed                               ')


