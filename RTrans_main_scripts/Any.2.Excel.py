__author__ = 'George'
import sys,math,os,math,shlex
import xlsxwriter
import argparse,glob
import statistics
from xlsxwriter.utility import xl_range
import numpy

from Excel_formats import Create_workbook_Formats, Create_workbook_Formats_whitespace

def ResizeArray(src_array, desired_length):
    initial_length = len(src_array)
    if initial_length == desired_length:
        return numpy.array(src_array)

    res_array = numpy.zeros(desired_length,dtype=float)
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

def Smoothed_Trimmed_Mean(array, from_percentile, to_percentile):
    if from_percentile > to_percentile:
        temp = from_percentile
        from_percentile = to_percentile
        to_percentile = temp
    from_percentile = int(max(0,from_percentile))
    to_percentile = int(min(100,to_percentile))
    weights = ResizeArray([0]*from_percentile + [1]*(to_percentile - from_percentile) + [0]*(100 - to_percentile),len(array))
    sorted_array = sorted(array)
    # return sum([weights[n]*sorted_array[n] for n in range(len(sorted_array))])/sum(weights)
    return numpy.average(sorted_array, weights=weights)

def Smoothed_Trimmed_StDev(array, from_percentile, to_percentile, ForcedAverage = None):
    if from_percentile > to_percentile:
        temp = from_percentile
        from_percentile = to_percentile
        to_percentile = temp
    from_percentile = int(max(0,from_percentile))
    to_percentile = int(min(100,to_percentile))
    weights = ResizeArray([0]*from_percentile + [1]*(to_percentile - from_percentile) + [0]*(100 - to_percentile),len(array))
    sorted_array = sorted(array)
    # return sum([weights[n]*sorted_array[n] for n in range(len(sorted_array))])/sum(weights)
    # print(sorted_array)
    # print(weights)
    return weighted_std(sorted_array, weights=weights,ForcedAverage = ForcedAverage)

def Smoothed_Trimmed_RelStDev(array,from_percentile,to_percentile,ForcedAverage = None):
    if from_percentile > to_percentile:
        temp = from_percentile
        from_percentile = to_percentile
        to_percentile = temp
    from_percentile = int(max(0,from_percentile))
    to_percentile = int(min(100,to_percentile))
    weights = ResizeArray([0]*from_percentile + [1]*(to_percentile - from_percentile) + [0]*(100 - to_percentile),len(array))
    sorted_array = sorted(array)
    # return sum([weights[n]*sorted_array[n] for n in range(len(sorted_array))])/sum(weights)
    # print(sorted_array)
    # print(weights)
    return weighted_std(sorted_array, weights=weights, ForcedAverage = ForcedAverage)/numpy.average(sorted_array, weights=weights)


def geo_mean_overflow(iterable):
    a = numpy.log(iterable)
    return numpy.exp(a.sum()/len(a))

def Smoothed_Trimmed_GeoMean(array, from_percentile, to_percentile):
    if from_percentile > to_percentile:
        temp = from_percentile
        from_percentile = to_percentile
        to_percentile = temp
    from_percentile = int(max(0,from_percentile))
    to_percentile = int(min(100,to_percentile))
    weights = ResizeArray([0]*from_percentile + [1]*(to_percentile - from_percentile) + [0]*(100 - to_percentile),len(array))
    sorted_array = sorted(array)
    # return sum([weights[n]*sorted_array[n] for n in range(len(sorted_array))])/sum(weights)
    return numpy.exp(numpy.average(numpy.log(sorted_array), weights=weights))


ScoreColorPoints = [(225,225,225),(0,0,0)]

def ConvertableToNumeric(text):
    try:
        d = float(text)
        return not (math.isnan(d) or math.isinf(d))
    except ValueError:
        return False

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

TwoD_ColorPointsSet5 = [
    ['4771dd','ffffff', 'fa8323'],
    ['4771dd','ffffff', 'fa8323']
    ]

Color_gradient_schemas = [TwoD_ColorPointsSet1, TwoD_ColorPointsSet2, TwoD_ColorPointsSet3, TwoD_ColorPointsSet4, TwoD_ColorPointsSet5]


# Sparklines_CPM_Groups_Gradient = ['47a7d0', '418fde']
Sparklines_CPM_Groups_Gradient = ['418fde', '418fde']

Sparklines_Preds_Gradient = ['9a9d94', '88be2a', 'b3bd3f', 'c09b19']
Sparklines_Preds_Gradient = ['c1bf18', '88be2a', 'c08b19']
Sparklines_Preds_Gradient = ['f7a24c', 'cfcb47', '6dd681', '5eddc9', '56b7f7']
Sparklines_Preds_Gradient = ['f59342', 'f0bf0d', 'b5e54a', '79dcf1', '81aeea']

HM_Preds_Gradient = ['f59342', 'f0bf0d', 'b5e54a', '79dcf1', '81aeea']


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

#
#
# def write_number__mod(sheet,rowN,colN,text,format=None):
#     try:
#         d = float(text)
#         if format != None:
#             if math.isnan(d) or math.isinf(d): sheet.write(rowN,colN,text,format)
#             else: sheet.write_number(rowN,colN,d,format)
#         else:
#             if math.isnan(d) or math.isinf(d): sheet.write(rowN,colN,text)
#             else: sheet.write_number(rowN,colN,d)
#     except:
#         if format != None:
#             sheet.write(rowN,colN,text,format)
#         else:
#             sheet.write(rowN,colN,text)


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

def ToBool(string, allow_None_as_True = True):
    if allow_None_as_True and string is None:
        return True
    string = str(string).casefold()
    if string in ['yes','y','on','enable','true','da','t']:
        return True
    elif string in ['no','n','off','disable','false','net','f']:
        return False
    print('Parameter "%s" should be either "yes" or "no". Please correct'%string)
    exit()

# def ToBool(word):
#     word = word.casefold()
#     if word in ['yes','y','on']: return  True
#     elif word in ['no','n','off']: return  False
#     print('Incorrect input "%s"'%word)
#     exit(127)

def to_float(value):
  try:    return float(value)
  except ValueError:   return float('NaN')

def WSparkLine_Offset(args, ColN, for_CPM_heatmap_cells = False, for_LogFC_heatmap_cells = False):
    if args.CPM_sparklines_start_col is None and args.LogFC_sparklines_start_col is None:
        return ColN
    
    elif args.CPM_sparklines_start_col is not None and args.LogFC_sparklines_start_col is None:
        if ColN >= args.CPM_sparklines_start_col - 1:
            adj_ColN = ColN + len(args.CPM_sparklines_groups) + args.insert_empty_col_before_CPM_sparklines
            if args.include_CPM_heatmaps and not for_CPM_heatmap_cells:
                adj_ColN += len(args.CPM_sparklines_groups) + 1 + sum([len(x) for x in args.CPM_sparklines_groups.values()])
            return adj_ColN
        else:
            return ColN

    elif args.CPM_sparklines_start_col is None and args.LogFC_sparklines_start_col is not None:
        if ColN >= args.LogFC_sparklines_start_col - 1:
            adj_ColN = ColN + len(args.LogFC_sparklines_groups) + args.insert_empty_col_before_LogFC_sparklines
            if args.include_LogFC_heatmaps and not for_LogFC_heatmap_cells:
                adj_ColN += len(args.LogFC_sparklines_groups) + 1 + sum([len(x) for x in args.LogFC_sparklines_groups.values()])
            return adj_ColN
        else:
            return ColN

    else:
        if ColN >= args.CPM_sparklines_start_col - 1 and ColN >= args.LogFC_sparklines_start_col - 1:
            adj_ColN = ColN + len(args.CPM_sparklines_groups) + args.insert_empty_col_before_CPM_sparklines + len(args.LogFC_sparklines_groups) + args.insert_empty_col_before_LogFC_sparklines
            if args.include_CPM_heatmaps and not for_CPM_heatmap_cells:
                adj_ColN += len(args.CPM_sparklines_groups) + 1 + sum([len(x) for x in args.CPM_sparklines_groups.values()])
            if args.include_LogFC_heatmaps and not for_LogFC_heatmap_cells:
                adj_ColN += len(args.LogFC_sparklines_groups) + 1 + sum([len(x) for x in args.LogFC_sparklines_groups.values()])
            return adj_ColN
        
        elif ColN >= args.CPM_sparklines_start_col - 1:
            adj_ColN = ColN + len(args.CPM_sparklines_groups) + args.insert_empty_col_before_CPM_sparklines
            if args.include_CPM_heatmaps and not for_CPM_heatmap_cells:
                adj_ColN += len(args.CPM_sparklines_groups) + 1 + sum([len(x) for x in args.CPM_sparklines_groups.values()])
            return adj_ColN
                
        elif ColN >= args.LogFC_sparklines_start_col - 1:
            adj_ColN = ColN + len(args.LogFC_sparklines_groups) + args.insert_empty_col_before_LogFC_sparklines
            if args.include_LogFC_heatmaps and not for_LogFC_heatmap_cells:
                adj_ColN += len(args.LogFC_sparklines_groups) + 1 + sum([len(x) for x in args.LogFC_sparklines_groups.values()])
            return adj_ColN

        else:
            return ColN



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creating Excel workbooks from txt|tsv LogFC, logCMP tables')
    parser.add_argument('-i','--in', dest='Input_FileNames', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('-d','--dir', dest='Input_Dir', nargs='?',action='store', required=False,default=None, help='')
    parser.add_argument('-ext','--extensions', dest='file_extensions', nargs='?',action='store', required=False,default=['txt','tsv'], help='')
    parser.add_argument('-m','--max-depth', dest='max_depth', nargs='?',action='store', required=False,default='10000', help='')
    parser.add_argument('-r','--recursive', dest='Recursive', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('-o','--out-excel', dest='Workbook_FileName', nargs='?',action='store', required=False,default=None, help='')
    parser.add_argument('--one-book', dest='Generate_SingleBook', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--max-formatted-cells', dest='max_formatted_cells', nargs='?',action='store', required=False,default='35000', help='')
    parser.add_argument('--sheet-names', dest='Sheet_names', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--metagenome-mode', dest='Metagenome_mode', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--simple-logfc-format', dest='Simple_LogFC_format', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--bidirectional-bars-in-logfc-cells', dest='Bidirectional_bars_in_LogFC_cells', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--add-cpm-cells-formatting', dest='Add_cpm_cells_formatting', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--cpm-cells-formatting--logfc-max', dest='cpm_cells_formatting__LogFC_max', nargs='?',action='store', required=False,default='1.5', help='')
    parser.add_argument('--cpm-cells-formatting--cpm-max', dest='cpm_cells_formatting__cpm_max', nargs='?',action='store', required=False,default='500', help='')
    parser.add_argument('--cpm-cells-formatting--cpm-degree', dest='cpm_cells_formatting__cpm_degree', nargs='?',action='store', required=False, default='0.7', help='')
    parser.add_argument('--cpm-cells-formatting--const-add', dest='cpm_cells_formatting__const_add', nargs='?',action='store', required=False, default='3', help='')
    parser.add_argument('--cpm-cells-formatting--trimming-low', dest='cpm_cells_formatting__trimming_low', nargs='?',action='store', required=False, default='0', help='')
    parser.add_argument('--cpm-cells-formatting--trimming-high', dest='cpm_cells_formatting__trimming_high', nargs='?',action='store', required=False, default='0', help='')
    parser.add_argument('--cpm-cells-formatting--gradient-color-schema', dest='cpm_cells_formatting__gradient_color_schema', nargs='?',action='store', required=False,default='1', help='')
    parser.add_argument('--logfc-cells-formatting--gradient-color-schema', dest='logfc_cells_formatting__gradient_color_schema', nargs='?',action='store', required=False, default='1', help='')
    parser.add_argument('--logfc-cells-formatting--gradient-range', dest='logfc_cells_formatting__gradient_range', nargs='?',action='store', required=False, default='3', help='')
    parser.add_argument('-v', '--verbose', dest='verbose', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--forced-cols', dest='Forced_ColTypes_by_ColNumber', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--forced-col-types-by-names', dest='Forced_ColTypes_by_ColName', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--forced-col-types-by-names-quote-symbol', dest='Forced_ColTypes_by_ColName__quote_symbol', nargs='?',action='store', required=False,default='&', help='')
    parser.add_argument('--forced-col-widths-by-names', dest='Forced_ColWidths_by_ColName', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--forced-col-widths-by-names-quote-symbol', dest='Forced_ColWidths_by_ColName__quote_symbol', nargs='?',action='store', required=False,default='&', help='')
    parser.add_argument('--forced-col-widths', dest='Forced_ColWidths', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--col-full-names', dest='ColFullNames', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--disable-coltype-autoassign', dest='Disable_ColType_autoassign', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--predictor-rows-count', dest='Predictor_rows_count', nargs='?',action='store', required=False,default='0', help='')
    parser.add_argument('--group-predictors', dest='Group_predictors', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--predictor-light-font', dest='Predictor_light_font', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--predictor-blue-yellow-formatting', dest='Predictor_blue_yellow_formatting', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--cpm-digits', dest='CPM_digits', nargs='?',action='store', required=False,default='1', help='')
    parser.add_argument('--cpm-narrow-digits', dest='CPM_narrow_digits', nargs='?',action='store', required=False, default='0', help='')
    parser.add_argument('--freeze-cols', dest='freeze_cols_count', nargs='?',action='store', required=False, default='0', help='')
    parser.add_argument('--freeze-rows', dest='freeze_rows_count', nargs='?',action='store', required=False, default='0', help='')
    parser.add_argument('--whitespace-mode', dest='create_White_background', nargs='?',action='store', required=False, default='no', help='')

    parser.add_argument('--logfc-abs-max', dest='LogFC_cell_abs_max', nargs='?',action='store', required=False,default='5.8', help='')
    parser.add_argument('--logfc-array-abs-max', dest='LogFC_array_abs_max', nargs='?',action='store', required=False,default='2.5', help='')
    parser.add_argument('--entry-type', dest='entry_type', nargs='?',action='store', required=False,default=None, help='')
    parser.add_argument('--first-col-italic', dest='first_col_italic', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--first-col-width', dest='first_col_width', nargs='?',action='store', required=False,default=None, help='')
    parser.add_argument('--first-col-is-formatted', dest='first_col_is_formatted', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--add-predictor-sparklines', dest='add_predictor_sparklines', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--insert-spacer-between-rows', dest='Insert_spasers_between_rows', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--max-strings', dest='max_strings', nargs='?',action='store', required=False, default='unlimited', help='')
    parser.add_argument('--max-chars-in-cell', dest='max_char_in_cell', nargs='?',action='store', required=False, default='unlimited', help='')
    parser.add_argument('--skip-absent-files', dest='Skip_absent_files', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--check-for-na-in-cpm-cols', dest='check_for_NA_in_CPM_cols', nargs='?',action='store', required=False, default='no', help='')
    
    parser.add_argument('--cpm-sparklines-start-col', dest='CPM_sparklines_start_col', nargs='?',action='store', required=False,default=None, help='')
    parser.add_argument('--cpm-sparklines-groups', dest='CPM_sparklines_groups', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--cpm-sparklines-start-col--by-sheet', dest='CPM_sparklines_start_col__by_sheet', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--cpm-sparklines-groups--by-sheet', dest='CPM_sparklines_groups__by_sheet', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--cpm-sparklines-x--by-sheet--quote-symbol', dest='CPM_sparklines__quote_symbol', nargs='?',action='store', required=False,default='&', help='')
    parser.add_argument('--cpm-sparklines-mode', dest='CPM_sparklines_mode', nargs='?',action='store', required=False,default='linear', help='')
    parser.add_argument('--cpm-sparklines-rel-scale-cpm-const-add', dest='CPM_sparklines_rel_scale_cpm_const_add', nargs='?',action='store', required=False,default='1', help='')
    parser.add_argument('--cpm-sparklines-rel-scale-cpm-trimming-low', dest='CPM_sparklines_rel_scale_cpm_trimming_low', nargs='?',action='store', required=False, default='0', help='')
    parser.add_argument('--cpm-sparklines-rel-scale-cpm-trimming-high', dest='CPM_sparklines_rel_scale_cpm_trimming_high', nargs='?',action='store', required=False, default='0', help='')
    parser.add_argument('--cpm-sparklines-col-width-factor', dest='CPM_sparklines_col_width_factor', nargs='?',action='store', required=False,default='0.5', help='')
    parser.add_argument('--cpm-sparklines-equal-barwidths', dest='CPM_sparklines_equal_barwidths', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--cpm-sparklines-logfc-mode-axis-limit', dest='CPM_sparklines_logfc_mode_axis_limit', nargs='?',action='store', required=False,default='2', help='')
    parser.add_argument('--insert-empty-col-before-CPM-sparklines', dest='insert_empty_col_before_CPM_sparklines', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--include-cpm-heatmaps', dest='include_CPM_heatmaps', nargs='?',action='store', required=False, default='no', help='')

    parser.add_argument('--logfc-sparklines-start-col', dest='LogFC_sparklines_start_col', nargs='?',action='store', required=False,default=None, help='')
    parser.add_argument('--logfc-sparklines-groups', dest='LogFC_sparklines_groups', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--logfc-sparklines-start-col--by-sheet', dest='LogFC_sparklines_start_col__by_sheet', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--logfc-sparklines-groups--by-sheet', dest='LogFC_sparklines_groups__by_sheet', nargs='*',action='store', required=False,default=None, help='')
    parser.add_argument('--logfc-sparklines-x--by-sheet--quote-symbol', dest='LogFC_sparklines__quote_symbol', nargs='?',action='store', required=False,default='&', help='')
    parser.add_argument('--logfc-sparklines-col-width-factor', dest='LogFC_sparklines_col_width_factor', nargs='?',action='store', required=False,default='0.5', help='')
    parser.add_argument('--logfc-sparklines-equal-barwidths', dest='LogFC_sparklines_equal_barwidths', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--logfc-sparklines-axis-limit', dest='LogFC_sparklines_axis_limit', nargs='?',action='store', required=False,default='3', help='')
    parser.add_argument('--insert-empty-col-before-logfc-sparklines', dest='insert_empty_col_before_LogFC_sparklines', nargs='?',action='store', required=False,default='yes', help='')
    parser.add_argument('--include-logfc-heatmaps', dest='include_LogFC_heatmaps', nargs='?',action='store', required=False, default='no', help='')

    parser.add_argument('--dist-matrix-mode', dest='dist_matrix_mode', nargs='?',action='store', required=False,default='no', help='')
    parser.add_argument('--dist-matrix--digits', dest='dist_matrix__digits', nargs='?',action='store', required=False,default='1', help='')
    parser.add_argument('--narrow-cols-width', dest='narrow_cols_width', nargs='?', action='store', required=False, default='3.5', help='')

    parser.add_argument('--red-bars-color', dest='red_bars_color', nargs='?', action='store', required=False, default='#fa8323', help='')
    parser.add_argument('--red-bars-max-value', dest='red_bars_max_value', nargs='?', action='store', required=False, default='99', help='')
    parser.add_argument('--red-bars-max-type', dest='red_bars_max_type', nargs='?', action='store', required=False, default='percentile', help='')
    parser.add_argument('--red-bars-min-value', dest='red_bars_min_value', nargs='?', action='store', required=False, default='0', help='')
    parser.add_argument('--red-bars-min-type', dest='red_bars_min_type', nargs='?', action='store', required=False, default='num', help='')
    parser.add_argument('--blue-bars-color', dest='blue_bars_color', nargs='?', action='store', required=False, default='#4771dd', help='')
    parser.add_argument('--blue-bars-max-value', dest='blue_bars_max_value', nargs='?', action='store', required=False, default='99', help='')
    parser.add_argument('--blue-bars-max-type', dest='blue_bars_max_type', nargs='?', action='store', required=False, default='percentile', help='')
    parser.add_argument('--blue-bars-min-value', dest='blue_bars_min_value', nargs='?', action='store', required=False, default='0', help='')
    parser.add_argument('--blue-bars-min-type', dest='blue_bars_min_type', nargs='?', action='store', required=False, default='num', help='')

    parser.add_argument('--red-gradient-color', dest='red_gradient_color', nargs='?', action='store', required=False, default='#fa8323', help='')
    parser.add_argument('--red-gradient-max-value', dest='red_gradient_max_value', nargs='?', action='store', required=False, default='99', help='')
    parser.add_argument('--red-gradient-max-type', dest='red_gradient_max_type', nargs='?', action='store', required=False, default='percentile', help='')
    parser.add_argument('--red-gradient-min-value', dest='red_gradient_min_value', nargs='?', action='store', required=False, default='0', help='')
    parser.add_argument('--red-gradient-min-type', dest='red_gradient_min_type', nargs='?', action='store', required=False, default='num', help='')
    parser.add_argument('--blue-gradient-color', dest='blue_gradient_color', nargs='?', action='store', required=False, default='#4771dd', help='')
    parser.add_argument('--blue-gradient-max-value', dest='blue_gradient_max_value', nargs='?', action='store', required=False, default='99', help='')
    parser.add_argument('--blue-gradient-max-type', dest='blue_gradient_max_type', nargs='?', action='store', required=False, default='percentile', help='')
    parser.add_argument('--blue-gradient-min-value', dest='blue_gradient_min_value', nargs='?', action='store', required=False, default='0', help='')
    parser.add_argument('--blue-gradient-min-type', dest='blue_gradient_min_type', nargs='?', action='store', required=False, default='num', help='')

    parser.add_argument('--index-cols-min-type', dest='index_cols_min_type', nargs='?', action='store', required=False, default='percentile', help='')
    parser.add_argument('--index-cols-min-value', dest='index_cols_min_value', nargs='?', action='store', required=False, default='0', help='')
    parser.add_argument('--index-cols-mid-type', dest='index_cols_mid_type', nargs='?', action='store', required=False, default='percentile', help='')
    parser.add_argument('--index-cols-mid-value', dest='index_cols_mid_value', nargs='?', action='store', required=False, default='50', help='')
    parser.add_argument('--index-cols-max-type', dest='index_cols_max_type', nargs='?', action='store', required=False, default='percentile', help='')
    parser.add_argument('--index-cols-max-value', dest='index_cols_max_value', nargs='?', action='store', required=False, default='98', help='')

    parser.add_argument('--cov100norm-range', dest='cov100norm_range', nargs='*',action='store', required=False, default=['40', '100', '160'], help='')
    parser.add_argument('--cov100norm-colors', dest='cov100norm_colors', nargs='*',action='store', required=False, default=['#4771dd', '#ffffff', '#fa8323'], help='')

    parser.add_argument('--additional-text-file', dest='additional_text_file', nargs='?',action='store', required=False, default=None, help='')
    
    # parser.add_argument('--exclude-cols', dest='exclude_Cols', nargs='*',action='store', required=False,default=None, help='')
    args = parser.parse_args()


    for x in ['Generate_SingleBook','Recursive','Metagenome_mode','Simple_LogFC_format','Add_cpm_cells_formatting','verbose',
              'Group_predictors', 'Predictor_light_font', 'Predictor_blue_yellow_formatting', 'Disable_ColType_autoassign',
              'first_col_italic', 'insert_empty_col_before_CPM_sparklines', 'insert_empty_col_before_LogFC_sparklines',
              'first_col_is_formatted', 'Insert_spasers_between_rows', 'create_White_background', 'Skip_absent_files',
              'Bidirectional_bars_in_LogFC_cells', 'dist_matrix_mode', 'check_for_NA_in_CPM_cols',
              'include_CPM_heatmaps', 'include_LogFC_heatmaps']:
        #print(x)
        setattr(args,x, ToBool(getattr(args,x)))

    if args.Simple_LogFC_format and args.Bidirectional_bars_in_LogFC_cells:
        print('Please turn on either args.Simple_LogFC_format or args.Bidirectional_bars_in_LogFC_cells. Not both')
        exit(127)

    if args.create_White_background:
        Create_workbook_Formats_uni = Create_workbook_Formats_whitespace
    else:
        Create_workbook_Formats_uni = Create_workbook_Formats

    args.max_formatted_cells = int(args.max_formatted_cells)
    args.max_depth = int(args.max_depth)

    args.Predictor_rows_count = int(args.Predictor_rows_count)

    args.cpm_cells_formatting__LogFC_max = float(args.cpm_cells_formatting__LogFC_max)
    args.cpm_cells_formatting__cpm_max = float(args.cpm_cells_formatting__cpm_max)
    args.cpm_cells_formatting__cpm_degree = float(args.cpm_cells_formatting__cpm_degree)
    args.cpm_cells_formatting__const_add = float(args.cpm_cells_formatting__const_add)
    args.cpm_cells_formatting__trimming_low = float(args.cpm_cells_formatting__trimming_low)
    if not (0 <= args.cpm_cells_formatting__trimming_low < 100):
        print('args.cpm_cells_formatting__trimming_low should be from >= 0 and < 100. Instead, it is ' + str(args.cpm_cells_formatting__trimming_low))
        exit(127)
    args.cpm_cells_formatting__trimming_high = float(args.cpm_cells_formatting__trimming_high)
    if not (0 <= args.cpm_cells_formatting__trimming_high < 100):
        print('args.cpm_cells_formatting__trimming_high should be from >= 0 and < 100. Instead, it is ' + str(args.cpm_cells_formatting__trimming_high))
        exit(127)
    if args.cpm_cells_formatting__trimming_low + args.cpm_cells_formatting__trimming_high >= 100:
        print('The sum of args.cpm_cells_formatting__trimming_low + args.cpm_cells_formatting__trimming_high should be less than 100. Instead it is %g + %g = %g'%(
            args.cpm_cells_formatting__trimming_low, args.cpm_cells_formatting__trimming_high, args.cpm_cells_formatting__trimming_low + args.cpm_cells_formatting__trimming_high))
    args.cpm_cells_formatting__gradient_color_schema = int(args.cpm_cells_formatting__gradient_color_schema)
    args.logfc_cells_formatting__gradient_color_schema = int(args.logfc_cells_formatting__gradient_color_schema)
    args.logfc_cells_formatting__gradient_range = float(args.logfc_cells_formatting__gradient_range)
    args.CPM_digits = int(args.CPM_digits)
    args.CPM_narrow_digits = int(args.CPM_narrow_digits)
    args.freeze_cols_count = int(args.freeze_cols_count)
    args.freeze_rows_count = int(args.freeze_rows_count)
    args.LogFC_cell_abs_max = float(args.LogFC_cell_abs_max)
    args.LogFC_array_abs_max = float(args.LogFC_array_abs_max)
    args.dist_matrix__digits = int(args.dist_matrix__digits)
    args.narrow_cols_width = float(args.narrow_cols_width)

    def Check_for_cond_formatting_type_and_value_consistency(a_type, a_value, type_argument_name = None, value_argument_name = None):
        if (value_argument_name is None) != (type_argument_name is None):
            raise RuntimeError('Both type_argument_name and value_argument_name should be either specified or nor specified')

        if type(a_type) is not str:
            raise RuntimeError('type(a_type) should be str. Instead, it is %s. a_type = %s'%(type(a_type), a_type))

        if a_type in ('min', 'max'):  return a_value

        try:
            a_value = float(a_value)
        except ValueError:
            if value_argument_name is not None:
                print('%s should be float. Instead, it is "%s"'%(value_argument_name, str(a_value)))
                exit(127)
            else:
                raise('a_value should be float. Instead, it is "%s"'%(str(a_value)))

        if a_type == 'percentile':
            if not (0 <= a_value < 100):
                if value_argument_name is not None:
                    print('%s should be set from 0 to 100, if %s == "percentile". Instead, it is set as %s'%(value_argument_name, type_argument_name, str(a_value)))
                    exit(127)
                else:
                    raise('a_value should be set from 0 to 100, if a_type == "percentile". Instead, it is set as %s'%(str(a_value)))

        elif a_type not in ('num'):
            if type_argument_name is not None:
                print('%s should be either min, max, num, percentile. Instead, it is set as %s'%(type_argument_name, str(a_type)))
                exit(127)
            else:
                raise('a_type should be either min, max, num, percentile. Instead, it is set as %s'%(str(a_type)))

        return a_value


    args.red_bars_max_value = Check_for_cond_formatting_type_and_value_consistency(
        args.red_bars_max_type, args.red_bars_max_value, 'args.red_bars_max_type', 'args.red_bars_max_value')
    args.blue_bars_max_value = Check_for_cond_formatting_type_and_value_consistency(
        args.blue_bars_max_type, args.blue_bars_max_value, 'args.blue_bars_max_type', 'args.blue_bars_max_value')
    args.red_bars_min_value = Check_for_cond_formatting_type_and_value_consistency(
        args.red_bars_min_type, args.red_bars_min_value, 'args.red_bars_min_type', 'args.red_bars_min_value')
    args.blue_bars_min_value = Check_for_cond_formatting_type_and_value_consistency(
        args.blue_bars_min_type, args.blue_bars_min_value, 'args.blue_bars_min_type', 'args.blue_bars_min_value')

    args.red_gradient_max_value = Check_for_cond_formatting_type_and_value_consistency(
        args.red_gradient_max_type, args.red_gradient_max_value, 'args.red_gradient_max_type', 'args.red_gradient_max_value')
    args.blue_gradient_max_value = Check_for_cond_formatting_type_and_value_consistency(
        args.blue_gradient_max_type, args.blue_gradient_max_value, 'args.blue_gradient_max_type', 'args.blue_gradient_max_value')
    args.red_gradient_min_value = Check_for_cond_formatting_type_and_value_consistency(
        args.red_gradient_min_type, args.red_gradient_min_value, 'args.red_gradient_min_type', 'args.red_gradient_min_value')
    args.blue_gradient_min_value = Check_for_cond_formatting_type_and_value_consistency(
        args.blue_gradient_min_type, args.blue_gradient_min_value, 'args.blue_gradient_min_type', 'args.blue_gradient_min_value')

    args.index_cols_min_value = Check_for_cond_formatting_type_and_value_consistency(
        args.index_cols_min_type, args.index_cols_min_value, 'args.index_cols_min_type', 'args.index_cols_min_value')
    args.index_cols_mid_value = Check_for_cond_formatting_type_and_value_consistency(
        args.index_cols_mid_type, args.index_cols_mid_value, 'args.index_cols_mid_type', 'args.index_cols_mid_value')
    args.index_cols_max_value = Check_for_cond_formatting_type_and_value_consistency(
        args.index_cols_max_type, args.index_cols_max_value, 'args.index_cols_max_type', 'args.index_cols_max_value')

    
    if len(args.cov100norm_range) != 3:
        print('Length of cov100norm_range arguments should be 3')
        exit(127)

    args.cov100norm_range = [int(x) for x in args.cov100norm_range]

    if args.max_strings not in ['auto', 'unlimited', '{auto}']:
        try:
            args.max_strings = int(args.max_strings)
        except:
            print('Incorrect args.max_strings argument')
            exit(127)

    if args.max_char_in_cell == 'unlimited':  args.max_char_in_cell = None
    if args.max_char_in_cell is not None:
        args.max_char_in_cell = int(args.max_char_in_cell)

    if args.dist_matrix_mode:
        args.Disable_ColType_autoassign = True
        args.first_col_italic = False
        args.entry_type = ''


    if not args.cpm_cells_formatting__gradient_color_schema in range(1, len(Color_gradient_schemas) + 1):
        print('Incorrect args.cpm_cells_formatting__gradient_color_schema argument. Should be from 1 to %d'%(len(Color_gradient_schemas) + 1))
        exit(127)

    if args.first_col_width != None:   args.first_col_width = int(args.first_col_width)
    
    if args.CPM_sparklines_groups == []:
        args.CPM_sparklines_groups = None
        args.CPM_sparklines_start_col = None
        
    if args.CPM_sparklines_start_col != None:
        args.CPM_sparklines_start_col = int(args.CPM_sparklines_start_col)

    if (not args.CPM_sparklines_start_col is None) and (not args.CPM_sparklines_start_col__by_sheet is None):
        print('both args.CPM_sparklines_start_col and args.CPM_sparklines_start_col__by_sheet cannot be supplied')
        exit(127)

    if (not args.CPM_sparklines_groups is None) and (not args.CPM_sparklines_groups__by_sheet is None):
        print('both args.CPM_sparklines_groups and args.CPM_sparklines_groups__by_sheet cannot be supplied')
        exit(127)

    if (args.CPM_sparklines_start_col is None or args.CPM_sparklines_start_col__by_sheet is None) != (args.CPM_sparklines_groups is None or args.CPM_sparklines_groups__by_sheet is None):
        print('args.CPM_sparklines_groups (or by-sheet) and args.CPM_sparklines_start_col (or by-sheet) should be either defined or not defined')
        exit(127)

    # if args.CPM_sparklines_start_col__by_sheet is not None:
    #     if len(args.CPM_sparklines_start_col__by_sheet) != len(args.CPM_sparklines_groups__by_sheet):
    #         print('args.CPM_sparklines_start_col__by_sheet and args.CPM_sparklines_groups__by_sheet should have equal length:')
    #         print(args.CPM_sparklines_start_col__by_sheet)
    #         print(args.CPM_sparklines_groups__by_sheet)
    #         exit(127)


    if args.CPM_sparklines_mode.casefold() not in ['relative', 'logfc', 'log', 'cpm', 'linear', 'logcpm']:
        print('Incorrect args.CPM_sparklines_mode argument (%s). Should be either linear or log'%args.CPM_sparklines_mode)
        exit(127)

    if args.CPM_sparklines_mode.casefold() in ['relative', 'linear', 'cpm']:
        args.CPM_sparklines_mode = 'linear'
    if args.CPM_sparklines_mode.casefold() in ['logfc', 'log', 'logcpm']:
        args.CPM_sparklines_mode = 'log'


    args.CPM_sparklines_col_width_factor = float(args.CPM_sparklines_col_width_factor)
    args.CPM_sparklines_rel_scale_cpm_const_add = float(args.CPM_sparklines_rel_scale_cpm_const_add)
    args.CPM_sparklines_rel_scale_cpm_trimming_low = float(args.CPM_sparklines_rel_scale_cpm_trimming_low)
    if not (0 <= args.CPM_sparklines_rel_scale_cpm_trimming_low < 100):
        print('args.CPM_sparklines_rel_scale_cpm_trimming_low should be from >= 0 and < 100. Instead, it is ' + str(args.CPM_sparklines_rel_scale_cpm_trimming_low))
        exit(127)
    args.CPM_sparklines_rel_scale_cpm_trimming_high = float(args.CPM_sparklines_rel_scale_cpm_trimming_high)
    if not (0 <= args.CPM_sparklines_rel_scale_cpm_trimming_high < 100):
        print('args.CPM_sparklines_rel_scale_cpm_trimming_high should be from >= 0 and < 100. Instead, it is ' + str(args.CPM_sparklines_rel_scale_cpm_trimming_high))
        exit(127)
    if args.CPM_sparklines_rel_scale_cpm_trimming_low + args.CPM_sparklines_rel_scale_cpm_trimming_high >= 100:
        print('The sum of args.CPM_sparklines_rel_scale_cpm_trimming_low + args.CPM_sparklines_rel_scale_cpm_trimming_high should be less than 100. Instead it is %g + %g = %g'%(
            args.CPM_sparklines_rel_scale_cpm_trimming_low, args.CPM_sparklines_rel_scale_cpm_trimming_high, args.CPM_sparklines_rel_scale_cpm_trimming_low + args.CPM_sparklines_rel_scale_cpm_trimming_high))
    args.CPM_sparklines_logfc_mode_axis_limit = float(args.CPM_sparklines_logfc_mode_axis_limit)
    args.CPM_sparklines_equal_barwidths = ToBool(args.CPM_sparklines_equal_barwidths)
    if args.CPM_sparklines_equal_barwidths:
        args.CPM_sparklines_col_width_factor = 1

    # if args.Forced_ColTypes_by_ColName != None and args.Forced_ColTypes_by_ColNumber != None:
    #     print('either Forced_ColTypes_by_ColNumber or Forced_ColTypes_by_ColName should be specified')
    #     exit(127)

    if args.Forced_ColTypes_by_ColName != None:
        args.Forced_ColTypes_by_ColName = [x.replace(args.Forced_ColTypes_by_ColName__quote_symbol, '"') for x in args.Forced_ColTypes_by_ColName]
    
    if args.Forced_ColWidths_by_ColName != None:
        args.Forced_ColWidths_by_ColName = [x.replace(args.Forced_ColWidths_by_ColName__quote_symbol, '"') for x in args.Forced_ColWidths_by_ColName]

    if args.CPM_sparklines_groups__by_sheet != None:
        args.CPM_sparklines_groups__by_sheet = [x.replace(args.CPM_sparklines__quote_symbol, '"') for x in args.CPM_sparklines_groups__by_sheet]

    if args.CPM_sparklines_start_col__by_sheet != None:
        args.CPM_sparklines_start_col__by_sheet = [x.replace(args.CPM_sparklines__quote_symbol, '"') for x in args.CPM_sparklines_start_col__by_sheet]
    
    # if args.Forced_ColWidths_by_ColName != None and args.Forced_ColWidths != None:
    #     print('both args.Forced_ColWidths and args.Forced_ColWidths_by_ColName cannot be specified')
    #     exit(127)

    if args.LogFC_sparklines_groups == []:
        args.LogFC_sparklines_groups = None
        args.LogFC_sparklines_start_col = None

    if args.LogFC_sparklines_start_col != None:
        args.LogFC_sparklines_start_col = int(args.LogFC_sparklines_start_col)


    if (not args.LogFC_sparklines_start_col is None) and (not args.LogFC_sparklines_start_col__by_sheet is None):
        print('both args.LogFC_sparklines_start_col and args.LogFC_sparklines_start_col__by_sheet cannot be supplied')
        exit(127)

    if (not args.LogFC_sparklines_groups is None) and (not args.LogFC_sparklines_groups__by_sheet is None):
        print('both args.LogFC_sparklines_groups and args.LogFC_sparklines_groups__by_sheet cannot be supplied')
        exit(127)

    if (args.LogFC_sparklines_start_col is None or args.LogFC_sparklines_start_col__by_sheet is None) != (args.LogFC_sparklines_groups is None or args.LogFC_sparklines_groups__by_sheet is None):
        print('args.LogFC_sparklines_groups (or by-sheet) and args.LogFC_sparklines_start_col (or by-sheet) should be either defined or not defined')
        exit(127)


    args.LogFC_sparklines_col_width_factor = float(args.LogFC_sparklines_col_width_factor)
    args.LogFC_sparklines_axis_limit = float(args.LogFC_sparklines_axis_limit)
    args.LogFC_sparklines_equal_barwidths = ToBool(args.LogFC_sparklines_equal_barwidths)
    if args.LogFC_sparklines_equal_barwidths:
        args.LogFC_sparklines_col_width_factor = 1

    # if args.Forced_ColTypes_by_ColName != None and args.Forced_ColTypes_by_ColNumber != None:
    #     print('either Forced_ColTypes_by_ColNumber or Forced_ColTypes_by_ColName should be specified')
    #     exit(127)


    if args.CPM_sparklines_groups__by_sheet != None:
        args.CPM_sparklines_groups__by_sheet = [x.replace(args.CPM_sparklines__quote_symbol, '"') for x in args.CPM_sparklines_groups__by_sheet]

    if args.CPM_sparklines_start_col__by_sheet != None:
        args.CPM_sparklines_start_col__by_sheet = [x.replace(args.CPM_sparklines__quote_symbol, '"') for x in args.CPM_sparklines_start_col__by_sheet]

    if args.LogFC_sparklines_groups__by_sheet != None:
        args.LogFC_sparklines_groups__by_sheet = [x.replace(args.LogFC_sparklines__quote_symbol, '"') for x in args.LogFC_sparklines_groups__by_sheet]

    if args.LogFC_sparklines_start_col__by_sheet != None:
        args.LogFC_sparklines_start_col__by_sheet = [x.replace(args.LogFC_sparklines__quote_symbol, '"') for x in args.LogFC_sparklines_start_col__by_sheet]



    allowed_forced_col_types = ['cpm', 'cpm2', 'cpm_narrow', 'cpm_array', 'logfc_array', 'logfc_array2', 'logfc', 'logfc_narrow', 'logcpm', 'cov100norm',
     'p', 'fdr', 'filter', 'score', 'spacer', 'corr', 'means', 'raw_means', 'blue_bars_custom', 'red_bars_custom',
     'blue_gradient_custom', 'red_gradient_custom',
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



    CPM_sparklines_groups = None
    CPM_sparklines_groups_names = None

    CPM_sparklines_groups__by_sheet = None
    CPM_sparklines_groups_names__by_sheet = None
    CPM_sparklines_start_col__by_sheet = None

    if args.CPM_sparklines_groups != None:
        CPM_sparklines_groups = dict()
        CPM_sparklines_groups_names = []
        # Forced_ColNumbers_by_ColType_dict = dict()
        for expr in args.CPM_sparklines_groups:
            if not ':' in expr:
                print('[13] Incorrect args.CPM_sparklines_groups argument (%s). Should looks like Target:1,2,3,5  Control:4,8,9   ...'%args.CPM_sparklines_groups)
                exit(127)
            expr = expr.split(':')
            if len(expr) != 2:
                print('[10] Incorrect args.CPM_sparklines_groups argument (%s). Should looks like Target:1,2,3,5  Control:4,8,9   ...'%args.CPM_sparklines_groups)
                exit(127)

            try:
                cols = [int(x) for x in expr[1].split(',')]
            except ValueError:
                print('[14] Incorrect args.CPM_sparklines_groups argument (%s). Should looks like Target:1,2,3,5  Control:4,8,9   ...'%args.CPM_sparklines_groups)
                exit(127)

            CPM_sparklines_groups[expr[0]] = cols
            CPM_sparklines_groups_names.append(expr[0])
        args.CPM_sparklines_groups = CPM_sparklines_groups

    elif args.CPM_sparklines_groups__by_sheet != None:
        CPM_sparklines_groups__by_sheet = dict()
        CPM_sparklines_groups_names__by_sheet = dict()
        CPM_sparklines_start_col__by_sheet = dict()

        for expr in args.CPM_sparklines_groups__by_sheet:

            'Should looks like "sheet1":"Target":1,2,3,5  "sheet1":"Control":1,2,3,5 ... '

            my_splitter = shlex.shlex(expr, posix=True)
            my_splitter.whitespace = ':'
            my_splitter.whitespace_split = True
            try:
                split_data = list(my_splitter)
            except:
                print('[15] Incorrect args.CPM_sparklines_groups__by_sheet argument (%s), part %s. Should looks like "sheet1":"Target":1,2,3,5  "sheet1":"Control":1,2,3,5 ... '%(args.CPM_sparklines_groups__by_sheet, expr))
                exit(127)

            if len(split_data) != 3:
                print('[16] Incorrect args.CPM_sparklines_groups__by_sheet argument (%s), field %s. Should looks like "sheet1":"Target":1,2,3,5  "sheet1":"Control":1,2,3,5 ... '%(args.CPM_sparklines_groups__by_sheet, expr))
                exit(127)

            sheet_name = split_data[0]
            if sheet_name not in CPM_sparklines_groups__by_sheet:
                CPM_sparklines_groups__by_sheet[sheet_name] = dict()
                CPM_sparklines_groups_names__by_sheet[sheet_name] = []

            group_name = split_data[1]
            cols_data = split_data[2]

            try:
                cols = [int(x) for x in cols_data.split(',')]
            except ValueError:
                print('[17] Incorrect args.CPM_sparklines_groups__by_sheet argument (%s), field %s. Should looks like "sheet1":"Target":1,2,3,5  "sheet1":"Control":1,2,3,5 ... '%(args.CPM_sparklines_groups__by_sheet, expr))
                exit(127)

            CPM_sparklines_groups__by_sheet[sheet_name][group_name] = cols
            CPM_sparklines_groups_names__by_sheet[sheet_name].append(group_name)

        for expr in args.CPM_sparklines_start_col__by_sheet:

            'Should looks like "sheet1":15  "sheet2":9 ... '

            my_splitter = shlex.shlex(expr, posix=True)
            my_splitter.whitespace = ':'
            my_splitter.whitespace_split = True
            try:
                split_data = list(my_splitter)
            except:
                print('[18] Incorrect args.CPM_sparklines_start_col__by_sheet argument (%s), part %s. Should looks like "sheet1":15  "sheet2":9 ... '%(args.CPM_sparklines_start_col__by_sheet, expr))
                exit(127)

            if len(split_data) != 2:
                print('[19] Incorrect args.CPM_sparklines_start_col__by_sheet argument (%s), field %s. Should looks like "sheet1":15  "sheet2":9 ... '%(args.CPM_sparklines_start_col__by_sheet, expr))
                exit(127)

            sheet_name = split_data[0]

            try:
                col = int(split_data[1])
            except ValueError:
                print('[20] Incorrect args.CPM_sparklines_start_col__by_sheet argument (%s), field %s. Should looks like "sheet1":15  "sheet2":9 ... '%(args.CPM_sparklines_start_col__by_sheet, expr))
                exit(127)

            if sheet_name in CPM_sparklines_start_col__by_sheet:
                print('[21] Duplicated values in args.CPM_sparklines_start_col__by_sheet argument (%s), field %s. Should looks like "sheet1":15  "sheet2":9 ... '%(args.CPM_sparklines_start_col__by_sheet, expr))
                exit(127)
            
            CPM_sparklines_start_col__by_sheet[sheet_name] = col




    LogFC_sparklines_groups = None
    LogFC_sparklines_groups_names = None

    LogFC_sparklines_groups__by_sheet = None
    LogFC_sparklines_groups_names__by_sheet = None
    LogFC_sparklines_start_col__by_sheet = None

    if args.LogFC_sparklines_groups != None:
        LogFC_sparklines_groups = dict()
        LogFC_sparklines_groups_names = []
        # Forced_ColNumbers_by_ColType_dict = dict()
        for expr in args.LogFC_sparklines_groups:
            if not ':' in expr:
                print('[27] Incorrect args.LogFC_sparklines_groups argument (%s). Should looks like Target:1,2,3,5  Control:4,8,9   ...'%args.LogFC_sparklines_groups)
                exit(127)
            expr = expr.split(':')
            if len(expr) != 2:
                print('[28] Incorrect args.LogFC_sparklines_groups argument (%s). Should looks like Target:1,2,3,5  Control:4,8,9   ...'%args.LogFC_sparklines_groups)
                exit(127)

            try:
                cols = [int(x) for x in expr[1].split(',')]
            except ValueError:
                print('[29] Incorrect args.LogFC_sparklines_groups argument (%s). Should looks like Target:1,2,3,5  Control:4,8,9   ...'%args.LogFC_sparklines_groups)
                exit(127)

            LogFC_sparklines_groups[expr[0]] = cols
            LogFC_sparklines_groups_names.append(expr[0])
        args.LogFC_sparklines_groups = LogFC_sparklines_groups

    elif args.LogFC_sparklines_groups__by_sheet != None:
        LogFC_sparklines_groups__by_sheet = dict()
        LogFC_sparklines_groups_names__by_sheet = dict()
        LogFC_sparklines_start_col__by_sheet = dict()

        for expr in args.LogFC_sparklines_groups__by_sheet:

            'Should looks like "sheet1":"Target":1,2,3,5  "sheet1":"Control":1,2,3,5 ... '

            my_splitter = shlex.shlex(expr, posix=True)
            my_splitter.whitespace = ':'
            my_splitter.whitespace_split = True
            try:
                split_data = list(my_splitter)
            except:
                print('[30] Incorrect args.LogFC_sparklines_groups__by_sheet argument (%s), part %s. Should looks like "sheet1":"Target":1,2,3,5  "sheet1":"Control":1,2,3,5 ... '%(args.LogFC_sparklines_groups__by_sheet, expr))
                exit(127)

            if len(split_data) != 3:
                print('[31] Incorrect args.LogFC_sparklines_groups__by_sheet argument (%s), field %s. Should looks like "sheet1":"Target":1,2,3,5  "sheet1":"Control":1,2,3,5 ... '%(args.LogFC_sparklines_groups__by_sheet, expr))
                exit(127)

            sheet_name = split_data[0]
            if sheet_name not in LogFC_sparklines_groups__by_sheet:
                LogFC_sparklines_groups__by_sheet[sheet_name] = dict()
                LogFC_sparklines_groups_names__by_sheet[sheet_name] = []

            group_name = split_data[1]
            cols_data = split_data[2]

            try:
                cols = [int(x) for x in cols_data.split(',')]
            except ValueError:
                print('[32] Incorrect args.LogFC_sparklines_groups__by_sheet argument (%s), field %s. Should looks like "sheet1":"Target":1,2,3,5  "sheet1":"Control":1,2,3,5 ... '%(args.LogFC_sparklines_groups__by_sheet, expr))
                exit(127)

            LogFC_sparklines_groups__by_sheet[sheet_name][group_name] = cols
            LogFC_sparklines_groups_names__by_sheet[sheet_name].append(group_name)

        for expr in args.LogFC_sparklines_start_col__by_sheet:

            'Should looks like "sheet1":15  "sheet2":9 ... '

            my_splitter = shlex.shlex(expr, posix=True)
            my_splitter.whitespace = ':'
            my_splitter.whitespace_split = True
            try:
                split_data = list(my_splitter)
            except:
                print('[33] Incorrect args.LogFC_sparklines_start_col__by_sheet argument (%s), part %s. Should looks like "sheet1":15  "sheet2":9 ... '%(args.LogFC_sparklines_start_col__by_sheet, expr))
                exit(127)

            if len(split_data) != 2:
                print('[34] Incorrect args.LogFC_sparklines_start_col__by_sheet argument (%s), field %s. Should looks like "sheet1":15  "sheet2":9 ... '%(args.LogFC_sparklines_start_col__by_sheet, expr))
                exit(127)

            sheet_name = split_data[0]

            try:
                col = int(split_data[1])
            except ValueError:
                print('[35] Incorrect args.LogFC_sparklines_start_col__by_sheet argument (%s), field %s. Should looks like "sheet1":15  "sheet2":9 ... '%(args.LogFC_sparklines_start_col__by_sheet, expr))
                exit(127)

            if sheet_name in LogFC_sparklines_start_col__by_sheet:
                print('[36] Duplicated values in args.LogFC_sparklines_start_col__by_sheet argument (%s), field %s. Should looks like "sheet1":15  "sheet2":9 ... '%(args.LogFC_sparklines_start_col__by_sheet, expr))
                exit(127)
            
            LogFC_sparklines_start_col__by_sheet[sheet_name] = col




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

    if args.verbose:  print('Total %d src-files were found'%(len(all_src_FileNames)))

    if args.Sheet_names is not None:
        if len(args.Sheet_names) != len(all_src_FileNames):
            print('Count of the forced sheet names (%d) does not correspond to the number of input files (%d)'%(len(args.Sheet_names), len(all_src_FileNames)))
            exit(127)

    if args.Generate_SingleBook:
        if args.Workbook_FileName is None:
            print('Excel workbook filename should be specified when Generate_SingleBook==True')
            exit(127)

        Workbook = xlsxwriter.Workbook(args.Workbook_FileName)
        existing_CPM_formats = dict()
        # print('HHHHHHHH')
        (bold, bold_light, bold_ww, italic, bold_italic, bold_italic_ww, format0, format1, format2, format1_narrow, format2_narrow, format1_center, format2_center,
            format2_center_bold, format_center, format0_center, format_perc, default_Rank_format, headers_format, small_headers_format,
            secondary_headers_format, merge_format, rotated_format, rotated_format_underline, rotated_format_center, whitespace_format,
            whitespace_format__right_border, whitespace_format__bottom_border,
            whitespace_format1, percentage_format0, whitespace_format__header, empty_logfc_format, sparklines_cell_format, sparklines_cell_format_no_bg,
            sparklines_cell_format_white_bg, pvalue_formats, format_quasi_invisible,
            format_quasi_invisible_left_border, format_quasi_invisible_right_border, format_quasi_invisible_left_right_borders,
            heatmap_cell_format_left_border, heatmap_cell_format_right_border, heatmap_cell_format_left_right_borders) = Create_workbook_Formats_uni(Workbook)



    file_N = 0
    trunc_sheet_number = 1
    # print('\n')

    for FN in all_src_FileNames:
        file_N += 1
        if args.verbose:  sys.stdout.write('\rProcessing file %d of %d...'%(file_N,len(all_src_FileNames)))
        if not os.path.exists(FN):
            print('File %s DOES NOT exist'%FN)
            if args.Skip_absent_files:  continue
            else:   exit(127)
            
        f = open(FN,'r')
        header_cells_count = f.readline().count('\t') + 1
        real_cells_count = f.readline().count('\t') + 1
        if header_cells_count == real_cells_count - 1:  row_names_is_present = True
        elif header_cells_count == real_cells_count:  row_names_is_present = False
        else:
            if real_cells_count == 1 or real_cells_count == 0:  continue
            print('file %s: Incorrect header and body cells count (%d and %d)'%(FN, header_cells_count, real_cells_count));  exit(127)
        f.close()

        if args.max_strings == 'auto' or args.max_strings == '{auto}':
            current_max_strings = 200 * 16000 / header_cells_count
        elif args.max_strings == 'unlimited':
            current_max_strings = 100000000000000
        else:
            current_max_strings = args.max_strings

        f = open(FN,'r')
        header_cells = f.readline().replace('\r','').replace('\n','').split('\t')
        header_cells = [C_el[1:-1] if ((C_el.startswith('"') and C_el.endswith('"')) or (C_el.startswith("'") and C_el.endswith("'"))) else C_el for C_el in header_cells]

        if row_names_is_present:
            if args.entry_type == None:  header_cells = ['entry'] + header_cells
            else:  header_cells = [args.entry_type] + header_cells

        Cols_count = len(header_cells)

        header_cells_cf = [x.casefold().replace(',',' ').replace(';',' ').replace('.',' ') for x in header_cells]

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

        if 'cov100norm' in Forced_ColNumbers_by_ColType_dict:   cov100norm_cols = Forced_ColNumbers_by_ColType_dict['cov100norm']
        elif 'cov100norm' in Forced_ColNames_by_ColTypes_dict: cov100norm_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['cov100norm'])]
        elif not args.Disable_ColType_autoassign:    cov100norm_cols = [x for x in range(len(header_cells)) if any([y == 'cov100norm' for y in header_cells_cf[x].split(' ')]) and x not in all_forced_cols]
        else:  cov100norm_cols = []

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

        if 'blue_bars_custom' in Forced_ColNumbers_by_ColType_dict:   blue_bars_custom_cols = Forced_ColNumbers_by_ColType_dict['blue_bars_custom']
        elif 'blue_bars_custom' in Forced_ColNames_by_ColTypes_dict: blue_bars_custom_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['blue_bars_custom'])]
        else:  blue_bars_custom_cols = []        

        if 'red_bars_custom' in Forced_ColNumbers_by_ColType_dict:   red_bars_custom_cols = Forced_ColNumbers_by_ColType_dict['red_bars_custom']
        elif 'red_bars_custom' in Forced_ColNames_by_ColTypes_dict: red_bars_custom_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['red_bars_custom'])]
        else:  red_bars_custom_cols = []        

        if 'blue_gradient_custom' in Forced_ColNumbers_by_ColType_dict:   blue_gradient_custom_cols = Forced_ColNumbers_by_ColType_dict['blue_gradient_custom']
        elif 'blue_gradient_custom' in Forced_ColNames_by_ColTypes_dict: blue_gradient_custom_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['blue_gradient_custom'])]
        else:  blue_gradient_custom_cols = []        

        if 'red_gradient_custom' in Forced_ColNumbers_by_ColType_dict:   red_gradient_custom_cols = Forced_ColNumbers_by_ColType_dict['red_gradient_custom']
        elif 'red_gradient_custom' in Forced_ColNames_by_ColTypes_dict: red_gradient_custom_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['red_gradient_custom'])]
        else:  red_gradient_custom_cols = []        

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

        if 'cpm_narrow' in Forced_ColNumbers_by_ColType_dict:  CPM_narrow_cols = Forced_ColNumbers_by_ColType_dict['cpm_narrow']
        elif 'cpm_narrow' in Forced_ColNames_by_ColTypes_dict:  CPM_narrow_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['cpm_narrow'])]
        else:  CPM_narrow_cols = []        

        if 'logfc_array' in Forced_ColNumbers_by_ColType_dict:    LogFC_array_cols = Forced_ColNumbers_by_ColType_dict['logfc_array']
        elif 'logfc_array' in Forced_ColNames_by_ColTypes_dict:  LogFC_array_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['logfc_array'])]
        else:  LogFC_array_cols = []

        if 'logfc_array2' in Forced_ColNumbers_by_ColType_dict:    LogFC_array2_cols = Forced_ColNumbers_by_ColType_dict['logfc_array2']
        elif 'logfc_array2' in Forced_ColNames_by_ColTypes_dict:  LogFC_array2_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['logfc_array2'])]
        else:  LogFC_array2_cols = []


        if args.Metagenome_mode:
            if 'means' in Forced_ColNumbers_by_ColType_dict:   means_CPM_cols = Forced_ColNumbers_by_ColType_dict['means']
            elif 'means' in Forced_ColNames_by_ColTypes_dict: means_CPM_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['means'])]
            elif not args.Disable_ColType_autoassign:   means_CPM_cols = [x for x in range(len(header_cells)) if x not in all_forced_cols and
                        ('group'in header_cells_cf[x].split(' ') and 'mean'in header_cells_cf[x].split(' ') and 'raw mean' not in header_cells_cf[x] ) and not any(
                            [y == 'logcpm' for y in header_cells_cf[x].split(' ')])]
            else:  means_CPM_cols = []

            if 'raw_means' in Forced_ColNumbers_by_ColType_dict:    raw_means_CPM_cols = Forced_ColNumbers_by_ColType_dict['raw_means']
            elif 'raw_means' in Forced_ColNames_by_ColTypes_dict: raw_means_CPM_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['raw_means'])]
            elif not args.Disable_ColType_autoassign:   raw_means_CPM_cols = [x for x in range(len(header_cells)) if x not in all_forced_cols and
                        ('group'in header_cells_cf[x].split(' ') and 'raw mean'in header_cells_cf[x]) and not any(
                            [y == 'logcpm' for y in header_cells_cf[x].split(' ')])]
            else:  raw_means_CPM_cols = []

            if len(means_CPM_cols) != len(raw_means_CPM_cols) and len(raw_means_CPM_cols) > 0:
                print('File %s: the number of "group means" and "group raw means" columns is not equal: %d and %d'%(FN, len(means_CPM_cols), len(raw_means_CPM_cols)))
                exit(131)

            last_stat_col = max([0] + Correlation_r_cols + FDR_cols + P_value_cols + logCPM_cols + Score_cols)

            if 'cpm_array' in Forced_ColNumbers_by_ColType_dict:  CPM_array_cols = Forced_ColNumbers_by_ColType_dict['cpm_array']
            elif 'cpm_array' in Forced_ColNames_by_ColTypes_dict: CPM_array_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['cpm_array'])]
            elif not args.Disable_ColType_autoassign:  CPM_array_cols = list(range(last_stat_col + 1, len(header_cells)))
            else:  CPM_array_cols = []

            if 'cpm' in Forced_ColNumbers_by_ColType_dict:    CPM_cols = Forced_ColNumbers_by_ColType_dict['cpm']
            elif 'cpm' in Forced_ColNames_by_ColTypes_dict: CPM_cols = [x for x in range(len(header_cells)) if (header_cells[x] in Forced_ColNames_by_ColTypes_dict['cpm'])]
            else:   CPM_cols = means_CPM_cols + CPM_array_cols

        else:
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


        if not args.Generate_SingleBook:
            Workbook = xlsxwriter.Workbook(args.Workbook_FileName)
            existing_CPM_formats = dict()
            # GO_Sheet = Workbook.add_worksheet('GO-centric DE profiles')
            (bold, bold_light, bold_ww, italic, bold_italic, bold_italic_ww, format0, format1, format2, format1_narrow, format2_narrow, format1_center, format2_center,
                format2_center_bold, format_center, format0_center, format_perc, default_Rank_format, headers_format, small_headers_format,
                secondary_headers_format, merge_format, rotated_format, rotated_format_underline, rotated_format_center, whitespace_format,
                whitespace_format__right_border, whitespace_format__bottom_border,
                whitespace_format1, percentage_format0, whitespace_format__header, empty_logfc_format, sparklines_cell_format, sparklines_cell_format_no_bg,
                sparklines_cell_format_white_bg, pvalue_formats, format_quasi_invisible,
                format_quasi_invisible_left_border, format_quasi_invisible_right_border, format_quasi_invisible_left_right_borders,
                heatmap_cell_format_left_border, heatmap_cell_format_right_border, heatmap_cell_format_left_right_borders) = Create_workbook_Formats_uni(Workbook)


        if args.Sheet_names is None:
            sheet_name_src = os.path.split(FN)[-1]
        else:
            sheet_name_src = args.Sheet_names[file_N - 1]

        sheet_name = sheet_name_src.replace('[', '(').replace(']', ')').replace('?', '_').replace(':', '(d)').replace('*','(m)').replace('\\', '_').replace('/', '_')

        CPM_sparklines_groups = None
        if args.CPM_sparklines_groups__by_sheet is not None:
            if sheet_name in CPM_sparklines_groups__by_sheet:
                CPM_sparklines_groups = CPM_sparklines_groups__by_sheet[sheet_name]
                CPM_sparklines_groups_names = CPM_sparklines_groups_names__by_sheet[sheet_name]
            elif sheet_name_src in CPM_sparklines_groups__by_sheet:
                CPM_sparklines_groups = CPM_sparklines_groups__by_sheet[sheet_name_src]
                CPM_sparklines_groups_names = CPM_sparklines_groups_names__by_sheet[sheet_name_src]
            elif os.path.split(FN)[-1] in CPM_sparklines_groups__by_sheet:
                CPM_sparklines_groups = CPM_sparklines_groups__by_sheet[os.path.split(FN)[-1]]
                CPM_sparklines_groups_names = CPM_sparklines_groups_names__by_sheet[os.path.split(FN)[-1]]
            elif FN in CPM_sparklines_groups__by_sheet:
                CPM_sparklines_groups = CPM_sparklines_groups__by_sheet[FN]
                CPM_sparklines_groups_names = CPM_sparklines_groups_names__by_sheet[FN]
            else:
                print('[22] Cannot find neither sheet name "%s" nor file name "%s" in argument CPM_sparklines_groups__by_sheet. Available sheet names: %s'%(sheet_name, FN, sorted(list(CPM_sparklines_groups__by_sheet.keys()))))
                exit(127)
            args.CPM_sparklines_groups = CPM_sparklines_groups
        
        elif args.CPM_sparklines_groups is not None:
            CPM_sparklines_groups = args.CPM_sparklines_groups



        CPM_sparklines_start_col = None
        if args.CPM_sparklines_start_col__by_sheet is not None:
            if sheet_name in CPM_sparklines_start_col__by_sheet:
                CPM_sparklines_start_col = CPM_sparklines_start_col__by_sheet[sheet_name]
            elif sheet_name_src in CPM_sparklines_start_col__by_sheet:
                CPM_sparklines_start_col = CPM_sparklines_start_col__by_sheet[sheet_name_src]
            elif os.path.split(FN)[-1] in CPM_sparklines_start_col__by_sheet:
                CPM_sparklines_start_col = CPM_sparklines_start_col__by_sheet[os.path.split(FN)[-1]]
            elif FN in CPM_sparklines_start_col__by_sheet:
                CPM_sparklines_start_col = CPM_sparklines_start_col__by_sheet[FN]
            else:
                print('[23] Cannot find neither sheet name "%s" nor file name "%s" in argument CPM_sparklines_start_col__by_sheet. Available sheet names: %s'%(sheet_name, FN, sorted(list(CPM_sparklines_start_col__by_sheet.keys()))))
                exit(127)
            args.CPM_sparklines_start_col = CPM_sparklines_start_col
        
        elif args.CPM_sparklines_start_col is not None:
            CPM_sparklines_start_col = args.CPM_sparklines_start_col

        if CPM_sparklines_start_col is not None:
            if args.LogFC_sparklines_groups is not None:
                CPM_sparklines_start_col += len(args.LogFC_sparklines_groups)
                if args.insert_empty_col_before_LogFC_sparklines:   CPM_sparklines_start_col += 1


        LogFC_sparklines_groups = None
        if args.LogFC_sparklines_groups__by_sheet is not None:
            if sheet_name in LogFC_sparklines_groups__by_sheet:
                LogFC_sparklines_groups = LogFC_sparklines_groups__by_sheet[sheet_name]
                LogFC_sparklines_groups_names = LogFC_sparklines_groups_names__by_sheet[sheet_name]
            elif sheet_name_src in LogFC_sparklines_groups__by_sheet:
                LogFC_sparklines_groups = LogFC_sparklines_groups__by_sheet[sheet_name_src]
                LogFC_sparklines_groups_names = LogFC_sparklines_groups_names__by_sheet[sheet_name_src]
            elif os.path.split(FN)[-1] in LogFC_sparklines_groups__by_sheet:
                LogFC_sparklines_groups = LogFC_sparklines_groups__by_sheet[os.path.split(FN)[-1]]
                LogFC_sparklines_groups_names = LogFC_sparklines_groups_names__by_sheet[os.path.split(FN)[-1]]
            elif FN in LogFC_sparklines_groups__by_sheet:
                LogFC_sparklines_groups = LogFC_sparklines_groups__by_sheet[FN]
                LogFC_sparklines_groups_names = LogFC_sparklines_groups_names__by_sheet[FN]
            else:
                print('[37] Cannot find neither sheet name "%s" nor file name "%s" in argument LogFC_sparklines_groups__by_sheet. Available sheet names: %s'%(sheet_name, FN, sorted(list(LogFC_sparklines_groups__by_sheet.keys()))))
                exit(127)
            args.LogFC_sparklines_groups = LogFC_sparklines_groups
        
        elif args.LogFC_sparklines_groups is not None:
            LogFC_sparklines_groups = args.LogFC_sparklines_groups



        LogFC_sparklines_start_col = None
        if args.LogFC_sparklines_start_col__by_sheet is not None:
            if sheet_name in LogFC_sparklines_start_col__by_sheet:
                LogFC_sparklines_start_col = LogFC_sparklines_start_col__by_sheet[sheet_name]
            elif sheet_name_src in LogFC_sparklines_start_col__by_sheet:
                LogFC_sparklines_start_col = LogFC_sparklines_start_col__by_sheet[sheet_name_src]
            elif os.path.split(FN)[-1] in LogFC_sparklines_start_col__by_sheet:
                LogFC_sparklines_start_col = LogFC_sparklines_start_col__by_sheet[os.path.split(FN)[-1]]
            elif FN in LogFC_sparklines_start_col__by_sheet:
                LogFC_sparklines_start_col = LogFC_sparklines_start_col__by_sheet[FN]
            else:
                print('[38] Cannot find neither sheet name "%s" nor file name "%s" in argument LogFC_sparklines_start_col__by_sheet. Available sheet names: %s'%(sheet_name, FN, sorted(list(LogFC_sparklines_start_col__by_sheet.keys()))))
                exit(127)
            args.LogFC_sparklines_start_col = LogFC_sparklines_start_col
        
        elif args.LogFC_sparklines_start_col is not None:
            LogFC_sparklines_start_col = args.LogFC_sparklines_start_col




        if len(sheet_name) > 30:
            sheet_name = sheet_name[:25] + '...%d'%trunc_sheet_number
            if args.Generate_SingleBook:
                trunc_sheet_number += 1

        all_Lines = f.readlines()

        sheet = Workbook.add_worksheet(sheet_name)

        # if args.create_White_background:
        #     for RowN in range(args.Predictor_rows_count + 2 + (len(all_Lines) - args.Predictor_rows_count + 1) * (1 + args.Insert_spasers_between_rows) ):
        #         sheet.set_row(RowN, 14.3)
                # for ColN in range(len(header_cells)*2 + 1):
                #     sheet.write(RowN, ColN, '', whitespace_format)

        for x in range(len(header_cells)):
            if header_cells[x] in ColFullNames_dict:
                cell_full_name = ColFullNames_dict[header_cells[x]]
            else:
                cell_full_name = header_cells[x]

            if x in Spacer_cols: continue
            elif args.dist_matrix_mode and x > 0:
                sheet.write(0, WSparkLine_Offset(args, x), cell_full_name, rotated_format)
                sheet.set_column(WSparkLine_Offset(args, x), WSparkLine_Offset(args, x), 3.5)
            elif x in logFC_narrow_cols or x in CPM_narrow_cols:
                sheet.write(0, WSparkLine_Offset(args, x), cell_full_name, rotated_format)
                sheet.set_column(WSparkLine_Offset(args, x), WSparkLine_Offset(args, x), args.narrow_cols_width)
            elif x in raw_means_CPM_cols:
                sheet.write(0,WSparkLine_Offset(args, x),'raw mean',secondary_headers_format)
            elif x in Rank_1_cols or x in Rank_2_cols:
                sheet.write(0,WSparkLine_Offset(args, x), cell_full_name, small_headers_format)
            else:
                sheet.write(0,WSparkLine_Offset(args, x), cell_full_name, headers_format)
        
        RowN = 1


        max_samples_count = 0
        for L in all_Lines[:args.Predictor_rows_count]:
            C = L.replace('\n','').replace('\r','').split('\t')
            C = [C_el[1:-1] if ((C_el.startswith('"') and C_el.endswith('"')) or (C_el.startswith("'") and C_el.endswith("'"))) else C_el for C_el in C]
            if args.Predictor_light_font:
                sheet.write(RowN, WSparkLine_Offset(args, 0), C[0], bold_light)
            else:
                sheet.write(RowN, WSparkLine_Offset(args, 0), C[0], bold)
            for x in range(1, len(C)):
                write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], bypass_if_NaN = False)
            max_samples_count = max(max_samples_count, len(C) - 1)
            RowN += 1


        #             sheet.add_sparkline(1 + pred_N, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines,
        #                                 {'location': xl_range(1 + pred_N, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines,
        #                                                       1 + pred_N, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines),
        #                                  'range': xl_range(1 + pred_N, , 1 + pred_N, ),
        #                                  'type': 'column', 'max': 100, 'min': 0,
        #                                  'series_color': Spark_color})

        CPM_adj_starts_ColN = []

        if CPM_sparklines_groups is not None:
            if args.include_CPM_heatmaps:
                # print(args.CPM_sparklines_start_col)
                last_CPM_adj_end_col = WSparkLine_Offset(args, args.CPM_sparklines_start_col, for_CPM_heatmap_cells = True) - 2
            else:
                last_CPM_adj_end_col = WSparkLine_Offset(args, Cols_count)
            
            for x in range(len(CPM_sparklines_groups_names)):
                Gn = CPM_sparklines_groups_names[x]
                CPM_profile_cols = CPM_sparklines_groups[Gn]
                sheet.write(0, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines, Gn, headers_format)
                CPM_adj_starts_ColN.append(last_CPM_adj_end_col + 2)

                out_of_range_cols = [x for x in CPM_profile_cols if x >= len(header_cells)]
                if len(out_of_range_cols) > 0:
                    print('Error. Sheet %s. Sparkline CPM group %s. The following columns are out of range: '%(sheet_name, Gn) + str(out_of_range_cols))
                    exit(127)

                if args.include_CPM_heatmaps:
                    sheet.merge_range(0, CPM_adj_starts_ColN[-1], 0, CPM_adj_starts_ColN[-1] + len(CPM_profile_cols) - 1, CPM_sparklines_groups_names[x], headers_format)
                else:
                    for col_n in range(len(CPM_profile_cols)):
                        if args.CPM_sparklines_mode == 'linear':
                            sheet.write(0, CPM_adj_starts_ColN[-1] + col_n, '(n) %s'%(header_cells[CPM_profile_cols[col_n]]), secondary_headers_format)
                        elif args.CPM_sparklines_mode == 'log':
                            sheet.write(0, CPM_adj_starts_ColN[-1] + col_n, '(log-R) %s'%(header_cells[CPM_profile_cols[col_n]]), secondary_headers_format)

                last_CPM_adj_end_col = CPM_adj_starts_ColN[-1] + len(CPM_profile_cols) - 1

            Cols_count = last_CPM_adj_end_col + 1


        LogFC_adj_starts_ColN = []
        if LogFC_sparklines_groups is not None:
            if args.include_LogFC_heatmaps:
                last_LogFC_adj_end_col = WSparkLine_Offset(args, args.LogFC_sparklines_start_col, for_LogFC_heatmap_cells = True) - 2
            else:
                last_LogFC_adj_end_col = WSparkLine_Offset(args, Cols_count)

            for x in range(len(LogFC_sparklines_groups_names)):
                Gn = LogFC_sparklines_groups_names[x]
                LogFC_profile_cols = LogFC_sparklines_groups[Gn]
                sheet.write(0, LogFC_sparklines_start_col - 1 + x + args.insert_empty_col_before_LogFC_sparklines, Gn, headers_format)
                LogFC_adj_starts_ColN.append(last_LogFC_adj_end_col + 2)

                out_of_range_cols = [x for x in LogFC_profile_cols if x >= len(header_cells)]
                if len(out_of_range_cols) > 0:
                    print('Error. Sheet %s. Sparkline LogFC group %s. The following columns are out of range: '%(sheet_name, Gn) + str(out_of_range_cols))
                    exit(127)

                if args.include_LogFC_heatmaps:
                    sheet.merge_range(0, LogFC_adj_starts_ColN[-1], 0, LogFC_adj_starts_ColN[-1] + len(LogFC_profile_cols) - 1, LogFC_sparklines_groups_names[x], headers_format)
                else:
                    for col_n in range(len(LogFC_profile_cols)):
                        sheet.write(0, LogFC_adj_starts_ColN[-1] + col_n, '(c) %s'%(header_cells[LogFC_profile_cols[col_n]]), secondary_headers_format)

                last_LogFC_adj_end_col = LogFC_adj_starts_ColN[-1] + len(LogFC_profile_cols) - 1

            Cols_count += last_LogFC_adj_end_col + 1


        ### sparklines for predictors - CPM sparklines
        if args.Predictor_rows_count > 0 and args.add_predictor_sparklines and max_samples_count > 0 and CPM_sparklines_groups_names != None:
            absent_predictor_marker_rel_size = 0.25
            for pred_N in range(args.Predictor_rows_count):
                C = all_Lines[pred_N].replace('\n','').replace('\r','').split('\t')
                C = [C_el[1:-1] if ((C_el.startswith('"') and C_el.endswith('"')) or (C_el.startswith("'") and C_el.endswith("'"))) else C_el for C_el in C]

                all_values = []
                all_values__groups = []
                for x in range(len(CPM_sparklines_groups_names)):
                    Gn = CPM_sparklines_groups_names[x]
                    CPM_profile_cols = CPM_sparklines_groups[Gn]
                    #sheet.write(0, CPM_sparklines_start_col - 1 + x, Gn, headers_format)

                    all_values__groups.append([len(all_values), len(all_values) + len(CPM_profile_cols) - 1])
                    for col_n in range(len(CPM_profile_cols)):
                        try:  all_values.append(float(C[CPM_profile_cols[col_n]]))
                        except ValueError:  all_values.append(float('NaN'))

                all_values_without_NaN = [x for x in all_values if not math.isnan(x)]
                if len(all_values_without_NaN) == 0:  continue
                max_value = max(all_values_without_NaN)
                negative_values_present = any([x < 0 for x in all_values])
                if negative_values_present:
                    min_value = min(all_values_without_NaN)
                    if max_value == min_value:  continue
                    all_values_with_markers = [x if not math.isnan(x) else 0 for x in all_values]
                else:
                    if max_value == 0:  max_value = 1
                    min_value = (-1) * max_value * absent_predictor_marker_rel_size
                    all_values_with_markers = [x if not math.isnan(x) else min_value for x in all_values]

                
                for x in range(len(CPM_sparklines_groups_names)):
                    Gn = CPM_sparklines_groups_names[x]
                    CPM_profile_cols = CPM_sparklines_groups[Gn]

                    if len(CPM_profile_cols) != all_values__groups[x][1] - all_values__groups[x][0] + 1:
                        print('Smth strange (3) %d != %d - %d + 1....'%(len(CPM_profile_cols), all_values__groups[x][0], all_values__groups[x][1]))
                        exit(127)

                    for col_n in range(len(CPM_profile_cols)):
                        if args.include_CPM_heatmaps:
                            # heatmap_cell_format_left_border, heatmap_cell_format_right_border, heatmap_cell_format_left_right_borders
                            if len(CPM_profile_cols) == 1:
                                write_number__mod(sheet, 1 + pred_N, CPM_adj_starts_ColN[x] + col_n, all_values_with_markers[all_values__groups[x][0] + col_n], format_quasi_invisible_left_right_borders)
                            elif col_n == 0:
                                write_number__mod(sheet, 1 + pred_N, CPM_adj_starts_ColN[x] + col_n, all_values_with_markers[all_values__groups[x][0] + col_n], format_quasi_invisible_left_border)
                            elif col_n == len(CPM_profile_cols) - 1:
                                write_number__mod(sheet, 1 + pred_N, CPM_adj_starts_ColN[x] + col_n, all_values_with_markers[all_values__groups[x][0] + col_n], format_quasi_invisible_right_border)
                            else:
                                write_number__mod(sheet, 1 + pred_N, CPM_adj_starts_ColN[x] + col_n, all_values_with_markers[all_values__groups[x][0] + col_n], format_quasi_invisible)
                        else:
                            write_number__mod(sheet, 1 + pred_N, CPM_adj_starts_ColN[x] + col_n, all_values_with_markers[all_values__groups[x][0] + col_n], format2)

                if args.Predictor_rows_count > 1:
                    HM_max_color = color(Gradient(HM_Preds_Gradient, pred_N / (args.Predictor_rows_count - 1)))
                else:  HM_max_color = color(Gradient(HM_Preds_Gradient, 1))

                NaNs_present = any([math.isnan(x) for x in all_values])
                if NaNs_present and not negative_values_present:
                    sheet.conditional_format(1 + pred_N, CPM_adj_starts_ColN[x], 1 + pred_N, CPM_adj_starts_ColN[-1] + len(CPM_sparklines_groups[CPM_sparklines_groups_names[-1]]) - 1,
                       {'type': '3_color_scale',
                        'min_color': "#aaaaaa ", 'mid_color': "#ffffff",
                        'max_color': HM_max_color,
                        'min_type': 'num', 'mid_type': 'num',
                        'max_type': 'num',
                        'min_value': min_value, 'mid_value': min(all_values_without_NaN),
                        'max_value': max(all_values_without_NaN)})
                else:
                    sheet.conditional_format(1 + pred_N, CPM_adj_starts_ColN[x], 1 + pred_N, CPM_adj_starts_ColN[-1] + len(CPM_sparklines_groups[CPM_sparklines_groups_names[-1]]) - 1,
                       {'type': '3_color_scale',
                        'min_color': "#ffffff",
                        'max_color': HM_max_color,
                        'min_type': 'num',
                        'max_type': 'num',
                        'min_value': min(all_values_without_NaN),
                        'max_value': max(all_values_without_NaN)})


                if args.insert_empty_col_before_CPM_sparklines:
                    sheet.write(1 + pred_N, CPM_sparklines_start_col - 1, ' ')


                for x in range(len(CPM_sparklines_groups_names)):
                    Gn = CPM_sparklines_groups_names[x]
                    CPM_profile_cols = CPM_sparklines_groups[Gn]

                    start_ColN = CPM_adj_starts_ColN[x]
                    end_ColN = CPM_adj_starts_ColN[x] + len(CPM_sparklines_groups[Gn]) - 1
                    # all_SL_Locations = [xl_range(c, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines,
                    #                              c, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines) for c in range(args.Predictor_rows_count + 2, Last_Row)]
                    # all_SL_Ranges = [xl_range(c, start_ColN, c, end_ColN) for c in range(args.Predictor_rows_count + 2, Last_Row)]

                    if args.Predictor_rows_count > 1:
                        Spark_color = color(Gradient(Sparklines_Preds_Gradient, pred_N / (args.Predictor_rows_count - 1)))
                    else:  Spark_color = color(Gradient(Sparklines_Preds_Gradient, 1))

                    sheet.add_sparkline(1 + pred_N, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines,
                                        {'location': xl_range(1 + pred_N, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines,
                                                              1 + pred_N, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines),
                                         'range': xl_range(1 + pred_N, start_ColN, 1 + pred_N, end_ColN),
                                         'type': 'column', 'max': max_value, 'min': min_value,
                                         'negative_points': True,
                                         'series_color': Spark_color,
                                         'negative_color': ('#919191' if negative_values_present else '#d1d1d1')})




        ### sparklines for predictors - LogFC sparklines
        if args.Predictor_rows_count > 0 and args.add_predictor_sparklines and max_samples_count > 0 and LogFC_sparklines_groups_names != None:
            absent_predictor_marker_rel_size = 0.25
            for pred_N in range(args.Predictor_rows_count):
                C = all_Lines[pred_N].replace('\n','').replace('\r','').split('\t')
                C = [C_el[1:-1] if ((C_el.startswith('"') and C_el.endswith('"')) or (C_el.startswith("'") and C_el.endswith("'"))) else C_el for C_el in C]

                all_values = []
                all_values__groups = []
                for x in range(len(LogFC_sparklines_groups_names)):
                    Gn = LogFC_sparklines_groups_names[x]
                    LogFC_profile_cols = LogFC_sparklines_groups[Gn]
                    #sheet.write(0, LogFC_sparklines_start_col - 1 + x, Gn, headers_format)

                    all_values__groups.append([len(all_values), len(all_values) + len(LogFC_profile_cols) - 1])
                    for col_n in range(len(LogFC_profile_cols)):
                        try:  all_values.append(float(C[LogFC_profile_cols[col_n]]))
                        except ValueError:  all_values.append(float('NaN'))

                all_values_without_NaN = [x for x in all_values if not math.isnan(x)]
                if len(all_values_without_NaN) == 0:  continue
                max_value = max(all_values_without_NaN)
                negative_values_present = any([x < 0 for x in all_values])
                if negative_values_present:
                    min_value = min(all_values_without_NaN)
                    if max_value == min_value:  continue
                    all_values_with_markers = [x if not math.isnan(x) else 0 for x in all_values]
                else:
                    if max_value == 0:  max_value = 1
                    min_value = (-1) * max_value * absent_predictor_marker_rel_size
                    all_values_with_markers = [x if not math.isnan(x) else min_value for x in all_values]

                for x in range(len(LogFC_sparklines_groups_names)):
                    Gn = LogFC_sparklines_groups_names[x]
                    LogFC_profile_cols = LogFC_sparklines_groups[Gn]

                    if len(LogFC_profile_cols) != all_values__groups[x][1] - all_values__groups[x][0] + 1:
                        print('Smth strange (3) %d != %d - %d + 1....'%(len(LogFC_profile_cols), all_values__groups[x][0], all_values__groups[x][1]))
                        exit(127)

                    for col_n in range(len(LogFC_profile_cols)):
                        if args.include_LogFC_heatmaps:
                            if len(LogFC_profile_cols) == 1:
                                write_number__mod(sheet, 1 + pred_N, LogFC_adj_starts_ColN[x] + col_n, all_values_with_markers[all_values__groups[x][0] + col_n], format_quasi_invisible_left_right_borders)
                            elif col_n == 0:
                                write_number__mod(sheet, 1 + pred_N, LogFC_adj_starts_ColN[x] + col_n, all_values_with_markers[all_values__groups[x][0] + col_n], format_quasi_invisible_left_border)
                            elif col_n == len(LogFC_profile_cols) - 1:
                                write_number__mod(sheet, 1 + pred_N, LogFC_adj_starts_ColN[x] + col_n, all_values_with_markers[all_values__groups[x][0] + col_n], format_quasi_invisible_right_border)
                            else:
                                write_number__mod(sheet, 1 + pred_N, LogFC_adj_starts_ColN[x] + col_n, all_values_with_markers[all_values__groups[x][0] + col_n], format_quasi_invisible)
                        else:
                            write_number__mod(sheet, 1 + pred_N, LogFC_adj_starts_ColN[x] + col_n, all_values_with_markers[all_values__groups[x][0] + col_n], format2)

                if args.Predictor_rows_count > 1:
                    HM_max_color = color(Gradient(HM_Preds_Gradient, pred_N / (args.Predictor_rows_count - 1)))
                else:  HM_max_color = color(Gradient(HM_Preds_Gradient, 1))

                NaNs_present = any([math.isnan(x) for x in all_values])
                if NaNs_present and not negative_values_present:
                    sheet.conditional_format(1 + pred_N, LogFC_adj_starts_ColN[x], 1 + pred_N, LogFC_adj_starts_ColN[-1] + len(LogFC_sparklines_groups[LogFC_sparklines_groups_names[-1]]) - 1,
                       {'type': '3_color_scale',
                        'min_color': "#aaaaaa ", 'mid_color': "#ffffff",
                        'max_color': HM_max_color,
                        'min_type': 'num', 'mid_type': 'num',
                        'max_type': 'num',
                        'min_value': min_value, 'mid_value': min(all_values_without_NaN),
                        'max_value': max(all_values_without_NaN)})
                else:
                    sheet.conditional_format(1 + pred_N, LogFC_adj_starts_ColN[x], 1 + pred_N, LogFC_adj_starts_ColN[-1] + len(LogFC_sparklines_groups[LogFC_sparklines_groups_names[-1]]) - 1,
                       {'type': '3_color_scale',
                        'min_color': "#ffffff",
                        'max_color': HM_max_color,
                        'min_type': 'num',
                        'max_type': 'num',
                        'min_value': min(all_values_without_NaN),
                        'max_value': max(all_values_without_NaN)})

                if args.insert_empty_col_before_LogFC_sparklines:
                    sheet.write(1 + pred_N, LogFC_sparklines_start_col - 1, ' ')


                for x in range(len(LogFC_sparklines_groups_names)):
                    Gn = LogFC_sparklines_groups_names[x]
                    LogFC_profile_cols = LogFC_sparklines_groups[Gn]

                    start_ColN = LogFC_adj_starts_ColN[x]
                    end_ColN = LogFC_adj_starts_ColN[x] + len(LogFC_sparklines_groups[Gn]) - 1
                    # all_SL_Locations = [xl_range(c, LogFC_sparklines_start_col - 1 + x + args.insert_empty_col_before_LogFC_sparklines,
                    #                              c, LogFC_sparklines_start_col - 1 + x + args.insert_empty_col_before_LogFC_sparklines) for c in range(args.Predictor_rows_count + 2, Last_Row)]
                    # all_SL_Ranges = [xl_range(c, start_ColN, c, end_ColN) for c in range(args.Predictor_rows_count + 2, Last_Row)]

                    if args.Predictor_rows_count > 1:
                        Spark_color = color(Gradient(Sparklines_Preds_Gradient, pred_N / (args.Predictor_rows_count - 1)))
                    else:  Spark_color = color(Gradient(Sparklines_Preds_Gradient, 1))

                    sheet.add_sparkline(1 + pred_N, LogFC_sparklines_start_col - 1 + x + args.insert_empty_col_before_LogFC_sparklines,
                                        {'location': xl_range(1 + pred_N, LogFC_sparklines_start_col - 1 + x + args.insert_empty_col_before_LogFC_sparklines,
                                                              1 + pred_N, LogFC_sparklines_start_col - 1 + x + args.insert_empty_col_before_LogFC_sparklines),
                                         'range': xl_range(1 + pred_N, start_ColN, 1 + pred_N, end_ColN),
                                         'type': 'column', 'max': max_value, 'min': min_value,
                                         'negative_points': True,
                                         'series_color': Spark_color,
                                         'negative_color': ('#919191' if negative_values_present else '#d1d1d1')})



        

        ## conditional formatting for predictors
        if max_samples_count > 0 and args.Predictor_rows_count > 0:
            if args.Group_predictors:
                if args.Predictor_blue_yellow_formatting:
                    sheet.conditional_format(1, WSparkLine_Offset(args, 1), args.Predictor_rows_count, WSparkLine_Offset(args, max_samples_count),
                       {'type': '3_color_scale',
                           'min_color': "#0cb6f4", 'mid_color': "#ffffff",
                           'max_color': "#fa8525",
                           'min_type': 'num', 'mid_type': 'percentile',
                           'max_type': 'percentile',
                           'min_value': 0, 'mid_value': 50,
                           'max_value': 100})
                else:
                    sheet.conditional_format(1, WSparkLine_Offset(args, 1), args.Predictor_rows_count, WSparkLine_Offset(args, max_samples_count),
                       {'type': '3_color_scale',
                           'min_color': "#ffffff ", 'mid_color': "#ffe389",
                           'max_color': "#9ece49",
                           'min_type': 'num', 'mid_type': 'percentile',
                           'max_type': 'percentile',
                           'min_value': 0, 'mid_value': 50,
                           'max_value': 100})

            else:
                for rowN in range(1, 1 + args.Predictor_rows_count):
                    if args.Predictor_blue_yellow_formatting:
                        sheet.conditional_format(rowN, WSparkLine_Offset(args, 1), rowN, WSparkLine_Offset(args, max_samples_count),
                           {'type': '3_color_scale',
                               'min_color': "#0cb6f4", 'mid_color': "#ffffff",
                               'max_color': "#fa8525",
                               'min_type': 'num', 'mid_type': 'percentile',
                               'max_type': 'percentile',
                               'min_value': 3, 'mid_value': 50,
                               'max_value': 97})
                    else:
                        sheet.conditional_format(rowN, WSparkLine_Offset(args, 1), rowN, WSparkLine_Offset(args, max_samples_count),
                           {'type': '3_color_scale',
                               'min_color': "#ffffff ", 'mid_color': "#ffe389",
                               'max_color': "#9ece49",
                               'min_type': 'num', 'mid_type': 'percentile',
                               'max_type': 'percentile',
                               'min_value': 0, 'mid_value': 50,
                               'max_value': 100})

        RowN += 1
        Max_ColN = 0
        for L in all_Lines[args.Predictor_rows_count:]:
            if RowN > current_max_strings:  break
            C = L.replace('\n','').replace('\r','').split('\t')
            C = [C_el[1:-1] if ((C_el.startswith('"') and C_el.endswith('"')) or (C_el.startswith("'") and C_el.endswith("'"))) else C_el for C_el in C]
            
            CPM_values = []
            CPM_values2 = []
            CPM_narrow_values = []

            CPM_values__clean = []
            CPM_values2__clean = []
            CPM_narrow_values__clean = []

            Max_ColN = max(Max_ColN, len(C))
            for x in range(len(C)):
                if x in Spacer_cols: continue

                if x in CPM_cols:
                    try:  x_float = float(C[x]);    CPM_values += [x_float];           CPM_values__clean += [x_float]
                    except ValueError:   CPM_values += [None]
                
                if x in CPM_cols2:
                    try:   x_float = float(C[x]);   CPM_values2 += [x_float];          CPM_values2__clean += [x_float]
                    except ValueError:   CPM_values2 += [None]

                if x in CPM_narrow_cols:
                    try:   x_float = float(C[x]);   CPM_narrow_values += [x_float];    CPM_narrow_values__clean += [x_float]
                    except ValueError:   CPM_narrow_values += [None]

                if C[x].casefold() in ['na','nan']:
                    if x in logFC_narrow_cols:
                        sheet.write(RowN, WSparkLine_Offset(args, x), '', empty_logfc_format)
                    continue

                if args.max_char_in_cell is not None:
                    if len(C[x]) > args.max_char_in_cell:
                        C[x] = '%s .... [truncated]'%C[x][ : max(1, args.max_char_in_cell - 17)]

                if x == 0 and not args.first_col_is_formatted:
                    if args.first_col_italic:  write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], italic)
                    else:  write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x])
                elif x > 0 and args.dist_matrix_mode:
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], format1_center)
                elif (x in P_value_cols) or (x in FDR_cols):
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], pvalue_mode = True, pvalue_formats = pvalue_formats)
                elif (x in logFC_cols) or (x in CPM_cols) or (x in Correlation_r_cols) or (x in CPM_cols2):
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], format2)
                elif x in logFC_narrow_cols:
                    if ConvertableToNumeric(C[x]):
                        write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], format2_narrow)
                    else:
                        sheet.write(RowN, WSparkLine_Offset(args, x), C[x], empty_logfc_format)
                elif x in CPM_narrow_cols:
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], format1_narrow)
                elif (x in Score_cols) or (x in raw_means_CPM_cols) or (x in Index1_cols) or (x in Index1Err_cols) or (x in LogFC_array_cols) or (x in LogFC_array2_cols):
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], format1_center)
                elif (x in Index2_cols) or (x in Index2Err_cols):
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], format2_center)
                elif x in Index_cols or x in gradual_30_cols or x in FFPM_cols:
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], format_center)
                elif x in logCPM_cols or x in percentage_cols:
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], format1)
                elif x in gray_gradient_cols + warm_gradient_cols + cold_gradient_cols + cov100norm_cols:
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], format0_center)
                elif x in delta_logfc_cols:
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], whitespace_format1)
                elif x in percents_updown_cols:
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], percentage_format0)
                elif x in red_bars_custom_cols or x in blue_bars_custom_cols:
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], format1)
                elif x in red_gradient_custom_cols or x in blue_gradient_custom_cols:
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], format0_center)
                elif x in Taxon_cols:
                    sheet.write(RowN, WSparkLine_Offset(args, x), C[x], italic)
                elif x in Color_cols:
                    if not C[x].startswith('#') or len(C[x]) != 7:
                        print('Incorrect color format %s'%C[x])
                        exit(127)
                    CurrentFormat = Workbook.add_format()
                    CurrentFormat.set_bg_color(C[x])
                    sheet.write(RowN, WSparkLine_Offset(args, x), ' ', CurrentFormat)
                elif x in Filter_cols:
                    cf = C[x].casefold()
                    if 'true' == cf or 'ok' == cf or 'passed' == cf:
                        sheet.write(RowN, WSparkLine_Offset(args, x),'yes', format_center)
                elif x in Counts_cols:
                    write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x], format0)
                elif x in Rank_1_cols or x in Rank_2_cols:
                    WriteFormatted_Rank_Cell(Workbook, sheet, RowN, WSparkLine_Offset(args, x), C[x],
                        upper_part = 17, lower_part = 17, extent = 1.7, default_format = default_Rank_format)
                else:  write_number__mod(sheet, RowN, WSparkLine_Offset(args, x), C[x])

                if not args.Simple_LogFC_format and not args.Bidirectional_bars_in_LogFC_cells and x in logFC_cols:
                    FormatLogFC_Cell(sheet, 2**float(C[x]), RowN, WSparkLine_Offset(args, x), coeff1 = 7/args.LogFC_cell_abs_max, coeff2 = 7/args.LogFC_cell_abs_max)


            CPM_values_pack = [CPM_values, CPM_values2, CPM_narrow_values]
            CPM_values__clean_pack = [CPM_values__clean, CPM_values2__clean, CPM_narrow_values__clean]
            CPM_cols_pack = [CPM_cols, CPM_cols2, CPM_narrow_cols]
            is_narrow = [False, False, True]

            for current_CPM_values, current_CPM_values__clean, current_CPM_cols, current_is_narrow in zip(CPM_values_pack, CPM_values__clean_pack, CPM_cols_pack, is_narrow):
                if args.check_for_NA_in_CPM_cols:
                    if(len(current_CPM_values) != len(current_CPM_cols)):
                        print('Error 342. Line %d: len(current_CPM_values) = %d, len(current_CPM_cols) = %d. They should be equal. Please check for missing or text values in CPM columns, which cannot be converted to float'%(
                            RowN, len(current_CPM_values), len(current_CPM_cols)))
                        exit(127)

                if args.Add_cpm_cells_formatting and len(current_CPM_values__clean) > 0:
                    current_CPM_values__clean_mod = [x + args.cpm_cells_formatting__const_add for x in current_CPM_values__clean]
                    if args.cpm_cells_formatting__trimming_low == 0 and args.cpm_cells_formatting__trimming_high == 0:
                        # avg_current_CPM = sum(current_CPM_values__clean) / len(current_CPM_values__clean)
                        gmean_current_CPM = geo_mean_overflow(current_CPM_values__clean_mod)
                    else:
                        # avg_current_CPM = Smoothed_Trimmed_Mean(current_CPM_values__clean, from_percentile = args.cpm_cells_formatting__trimming_low, to_percentile = 100 - args.cpm_cells_formatting__trimming_high)
                        gmean_current_CPM = Smoothed_Trimmed_GeoMean(current_CPM_values__clean_mod, from_percentile = args.cpm_cells_formatting__trimming_low, to_percentile = 100 - args.cpm_cells_formatting__trimming_high)

                    # gmean_current_CPM = 1
                    # for x in current_CPM_values__clean_mod:
                    #     gmean_current_CPM = gmean_current_CPM*(x/(avg_current_CPM + args.cpm_cells_formatting__const_add))
                    # gmean_current_CPM = (gmean_current_CPM**(1/len(current_CPM_values__clean_mod)))*(avg_current_CPM  + args.cpm_cells_formatting__const_add)
                    if gmean_current_CPM == 0:   continue

                    current_LogFCs = [None if x is None else math.log2((x + args.cpm_cells_formatting__const_add) / (gmean_current_CPM)) for x in current_CPM_values]
                    current_LogFCs = [None if x is None else min(max(x, (-1)*args.cpm_cells_formatting__LogFC_max), args.cpm_cells_formatting__LogFC_max) for x in current_LogFCs]
                    current_LogFCs_color_coords = [None if x is None else 0.5 + x / args.cpm_cells_formatting__LogFC_max / 2 for x in current_LogFCs]

                    avg_current_CPM = sum(current_CPM_values__clean) / len(current_CPM_values__clean)
                    cpm_color_coord = (avg_current_CPM / args.cpm_cells_formatting__cpm_max)**args.cpm_cells_formatting__cpm_degree

                    # print(CPM_cols)
                    for n in range(len(current_CPM_cols)):
                        col_n = current_CPM_cols[n]
                        try:  val = float(C[col_n])
                        except:  continue

                        if current_LogFCs_color_coords[n] is None:  continue

                        if math.isnan(val):   BgC = '#ffffff'
                        else:
                            try:   BgC = color(TwoD_Gradient(Color_gradient_schemas[args.cpm_cells_formatting__gradient_color_schema - 1], cpm_color_coord, current_LogFCs_color_coords[n]))
                            except:
                                print('[Error 236]  Incorrect CPM columns. Please check (line %d)'%RowN)
                                print('n = %d, len(current_LogFCs_color_coords) = %d, len(current_CPM_cols) = %d'%(n, len(current_LogFCs_color_coords), len(current_CPM_cols)))
                                try:  print('val = %s, current_LogFCs_color_coords[n] = %g, cpm_color_coord = %g'%(str(val), current_LogFCs_color_coords[n], cpm_color_coord))
                                except:  pass
                                exit(127)

                        if (BgC + str(current_is_narrow)) in existing_CPM_formats:  CurrentFormat = existing_CPM_formats[BgC + str(current_is_narrow)]
                        else:
                            CurrentFormat = Workbook.add_format()
                            CurrentFormat.set_align('center')
                            CurrentFormat.set_bg_color(BgC)
                            CurrentFormat.set_border(0)
                            if current_is_narrow:
                                CurrentFormat.set_font_size(8)
                                if args.CPM_narrow_digits > 0:   CurrentFormat.set_num_format('0.' + '0'*args.CPM_narrow_digits)
                                else:                            CurrentFormat.set_num_format('0')
                            else:                                CurrentFormat.set_num_format('0.' + '0'*args.CPM_digits)
                            existing_CPM_formats[BgC + str(current_is_narrow)] = CurrentFormat

                        write_number__mod(sheet, RowN, WSparkLine_Offset(args, col_n), val, format = CurrentFormat)

                    # args.cpm_cells_formatting__const_add
                    # args.cpm_cells_formatting__LogFC_max
                    # args.cpm_cells_formatting__cpm_max
                    # args.cpm_cells_formatting__cpm_degree
                    # args.cpm_cells_formatting__gradient_color_schema

            if CPM_sparklines_groups is not None:
                all_values = []
                all_values__groups = []
                for x in range(len(CPM_sparklines_groups_names)):
                    Gn = CPM_sparklines_groups_names[x]
                    CPM_profile_cols = CPM_sparklines_groups[Gn]
                    #sheet.write(0, CPM_sparklines_start_col - 1 + x, Gn, headers_format)

                    all_values__groups.append([len(all_values), len(all_values) + len(CPM_profile_cols) - 1])
                    for col_n in range(len(CPM_profile_cols)):
                        try:  all_values.append(float(C[CPM_profile_cols[col_n]]))
                        except ValueError:  all_values.append(None)

                if None not in all_values:
                    max_value = max(all_values)
                    if args.CPM_sparklines_mode == 'linear' and max_value > 0:
                        all_values = [x/max_value*100 for x in all_values]

                    elif args.CPM_sparklines_mode == 'log' and max_value > 0:
                        # all_values
                        # avg_current_CPM = sum(all_values) / len(all_values)
                        # median_current_CPM = statistics.median(all_values)
                        # current_CPM_values_mod = [x + 1 + args.CPM_sparklines_rel_scale_cpm_const_add for x in all_values]
                        current_CPM_values_mod = [x + args.CPM_sparklines_rel_scale_cpm_const_add for x in all_values]
                        
                        # gmean_current_CPM = 1
                        # for x in current_CPM_values_mod:
                        #     gmean_current_CPM = gmean_current_CPM * (x / (avg_current_CPM + 1 + args.CPM_sparklines_rel_scale_cpm_const_add))
                        # gmean_current_CPM = (gmean_current_CPM ** (1 / len(current_CPM_values_mod))) * (avg_current_CPM + 1 + args.CPM_sparklines_rel_scale_cpm_const_add)

                        if args.CPM_sparklines_rel_scale_cpm_trimming_low == 0 and args.CPM_sparklines_rel_scale_cpm_trimming_high == 0:
                            gmean_current_CPM = geo_mean_overflow(current_CPM_values_mod)
                        else:
                            gmean_current_CPM = Smoothed_Trimmed_GeoMean(current_CPM_values_mod, from_percentile = args.CPM_sparklines_rel_scale_cpm_trimming_low, to_percentile = 100 - args.CPM_sparklines_rel_scale_cpm_trimming_high)

                        current_LogFCs = [math.log2(x / (gmean_current_CPM)) for x in current_CPM_values_mod]
                        all_values = current_LogFCs
                        # current_LogFCs = [min(max(x, (-1) * args.cpm_cells_formatting__LogFC_max), args.cpm_cells_formatting__LogFC_max) for x in current_LogFCs]

                    for x in range(len(CPM_sparklines_groups_names)):
                        Gn = CPM_sparklines_groups_names[x]
                        CPM_profile_cols = CPM_sparklines_groups[Gn]

                        if len(CPM_profile_cols) != all_values__groups[x][1] - all_values__groups[x][0] + 1:
                            print('Smth strange (2) %d != %d - %d + 1....'%(len(CPM_profile_cols), all_values__groups[x][0], all_values__groups[x][1]))
                            exit(127)

                        for col_n in range(len(CPM_profile_cols)):
                            if args.include_CPM_heatmaps:
                                if len(CPM_profile_cols) == 1:
                                    write_number__mod(sheet, RowN, CPM_adj_starts_ColN[x] + col_n, all_values[all_values__groups[x][0] + col_n], format_quasi_invisible_left_right_borders)
                                elif col_n == 0:
                                    write_number__mod(sheet, RowN, CPM_adj_starts_ColN[x] + col_n, all_values[all_values__groups[x][0] + col_n], format_quasi_invisible_left_border)
                                elif col_n == len(CPM_profile_cols) - 1:
                                    write_number__mod(sheet, RowN, CPM_adj_starts_ColN[x] + col_n, all_values[all_values__groups[x][0] + col_n], format_quasi_invisible_right_border)
                                else:
                                    write_number__mod(sheet, RowN, CPM_adj_starts_ColN[x] + col_n, all_values[all_values__groups[x][0] + col_n], format_quasi_invisible)
                            else:
                                write_number__mod(sheet, RowN, CPM_adj_starts_ColN[x] + col_n, all_values[all_values__groups[x][0] + col_n], format2)

                if args.insert_empty_col_before_CPM_sparklines:
                    sheet.write(RowN, CPM_sparklines_start_col - 1, ' ')


            if LogFC_sparklines_groups is not None:
                all_values = []
                all_values__groups = []
                for x in range(len(LogFC_sparklines_groups_names)):
                    Gn = LogFC_sparklines_groups_names[x]
                    LogFC_profile_cols = LogFC_sparklines_groups[Gn]
                    #sheet.write(0, LogFC_sparklines_start_col - 1 + x, Gn, headers_format)

                    all_values__groups.append([len(all_values), len(all_values) + len(LogFC_profile_cols) - 1])
                    for col_n in range(len(LogFC_profile_cols)):
                        try:  all_values.append(float(C[LogFC_profile_cols[col_n]]))
                        except ValueError:  all_values.append(None)

                # if None not in all_values:
                for x in range(len(LogFC_sparklines_groups_names)):
                    Gn = LogFC_sparklines_groups_names[x]
                    LogFC_profile_cols = LogFC_sparklines_groups[Gn]

                    if len(LogFC_profile_cols) != all_values__groups[x][1] - all_values__groups[x][0] + 1:
                        print('Smth strange (5) %d != %d - %d + 1....'%(len(LogFC_profile_cols), all_values__groups[x][0], all_values__groups[x][1]))
                        exit(127)

                    for col_n in range(len(LogFC_profile_cols)):
                        val = all_values[all_values__groups[x][0] + col_n]
                        if val is not None:
                            if args.include_LogFC_heatmaps:
                                if len(LogFC_profile_cols) == 1:
                                    write_number__mod(sheet, RowN, LogFC_adj_starts_ColN[x] + col_n, val, format_quasi_invisible_left_right_borders)
                                elif col_n == 0:
                                    write_number__mod(sheet, RowN, LogFC_adj_starts_ColN[x] + col_n, val, format_quasi_invisible_left_border)
                                elif col_n == len(LogFC_profile_cols) - 1:
                                    write_number__mod(sheet, RowN, LogFC_adj_starts_ColN[x] + col_n, val, format_quasi_invisible_right_border)
                                else:
                                    write_number__mod(sheet, RowN, LogFC_adj_starts_ColN[x] + col_n, val, format_quasi_invisible)
                            else:
                                write_number__mod(sheet, RowN, LogFC_adj_starts_ColN[x] + col_n, val, format2)

                if args.insert_empty_col_before_LogFC_sparklines:
                    sheet.write(RowN, LogFC_sparklines_start_col - 1, ' ')


            if args.Insert_spasers_between_rows:
                RowN += 1
                sheet.set_row(RowN, 3.5)

            RowN += 1

        Last_Row = RowN

        if args.Simple_LogFC_format:
            for ColN in logFC_cols + logFC_narrow_cols:
                if args.logfc_cells_formatting__gradient_color_schema == 1:
                    sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                              'min_color': "#1170f1",'mid_color': "#ffffff",'max_color': "#ec4a18",
                                                              'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                              'min_value': -args.logfc_cells_formatting__gradient_range, 'mid_value': 0.0, 'max_value': args.logfc_cells_formatting__gradient_range})

                elif args.logfc_cells_formatting__gradient_color_schema == 2:
                    sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                              'min_color': "#06ade0",'mid_color': "#ffffff",'max_color': "#efb00e",
                                                              'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                              'min_value': -args.logfc_cells_formatting__gradient_range, 'mid_value': 0.0, 'max_value': args.logfc_cells_formatting__gradient_range})

                else:
                    print('Incorrect args.logfc_cells_formatting__gradient_color_schema. Should be either 1 or 2')
                    exit(127)


        if args.Bidirectional_bars_in_LogFC_cells:
            for ColN in logFC_cols + logFC_narrow_cols:
                sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': 'data_bar',
                                                          'bar_color': '#f28241',
                                                          'bar_border_color': '#b3541d',
                                                          'bar_negative_color':'#558cfa',
                                                          'bar_negative_border_color': '#3565c4',
                                                          'bar_axis_position': 'middle',
                                                          'bar_no_border': False,
                                                          'bar_negative_border_color_same': False,
                                                          'bar_negative_color_same': False,
                                                          'min_type':'num','max_type':'num',
                                                          'min_value':-3.1,'max_value': 3.1,
                                                          'min_length':0, 'max_length':100})
                
                if args.Insert_spasers_between_rows or args.create_White_background:
                    sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 12, whitespace_format)
                else:
                    sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 12)

        for ColN in Correlation_r_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                      'min_color': "#5c86f1",'mid_color': "#ffffff",'max_color': "#ee6c44",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': -0.98, 'mid_value': 0.0, 'max_value': 0.98})

        for ColN in Counts_cols + Index_cols + Index1_cols + Index2_cols:
            # sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
            #                                           'min_color': "#f2dc49",'mid_color': "#9af762",'max_color': "#7fa5f2",
            #                                           'min_type': 'percent','mid_type': 'percent','max_type': 'percent',
            #                                           'min_value': 0, 'mid_value': 50, 'max_value': 100})

            # sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
            #                                           'min_color': "#f96f45",'mid_color': "#ffffff",'max_color': "#71e27c",
            #                                           'min_type': 'percentile','mid_type': 'percentile','max_type': 'percentile',
            #                                           'min_value': 0, 'mid_value': 50, 'max_value': 98})

            # sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
            #                                           'min_color': "#5c86f1",'mid_color': "#ffffff",'max_color': "#ee6c44",
            #                                           'min_type': 'num','mid_type': 'num','max_type': 'num',
            #                                           'min_value': 0.0, 'mid_value': 50.0, 'max_value': 100.0})

            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                      'min_color': "#ffffff",'mid_color': "#afff91",'max_color': "#a1cdff",
                                                      'min_type': args.index_cols_min_type,'mid_type': args.index_cols_mid_type, 'max_type': args.index_cols_max_type,
                                                      'min_value': args.index_cols_min_value, 'mid_value': args.index_cols_mid_value, 'max_value': args.index_cols_max_value})

        

        for ColN in gradual_30_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                      'min_color': "#ffffff",'mid_color': "#ffe08d",'max_color': "#9ece49",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': 0, 'mid_value': 14, 'max_value': 30})

        for ColN in FFPM_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                      'min_color': "#ffffff",'mid_color': "#ffe08d",'max_color': "#9ece49",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': 0, 'mid_value': 2, 'max_value': 5})

        for ColN in delta_logfc_cols:
            #'min_color': "#06ade0",'mid_color': "#ffffff",'max_color': "#efb00e",
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': 'data_bar',
                                                      'bar_color': '#efb00e',
                                                      'bar_negative_color':'#06ade0',
                                                      'bar_border_color': '#bc8803',
                                                      'bar_negative_border_color': '#0387af',
                                                      'min_type':'num','max_type':'num',
                                                      'min_value':-3.1,'max_value': 3.1,
                                                      'min_length':0, 'max_length':100})
            if ColN in Hidden_cols:  continue
            if args.Insert_spasers_between_rows or args.create_White_background:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 15, whitespace_format)
            else:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 15)

        for ColN in percents_updown_cols:
            #'min_color': "#06ade0",'mid_color': "#ffffff",'max_color': "#efb00e",
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                      'min_color': "#317bf5",'mid_color': "#ffffff",'max_color': "#ec692f",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': -1, 'mid_value': 0, 'max_value': 4})

        for ColN in red_bars_custom_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': 'data_bar',
                                                      'bar_color': args.red_bars_color,#'#f25e61',
                                                      # 'bar_border_color': '#d13033',
                                                      'min_type':args.red_bars_min_type,'max_type':args.red_bars_max_type,
                                                      'min_value':args.red_bars_min_value, 'max_value': args.red_bars_max_value,
                                                      'min_length':0, 'max_length':100})
            
        for ColN in blue_bars_custom_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': 'data_bar',
                                                      'bar_color': args.blue_bars_color,#'#5b9beb',
                                                      # 'bar_border_color': '#185096',
                                                      'min_type':args.blue_bars_min_type,'max_type':args.blue_bars_max_type,
                                                      'min_value':args.blue_bars_min_value, 'max_value': args.blue_bars_max_value,
                                                      'min_length':0, 'max_length':100})

        for ColN in red_gradient_custom_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '2_color_scale',
                                                      'min_color': "#ffffff", 'max_color': args.red_gradient_color,
                                                      'min_type': args.red_gradient_min_type, 'max_type': args.red_gradient_max_type,
                                                      'min_value': args.red_gradient_min_value, 'max_value': args.red_gradient_max_value})

        for ColN in blue_gradient_custom_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '2_color_scale',
                                                      'min_color': "#ffffff", 'max_color': args.blue_gradient_color,
                                                      'min_type': args.blue_gradient_min_type, 'max_type': args.blue_gradient_max_type,
                                                      'min_value': args.blue_gradient_min_value, 'max_value': args.blue_gradient_max_value})

            
            # sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': 'data_bar',
            #                                           'bar_color': '#efb00e',
            #                                           'bar_negative_color':'#06ade0',
            #                                           'bar_border_color': '#bc8803',
            #                                           'bar_negative_border_color': '#0387af',
            #                                           'min_type':'num','max_type':'num',
            #                                           'min_value': -1.0,'max_value': 3.0,
            #                                           'min_length':0, 'max_length':100})

            # if ColN in Hidden_cols:  continue
            # sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 15)

        if args.Insert_spasers_between_rows or args.create_White_background:
            for ColN in range(0, Cols_count*3):
                if ColN in Hidden_cols:  continue
                sheet.set_column(ColN, ColN, 9, whitespace_format)

        for ColN in Rank_1_cols:
            if args.Insert_spasers_between_rows or args.create_White_background:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 2.5, whitespace_format)
            else:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 2.5)
        
        for ColN in Rank_2_cols:
            if args.Insert_spasers_between_rows or args.create_White_background:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 3.5, whitespace_format)
            else:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 3.5)

        for ColN in Taxon_cols:
            if args.Insert_spasers_between_rows or args.create_White_background:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 35, whitespace_format)
            else:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 35)


        for ColN in logCPM_cols:
            if args.Metagenome_mode:  max_value = 20
            else:  max_value = 10
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': 'data_bar','bar_color': '#d89c0d',
                                                                                  'min_type':'num','max_type':'num',
                                                                                  'min_value':0,'max_value': max_value})

        for ColN in cov100norm_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                      'min_color': args.cov100norm_colors[0], 'mid_color': args.cov100norm_colors[1], 'max_color': args.cov100norm_colors[2],
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': args.cov100norm_range[0], 'mid_value': args.cov100norm_range[1], 'max_value': args.cov100norm_range[2]})

        for ColN in percentage_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': 'data_bar','bar_color': '#1cbaff',
                                                                                  'min_type':'num','max_type':'num',
                                                                                  'min_value':0,'max_value': 100})

        for ColN in warm_gradient_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                      'min_color': "#ffffff",'mid_color': "#ffeb84",'max_color': "#f8696b",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': 0, 'mid_value': 31, 'max_value': 82})

        for ColN in cold_gradient_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                      'min_color': "#ffffff",'mid_color': "#85ebe6",'max_color': "#3c86e8",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': 0, 'mid_value': 32, 'max_value': 90})

        for ColN in gray_gradient_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                      'min_color': "#ffffff",'mid_color': "#d9d9d9",'max_color': "#808080",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': 0, 'mid_value': 32, 'max_value': 90})

        for ColN in P_value_cols + FDR_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                      'min_color': "#9ece49",'mid_color': "#ffe08d",'max_color': "#ffffff",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': 0, 'mid_value': 0.0002, 'max_value': 0.07})

        for ColN in LogFC_array_cols + LogFC_array2_cols:
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                      'min_color': "#06ade0",'mid_color': "#ffffff",'max_color': "#efb00e",
                                                      'min_type': 'num','mid_type': 'num','max_type': 'num',
                                                      'min_value': -args.LogFC_array_abs_max, 'mid_value': 0.0, 'max_value': args.LogFC_array_abs_max})

        if args.Metagenome_mode:
            for ColN in Score_cols:
                sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type': '3_color_scale',
                                                                   'min_color': "#ffffff", 'mid_color': "#ffe08d",
                                                                   'max_color': "#9ece49",
                                                                   'min_type': 'num', 'mid_type': 'num',
                                                                   'max_type': 'num',
                                                                   'min_value': 2, 'mid_value': 12,
                                                                   'max_value': 25})

        # Light red fill with dark red text.
        format_lincRNA = Workbook.add_format({'bg_color':   '#c3d594'})

        # Light yellow fill with dark yellow text.
        format_antisense = Workbook.add_format({'bg_color':   '#d9c188'})

        # Green fill with dark green text.
        format_pseudogene = Workbook.add_format({'bg_color':   '#a0b4c7'})

        for ColN in Biotype_cols:
            if ColN in Hidden_cols:  continue
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type':'text', 'criteria': 'containing', 'value': 'lincRNA', 'format':   format_lincRNA})
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type':'text', 'criteria': 'containing', 'value': 'lncRNA', 'format':   format_lincRNA})
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type':'text', 'criteria': 'containing', 'value': 'antisense', 'format':   format_antisense})
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type':'text', 'criteria': 'containing', 'value': 'AS', 'format':   format_antisense})
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type':'text', 'criteria': 'containing', 'value': 'pseudogene', 'format':   format_pseudogene})
            sheet.conditional_format(1, WSparkLine_Offset(args, ColN), Last_Row, WSparkLine_Offset(args, ColN), {'type':'text', 'criteria': 'containing', 'value': 'PG', 'format':   format_pseudogene})
            if args.Insert_spasers_between_rows or args.create_White_background:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 18, whitespace_format)
            else:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 18)

        if args.first_col_width != None:  fw = 10
        else:    fw = args.first_col_width
        if args.Insert_spasers_between_rows or args.create_White_background:
            sheet.set_column(WSparkLine_Offset(args, 0), WSparkLine_Offset(args, 0), fw, whitespace_format)
        else:
            sheet.set_column(WSparkLine_Offset(args, 0), WSparkLine_Offset(args, 0), fw)
            

        for ColN in Gene_Name_cols:
            if ColN in Hidden_cols:  continue
            if args.Insert_spasers_between_rows or args.create_White_background:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 35, whitespace_format)
            else:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 35)

        if args.Metagenome_mode:
            if args.Insert_spasers_between_rows or args.create_White_background:
                sheet.set_column(WSparkLine_Offset(args, 0), WSparkLine_Offset(args, 0), 19, whitespace_format)
            else:
                sheet.set_column(WSparkLine_Offset(args, 0), WSparkLine_Offset(args, 0), 19)

        for ColN in warm_gradient_cols + cold_gradient_cols + gray_gradient_cols:
            if ColN in Hidden_cols:  continue
            if args.Insert_spasers_between_rows or args.create_White_background:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 5.5, whitespace_format)
            else:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 5.5)

        for ColN in cov100norm_cols:
            if ColN in Hidden_cols:  continue
            if args.Insert_spasers_between_rows or args.create_White_background:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 5.5, whitespace_format)
            else:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 5.5)

        for ColN in Spacer_cols:
            if ColN in Hidden_cols:  continue
            if args.Insert_spasers_between_rows or args.create_White_background:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 2, whitespace_format)
            else:
                sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), 2)

        for ColN in Hidden_cols:
            sheet.set_column(WSparkLine_Offset(args, ColN), WSparkLine_Offset(args, ColN), None, None, {'hidden': True})

        all_CPM_sparklines_cols = []
        if CPM_sparklines_groups is not None:
            Sparklines_Elements_Count = []
            for x in range(len(CPM_sparklines_groups_names)):
                Gn = CPM_sparklines_groups_names[x]
                Sparklines_Elements_Count.append(len(CPM_sparklines_groups[Gn]))


            sum_d = sum([x**args.CPM_sparklines_col_width_factor for x in Sparklines_Elements_Count])
            standart_width = 5.5 * (((sum(Sparklines_Elements_Count) / len(Sparklines_Elements_Count)) /  3)**0.7)

            Sparklines_Column_Widths = [(x ** args.CPM_sparklines_col_width_factor) / sum_d * standart_width * len(CPM_sparklines_groups_names) for x in Sparklines_Elements_Count]
            if args.CPM_sparklines_equal_barwidths:
                Sparklines_Column_Widths = [x * 30 / max(Sparklines_Column_Widths + [30]) for x in Sparklines_Column_Widths]
            else:
                Sparklines_Column_Widths = [min(27, max(5, x)) for x in Sparklines_Column_Widths]

            for x in range(len(CPM_sparklines_groups_names)):
                Gn = CPM_sparklines_groups_names[x]
                CPM_profile_cols = CPM_sparklines_groups[Gn]

                start_ColN = CPM_adj_starts_ColN[x]
                end_ColN = CPM_adj_starts_ColN[x] + len(CPM_sparklines_groups[Gn]) - 1
                all_SL_Locations = [xl_range(c, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines,
                                             c, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines) for c in range(args.Predictor_rows_count + 2, Last_Row)]
                all_SL_Ranges = [xl_range(c, start_ColN, c, end_ColN) for c in range(args.Predictor_rows_count + 2, Last_Row)]
                if args.CPM_sparklines_mode == 'linear':
                    sheet.conditional_format(args.Predictor_rows_count + 1, start_ColN, Last_Row, end_ColN,
                                             {'type': '3_color_scale',
                                              'min_color': "#ffffff", 'mid_color': "#819ee8",
                                              'max_color': "#386cef",
                                              'min_type': 'num', 'mid_type': 'percent',
                                              'max_type': 'num',
                                              'min_value': 0, 'mid_value': 50,
                                              'max_value': 100})

                    if len(CPM_sparklines_groups) > 1:
                        Spark_color = color(Gradient(Sparklines_CPM_Groups_Gradient, x / (len(CPM_sparklines_groups) - 1)))
                    else:  Spark_color = color(Gradient(Sparklines_CPM_Groups_Gradient, 1))

                    sheet.add_sparkline(1 + args.Predictor_rows_count, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines,
                                        {'location': all_SL_Locations, 'range': all_SL_Ranges,
                                         'type': 'column', 'max': 100, 'min': 0,
                                         'series_color': Spark_color})


                elif args.CPM_sparklines_mode == 'log':
                    sheet.conditional_format(args.Predictor_rows_count + 1, start_ColN, Last_Row, end_ColN,
                                             {'type': '3_color_scale',
                                              'min_color': "#4771dd", 'mid_color': "#ffffff",
                                              'max_color': "#fa8323",
                                              'min_type': 'num', 'mid_type': 'num',
                                              'max_type': 'num',
                                              'min_value': -args.CPM_sparklines_logfc_mode_axis_limit - 0.01, 'mid_value': 0,
                                              'max_value': args.CPM_sparklines_logfc_mode_axis_limit + 0.01})

                    sheet.add_sparkline(1 + args.Predictor_rows_count, CPM_sparklines_start_col - 1 + x + args.insert_empty_col_before_CPM_sparklines,
                                        {'location': all_SL_Locations, 'range': all_SL_Ranges,
                                         'type': 'column', 'max': args.CPM_sparklines_logfc_mode_axis_limit + 0.01, 'min': -args.CPM_sparklines_logfc_mode_axis_limit - 0.01,
                                         'negative_points': True,
                                         'series_color': '#ff6e2e',
                                         'negative_color': '#327fff'})

                all_CPM_sparklines_cols += [CPM_sparklines_start_col + x - 1 + args.insert_empty_col_before_CPM_sparklines]
                if args.Insert_spasers_between_rows or args.create_White_background:
                    sheet.set_column(CPM_sparklines_start_col + x - 1 + args.insert_empty_col_before_CPM_sparklines,
                                     CPM_sparklines_start_col + x - 1 + args.insert_empty_col_before_CPM_sparklines, Sparklines_Column_Widths[x], sparklines_cell_format_white_bg)
                else:
                    sheet.set_column(CPM_sparklines_start_col + x - 1 + args.insert_empty_col_before_CPM_sparklines,
                                     CPM_sparklines_start_col + x - 1 + args.insert_empty_col_before_CPM_sparklines, Sparklines_Column_Widths[x], sparklines_cell_format_no_bg)

                if args.include_CPM_heatmaps:
                    if start_ColN > 0:
                        if args.Insert_spasers_between_rows or args.create_White_background:
                            sheet.set_column(start_ColN - 1, start_ColN - 1, 1, whitespace_format)
                        else:
                            sheet.set_column(start_ColN - 1, start_ColN - 1, 1)

                    # heatmap_cell_format_left_border, heatmap_cell_format_right_border, heatmap_cell_format_left_right_borders
                    if start_ColN == end_ColN:
                        sheet.set_column(start_ColN, start_ColN, Sparklines_Column_Widths[x] / Sparklines_Elements_Count[x], heatmap_cell_format_left_right_borders)
                    else:
                        sheet.set_column(start_ColN, start_ColN, Sparklines_Column_Widths[x] / Sparklines_Elements_Count[x], heatmap_cell_format_left_border)
                        for ColN in range(start_ColN + 1, end_ColN):
                            if args.Insert_spasers_between_rows or args.create_White_background:
                                sheet.set_column(ColN, ColN, Sparklines_Column_Widths[x] / Sparklines_Elements_Count[x], whitespace_format)
                            else:
                                sheet.set_column(ColN, ColN, Sparklines_Column_Widths[x] / Sparklines_Elements_Count[x])
                        sheet.set_column(end_ColN, end_ColN, Sparklines_Column_Widths[x] / Sparklines_Elements_Count[x], heatmap_cell_format_right_border)

                    if args.Insert_spasers_between_rows or args.create_White_background:
                        sheet.set_column(end_ColN + 1, end_ColN + 1, 1, whitespace_format)
                    else:
                        sheet.set_column(end_ColN + 1, end_ColN + 1, 1)

                    # worksheet1.add_sparkline(26, 1, {'location': ['A27', 'A28'],
                    #                                  'range': ['Sheet2!A5:J5', 'Sheet2!A6:J6'],
                    #                                  'markers': True})

                    #     monolyth = True
        #     for Gn in args.CPM_sparklines_groups:
        #         cols = args.CPM_sparklines_groups[Gn]
        #         col_prev = cols[0]
        #         for x in range(1,len(cols)):
        #             if cols[x] != col_prev + 1:
        #                 print('CPM sparkline Group %s (%s) seems to be not solid'%(Gn,str(cols)))
        #                 monolyth = False
        #             col_prev = cols[x]
        #
        #     if monolyth:
        #         for x in range(len(CPM_sparklines_groups_names)):
        #             Gn = CPM_sparklines_groups_names[x]
        #             cols = args.CPM_sparklines_groups[Gn]
        #             #cols[0]:cols[-1]
        #             sheet.write(0, CPM_sparklines_start_col + x, 'g."%s" CPMs'%Gn, headers_format)
        #             if len(args.CPM_sparklines_groups) > 1:
        #                 Spark_color = color(Gradient(Sparklines_CPM_Groups_Gradient, x / (len(args.CPM_sparklines_groups) - 1)))
        #             else:  Spark_color = color(Gradient(Sparklines_CPM_Groups_Gradient, 1))
        #
        #             for RowN in range(args.Predictor_rows_count + 1, Last_Row + 1):
        #                 sheet.add_sparkline(RowN, CPM_sparklines_start_col + x,
        #                                     {'range': xl_range(RowN, cols[0], RowN, cols[-1]),
        #                                      'type': 'column', 'min': 0,
        #                                      'series_color':Spark_color})
        #
        #     else:
        #         #Max_ColN
        #         pass

        all_LogFC_sparklines_cols = []
        if LogFC_sparklines_groups is not None:
            Sparklines_Elements_Count = []
            for x in range(len(LogFC_sparklines_groups_names)):
                Gn = LogFC_sparklines_groups_names[x]
                Sparklines_Elements_Count.append(len(LogFC_sparklines_groups[Gn]))


            sum_d = sum([x**args.LogFC_sparklines_col_width_factor for x in Sparklines_Elements_Count])
            standart_width = 5.5 * (((sum(Sparklines_Elements_Count) / len(Sparklines_Elements_Count)) /  3)**0.7)

            Sparklines_Column_Widths = [(x ** args.LogFC_sparklines_col_width_factor) / sum_d * standart_width * len(LogFC_sparklines_groups_names) for x in Sparklines_Elements_Count]
            if args.LogFC_sparklines_equal_barwidths:
                Sparklines_Column_Widths = [x * 30 / max(Sparklines_Column_Widths + [30]) for x in Sparklines_Column_Widths]
            else:
                Sparklines_Column_Widths = [min(27, max(5, x)) for x in Sparklines_Column_Widths]

            for x in range(len(LogFC_sparklines_groups_names)):
                Gn = LogFC_sparklines_groups_names[x]
                LogFC_profile_cols = LogFC_sparklines_groups[Gn]

                start_ColN = LogFC_adj_starts_ColN[x]
                end_ColN = LogFC_adj_starts_ColN[x] + len(LogFC_sparklines_groups[Gn]) - 1
                all_SL_Locations = [xl_range(c, LogFC_sparklines_start_col - 1 + x + args.insert_empty_col_before_LogFC_sparklines,
                                             c, LogFC_sparklines_start_col - 1 + x + args.insert_empty_col_before_LogFC_sparklines) for c in range(args.Predictor_rows_count + 2, Last_Row)]
                all_SL_Ranges = [xl_range(c, start_ColN, c, end_ColN) for c in range(args.Predictor_rows_count + 2, Last_Row)]

                sheet.conditional_format(args.Predictor_rows_count + 1, start_ColN, Last_Row, end_ColN,
                                         {'type': '3_color_scale',
                                          'min_color': "#06ade0", 'mid_color': "#ffffff",
                                          'max_color': "#efb00e",
                                          'min_type': 'num', 'mid_type': 'num',
                                          'max_type': 'num',
                                          'min_value': -args.LogFC_sparklines_axis_limit - 0.01, 'mid_value': 0,
                                          'max_value': args.LogFC_sparklines_axis_limit + 0.01})

                sheet.add_sparkline(1 + args.Predictor_rows_count, LogFC_sparklines_start_col - 1 + x + args.insert_empty_col_before_LogFC_sparklines,
                                    {'location': all_SL_Locations, 'range': all_SL_Ranges,
                                     'type': 'column', 'max': args.LogFC_sparklines_axis_limit + 0.01, 'min': -args.LogFC_sparklines_axis_limit - 0.01,
                                     'negative_points': True,
                                     'series_color': '#efb00e',
                                     'negative_color': '#1dc0f2'})

                all_LogFC_sparklines_cols += [LogFC_sparklines_start_col + x - 1 + args.insert_empty_col_before_LogFC_sparklines]
                if args.Insert_spasers_between_rows or args.create_White_background:
                    sheet.set_column(LogFC_sparklines_start_col + x - 1 + args.insert_empty_col_before_LogFC_sparklines,
                                     LogFC_sparklines_start_col + x - 1 + args.insert_empty_col_before_LogFC_sparklines, Sparklines_Column_Widths[x], sparklines_cell_format_white_bg)
                else:
                    sheet.set_column(LogFC_sparklines_start_col + x - 1 + args.insert_empty_col_before_LogFC_sparklines,
                                     LogFC_sparklines_start_col + x - 1 + args.insert_empty_col_before_LogFC_sparklines, Sparklines_Column_Widths[x], sparklines_cell_format_no_bg)


        if args.insert_empty_col_before_CPM_sparklines and args.CPM_sparklines_start_col is not None:
            if args.Insert_spasers_between_rows or args.create_White_background:
                sheet.set_column(CPM_sparklines_start_col - 1, CPM_sparklines_start_col - 1, 1, whitespace_format)
            else:
                sheet.set_column(CPM_sparklines_start_col - 1, CPM_sparklines_start_col - 1, 1)

        if args.insert_empty_col_before_LogFC_sparklines and args.LogFC_sparklines_start_col is not None:
            if args.Insert_spasers_between_rows or args.create_White_background:
                sheet.set_column(LogFC_sparklines_start_col - 1, LogFC_sparklines_start_col - 1, 1, whitespace_format)
            else:
                sheet.set_column(LogFC_sparklines_start_col - 1, LogFC_sparklines_start_col - 1, 1)

        for x in range(len(header_cells)):
            if header_cells[x] in Forced_ColWidths_by_ColNames:
                if args.Insert_spasers_between_rows or args.create_White_background:
                    sheet.set_column(WSparkLine_Offset(args, x), WSparkLine_Offset(args, x), Forced_ColWidths_by_ColNames[header_cells[x]], whitespace_format)
                else:
                    sheet.set_column(WSparkLine_Offset(args, x), WSparkLine_Offset(args, x), Forced_ColWidths_by_ColNames[header_cells[x]])

        if args.dist_matrix_mode:
            sheet.conditional_format(args.Predictor_rows_count + 1, 1, Last_Row, len(header_cells),
                                     {'type': '3_color_scale',
                                      'min_color': "#5fce7c", 'mid_color': "#ffe76d",
                                      'max_color': "#f8696b",
                                      'min_type': 'percentile', 'mid_type': 'percentile',
                                      'max_type': 'percentile',
                                      'min_value': 0, 'mid_value': 50,
                                      'max_value': 99})


        if args.freeze_cols_count > 0 or args.freeze_rows_count > 0:
            sheet.freeze_panes(args.freeze_rows_count, args.freeze_cols_count)

        # sheet.write(Last_Row + 1, 0, 'generated from: %s'%FN, italic)

        if args.additional_text_file is not None:
            if os.path.exists(args.additional_text_file):
                RowN = Last_Row + 1
                for L in open(args.additional_text_file, 'r').readlines():
                    cells = L.rstrip('\r\n\t ').split('\t')
                    for nC in range(len(cells)):
                        val = cells[nC]
                        try:
                            val = float(val)
                            sheet.write_number(RowN, WSparkLine_Offset(args, nC), val)
                        except ValueError:
                            sheet.write(RowN, WSparkLine_Offset(args, nC), cells[nC])

                    RowN += 1

        if not args.Generate_SingleBook:
            Workbook.close()
    if args.Generate_SingleBook:
        Workbook.close()

    if args.verbose:
        print('\rExcel workbook generation completed                               ')


