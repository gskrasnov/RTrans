def Create_workbook_Formats(Workbook):
    bold = Workbook.add_format({'bold': True, 'italic': False})
    bold_light = Workbook.add_format({'bold': True, 'italic': True, 'font_color': '#BBBBBB'})
    bold_ww = Workbook.add_format({'bold': True, 'italic': False})
    bold_ww.set_text_wrap()
    italic = Workbook.add_format({'bold': False, 'italic': True})
    bold_italic = Workbook.add_format({'bold': True, 'italic': True})
    bold_italic_ww = Workbook.add_format({'bold': True, 'italic': True})
    bold_italic_ww.set_text_wrap()
    format0 = Workbook.add_format({'num_format':'0', 'valign': 'vcenter'})
    format1 = Workbook.add_format({'num_format':'0.0', 'valign': 'vcenter'})
    format2 = Workbook.add_format({'num_format':'0.00', 'valign': 'vcenter'})
    format1_narrow = Workbook.add_format({'num_format': '0.0', 'font_size':8, 'valign': 'vcenter', 'align': 'center'})
    format2_narrow = Workbook.add_format({'num_format': '0.00', 'font_size':8, 'valign': 'vcenter', 'align': 'center'})
    format1_center = Workbook.add_format({'num_format':'0.0', 'align': 'center', 'valign': 'vcenter'})
    format2_center = Workbook.add_format({'align': 'center', 'num_format':'0.00', 'valign': 'vcenter'})
    format2_center_bold = Workbook.add_format({'align': 'center', 'bold': True, 'num_format': '0.00', 'valign': 'vcenter'})
    format_center = Workbook.add_format({'align': 'center', 'valign': 'vcenter'})
    format0_center = Workbook.add_format({'align': 'center', 'valign': 'vcenter', 'num_format':'0'})
    format_perc = Workbook.add_format({'num_format': '0.00%'})

    default_Rank_format = Workbook.add_format({'align':'center', 'num_format':'0.0', 'bg_color':'#ffffff', 'font_size':8})

    headers_format = Workbook.add_format({'bold': True, 'italic': False, 'border': True, 'align': 'center', 'valign': 'vcenter', 'text_wrap':True})
    small_headers_format = Workbook.add_format({'bold': False, 'italic': False, 'border': True, 'align': 'center', 'valign': 'vcenter', 'size': 8, 'text_wrap':True})
    secondary_headers_format = Workbook.add_format({'bold': False, 'italic': False, 'border': True, 'align': 'center', 'valign': 'vcenter', 'text_wrap':True})
    merge_format = Workbook.add_format({'bold': True, 'italic': True, 'border': True, 'align': 'center', 'valign': 'vcenter', 'text_wrap':True})

    rotated_format = Workbook.add_format({'bold': False, 'italic': False, 'font_size':10, 'text_wrap':False, 'rotation':90, 'align':'center'})
    rotated_format_underline = Workbook.add_format({'bold': False, 'italic': False, 'font_size':10, 'text_wrap':False, 'rotation':90, 'align':'center', 'bottom':2, 'bottom_color':'#000000'})
    rotated_format_center = Workbook.add_format({'bold': False, 'italic': False, 'align': 'center', 'valign': 'vcenter', 'font_size':9, 'text_wrap':True, 'rotation':90})

    whitespace_format = Workbook.add_format({'bg_color':'#ffffff'})
    whitespace_format__right_border = Workbook.add_format({'bg_color':'#ffffff', 'right':2})
    whitespace_format__bottom_border = Workbook.add_format({'bg_color':'#ffffff', 'bottom':2})
    whitespace_format1 = Workbook.add_format({'bg_color':'#ffffff', 'num_format': '0.0', 'valign': 'vcenter'})

    percentage_format0 = Workbook.add_format({'num_format': '0%'})

    whitespace_format__header = Workbook.add_format({'bg_color': '#ffffff', 'align': 'center', 'valign': 'bottom'})
    empty_logfc_format = Workbook.add_format({'bg_color':'#dddddd', 'num_format': '0.00', 'font_size':8})

    sparklines_cell_format = Workbook.add_format({'bg_color':'#ffffff', 'left':2, 'right':2})
    sparklines_cell_format_no_bg = Workbook.add_format({'left':2, 'right':2, 'left_color':'#000000', 'right_color':'#000000'})
    sparklines_cell_format_white_bg = Workbook.add_format({'bg_color':'#ffffff', 'left':2, 'right':2, 'left_color':'#000000', 'right_color':'#000000'})

    heatmap_cell_format_left_border = Workbook.add_format({'left':1, 'left_color':'#000000'})
    heatmap_cell_format_right_border = Workbook.add_format({'right':1, 'right_color':'#000000'})
    heatmap_cell_format_left_right_borders = Workbook.add_format({'left':1, 'right':1, 'left_color':'#000000', 'right_color':'#000000'})

    pvalue_formats = []
    for digits in range(5):
        pvalue_formats.append(Workbook.add_format())
        pvalue_formats[-1].set_num_format('0.%s'%('0'*digits))

    format_quasi_invisible = Workbook.add_format({'font_color':'#ffffff', 'font_size':1, 'align': 'center', 'valign': 'vcenter', 'num_format':' '})
    format_quasi_invisible_left_border = Workbook.add_format({'font_color':'#ffffff', 'font_size':1, 'align': 'center', 'valign': 'vcenter', 'num_format':' ', 'left':1, 'left_color':'#000000'})
    format_quasi_invisible_right_border = Workbook.add_format({'font_color':'#ffffff', 'font_size':1, 'align': 'center', 'valign': 'vcenter', 'num_format':' ', 'right':1, 'right_color':'#000000'})
    format_quasi_invisible_left_right_borders = Workbook.add_format({'font_color':'#ffffff', 'font_size':1, 'align': 'center', 'valign': 'vcenter', 'num_format':' ', 'left':1, 'right':1, 'left_color':'#000000', 'right_color':'#000000'})

    return (bold, bold_light, bold_ww, italic, bold_italic, bold_italic_ww, format0, format1, format2, format1_narrow, format2_narrow, format1_center, format2_center,
        format2_center_bold, format_center, format0_center, format_perc, default_Rank_format, headers_format, small_headers_format,
        secondary_headers_format, merge_format, rotated_format, rotated_format_underline, rotated_format_center, whitespace_format,
        whitespace_format__right_border, whitespace_format__bottom_border,
        whitespace_format1, percentage_format0, whitespace_format__header, empty_logfc_format, sparklines_cell_format, sparklines_cell_format_no_bg,
        sparklines_cell_format_white_bg, pvalue_formats, format_quasi_invisible,
        format_quasi_invisible_left_border, format_quasi_invisible_right_border, format_quasi_invisible_left_right_borders,
        heatmap_cell_format_left_border, heatmap_cell_format_right_border, heatmap_cell_format_left_right_borders)

def Create_workbook_Formats_whitespace(Workbook):
    bold = Workbook.add_format({'bold': True, 'italic': False, 'bg_color': '#ffffff'})
    bold_light = Workbook.add_format({'bold': True, 'italic': True, 'font_color': '#BBBBBB', 'bg_color': '#ffffff'})
    bold_ww = Workbook.add_format({'bold': True, 'italic': False, 'bg_color': '#ffffff'})
    bold_ww.set_text_wrap()
    italic = Workbook.add_format({'bold': False, 'italic': True, 'bg_color': '#ffffff'})
    bold_italic = Workbook.add_format({'bold': True, 'italic': True, 'bg_color': '#ffffff'})
    bold_italic_ww = Workbook.add_format({'bold': True, 'italic': True, 'bg_color': '#ffffff'})
    bold_italic_ww.set_text_wrap()
    format0 = Workbook.add_format({'num_format':'0', 'valign': 'vcenter', 'bg_color': '#ffffff'})
    format1 = Workbook.add_format({'num_format':'0.0', 'valign': 'vcenter', 'bg_color': '#ffffff'})
    format2 = Workbook.add_format({'num_format':'0.00', 'valign': 'vcenter', 'bg_color': '#ffffff'})
    format1_narrow = Workbook.add_format({'num_format': '0.0', 'font_size':8, 'valign': 'vcenter', 'align': 'center', 'bg_color': '#ffffff'})
    format2_narrow = Workbook.add_format({'num_format': '0.00', 'font_size':8, 'valign': 'vcenter', 'align': 'center', 'bg_color': '#ffffff'})
    format1_center = Workbook.add_format({'num_format':'0.0', 'align': 'center', 'valign': 'vcenter', 'bg_color': '#ffffff'})
    format2_center = Workbook.add_format({'align': 'center', 'num_format':'0.00', 'valign': 'vcenter', 'bg_color': '#ffffff'})
    format2_center_bold = Workbook.add_format({'align': 'center', 'bold': True, 'num_format': '0.00', 'valign': 'vcenter', 'bg_color': '#ffffff'})
    format_center = Workbook.add_format({'align': 'center', 'valign': 'vcenter', 'bg_color': '#ffffff'})
    format0_center = Workbook.add_format({'align': 'center', 'valign': 'vcenter', 'num_format':'0', 'bg_color': '#ffffff'})
    format_perc = Workbook.add_format({'num_format': '0.00%', 'bg_color': '#ffffff'})

    default_Rank_format = Workbook.add_format({'align':'center', 'num_format':'0.0', 'bg_color':'#ffffff', 'font_size':8})

    headers_format = Workbook.add_format({'bold': True, 'italic': False, 'border': True, 'align': 'center', 'valign': 'vcenter', 'text_wrap':True, 'bg_color': '#ffffff'})
    small_headers_format = Workbook.add_format({'bold': False, 'italic': False, 'border': True, 'align': 'center', 'valign': 'vcenter', 'size': 8, 'text_wrap':True, 'bg_color': '#ffffff'})
    secondary_headers_format = Workbook.add_format({'bold': False, 'italic': False, 'border': True, 'align': 'center', 'valign': 'vcenter', 'text_wrap':True, 'bg_color': '#ffffff'})
    merge_format = Workbook.add_format({'bold': True, 'italic': True, 'border': True, 'align': 'center', 'valign': 'vcenter', 'text_wrap':True, 'bg_color': '#ffffff'})

    rotated_format = Workbook.add_format({'bold': False, 'italic': False, 'font_size':10, 'text_wrap':False, 'rotation':90, 'align':'center', 'bg_color': '#ffffff'})
    rotated_format_underline = Workbook.add_format({'bold': False, 'italic': False, 'font_size':10, 'text_wrap':False, 'rotation':90, 'align':'center', 'bottom':2, 'bottom_color':'#000000', 'bg_color': '#ffffff'})
    rotated_format_center = Workbook.add_format({'bold': False, 'italic': False, 'align': 'center', 'valign': 'vcenter', 'font_size':9, 'text_wrap':True, 'rotation':90, 'bg_color': '#ffffff'})

    whitespace_format = Workbook.add_format({'bg_color':'#ffffff'})
    whitespace_format__right_border = Workbook.add_format({'bg_color':'#ffffff', 'right':2})
    whitespace_format__bottom_border = Workbook.add_format({'bg_color':'#ffffff', 'bottom':2})
    whitespace_format1 = Workbook.add_format({'bg_color':'#ffffff', 'num_format': '0.0', 'valign': 'vcenter'})

    percentage_format0 = Workbook.add_format({'num_format': '0%', 'bg_color': '#ffffff'})

    whitespace_format__header = Workbook.add_format({'bg_color': '#ffffff', 'align': 'center', 'valign': 'bottom'})
    empty_logfc_format = Workbook.add_format({'bg_color':'#dddddd', 'num_format': '0.00', 'font_size':8})

    sparklines_cell_format = Workbook.add_format({'bg_color':'#ffffff', 'left':2, 'right':2})
    sparklines_cell_format_no_bg = Workbook.add_format({'left':2, 'right':2, 'left_color':'#000000', 'right_color':'#000000', 'bg_color': '#ffffff'})
    sparklines_cell_format_white_bg = Workbook.add_format({'bg_color':'#ffffff', 'left':2, 'right':2, 'left_color':'#000000', 'right_color':'#000000'})
    
    heatmap_cell_format_left_border = Workbook.add_format({'bg_color':'#ffffff', 'left':1, 'left_color':'#000000'})
    heatmap_cell_format_right_border = Workbook.add_format({'bg_color':'#ffffff', 'right':1, 'right_color':'#000000'})
    heatmap_cell_format_left_right_borders = Workbook.add_format({'bg_color':'#ffffff', 'left':1, 'right':1, 'left_color':'#000000', 'right_color':'#000000'})

    pvalue_formats = []
    for digits in range(5):
        pvalue_formats.append(Workbook.add_format({'bg_color': '#ffffff'}))
        pvalue_formats[-1].set_num_format('0.%s'%('0'*digits))

    format_quasi_invisible = Workbook.add_format({'bg_color':'#ffffff', 'font_color':'#ffffff', 'font_size':1, 'align': 'center', 'valign': 'vcenter', 'num_format':' '})
    format_quasi_invisible_left_border = Workbook.add_format({'bg_color':'#ffffff', 'font_color':'#ffffff', 'font_size':1, 'align': 'center', 'valign': 'vcenter', 'num_format':' ', 'left':1, 'left_color':'#000000'})
    format_quasi_invisible_right_border = Workbook.add_format({'bg_color':'#ffffff', 'font_color':'#ffffff', 'font_size':1, 'align': 'center', 'valign': 'vcenter', 'num_format':' ', 'right':1, 'right_color':'#000000'})
    format_quasi_invisible_left_right_borders = Workbook.add_format({'bg_color':'#ffffff', 'font_color':'#ffffff', 'font_size':1, 'align': 'center', 'valign': 'vcenter', 'num_format':' ', 'left':1, 'right':1, 'left_color':'#000000', 'right_color':'#000000'})

    return (bold, bold_light, bold_ww, italic, bold_italic, bold_italic_ww, format0, format1, format2, format1_narrow, format2_narrow, format1_center, format2_center,
        format2_center_bold, format_center, format0_center, format_perc, default_Rank_format, headers_format, small_headers_format,
        secondary_headers_format, merge_format, rotated_format, rotated_format_underline, rotated_format_center, whitespace_format,
        whitespace_format__right_border, whitespace_format__bottom_border,
        whitespace_format1, percentage_format0, whitespace_format__header, empty_logfc_format, sparklines_cell_format, sparklines_cell_format_no_bg,
        sparklines_cell_format_white_bg, pvalue_formats, format_quasi_invisible,
        format_quasi_invisible_left_border, format_quasi_invisible_right_border, format_quasi_invisible_left_right_borders,
        heatmap_cell_format_left_border, heatmap_cell_format_right_border, heatmap_cell_format_left_right_borders)