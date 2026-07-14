function [] = pastePlotToExcel(cellToPaste, fileHandle, fig)
% This function will paste a plot to an OPEN excel file.
% The file must be open and will not close when the function terminates.
% cellToPaste will be the location for the upper-left corener of the plot
% fileHandel - of type actxserver('excel.application').
 print -dmeta;
 ActiveSheet  = fileHandle.ActiveSheet;
 ActiveSheetRange  = get(ActiveSheet,'Range',cellToPaste);
 ActiveSheetRange.Select;
 ActiveSheetRange.PasteSpecial; %.................Pasting the figure to the selected location
 end

