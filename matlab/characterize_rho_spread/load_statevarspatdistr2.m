function data = load_statevarspatdistr2(filename, startRow, endRow)
filename = 'Z:\damaskpdt\run\500g256x256x256\DAMASKPDT.SimID.250.Incr.250.EdgeDensitySpatDistr.csv';
delimiter = ';';
startRow = 4;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
data = [dataArray{1:end-1}];
%clearvars filename delimiter startRow formatSpec fileID dataArray ans;
