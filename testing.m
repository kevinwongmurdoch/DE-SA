%Step 1 : Read from file
fileexcel = 'C:\Users\32599642\Desktop\exceldata.xls';
% T = readtable(fileexcel);
B = xlsread(fileexcel);

%Step 2: Append new data to this table
newADDC = zeros(20,1);
newFMeasure = zeros(20,1);
size(T)
oldADDC= B(:,1)
oldFMeasure= B(:,2)

completeADDC = [oldADDC newADDC]
completeFMeasure = [oldFMeasure newFMeasure]


