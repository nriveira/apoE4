function make_inFile
% function make_inFile
% 

curDir = pwd;

tmpInfo = cell(5,1);

tmpInfo{1} = [curDir '\TTList.txt'];

for b = 1:4
tmpInfo{b+1} = [curDir '\begin' num2str(b)];
end %begin


tmpTable = cell2table(tmpInfo);
writetable(tmpTable, 'inFile.txt', 'WriteVariableNames', 0)

end %function