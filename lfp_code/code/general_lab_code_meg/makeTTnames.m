function TTnames = makeTTnames(tetNum)


TTnames = cell(20,1);
for u = 1:20
    TTnames{u} = ['TT' num2str(tetNum) '_' num2str(u) '.t'];
end %u


end %function