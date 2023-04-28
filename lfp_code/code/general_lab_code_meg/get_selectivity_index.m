function SI = get_selectivity_index(rm1, rm2)
% function get_selectivity_index(rm1, rm2)
% 
% PURPOSE:
%     Get the selectivity index of a cell between two environments.
%     Reference: Hwaun & Colgin 2019
% 
% INPUT:
%     rm1 = ratemap from environment 1
%     rm2 = ratemap from environment 2
% 
% OUTPUT:
%     SI = selectivity index, a score ranging from -1 to 1, with -1
%        indicating that the cell ws exclusively active in environment 1,
%        and 2 indicating that the cell was exclusively active in
%        environment 2.
% 
% MMD
% 08/2022
% Colgin Lab

% remove portions of the map that the rat did not visit in both environments                        
goodInds = ~isnan(rm1+rm2);
rm1 = rm1(goodInds);
rm2 = rm2(goodInds);

%get mean firing rates
mu1 = mean(rm1(:));
mu2 = mean(rm2(:));

%calculate SI
SI = (mu2 - mu1) ./ (mu2 + mu1);



end %function