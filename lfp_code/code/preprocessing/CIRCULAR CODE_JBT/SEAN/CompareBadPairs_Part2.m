clear all; close all; clc
load('I:\JBT_RERUN\bad_cell_list_JBT.mat')
load('I:\JBT_RERUN\Bad_cell_debug_SGT.mat')
for ii=1:size(jbt_badcell,1)
	for jj=1:size(bad_cell_list,1)
		if(sum(jbt_badcell(ii,:)==bad_cell_list(jj,:))==5)
			jbt_badcell(ii,:)=nan(1,5);
			bad_cell_list(jj,:)=nan(1,5);
			break;
		end
	end
end
jbt_idx=[];
for aa=1:size(jbt_badcell,1)

	if isnan(jbt_badcell(aa,1)) 
		jbt_idx=[jbt_idx, aa];
	end

end

sgt_idx=[];
for aa=1:size(bad_cell_list,1)
	if isnan(bad_cell_list(aa,1))
		sgt_idx=[sgt_idx,aa];
	end
end

jbt_badcell(jbt_idx,:)=[];
bad_cell_list(sgt_idx,:)=[];