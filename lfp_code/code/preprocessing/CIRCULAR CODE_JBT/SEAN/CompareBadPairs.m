fid0 = fopen('I:\JBT_RERUN\Discarded_MEC_Cell_Info.txt','r');
tmp = fgetl(fid0);
tmp = fgetl(fid0);

%% rat+session to aa_MAIN
	rsp = nan(6,7);
	rsp(1,1) = 23;
	rsp(1,2) = 24;
	rsp(1,3) = 25;
	rsp(1,4) = 26;
	rsp(1,5) = 27;
	rsp(1,6) = 28;
	rsp(1,7) = 29;
	rsp(2,1) = 1;
	rsp(2,2) = 2;
	rsp(2,3) = 3;
	rsp(2,4) = 4;
	rsp(2,5) = 5;
	rsp(3,1) = 6;
	rsp(3,2) = 7;
	rsp(3,3) = 8;
	rsp(3,4) = 9;
	rsp(4,1) = 10;
	rsp(5,1) = 11;
	rsp(5,2) = 12;
	rsp(5,3) = 13;
	rsp(5,4) = 14;
	rsp(5,5) = 15;
	rsp(5,6) = 16;
	rsp(5,7) = 17;
	rsp(6,1) = 18;
	rsp(6,2) = 19;
	rsp(6,3) = 20;
	rsp(6,4) = 21;
	rsp(6,5) = 22;

jbt_ctr=0;
    
while ~feof(fid0)
	tmp = strip(tmp);
	tmp = lower(tmp);
    fprintf('Checking line! [%s]\n',tmp);
	%% check to see if we're beginning a text block
	if(numel(tmp) > 8 && strcmp(tmp(1:8),'bad pair'))
		rat_sess_templine = strip(lower(fgetl(fid0))); %get the line with the rat and session number
		cell1_templine = strip(lower(fgetl(fid0)));  %get the line with cell 1
		cell2_templine = strip(lower(fgetl(fid0)));  %get the line with cell 2
		
		%pull the relevant coords
		tmp2 = strsplit(rat_sess_templine,{'rat ',',',' ses '});
		rsp_index_1 = str2num(tmp2{2});
		rsp_index_2 = str2num(tmp2{3});
		
		tmp3 = strsplit(cell1_templine,{': t',',',' u'});
		tt_1_n = str2num(tmp3{2});
		cell_1_n = str2num(tmp3{3});
		tmp3 = strsplit(cell2_templine,{': t',',',' u'});
		tt_2_n = str2num(tmp3{2});
		cell_2_n = str2num(tmp3{3});
		
		jbt_ctr=jbt_ctr+1;
		jbt_badcell(jbt_ctr,:) = [rsp(rsp_index_1,rsp_index_2), tt_1_n, cell_1_n, tt_2_n, cell_2_n];
	end
	tmp=fgetl(fid0);
end	
		[~,idx] = sort(jbt_badcell(:,1));
jbt_badcell = jbt_badcell(idx,:);
save('I:\JBT_RERUN\bad_cell_list_JBT.mat','jbt_badcell');