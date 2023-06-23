% Cycle through data struct to get plottable values
variables_of_interest = {'s2vs1_objA_discrim', 's2vs1_objB_discrim', 's2vs1s3_objA', 's2vs1s3_objB', 's2vs1s3_objA_avgDI', 's2vs1s3_objB_avgDI'};
for vof = 1:length(variables_of_interest)

var_of_interest = variables_of_interest{vof};
for g = 1:length(opa_sorted)
    % Discrimination index values for all rats (objectA/objectB)
    f_discrim = nan(length(opa_sorted(g).rat), 1);
    nl_discrim = nan(length(opa_sorted(g).rat), 1);
    no_discrim = nan(length(opa_sorted(g).rat), 1);
    nonl_discrim = nan(length(opa_sorted(g).rat), 1);

    for r = 1:length(opa_sorted(g).rat)
        familiar_sess = opa_sorted(g).rat(r).session(contains({opa_sorted(g).rat(r).session.name}, 'F'));
        
        f_sess_concat = [];

        for i = 1:length(familiar_sess)
            f_sess_concat = [f_sess_concat; mean([familiar_sess(i).(var_of_interest)])];
        end

        f_discrim(r) = mean(f_sess_concat, 'omitnan');
        nl_discrim(r) = opa_sorted(g).rat(r).session(strcmp({opa_sorted(g).rat(r).session.name}, 'NL')).(var_of_interest);
        no_discrim(r) = opa_sorted(g).rat(r).session(strcmp({opa_sorted(g).rat(r).session.name}, 'NO')).(var_of_interest);
        nonl_discrim(r) = opa_sorted(g).rat(r).session(strcmp({opa_sorted(g).rat(r).session.name}, 'NO_NL')).(var_of_interest);
    end

    opa_sorted(g).stats(vof).var_of_interest = var_of_interest;
    opa_sorted(g).stats(vof).f_discrim = f_discrim; 
    opa_sorted(g).stats(vof).f_discrim_mean = mean(f_discrim, 'omitnan');
    opa_sorted(g).stats(vof).f_discrim_sterr = std(f_discrim, 'omitnan') / sqrt(sum(~isnan(f_discrim)));

    opa_sorted(g).stats(vof).nl_discrim = nl_discrim;
    opa_sorted(g).stats(vof).nl_discrim_mean = mean(nl_discrim, 'omitnan');
    opa_sorted(g).stats(vof).nl_discrim_sterr = std(nl_discrim, 'omitnan') / sqrt(sum(~isnan(nl_discrim)));

    opa_sorted(g).stats(vof).no_discrim = no_discrim;
    opa_sorted(g).stats(vof).no_discrim_mean = mean(no_discrim, 'omitnan');
    opa_sorted(g).stats(vof).no_discrim_sterr = std(no_discrim, 'omitnan') / sqrt(sum(~isnan(no_discrim)));

    opa_sorted(g).stats(vof).nonl_discrim = nonl_discrim;
    opa_sorted(g).stats(vof).nonl_discrim_mean = mean(nonl_discrim, 'omitnan');
    opa_sorted(g).stats(vof).nonl_discrim_sterr = std(nonl_discrim, 'omitnan') / sqrt(sum(~isnan(nonl_discrim)));

    % Stats
    opa_sorted(g).stats(vof).discrim_stats = [opa_sorted(g).stats(vof).f_discrim, opa_sorted(g).stats(vof).nl_discrim, opa_sorted(g).stats(vof).no_discrim, opa_sorted(g).stats(vof).nonl_discrim];
    [opa_sorted(g).stats(vof).p, ~, opa_sorted(g).stats(vof).stats] = anova1(opa_sorted(g).stats(vof).discrim_stats);
    opa_sorted(g).stats(vof).results = multcompare(opa_sorted(g).stats(vof).stats);
    opa_sorted(g).stats(vof).tbl = array2table(opa_sorted(g).stats(vof).results,"VariableNames", ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    close all;
end



%
figure(1); clf; hold on;
x = 1:4;
discrim_index_e3 = [opa_sorted(1).stats(vof).f_discrim_mean, ...
                    opa_sorted(1).stats(vof).nl_discrim_mean, ...
                    opa_sorted(1).stats(vof).no_discrim_mean, ...
                    opa_sorted(1).stats(vof).nonl_discrim_mean];

discrim_error_e3 = [opa_sorted(1).stats(vof).f_discrim_sterr, ...
                     opa_sorted(1).stats(vof).nl_discrim_sterr, ...
                     opa_sorted(1).stats(vof).no_discrim_sterr, ...
                     opa_sorted(1).stats(vof).nonl_discrim_sterr];

discrim_value_e3 = [opa_sorted(1).stats(vof).f_discrim, ...
                     opa_sorted(1).stats(vof).nl_discrim, ...
                     opa_sorted(1).stats(vof).no_discrim, ...
                     opa_sorted(1).stats(vof).nonl_discrim];


discrim_index_e4 = [opa_sorted(2).stats(vof).f_discrim_mean, ...
                    opa_sorted(2).stats(vof).nl_discrim_mean, ...
                    opa_sorted(2).stats(vof).no_discrim_mean, ...
                    opa_sorted(2).stats(vof).nonl_discrim_mean];

discrim_error_e4 = [opa_sorted(2).stats(vof).f_discrim_sterr, ...
                     opa_sorted(2).stats(vof).nl_discrim_sterr, ...
                     opa_sorted(2).stats(vof).no_discrim_sterr, ...
                     opa_sorted(2).stats(vof).nonl_discrim_sterr];

discrim_value_e4 = [opa_sorted(2).stats(vof).f_discrim, ...
                     opa_sorted(2).stats(vof).nl_discrim, ...
                     opa_sorted(2).stats(vof).no_discrim, ...
                     opa_sorted(2).stats(vof).nonl_discrim];

bar(x, [discrim_index_e3; discrim_index_e4])
errorbar(x-0.15, discrim_index_e3, discrim_error_e3, 'k', "LineStyle","none");
scatter(x-0.15, discrim_value_e3, '.k');

errorbar(x+0.15, discrim_index_e4, discrim_error_e4, 'k', "LineStyle","none")
scatter(x+0.15, discrim_value_e4, '.k');


ylabel('Discrimination index')
yline(0.5, '--')
xticklabels({'F', '', 'NL', '', 'NO', '', 'NO+NL'})
title(strrep(var_of_interest, '_', ' '))
legend({'E3', 'E4'})
saveas(gcf, ['../figures/ADGrant/' var_of_interest '.eps'])

end % End of analysis

%% Run stats on groups combined
for s = 1:length(opa_sorted(1).stats)
    stats(s).var = opa_sorted(2).stats(s).var_of_interest;
    e4_stats = opa_sorted(2).stats(s).discrim_stats;
    e4_stats(6,:) = nan;

    stats(s).combined_stats = [opa_sorted(1).stats(s).discrim_stats, e4_stats];
    [stats(s).p, stats(s).tbl, stats(s).stats] = anova1(stats(s).combined_stats);
    stats(s).results = multcompare(stats(s).stats);
    tbl = array2table(stats(s).results,"VariableNames", ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
    close all;
end