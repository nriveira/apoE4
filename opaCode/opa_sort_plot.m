% Cycle through data struct to get plottable values

for vof = {'s2vs1_objA_discrim', 's2vs1_objB_discrim', 's2vs1s3_objA', 's2vs1s3_objB', 's2vs1s3_objA_avgDI', 's2vs1s3_objB_avgDI'}
var_of_interest = vof{1};
for g = 1:length(opa_sorted)
    % Discrimination index values for all rats (objectA/objectB)
    f_discrim = nan(length(opa_sorted(g).rat), 1);
    nl_discrim = nan(length(opa_sorted(g).rat), 1);
    no_discrim = nan(length(opa_sorted(g).rat), 1);
    nonl_discrim = nan(length(opa_sorted(g).rat), 1);

    for r = 1:length(opa_sorted(g).rat)
        f_discrim(r) = mean([opa_sorted(g).rat(r).session(contains({opa_sorted(g).rat(r).session.name}, 'F')).(var_of_interest)]);
        nl_discrim(r) = opa_sorted(g).rat(r).session(strcmp({opa_sorted(g).rat(r).session.name}, 'NL')).(var_of_interest);
        no_discrim(r) = opa_sorted(g).rat(r).session(strcmp({opa_sorted(g).rat(r).session.name}, 'NO')).(var_of_interest);
        nonl_discrim(r) = opa_sorted(g).rat(r).session(strcmp({opa_sorted(g).rat(r).session.name}, 'NO_NL')).(var_of_interest);
    end
    
    opa_sorted(g).f_discrim_mean = mean(f_discrim, 'omitnan');
    opa_sorted(g).f_discrim_sterr = std(f_discrim) / sqrt(sum(~isnan(f_discrim)));

    opa_sorted(g).nl_discrim_mean = mean(nl_discrim, 'omitnan');
    opa_sorted(g).nl_discrim_sterr = std(nl_discrim) / sqrt(sum(~isnan(nl_discrim)));

    opa_sorted(g).no_discrim_mean = mean(no_discrim, 'omitnan');
    opa_sorted(g).no_discrim_sterr = std(no_discrim) / sqrt(sum(~isnan(no_discrim)));

    opa_sorted(g).nonl_discrim_mean = mean(nonl_discrim, 'omitnan');
    opa_sorted(g).nonl_discrim_sterr = std(nonl_discrim, 'omitnan') / sqrt(sum(~isnan(nonl_discrim)));
end

figure(1); clf; hold on;
x = 1:4;
discrim_index_e3 = [opa_sorted(1).f_discrim_mean, ...
                    opa_sorted(1).nl_discrim_mean, ...
                    opa_sorted(1).no_discrim_mean, ...
                    opa_sorted(1).nonl_discrim_mean];

discrim_index_e4 = [opa_sorted(2).f_discrim_mean, ...
                    opa_sorted(2).nl_discrim_mean, ...
                    opa_sorted(2).no_discrim_mean, ...
                    opa_sorted(2).nonl_discrim_mean];

bar(x, [discrim_index_e3; discrim_index_e4])
ylabel('Discrimination index')
xticklabels({'F', '', 'NL', '', 'NO', '', 'NO+NL'})
title(strrep(var_of_interest, '_', ' '))
legend({'E3', 'E4'}, 'Location', 'southwest')
% saveas(gcf, ['../figures/discrimination_index/20230502_' var_of_interest '.png'])

end