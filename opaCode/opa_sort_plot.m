% Cycle through data struct to get plottable values

%for vof = {'s2vs1_objA_discrim', 's2vs1_objB_discrim', 's2vs1s3_objA', 's2vs1s3_objB', 's2vs1s3_objA_avgDI', 's2vs1s3_objB_avgDI'}
%var_of_interest = vof{1};
var_of_interest = 's2vs1_objA_discrim';
for g = 1:length(opa_sorted)
    % Discrimination index values for all rats (objectA/objectB)
    f_discrim = nan(length(opa_sorted(g).rat), 1);
    nl_discrim = nan(length(opa_sorted(g).rat), 1);
    no_discrim = nan(length(opa_sorted(g).rat), 1);
    nonl_discrim = nan(length(opa_sorted(g).rat), 1);

    f_prac = nan(length(opa_sorted(g).rat), 3);
    nl_prac = nan(length(opa_sorted(g).rat), 3);
    no_prac = nan(length(opa_sorted(g).rat), 3);
    nonl_prac = nan(length(opa_sorted(g).rat), 3);

    for r = 1:length(opa_sorted(g).rat)
        familiar_sess = opa_sorted(g).rat(r).session((contains({opa_sorted(g).rat(r).session.name}, 'F')) & (~{isempty({opa_sorted(g).rat(r).session.(var_of_interest)})}));
        f_sess_concat = [];
        for i = 1:length(familiar_sess)
            f_sess_concat = [f_sess_concat; familiar_sess(i).(var_of_interest)];
        end

        f_discrim(r) = mean(f_sess_concat);
        nl_discrim(r) = opa_sorted(g).rat(r).session(strcmp({opa_sorted(g).rat(r).session.name}, 'NL')).(var_of_interest);
        no_discrim(r) = opa_sorted(g).rat(r).session(strcmp({opa_sorted(g).rat(r).session.name}, 'NO')).(var_of_interest);
        nonl_discrim(r) = opa_sorted(g).rat(r).session(strcmp({opa_sorted(g).rat(r).session.name}, 'NO_NL')).(var_of_interest);

        f_prac_concat = [];
        for i = 1:length(familiar_sess)
            f_prac_concat = [f_prac_concat; familiar_sess(i).objA_discrim];
        end
        f_prac(r, :) = mean(f_prac_concat);
        nl_prac(r,:) = opa_sorted(g).rat(r).session(strcmp({opa_sorted(g).rat(r).session.name}, 'NL')).objA_discrim;
        no_prac(r,:) = opa_sorted(g).rat(r).session(strcmp({opa_sorted(g).rat(r).session.name}, 'NO')).objA_discrim;
        nonl_prac(r,:) = opa_sorted(g).rat(r).session(strcmp({opa_sorted(g).rat(r).session.name}, 'NO_NL')).objA_discrim;
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

%end