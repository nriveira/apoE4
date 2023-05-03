opa_build_struct;

%% Sort all of the data into comparable groups
opa_sorted = {};
groups = unique({opa_struct.group});
for g = 1:length(groups)
    opa_sorted(g).name = groups{g};
    rats = unique({opa_struct(strcmp({opa_struct.group}, groups{g})).name});
    for r = 1:length(rats)
        opa_sorted(g).rat(r).name = rats{r};
        opa_rat_struct = opa_struct(strcmp({opa_struct.name}, rats{r}) & strcmp({opa_struct.group}, groups{g}));
        sessions = unique({opa_rat_struct.sessionType});
        for s = 1:length(sessions)
            opa_sess_struct = opa_rat_struct(strcmp({opa_rat_struct.sessionType}, sessions{s}));
            opa_sorted(g).rat(r).session(s).name = sessions{s};
            opa_sorted(g).rat(r).session(s).session = opa_sess_struct;

            % Exclude practice sessions
            if(length(opa_sess_struct) > 2)
                s1_first3_objA = opa_sess_struct(strcmp({opa_sess_struct.sessionNumber}, 'S1')).first3_objA;
                s1_first3_objB = opa_sess_struct(strcmp({opa_sess_struct.sessionNumber}, 'S1')).first3_objB;
                s2_first3_objA = opa_sess_struct(strcmp({opa_sess_struct.sessionNumber}, 'S2')).first3_objA;
                s2_first3_objB = opa_sess_struct(strcmp({opa_sess_struct.sessionNumber}, 'S2')).first3_objB;
                s3_first3_objA = opa_sess_struct(strcmp({opa_sess_struct.sessionNumber}, 'S3')).first3_objA;
                s3_first3_objB = opa_sess_struct(strcmp({opa_sess_struct.sessionNumber}, 'S3')).first3_objB;
    
                % Discrimination index of S2 vs S1
                opa_sorted(g).rat(r).session(s).s2vs1_objA_discrim = s2_first3_objA / (s2_first3_objA + s1_first3_objA);
                opa_sorted(g).rat(r).session(s).s2vs1_objB_discrim = s2_first3_objB / (s2_first3_objB + s1_first3_objB);
                
                % Discrimination index of S2 vs S1+S3/2 (Average time)
                opa_sorted(g).rat(r).session(s).s2vs1s3_objA = s2_first3_objA / (0.5*(s1_first3_objA+s3_first3_objA) + s2_first3_objA);
                opa_sorted(g).rat(r).session(s).s2vs1s3_objB = s2_first3_objB / (0.5*(s1_first3_objB+s3_first3_objB) + s2_first3_objB);

                % Discrimination index of S2 vs S1+S3/2 (Average DI)
                di_s1_objA = s2_first3_objA / (s2_first3_objA + s1_first3_objA);
                di_s3_objA = s2_first3_objA / (s2_first3_objA + s3_first3_objA);

                di_s1_objB = s2_first3_objA / (s2_first3_objA + s1_first3_objA);
                di_s3_objB = s2_first3_objB / (s2_first3_objB + s3_first3_objB);

                opa_sorted(g).rat(r).session(s).s2vs1s3_objA_avgDI = 0.5*(di_s1_objA + di_s3_objA);
                opa_sorted(g).rat(r).session(s).s2vs1s3_objB_avgDI = 0.5*(di_s1_objB + di_s3_objB);
            end
        end
    end
end