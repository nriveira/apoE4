opa_build_struct;
firstXframes = 1800 * 3; % 30 fps -> 1800 frames 60 seconds

%% Sort all of the data into comparable groups
opa_sorted = {};
groups = unique({opa_struct.group});
for g = 1:length(groups) % Cycle through both groups
    opa_sorted(g).name = groups{g};
    rats = unique({opa_struct(strcmp({opa_struct.group}, groups{g})).name});
    for r = 1:length(rats) % Cycle through each rat
        opa_sorted(g).rat(r).name = rats{r};

        % Isolate opa struct sessions and practices
        temp_rat_struct = opa_struct(strcmp({opa_struct.name}, rats{r}) & strcmp({opa_struct.group}, groups{g}));
        opa_rat_pract = temp_rat_struct(contains({temp_rat_struct.sessionType}, 'P'));
        opa_rat_famil = temp_rat_struct(contains({temp_rat_struct.sessionType}, 'F'));
        opa_rat_novel = temp_rat_struct((~contains({temp_rat_struct.sessionType}, 'P') & ~contains({temp_rat_struct.sessionType}, 'F')));

        novel = unique({opa_rat_novel.sessionType});
        fam = unique({opa_rat_famil.sessionType});
        temp_rat_struct = [];
        
        % Analyze all of the novel conditions
        for s = 1:length(novel) % Cycle through each novel condition
            % Analyze between session
            opa_sess_struct = opa_rat_novel(strcmp({opa_rat_novel.sessionType}, novel{s}));
            opa_sorted(g).rat(r).session(s).name = novel{s};
            opa_sorted(g).rat(r).session(s).session = opa_sess_struct;
            if(length(opa_sess_struct) > 2) % Only do sessions with S1, S2, and S3
                s1_first3_objA = opa_sess_struct(strcmp({opa_sess_struct.sessionNumber}, 'S1')).first3_objA;
                s1_first3_objB = opa_sess_struct(strcmp({opa_sess_struct.sessionNumber}, 'S1')).first3_objB;
                s2_first3_objA = opa_sess_struct(strcmp({opa_sess_struct.sessionNumber}, 'S2')).first3_objA;
                s2_first3_objB = opa_sess_struct(strcmp({opa_sess_struct.sessionNumber}, 'S2')).first3_objB;
                s3_first3_objA = opa_sess_struct(strcmp({opa_sess_struct.sessionNumber}, 'S3')).first3_objA;
                s3_first3_objB = opa_sess_struct(strcmp({opa_sess_struct.sessionNumber}, 'S3')).first3_objB;
    
                % Discrimination index of S2 vs S1
                opa_sorted(g).rat(r).session(s).s2vs1_objA_discrim = s2_first3_objA / (s2_first3_objA + s1_first3_objA);
                opa_sorted(g).rat(r).session(s).s2vs1_objB_discrim = s2_first3_objB / (s2_first3_objB + s1_first3_objB);

                % Discrimination index of S2 vs S3
                opa_sorted(g).rat(r).session(s).s2vs3_objA_discrim = s2_first3_objA / (s2_first3_objA + s3_first3_objA);
                opa_sorted(g).rat(r).session(s).s2vs3_objB_discrim = s2_first3_objB / (s2_first3_objB + s3_first3_objB);

                % Discrimination index of S1 vs S3
                opa_sorted(g).rat(r).session(s).s2vs3_objA_discrim = s1_first3_objA / (s1_first3_objA + s3_first3_objA);
                opa_sorted(g).rat(r).session(s).s2vs3_objB_discrim = s1_first3_objB / (s1_first3_objB + s3_first3_objB);
                
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

        % Analyze the familiar conditions
        for f = 1:length(fam)
            opa_rat_temp = opa_rat_famil(strcmp({opa_rat_famil.sessionType}, fam{f}));
            opa_sorted(g).rat(r).session(f+length(novel)).name = opa_rat_temp.sessionType;
            opa_sorted(g).rat(r).session(f+length(novel)).session = opa_rat_temp;

            if(length(opa_rat_temp) > 2) % Only do sessions with S1, S2, and S3
                s1_first3_objA = opa_rat_temp(strcmp({opa_rat_temp.sessionNumber}, 'S1')).first3_objA;
                s1_first3_objB = opa_rat_temp(strcmp({opa_rat_temp.sessionNumber}, 'S1')).first3_objB;
                s2_first3_objA = opa_rat_temp(strcmp({opa_rat_temp.sessionNumber}, 'S2')).first3_objA;
                s2_first3_objB = opa_rat_temp(strcmp({opa_rat_temp.sessionNumber}, 'S2')).first3_objB;
                s3_first3_objA = opa_rat_temp(strcmp({opa_rat_temp.sessionNumber}, 'S3')).first3_objA;
                s3_first3_objB = opa_rat_temp(strcmp({opa_rat_temp.sessionNumber}, 'S3')).first3_objB;
    
                % Discrimination index of S2 vs S1
                opa_sorted(g).rat(r).session(f+length(novel)).s2vs1_objA_discrim = s2_first3_objA / (s2_first3_objA + s1_first3_objA);
                opa_sorted(g).rat(r).session(f+length(novel)).s2vs1_objB_discrim = s2_first3_objB / (s2_first3_objB + s1_first3_objB);

                % Discrimination index of S2 vs S3
                opa_sorted(g).rat(r).session(f+length(novel)).s2vs3_objA_discrim = s2_first3_objA / (s2_first3_objA + s3_first3_objA);
                opa_sorted(g).rat(r).session(f+length(novel)).s2vs3_objB_discrim = s2_first3_objB / (s2_first3_objB + s3_first3_objB);

                % Discrimination index of S1 vs S3
                opa_sorted(g).rat(r).session(f+length(novel)).s2vs3_objA_discrim = s1_first3_objA / (s1_first3_objA + s3_first3_objA);
                opa_sorted(g).rat(r).session(f+length(novel)).s2vs3_objB_discrim = s1_first3_objB / (s1_first3_objB + s3_first3_objB);
                
                % Discrimination index of S2 vs S1+S3/2 (Average time)
                opa_sorted(g).rat(r).session(f+length(novel)).s2vs1s3_objA = s2_first3_objA / (0.5*(s1_first3_objA+s3_first3_objA) + s2_first3_objA);
                opa_sorted(g).rat(r).session(f+length(novel)).s2vs1s3_objB = s2_first3_objB / (0.5*(s1_first3_objB+s3_first3_objB) + s2_first3_objB);

                % Discrimination index of S2 vs S1+S3/2 (Average DI)
                di_s1_objA = s2_first3_objA / (s2_first3_objA + s1_first3_objA);
                di_s3_objA = s2_first3_objA / (s2_first3_objA + s3_first3_objA);

                di_s1_objB = s2_first3_objA / (s2_first3_objA + s1_first3_objA);
                di_s3_objB = s2_first3_objB / (s2_first3_objB + s3_first3_objB);

                opa_sorted(g).rat(r).session(f+length(novel)).s2vs1s3_objA_avgDI = 0.5*(di_s1_objA + di_s3_objA);
                opa_sorted(g).rat(r).session(f+length(novel)).s2vs1s3_objB_avgDI = 0.5*(di_s1_objB + di_s3_objB);
            end
        end

        % Analyze the familiar compared to each practice sessions
        for p = 1:length(opa_rat_pract)
            temp_session_prac = opa_rat_pract(p);
            for f = 1:length(opa_rat_famil)
                opa_rat_temp = opa_rat_famil(f);

                objA_coords = [opa_rat_temp.session.objA_x_coord opa_rat_temp.session.objA_y_coord];
                objB_coords = [opa_rat_temp.session.objB_x_coord opa_rat_temp.session.objB_y_coord];
            
                distance_objA = sqrt((temp_session_prac.session.nose_x_coord-objA_coords(1)).^2 + (temp_session_prac.session.nose_y_coord-objA_coords(2)).^2)./parameters.est_conversion;
                distance_objB = sqrt((temp_session_prac.session.nose_x_coord-objB_coords(1)).^2 + (temp_session_prac.session.nose_y_coord-objB_coords(2)).^2)./parameters.est_conversion;
            
                prac_time_objA = sum(distance_objA(1:firstXframes) < parameters.cm_from_object)/parameters.fps;
                prac_time_objB = sum(distance_objB(1:firstXframes) < parameters.cm_from_object)/parameters.fps;
            
                sess_time_objA = opa_rat_temp.first3_objA;
                sess_time_objB = opa_rat_temp.first3_objB;    

                opa_sorted(g).rat(r).practice(p).familiar(f).sessionType = opa_rat_temp.sessionType;
                opa_sorted(g).rat(r).practice(p).familiar(f).sessionNumber = opa_rat_temp.sessionNumber;
                opa_sorted(g).rat(r).practice(p).familiar(f).pracNumber = temp_session_prac.sessionType;
                opa_sorted(g).rat(r).practice(p).familiar(f).objA_discrim = sess_time_objA / (prac_time_objA+sess_time_objA);
                opa_sorted(g).rat(r).practice(p).familiar(f).objB_discrim = sess_time_objB / (prac_time_objB+sess_time_objB);
            end
        end
    end
end