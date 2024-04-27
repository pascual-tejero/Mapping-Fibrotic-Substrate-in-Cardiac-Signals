function [matrix_AT,matrix_RT] = calculation_AT_RT(matrix,cut_simulation,a)

    sampleFrequency = 333;
    mode = 1;
    k = 1;
    lenghtOfWin = 1;

    % Obtention of the time window
    accpt = zeros(size(matrix,1),size(matrix,2));
    
    % matrix_eval = matrix(:,200:end);
    matrix_eval = matrix;
    
    % Location of the last time window in order to calculate the AT
    y_lins = linspace(cut_simulation(2),cut_simulation(1),size(matrix,2));
    
    for i = 1:size(matrix,1)
    %     [stepFuncFinal(:,i),segments{:,i}] = getActiveSegmentsFromNLEOuni(signalOut(:,i),sampleFrequency,lenghtOfWin,k);
        [accp,timelngth] = Hactseg(matrix_eval(i,:),sampleFrequency,a);
        accpt(i,1+(size(matrix,2)-length(accp)):end) = accp;

    %% Obtaining the AT of the last window
             
    %Calculation of the last time window
    pos_last1_accpt_lw = find(accpt(i,:),1,'last');

    if isempty(pos_last1_accpt_lw)
        matrix_AT(i) = NaN;
        matrix_RT(i) = NaN;
        continue
    end

    vector_values_lw = accpt(i,1:pos_last1_accpt_lw);
    pos_0_accpt_lw = find(vector_values_lw == 0);       

    pos_window_lw(i,1) = pos_0_accpt_lw(end); pos_window_lw(i,1) = pos_window_lw(i,1) + 1;
    pos_window_lw(i,2) = pos_last1_accpt_lw;

    window_B_lw = matrix(i,pos_window_lw(i,1):pos_window_lw(i,2));

    pos_window_lw_s(i,1) = pos_window_lw(i,1)/size(matrix,2)*2000;
    pos_window_lw_s(i,2) = pos_window_lw(i,2)/size(matrix,2)*2000;



    %Calculation RT
    vector_values_RT = accpt(i,1:pos_window_lw(i,1)); vector_values_RT(end) = 0;
    pos_last1_accpt_RT = find(vector_values_RT,1,'last');

    if isempty(pos_last1_accpt_RT)
        matrix_AT(i) = NaN;
        matrix_RT(i) = NaN;
        continue
    end

%         disp(i) %549
    pos_window_RT(i,1) = pos_last1_accpt_RT; pos_window_RT(i,1) = pos_window_RT(i,1) + 1;
    pos_window_RT(i,2) = pos_window_lw(i,1); pos_window_RT(i,2) = pos_window_RT(i,2) - 1;

    window_B_RT = matrix(i,pos_window_RT(i,1):pos_window_RT(i,2));

    pos_window_RT_s(i,1) = pos_window_RT(i,1)/size(matrix,2)*2000;
    pos_window_RT_s(i,2) = pos_window_RT(i,2)/size(matrix,2)*2000;


    %Calculation AT
    pos_last1_accpt_AT = pos_last1_accpt_RT;
    vector_values_AT = accpt(i,1:pos_last1_accpt_AT); 
    vector_values_find0_AT = find(vector_values_AT == 0); 
    pos_first1_accpt_AT = vector_values_find0_AT(end);

    pos_window_AT(i,1) = pos_first1_accpt_AT; pos_window_AT(i,1) = pos_window_AT(i,1) + 1;
    pos_window_AT(i,2) = pos_last1_accpt_AT; 

    window_B_AT = matrix(i,pos_window_AT(i,1):pos_window_AT(i,2));

    pos_window_RT_s(i,1) = pos_window_AT(i,1)/size(matrix,2)*2000;
    pos_window_RT_s(i,2) = pos_window_AT(i,2)/size(matrix,2)*2000;



%         value_window_lw_s = y_lins(pos_window_lw(i,1):pos_window_lw(i,2));
    value_window_RT_s = y_lins(pos_window_RT(i,1):pos_window_RT(i,2));
    value_window_AT_s = y_lins(pos_window_AT(i,1):pos_window_AT(i,2));

%         window_s_lw(i,1:length(value_window_lw_s)) = value_window_lw_s;
    window_s_RT(i,1:length(value_window_RT_s)) = value_window_RT_s;
    window_s_AT(i,1:length(value_window_AT_s)) = value_window_AT_s;

% disp(i);
%         Plot
%         figure; plot(y_lins,matrix(i,:)); hold on; 
%         plot(y_lins,accpt(i,:)); hold on; 
%         plot(value_window_AT_s, window_B_AT); hold on;
%         plot(value_window_RT_s,window_B_RT);
%         
%         legend("Electrogram","Detection function","Window AT","Window RT");
%         xlabel("Time (ms)"); ylabel("Amplitude voltage (mV)");
%         title(['Fibrotic tissue - Reentry case - Electrode: ', num2str(i)]); 
%         title(['Fibrotic tissue - Reentry case - Source point: ', num2str(i)]); 
% %         
    %Calculation of AT
    derivative_AT = gradient(window_B_AT);
    [~,at(i)] = max((derivative_AT),[],2);
    at(i) = at(i) + pos_window_AT(i,1);

    matrix_AT(i) = y_lins(at(i));

    %Calculation of RT
    if length(window_B_RT) <= 2
        pos_RT(i) = pos_window_RT(i,1);
        matrix_RT(i) = y_lins(floor(pos_RT(i)));
        continue;
    end
    
    [accp_RT,timelngth] = Hactseg(window_B_RT,6,a);
    [pks_accp_RT,loc_accp_RT] = findpeaks(accp_RT);

%     figure; plot(window_B_RT); hold on; plot(accp_RT);

    if (sum(accp_RT) == 0) || (sum(accp_RT)/length(accp_RT) == 1)
        value_accpRT_first1 = 0;
        value_accpRT_first0(1) = 0;
        accp_RT_first0 = accp_RT;
        
    elseif length(pks_accp_RT) == 1
        value_accpRT_first1 = find(accp_RT,1);
        accp_RT_first1 = accp_RT(value_accpRT_first1:end);
        
        value_accpRT_first0 = 0;
        accp_RT_first0 = accp_RT_first1;  
               
    else
        value_accpRT_first1 = find(accp_RT,1);
        accp_RT_first1 = accp_RT(value_accpRT_first1:end);
    
        value_accpRT_first0 = find(accp_RT_first1 == 0);
        accp_RT_first0 = accp_RT_first1(value_accpRT_first0(1):end);     
    end
    
%     window_B_RT_new = window_B_RT(value_accpRT_first0(1):end-15);
    
%     value_accpRT_first1_2 = find(accp_RT_first0,1);
%     accp_RT_first1_2 = accp_RT_first0(value_accpRT_first1_2:end);
%     
%     accp_RT_new = accp_RT_first1_2(1:end-15);
    
    if (sum(accp_RT_first0) == 0) || (sum(accp_RT_first0)/length(accp_RT_first0) == 1)
        window_B_RT_new = window_B_RT;
    elseif (value_accpRT_first1+value_accpRT_first0(1) + 20) >= length(window_B_RT)
        window_B_RT_new = window_B_RT((value_accpRT_first1+value_accpRT_first0(1)+1):end);
    else
        window_B_RT_new = window_B_RT((value_accpRT_first1+value_accpRT_first0(1)):end-20);
    end
%     figure; plot(window_B_RT);

    if length(window_B_RT_new) <= 2
        [~,pos_RT(i)] = max(window_B_RT_new);
        pos_RT(i) = pos_RT(i) + value_accpRT_first1 + value_accpRT_first0(1) + pos_window_RT(i,1);
        continue;
    end

        [pks,locs] = findpeaks(abs(window_B_RT_new),'MinPeakProminence',0.01);
        %     figure; plot(abs(window_B_RT_new));
        
   
    

    if length(pks) == 1
        pos_RT(i) = locs + value_accpRT_first1 + value_accpRT_first0(1) + pos_window_RT(i,1);
    elseif length(pks) > 1
        pos_RT(i) = mean(locs);
        pos_RT(i) = pos_RT(i) + value_accpRT_first1 + value_accpRT_first0(1) + pos_window_RT(i,1);
    else
        [value_RT, pos_RT(i)] = max(window_B_RT_new); %ver 23 en fuentes
        
        if (value_RT == window_B_RT_new(end)) 
            [~,pos_RT(i)] = max(gradient(window_B_RT_new));
            pos_RT(i) = pos_RT(i) + value_accpRT_first1 + value_accpRT_first0(1) + pos_window_RT(i,1);
        else
            pos_RT(i) = pos_RT(i) + value_accpRT_first1 + value_accpRT_first0(1) + pos_window_RT(i,1);
        end

    end
%
% 
% disp(i)
% 

    matrix_RT(i) = y_lins(floor(pos_RT(i)));
%         
%         figure;
%         plot(window_B_RT); hold on;
%         plot(derivative_RT); hold on;
%         plot(derivative_2_RT);
%         
%         legend('Signal','First derivative','Second derivative');
%         title("Sources - Control case");

%         for j = 1:length(derivative_2_RT)
%             if (j < length(derivative_2_RT))
%                 if((derivative_2_RT(j) > 0) & (derivative_2_RT(j+1) < 0))
%                     if (derivative_2_RT(j) > derivative_2_RT(j+1))
%                         pos_RT(i) = j+1+pos_window_RT(i,1);
%                         break;
%                     else 
%                         pos_RT(i) = j+pos_window_RT(i,1);
%                         break;
%                     end
%                 end
%             else
%                 if(derivative_2_RT > 0)
%                     [~,pos_RT(i)] = min(derivative_2_RT);
%                     pos_RT(i) = pos_RT(i)+pos_window_RT(i,1);
%                     break;
%                 else
%                     [~,pos_RT(i)] = max(derivative_2_RT);
%                     pos_RT(i) = pos_RT(i) + pos_window_RT(i,1);
%                     break;
%                 end
%             end
%         end







    end
end

    % [~,AT_pos(1,1)] = max(at); [~,AT_pos(1,2)] = min(at);
    % value = window_s(at(i));
    % 
    % AT_matrix_final(1,1) = window_s(AT_pos(1,1),at(AT_pos(1,1)));
    % AT_matrix_final(1,2) = window_s(AT_pos(1,2),at(AT_pos(1,2)));
%   find((value_at < 1850) & (value_at > 1830))
%     for i=1:length(aa)
%         hold on; 
%         plot(matrix(aa(i),:))
%     end



