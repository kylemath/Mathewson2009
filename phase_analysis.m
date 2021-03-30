% function [circ_grand_mean F] = fft_analysis_phase(channel)

clear all
close all

% If you need to load subjects

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
subnums = ['424' '431' '334' '423' '409' '433' '434' '435' '067' '238' '437']; % 
numsubs = length(subnums)/3;

for x = 1:3:(length(subnums)-2)

    subnum = subnums(1,x:(x+2));
    path =  ['C:\DATA\\omm\\eeg\\' subnum 'EEGLAB\\'];

    EEG = pop_loadset( 'filename', [subnum 'regular miss.set'], 'filepath', path);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = pop_loadset( 'filename', [subnum 'regular hit.set'], 'filepath', path);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

end

time_period_low = 131;   %200 ms pre is 131 to 170; 400ms pre is 91 to 170
time_period_high = 170;
n_chan = 23;


Fs = 200;                    % Sampling frequency
T = 1/Fs;                    % Time
L = (time_period_high-time_period_low+1);                     % Trial Length
f = Fs/2*linspace(0,1,L/2);   %divide the frequency into bins from 0 to half the sampling frequency
low_f = 3;


for i_sub = 1:11                                                  %for each subject
    for i_cond = 1:2                                      %for each condition
        for i_chan = 17:17
            for i_trial=1:size(ALLEEG(1,(2*(i_sub-1))+ i_cond).data,3)           %for each single trial

                data = ALLEEG(1,(2*(i_sub-1))+ i_cond).data(i_chan,time_period_low:time_period_high,i_trial);   %get the trial data

                y = fft(data);
                m = abs(y)/(L/2);             %get the magnitude of each fbin
                pow = m(1:L/2).^2;            %square to get the power, only of the first half, second is redundant

                
                p = angle(y(1:L/2));          %get the phase in radians of each bin, only first half
                p_d = p * (180/pi);           %Convert phase to degrees
             
                running_power(i_trial,((i_sub-1)*2) + i_cond,i_chan ) = log(pow(1,low_f));  %keep a list for each trial of the power in the bin
                
                if p_d(1,low_f) >=0
                    running_phase(i_trial,((i_sub-1)*2) + i_cond,i_chan ) = p_d(1,low_f);  %and again for the phase
                else
                    running_phase(i_trial,((i_sub-1)*2) + i_cond,i_chan ) = 360 - abs((p_d(1,low_f)));
                end

            end
       end
    end
end


for i_chan = 17:17
    for i_sub = 1:11
        for i_cond = 1:2
            [circ_mean,range,X,Y,cos_a,sin_a] = circle_mean(running_phase(:,((i_sub-1)*2)+ i_cond,i_chan));    % Computes the circular mean for each condition in each subject
            phase_out(i_sub,i_cond,i_chan) = circ_mean;                        %record the condition circular phase mean
            range_out(i_sub,i_cond,i_chan) = range;                            %record the range (concentration) of each of these means 
            X_out(i_sub,i_cond,i_chan) = X;                                    %record the cosine component of the mean
            Y_out(i_sub,i_cond,i_chan) = Y;                                    %record the sine component of the mean
            cos_out(i_sub,i_cond,i_chan) = cos_a;      
            sin_out(i_sub,i_cond,i_chan) = sin_a;
           
            
            pow_sum = 0;
            for i_pow_row = 1:size(running_power(:,((i_sub-1)*2) + i_cond,i_chan))
                if running_power(i_pow_row,((i_sub-1)*2) + i_cond,i_chan) == 0
                    break
                end
                if running_power(i_pow_row,((i_sub-1)*2) + i_cond,i_chan) ~= 0
                    pow_sum = pow_sum + running_power(i_pow_row,((i_sub-1)*2) + i_cond,i_chan);
                end
            end
            power_out(i_sub,i_cond,i_chan) = pow_sum/(i_pow_row-1);
            
        end
    end
end

for i_chan = 17:17
    [circ_grand_mean(i_chan,1), circ_grand_range(i_chan,1), X_bar(i_chan,1), Y_bar(i_chan,1), cos_bar(i_chan,1), sin_bar(i_chan,1)] = circle_grand_mean(phase_out(:,1,i_chan),range_out(:,1,i_chan));
    [circ_grand_mean(i_chan,2), circ_grand_range(i_chan,2), X_bar(i_chan,2), Y_bar(i_chan,2), cos_bar(i_chan,2), sin_bar(i_chan,2)] = circle_grand_mean(phase_out(:,2,i_chan),range_out(:,2,i_chan));
    circ_grand_mean(i_chan,3) = circ_grand_mean(i_chan,1) - circ_grand_mean(i_chan,2);

    if abs(circ_grand_mean(i_chan,3)) > 180
        if circ_grand_mean(i_chan,3) < 0
            circ_grand_mean(i_chan,3) = -360 - circ_grand_mean(i_chan,3);
        else
            circ_grand_mean(i_chan,3) = 360 - circ_grand_mean(i_chan,3);
        end
    end
    
    
    [circ_grand_mean(i_chan,4)] = circle_test(cos_out(:,:,i_chan),sin_out(:,:,i_chan),cos_bar(i_chan,:),sin_bar(i_chan,:));
    
    pow_grand_mean(i_chan,1) = mean(power_out(:,1,i_chan));
    pow_grand_mean(i_chan,2) = mean(power_out(:,2,i_chan));
    pow_grand_mean(i_chan,3) = pow_grand_mean(i_chan,1) - pow_grand_mean(i_chan,2);    
    
end


