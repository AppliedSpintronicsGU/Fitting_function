clear all
close all
clc

%% Procedure and output format
% 
% This code allows to find as many Lorentzians as the user sets and also
% finds low power and closely positioned peaks.
% 
% The fitted data is then organized in groups so that it will be easier to
% sort them out. Groups with only one peak are regarded as wrong fits and
% are discarded.
% 
% The output format is given mainly in the export_vector which form is:
% 
% | Angle | Current | Resistance | Peak_freq_G1 | Peak_power_G1 | Int_power_G1 | Linewidth_G1 | Peak_freq_G2 | etc
% 
%  Each column is a vector of length specified by the current.

%% How to
% 
% The variables below are the control parameters of the fitting.
% 
% When correction is set to zero, the code works as the previous version 
% while setting it to 1 would try to fit any small peaks that can find 
% (even close peaks). Setting correction to 2 is useful for noisy data
% where noise is fitted as peaks.
% 
% Max_peak is the maximum number of peaks that the fitting will try to
% find. This is to prevent irrelevant features to be fitted. Set it to the
% number of important peaks you expect to find or to a reasonable ammount
% to speed up the fitting.
% 
% The threshold only appears in the peak_finder_auto. By default is 2.
% Lowering might help you to find low power peaks but if too many are found
% the fitting might be wrong.
% 
% N_attemp is the number of times the correction algorithm will attempt to
% find low power peaks. The correction is basd on randomly assigning peak
% positions so the lower the peak, the greater the number of attempts
% needed to find it. The default value 5 is good for reasonably strong
% peaks.
% 
% Gain and refl_corr are experimental corrections for the amplifier gain
% and the reflections due to impedance mimatch
% 
% dF is a frequency step using to organize the data. It should correspond
% to the maximum step observed for a single "branch".
% 
% REGROUP is a flag. When set to 1 triggers a section of the code that
% tries to sort the fitted data based on frequency similarities. Due to
% this, if there are huge frequency steps, the grouping will sort them as
% two different peaks. Setting it to 2 will trigger another section of the
% code that tries to readjust the previous grouping. Both results can
% improve or neglect data so it is recommended to check the "_org2"
% variables to track down the changes.
% 
% LowFreq sets the minimum frequency (in GHz) to cut out from the measured
% spectra. Neglecting low frequency peaks leads to a good definition of the
% "0 power" and correct fitting of Lorentzian functions.
% 
% NoiseFreq is the window (in GHz) where the code looks for a noise
% reference. This means that, from the first spectra, the frequencies
% ranging from LowFreq to NoiseFreq are taken to calculate the noise level
% (mean value) and minimum power (standard deviation*3). So one must be
% sure that there is only noise in the given window.
% 
% MaxLinewidth is the maximum linewidth that the fitting routine will
% attempt to fit. Beyond this value, the linewith will be saturated.


%% Important Variables!
correction = 2;          % Allow power correction, 0 none, 1 low power,
                         % 2 high power (good for noisy data)
Max_peak = 2;            % Maximum peaks to fit
threshold = 2;           % Threshold for peak detection. NOTE this is not
                         % a threshold above which maxima are found, but
                         % a parameter for the peak_det function. A value
                         % of 2 is suggested.
N_attempt = 1;           % Allow several attempts to find the best fit. For
                         % low power peaks the default is 5, while for data
                         % with well defined peaks setting it to 1 is good
                         % enough.
gain = 35;               % Gain of the amplifier used during measurements
refl_corr = 5;           % Correction for the reflections due to impedance mismatch
dF = 150e6;              % Maximum frequency step in Hz for grouping. Default 50e6.
REGROUP = 1;             % Set to 1 to sort based on frequency similarity
                         % Set to 2 to readjust groups
                         % Roughly set 1 when the data is "good" and there
                         % are few peaks, set to 2 otherwise. To be sure
                         % that the grouping was succesfull, compare with
                         % the "_org2" files.
LowFreq = 5;             % Minimum frequency to start scan in GHz
WindowFreq = [5 30];     % Window to estimate the measurement noise in GHz.
MaxLinewidth = 1;        % Maximum linewidth before saturation in GHz.

% String where the text file will be saved.
Folder = 'C:\Users\Dumas\Desktop';
                         

%% Code
[log_name, pathname] = uigetfile({'*.log','Log Files (*.log)';...
    '*.*','All Files (*.*)'},'Choose a file...','Multiselect','off');

% If no file is selected, stops the program.
if pathname==0
    break;

% Otherwise, runs the main program
else
    % Load the .log file in a matrix, skipping the header. The location of the
    % file containing the spectra is the last string column.
    log_path=strcat(pathname,log_name);
    [current_log angle_log spectrum_path]=textread(log_path,...
        '%f %f %s',-1,'headerlines',7,'delimiter','\t'); % 7 is standard
end

%%
% Count time
tic;
% Reads sequentially the rows of the .log file, and points each time to a
% new file with the spectrum.
steps = length(angle_log);
current = zeros(steps,1);
angle = zeros(steps,1);
resistance = zeros(steps,1);
iteration = length(angle_log);

for i = 1 : length(angle_log)
    angle_log(i);
    i
    clear fid
    fid=spectrum_path{i};
    fid = regexp(fid, '\\', 'split');
    fid = fid(end);
    fid=char(strcat(pathname,fid));
    
    % Reads the measurement setting from the header of each spectrum.
    [temp1 current(i)]=textread(fid,'%s %f',1,'headerlines',5);
    [temp1 angle(i)]=textread(fid,'%s %f',1,'headerlines',6);
    [temp1 resistance(i)]=textread(fid,'%s %f',1,'headerlines',8);
    [temp1 temp2 vbw]=textread(fid,'%s %s %f',1,'headerlines',9);
    [temp1 temp2 rbw]=textread(fid,'%s %s %f',1,'headerlines',10);
    [temp1 temp2 points]=textread(fid,'%s %s %f',1,'headerlines',11);

    % Recalculate the actual resistance of the STO. Uncomment next line
    % if the correction is needed.
%     resistance = resistance - 6.7; % 6.7 Ohm is the measured lead
                                   % resistance including bias tee and
                                   % circulator.
    
    % Create a matrix with the spectrum, whose length is the number of data
    % points acquired by the spectrum analyzer.
    [frequency_pre spectrum_pre clean_spectrum_pre]=textread(fid,'%f %f %f',points,...
        'headerlines',18,'delimiter','\t');
    
    % Identifies singular points created by Labview communication errors,
    % and set the value of the clean spectrum to be 0 at the single point.
    error_ind=find(abs(clean_spectrum_pre)>40);
    if isempty(error_ind)==0
        clean_spectrum_pre(error_ind)=0;
    end
    
    error2_ind=find(clean_spectrum_pre<-2);
    if isempty(error2_ind)==0
        clean_spectrum_pre(error2_ind)=0;
    end
    
    error3_ind=find(isnan(clean_spectrum_pre)==1);
    if isempty(error3_ind)==0
        clean_spectrum_pre(error3_ind)=0;
    end
    
    % Deteles the low frequencies of the spectrum
    if i == 1
        low_point = find(frequency_pre <= LowFreq*1e9);
    end
    clean_spectrum = clean_spectrum_pre(low_point(end):end);
    spectrum = spectrum_pre(low_point(end):end);
    frequency = frequency_pre(low_point(end):end);
    
    % Remove NaN values from the raw data. Uses interpolation to create a data point.
    clean_spectrum=naninterp(clean_spectrum);
    spectrum=naninterp(spectrum);
    
    % Normalize spectrum and calculate noise floor.
    if i == 1
        noise_dBW = max(spectrum) - max(clean_spectrum) - 30; % from dBm to dBW'
        clean_spectrum_dBw = clean_spectrum + noise_dBW - gain;
        real_dummy = (1/rbw)*refl_corr*10.^(clean_spectrum_dBw/10);
        noise_select_down = find(frequency <= WindowFreq(1) * 1E9);
        noise_select_up = find(frequency >= WindowFreq(2) * 1E9);
        noise_level = mean(real_dummy(noise_select_down(end):noise_select_up(1)));
        min_power = 3 * std(real_dummy - noise_level);
    else
        clean_spectrum_dBw = clean_spectrum + noise_dBW - gain;
    end

    %% Fitting rountine
    % Check the best of N attempts
    s = 1;
    
    % Allocate dummy variables
    for ii = 1:N_attempt
        peak_power_dummy{ii} = {NaN};
        int_power_dummy{ii} = {NaN};
        linewidth_dummy{ii} = {NaN};
        peak_freq_dummy{ii} = {NaN};
        PSDreal_dummy{ii} = {NaN};
        yf{ii} = {NaN};
    end
    error = zeros(1,N_attempt);
    
    % Perform N_attempt fittings
    for kk = 1 : N_attempt
        [peak_power_dummy{kk},int_power_dummy{kk},linewidth_dummy{kk},...
            peak_freq_dummy{kk},PSDreal_dummy{kk},error(kk),yf{kk}]=...
            peak_fit_routine(frequency,spectrum,clean_spectrum_dBw,...
            gain,rbw,refl_corr,s,correction,Max_peak,threshold,...
            noise_level,min_power,MaxLinewidth);
    end
    
    % Finds the fit with minimum error
    [dummy,select] = min(error);
    if dummy~=Inf
        % Stores the best fit information
        peak_power(i,1:length(peak_power_dummy{select})) = peak_power_dummy{select};
        int_power(i,1:length(peak_power_dummy{select})) = int_power_dummy{select};
        linewidth(i,1:length(peak_power_dummy{select})) = linewidth_dummy{select};
        peak_freq(i,1:length(peak_power_dummy{select})) = peak_freq_dummy{select};
        PSD_real(i,1:length(PSDreal_dummy{select})) = PSDreal_dummy{select};
        PSD_fit(i,1:length(yf{select})) = yf{select};
    else
        % Stores zeros
        peak_power(i,:) = 0;
        int_power(i,:) = 0;
        linewidth(i,:) = 0;
        peak_freq(i,:) = 0;
        PSD_real(i,:) = 0;
        PSD_fit(i,:) = 0;
    end
    f_I_matrix(i,:) = clean_spectrum;
end
% Removes zeroes for NaN
[long,wide] = size(peak_freq);
for ii=1:long
    for jj=1:wide
        if peak_freq(ii,jj) == 0
            int_power(ii,jj) = NaN;
            linewidth(ii,jj) = NaN;
            peak_power(ii,jj) = NaN;
            peak_freq(ii,jj) = NaN;
        end
    end
end

%% Peak grouping
% Define variables
[loop,Max_peaks] = size(peak_freq);
group_no = 0;
output{loop+1,group_no+1} = {0};
flag = 0;
l = 0;
number = 1;
flag_break = 0;

% Build a cell matrix with the indices of the sorted peak frequencies
for i = 1 : loop
    for j = 1 : Max_peaks
        % Check for valid peak frequencies
        if (peak_freq(i,j) ~= NaN) && (peak_freq(i,j) > frequency(1)) && (peak_freq(i,j) < frequency(end))
            % If it is the first, save without comparing
            if i == 1
                group_no = group_no + 1;
                output{i,group_no} = [i,j];
            % If the first has been already stored, check for similarities
            else
                % Check if the peak is similar to previously existing peaks
                for k = 1: group_no
                    % Search for last stored peak of every group
                    while isempty(output{i-l,k}) == 1
                        l = l + 1;
                        % Break the while if l = 1
                        if (l == i) || (l == 10)
                            flag_break = 1;
                            break;
                        end
                    end
                    % Compare
                    if flag_break == 0
                        if abs(peak_freq(i,j) - peak_freq(output{i-l,k}(1),output{i-l,k}(2))) < dF
                            flag = 1;
                            similar(number) = k;
                            difference(number) = abs(peak_freq(i,j) - peak_freq(output{i-l,k}(1),output{i-l,k}(2)));
                            number = number + 1;
                        end
                    end
                    l = 1;
                    flag_break = 0;
                end
                number = 1;
                % Choose closer peak and put it in the similar group
                if flag == 1
                    [dummy,select] = find(difference - min(difference) == 0);
                    output{i,similar(select(1))} = [i,j];
                    flag = 0;
                    clear similar
                    clear difference
                % If not, create a separate group
                elseif flag == 0
                    % Check also if the peak is too close to the current peaks (10MHz)
                    for k = 1 : group_no
                        if isempty(output{i,k}) ~= 1
                            if abs(peak_freq(i,j) - peak_freq(output{i,k}(1),output{i,k}(2))) < 10e6
                                flag = 2;
                            end
                        end
                    end
                    if flag == 0
                        group_no = group_no + 1;
                        output{i,group_no} = [i,j];
                    end
                end
            end
        end
    end
end

% Create organized data
peak_freq_org(loop,group_no) = NaN;
peak_power_org(loop,group_no) = NaN;
int_power_org(loop,group_no) = NaN;
linewidth_org(loop,group_no) = NaN;
for i = 1 : loop
    for j = 1 : group_no
        if isempty(output{i,j}) ~= 1
            if output{i,j}(1) ~= 0
                peak_freq_org(i,j) = peak_freq(output{i,j}(1),output{i,j}(2));
                peak_power_org(i,j) = peak_power(output{i,j}(1),output{i,j}(2));
                int_power_org(i,j) = int_power(output{i,j}(1),output{i,j}(2));
                linewidth_org(i,j) = linewidth(output{i,j}(1),output{i,j}(2));
            end
        end
    end
end
clear output

% Remove single points
j = 1;
for i = 1 : group_no
    [dummy]=find(peak_freq_org(:,i)~=0);
    if length(dummy) > 5
        peak_freq_org2(:,j) = peak_freq_org(:,i);
        peak_power_org2(:,j) = peak_power_org(:,i);
        int_power_org2(:,j) = int_power_org(:,i);
        linewidth_org2(:,j) = linewidth_org(:,i);
        j = j + 1;
    end
end

%% Regroup script
switch REGROUP
    case 1
        % Reduces the number of groups based on the desired Max_peak
        regroup = Max_peaks;
        peak_freq_org3(:,1:Max_peaks) = peak_freq_org2(:,1:Max_peaks);
        peak_power_org3(:,1:Max_peaks) = peak_power_org2(:,1:Max_peaks);
        int_power_org3(:,1:Max_peaks) = int_power_org2(:,1:Max_peaks);
        linewidth_org3(:,1:Max_peaks) = linewidth_org2(:,1:Max_peaks);
        for i = Max_peaks + 1 : j - 1
            % Find the first value of the extra group
            found = find(peak_freq_org2(:,i)>0);

            % Chooses if can be added to another group
            found_pre = find(peak_freq_org2(found(1),1:Max_peaks) == 0);
            if ~isempty(found_pre)
                if length(found_pre) == 1
                    peak_freq_org3(:,found_pre) = peak_freq_org3(:,found_pre) + peak_freq_org2(:,i);
                    peak_power_org3(:,found_pre) = peak_power_org3(:,found_pre) + peak_power_org2(:,i);
                    int_power_org3(:,found_pre) = int_power_org3(:,found_pre) + int_power_org2(:,i);
                    linewidth_org3(:,found_pre) = linewidth_org3(:,found_pre) + linewidth_org2(:,i);
                else
                    freq_comp = abs(peak_freq_org2(found(1),i)*ones(1,length(found_pre)) - peak_freq_org2(found(1),found_pre));
                    select = find((freq_comp) == min(freq_comp));
                    peak_freq_org3(:,select(1)) = peak_freq_org3(:,select(1)) + peak_freq_org2(:,i);
                    peak_power_org3(:,select(1)) = peak_power_org3(:,select(1)) + peak_power_org2(:,i);
                    int_power_org3(:,select(1)) = int_power_org3(:,select(1)) + int_power_org2(:,i);
                    linewidth_org3(:,select(1)) = linewidth_org3(:,select(1)) + linewidth_org2(:,i);
                end
            else
            regroup = regroup + 1;
            peak_freq_org3(:,regroup) = peak_freq_org2(:,i);
            peak_power_org3(:,regroup) = peak_power_org2(:,i);
            int_power_org3(:,regroup) = int_power_org2(:,i);
            linewidth_org3(:,regroup) = linewidth_org2(:,i);
            end
        end
        for i = 1 : loop
            for k = 1 : regroup
                if peak_freq_org3(i,k) == 0
                    peak_freq_org3(i,k) = NaN;
                    peak_power_org3(i,k) = NaN;
                    int_power_org3(i,k) = NaN;
                    linewidth_org3(i,k) = NaN;
                end
            end
        end
        
    case 2
        % Re-sort the data in the minimum ammount of peaks based on dF as a
        % comparison parameter

        % Allocate
        regroup = 1;
        log_regroup = 1;
        peak_freq_org3(:,regroup) = zeros(iteration,1);
        peak_power_org3(:,regroup) = zeros(iteration,1);
        int_power_org3(:,regroup) = zeros(iteration,1);
        linewidth_org3(:,regroup) = zeros(iteration,1);

        % Expand data in chunks
        count = 0;
        chunk = 1;
        for i = 1 : j - 1
            for k = 1 : iteration
                if peak_freq_org2(k,i) == 0
                    if count >= 5
                        output{chunk} = [i,k-count:k-1];
                        chunk = chunk + 1;
                        count = 0;
                    else
                        count = 0;
                    end
                else
                    count = count + 1;
                end
            end
        end

        % Assign new groups and re-group if a similar group is found

        % First group
        peak_freq_org3(output{1}(2:end),regroup) = peak_freq_org2(output{1}(2:end),output{1}(1));
        peak_power_org3(output{1}(2:end),regroup) = peak_power_org2(output{1}(2:end),output{1}(1));
        int_power_org3(output{1}(2:end),regroup) = int_power_org2(output{1}(2:end),output{1}(1));
        linewidth_org3(output{1}(2:end),regroup) = linewidth_org2(output{1}(2:end),output{1}(1));

        % Assign and compare following groups
        for i = 2 : chunk - 1
            % Allocate
            check = 0;
            dummy_index = 1;
            clear similar
            clear select
            clear select2

            % Check if the current chunk could fit another group
            for k = 1 : regroup
                % Find the last chunk of each group
                clear select_chunk
                select_chunk = find(log_regroup == k);

                % Save similar group index
                if (output{i}(2) > output{select_chunk(end)}(end)) || (output{i}(end) < output{select_chunk(end)}(2))
                    check(dummy_index,1) = select_chunk(end);

                    % Save flag of comparison type
                    if output{i}(2) > output{select_chunk(end)}(end)
                        check(dummy_index,2) = 1;
                    else
                        check(dummy_index,2) = 2;
                    end

                    dummy_index = dummy_index + 1;
                end
            end

            % If one or more fitting groups are found, check if the frequency
            % difference is small enough
            if check ~= 0
                for k = 1 : dummy_index - 1
                    % Check difference of each possible group
                    if check(k,2) == 1
                        similar(k) = abs(peak_freq_org3(output{check(k,1)}(end),log_regroup(check(k,1))) - peak_freq_org2(output{i}(2),output{i}(1)));
                    else
                        similar(k) = abs(peak_freq_org3(output{check(k,1)}(2),log_regroup(check(k,1))) - peak_freq_org2(output{i}(end),output{i}(1)));
                    end
                end

                % Check if any group satisfy the minimum frequency condition
                select = find(similar <= dF);

                % If one or more do, save the one with minimum frequency and
                % save
                if isempty(select) == 0
                    select2 = find(similar == min(similar));

                    % Find in which "regroup" that chunk is
                    log_regroup(i) = log_regroup(check(select2(1),1));
                    peak_freq_org3(output{i}(2:end),log_regroup(i)) = peak_freq_org2(output{i}(2:end),output{i}(1));
                    peak_power_org3(output{i}(2:end),log_regroup(i)) = peak_power_org2(output{i}(2:end),output{i}(1));
                    int_power_org3(output{i}(2:end),log_regroup(i)) = int_power_org2(output{i}(2:end),output{i}(1));
                    linewidth_org3(output{i}(2:end),log_regroup(i)) = linewidth_org2(output{i}(2:end),output{i}(1));
                else
                    regroup = regroup + 1;
                    log_regroup(i) = regroup;
                    peak_freq_org3(output{i}(2:end),regroup) = peak_freq_org2(output{i}(2:end),output{i}(1));
                    peak_power_org3(output{i}(2:end),regroup) = peak_power_org2(output{i}(2:end),output{i}(1));
                    int_power_org3(output{i}(2:end),regroup) = int_power_org2(output{i}(2:end),output{i}(1));
                    linewidth_org3(output{i}(2:end),regroup) = linewidth_org2(output{i}(2:end),output{i}(1));
                end

            % Otherwise save in a new group
            else
                regroup = regroup + 1;
                log_regroup(i) = regroup;
                peak_freq_org3(output{i}(2:end),regroup) = peak_freq_org2(output{i}(2:end),output{i}(1));
                peak_power_org3(output{i}(2:end),regroup) = peak_power_org2(output{i}(2:end),output{i}(1));
                int_power_org3(output{i}(2:end),regroup) = int_power_org2(output{i}(2:end),output{i}(1));
                linewidth_org3(output{i}(2:end),regroup) = linewidth_org2(output{i}(2:end),output{i}(1));
            end
        end

        % Substitute and correct
        [dummy, group_test] = size(peak_freq_org3);
        if group_test > regroup
            for i = 1 : iteration
                for k = regroup + 1 : group_test
                    peak_freq_org3(i,k) = NaN;
                    peak_power_org3(i,k) = NaN;
                    int_power_org3(i,k) = NaN;
                    linewidth_org3(i,k) = NaN;
                end
            end
        end
        for i = 1 : iteration
            for k = 1 : regroup
                if peak_freq_org3(i,k) == 0
                    peak_freq_org3(i,k) = NaN;
                    peak_power_org3(i,k) = NaN;
                    int_power_org3(i,k) = NaN;
                    linewidth_org3(i,k) = NaN;
                end
            end
        end
   
    case 0
        % Substitute zeros for NaN
        for i = 1 : loop
            for k = 1 : j - 1
                if peak_freq_org(i,k) == 0
                    peak_freq_org(i,k) = NaN;
                    peak_power_org(i,k) = NaN;
                    int_power_org(i,k) = NaN;
                    linewidth_org(i,k) = NaN;
                end
            end
        end
        peak_freq_org3 = peak_freq_org2;
        peak_power_org3 = peak_power_org2;
        int_power_org3 = int_power_org2;
        linewidth_org3 = linewidth_org2;
end

% Final number of groups
[dummy GROUP] = size(peak_freq_org3);


%% Export variables
export_vector = [angle current resistance];
for ii = 1 : GROUP
    export_vector(:,ii*4) = peak_freq_org3(:,ii);
    export_vector(:,ii*4 + 1) = peak_power_org3(:,ii);
    export_vector(:,ii*4 + 2) = int_power_org3(:,ii);
    export_vector(:,ii*4 + 3) = linewidth_org3(:,ii);
end

% Write text file
name = strcat(Folder,'\Data_',log_name,'.txt');
dlmwrite(name, export_vector,  'delimiter', '\t')

%% Shows final time
toc

%% Plot
% figure
% freq_axis=12:16/(points-1):28;
% surf(1E3*abs(current_log),freq_axis,f_I_matrix')
% box on
% shading interp
% xlabel('\it I \rm (mA)','FontSize',12)
% ylabel('f (GHz)','FontSize',12)
% caxis([0 10])
% axis([45 90 19 23 -70 30])
% %set(gca,'YTick',4:2:16,'FontSize',12)
% 
% %%
% figure
% hold on
% plot(1E3*abs(current_log),peak_freq,'.k','MarkerFaceColor','k')
% % plot(1E3*abs(current_export2),f_peak2,'.r','MarkerFaceColor','r')
% xlabel('I (mA)','FontSize',12)
% ylabel('Frequency (GHz)','FontSize',12)
% axis([45 90 19 23])
% set(gca,'Linewidth',2,'FontSize',12)
% box on
% 
% 
% R_device =2.8;
% fg = 1/(1 - abs((R_device-50)/(50+R_device)));
% figure
% hold on
% plot(1E3*abs(current_log),1E12*int_power*fg,'.k',...
%     'MarkerFaceColor','k')
% % plot(1E3*abs(current_export2),1E12*int_power_export2*fg,'.r',...
% %     'MarkerFaceColor','r')
% xlabel('I (mA)','FontSize',12)
% ylabel('Integrated power (pW)','FontSize',12)
% axis([45 90 0 120])
% set(gca,'Linewidth',2,'FontSize',12)
% box on
% 
% figure
% hold on
% plot(1E3*abs(current_log),linewidth*1E-6,'.k',...
%     'MarkerFaceColor','k')
% % plot(1E3*abs(current_export2),linewidth_export2*1E-6,'.r',...
% %     'MarkerFaceColor','r')
% xlabel('I (mA)','FontSize',12)
% ylabel('Linewidth (MHz)','FontSize',12)
% axis([45 90 0 350])
% set(gca,'Linewidth',2,'FontSize',12)
% box on


% Plot extracted values of peak frequency, integrated power and linewidth

figure; clf;
plot(abs(angle_log)*1e3, peak_freq_org3/1e9, '.');

figure; clf;
plot(abs(angle_log)*1e3, int_power_org3/1e-12, '.');

figure; clf;
plot(abs(angle_log)*1e3, linewidth_org3/1e6, '.');