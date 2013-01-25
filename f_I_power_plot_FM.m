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
% | Current[mA] | Resistance | Field | Peak_freq_G1 | Peak_power_G1 | Int_power_G1 | Linewidth_G1 | Peak_freq_G2 | etc
% 
%  Each column is a vector of length specified by the current.

%% How to
% 
% The variables below are the control parameters of the fitting.
% 
% When correction is set to zero, the code works as the previous version 
% while setting it to 1 would try to fit any small peaks that can find 
% (even close peaks).
% 
% Max_peak is the maximum number of peaks that the fitting will try to
% find. This is to prevent irrelevant features to be fitted. Set it to the
% number of important peaks you expect to find or to a reasonable ammount
% to speed up the fitting.
% 
% N_attemp is the number of times the correction algorithm will attempt to
% find low power peaks. The correction is basd on randomly assigning peak
% positions so the lower the peak, the greater the number of attempts
% needed to find it. The default value 5 is good for reasonably strong
% peaks.
% 
% The threshold only appears in the peak_finder_auto. By default is 2.
% Lowering might help you to find low power peaks but if too many are found
% the fitting might be wrong.

%% Important Variables!
correction = 1;          % Allow low-power correction
Max_peak = 40;           % Maximum peaks to fit
N_attempt = 5;           % Allow several attempts to find the best fit
threshold = 1.7;         % Threshold for initial peak detection
gain = 22;               % Gain of the amplifier used during measurements
REGROUP = 1;             % Set to 1 for the groups to be "regrouped"
FM = 1;                  % Flag to regoup when fitting FM data
bp_trsh = 30;            % Set the threshold value for the bad point detection.
                         % Bad points are caused by miscommunication in LabView
                         % and are usually present at the beginning of spectra
                         % Adjust the "bp_trsh" value in such a way so that
                         % it is always higher then the maximum peak level,
                         % otherwise useful data might be lost.

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
    [RFP1_log RFF1_log RFP2_log RFF2_log current_log spectrum_path]=textread(log_path,...
        '%f %f %f %f %f %s',-1,'headerlines',7,'delimiter','\t'); % 7 is standard
end

%%

seqfig1_1=RFF1_log/1e6;
seqfig1_1name=' F_m_1 = ';
seqfig1_1units=' MHz';
seqfig1_2=RFF2_log/1e6;
seqfig1_2name=' F_m_2 = ';
seqfig1_2units=' MHz';
seqfig2=current_log*1e3;
seqfig2name=' I_d_c = ';
seqfig2units=' mA';
seqfig3=RFP2_log*1e3;
seqfig3name=' RFP_2 =';
seqfig3units=' mV';
seqfig3name_Iac2=' I_m_2 =';
seqfig3units_Iac2=' mA';
seqaxis=RFP1_log*1e3;
axisname='RF Power 1 [mV]';
% seqaxis=RFF_log/1e6;
% axisname='RF Frequency [MHz]';

seqfig1_1_unique=unique(seqfig1_1);
seqfig1_2_unique=unique(seqfig1_2);
seqfig2_unique=unique(seqfig2);
seqfig3_unique=unique(seqfig3);
seqaxis_unique=unique(seqaxis);

% Reads sequentially the rows of the .log file, and points each time to a
% new file with the spectrum.
for lll=1:length(seqfig2_unique)
% for lll=1

for mmm=1:length(seqfig1_2_unique)
% for mmm=1
    
for jjj=1:length(seqfig1_1_unique)
% for jjj=4

for nnn=1:length(seqfig3_unique)
% for nnn=27
    
    logindex=find((abs(seqfig1_1-seqfig1_1_unique(jjj))<1e-1)&(abs(seqfig1_2-seqfig1_2_unique(mmm))<1e-1)&(abs(seqfig3-seqfig3_unique(nnn))<1e-1)&(abs(seqfig2-seqfig2_unique(lll))<1e-1));

for z=1:length(logindex)
%     current_log(i);
    z
    kkk=logindex(z);
    clear fid
    fid=spectrum_path{kkk};
    fid = regexp(fid, '\\', 'split');
    fid = fid(end);
    fid=char(strcat(pathname,fid));
    
%     steps = length(logindex);
%     rfp1 = zeros(steps,1);
%     rff1 = zeros(steps,1);
%     rfp2 = zeros(steps,1);
%     rff2 = zeros(steps,1);
%     current = zeros(steps,1);
%     field = zeros(steps,1);
%     resistance = zeros(steps,1);

    % Reads the measurement setting from the header of each spectrum.
    [temp1 rfp1(z)]=textread(fid,'%s %f',1,'headerlines',5);
    [temp1 rff1(z)]=textread(fid,'%s %f',1,'headerlines',6);
	[temp1 rfp2(z)]=textread(fid,'%s %f',1,'headerlines',7);
	[temp1 rff2(z)]=textread(fid,'%s %f',1,'headerlines',8);
    [temp1 current(z)]=textread(fid,'%s %f',1,'headerlines',9);
    [temp1 temp2 field(z)]=textread(fid,'%s %s %f',1,'headerlines',10);
    [temp1 resistance(z)]=textread(fid,'%s %f',1,'headerlines',11);
    [temp1 temp2 vbw]=textread(fid,'%s %s %f',1,'headerlines',12);
    [temp1 temp2 rbw]=textread(fid,'%s %s %f',1,'headerlines',13);
    [temp1 temp2 points]=textread(fid,'%s %s %f',1,'headerlines',14);
    [temp1 temp2 startf]=textread(fid,'%s %s %f',1,'headerlines',15);
    [temp1 temp2 stopf]=textread(fid,'%s %s %f',1,'headerlines',16);

    % Recalculate the actual resistance of the STO. Uncomment next line
    % if the correction is needed.
%     resistance = resistance - 6.7; % 6.7 Ohm is the measured lead
                                   % resistance including bias tee and
                                   % circulator.
    
    % Create a matrix with the spectrum, whose length is the number of data
    % points acquired by the spectrum analyzer.
    [frequency spectrum clean_spectrum]=textread(fid,'%f %f %f',points,...
        'headerlines',21,'delimiter','\t');
    
%     numberOfFrequencyPoints = length(frequency);
%     f_I_matrix = NaN(length(logindex),numberOfFrequencyPoints);
    f_I_matrix = NaN(length(logindex),points);
    
    % Identifies singular points created by Labview communication errors,
    % and set the value of the clean spectrum to be 0 at the single point.
    error_ind=find(abs(clean_spectrum)>40);
    if isempty(error_ind)==0
        clean_spectrum(error_ind)=0;
    end
    
    error2_ind=find(clean_spectrum<-2);
    if isempty(error2_ind)==0
        clean_spectrum(error2_ind)=0;
    end
    
    error3_ind=find(isnan(clean_spectrum)==1);
    if isempty(error3_ind)==0
        clean_spectrum(error3_ind)=0;
    end
    
    % Remove NaN values from the raw data. Uses interpolation to create a data point.
    clean_spectrum=naninterp(clean_spectrum);
    spectrum=naninterp(spectrum);
    
    % Normalize spectrum
    if z == 1
        noise_dBW = max(spectrum) - max(clean_spectrum) - 30; % from dBm to dBW'
        clean_spectrum = clean_spectrum + noise_dBW - gain;
        real_dummy = (1/rbw)*10.^(clean_spectrum/10);
        noise_level = mean(real_dummy(1:100));
    else
%         noise_level = 1.7612e-20;
%         noise_dBW = -115.9744;
        clean_spectrum = clean_spectrum + noise_dBW - gain;
    end
    
    % Cutting data
%     frequency_dummy = frequency(701:1102);
%     spectrum_dummy = spectrum(701:1102);
%     clean_spectrum_dummy = clean_spectrum(701:1102);
%     clear frequency spectrum clean_spectrum
%     frequency = frequency_dummy;
%     spectrum = spectrum_dummy;
%     clean_spectrum = clean_spectrum_dummy;

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
    for kk=1:N_attempt
        [peak_power_dummy{kk},int_power_dummy{kk},linewidth_dummy{kk},...
            peak_freq_dummy{kk},PSDreal_dummy{kk},error(kk),yf{kk}]=...
            peak_fit_routine(frequency,spectrum,clean_spectrum,...
            gain,rbw,s,correction,Max_peak,threshold,noise_level);
    end
    
    % Finds the fit with minimum error
    [dummy,select] = min(error);
    if dummy~=Inf
        % Stores the best fit information
        peak_power(z,1:length(peak_power_dummy{select})) = peak_power_dummy{select};
        int_power(z,1:length(peak_power_dummy{select})) = int_power_dummy{select};
        linewidth(z,1:length(peak_power_dummy{select})) = linewidth_dummy{select};
        peak_freq(z,1:length(peak_power_dummy{select})) = peak_freq_dummy{select};
        PSD_real(z,1:length(PSDreal_dummy{select})) = PSDreal_dummy{select};
        PSD_fit(z,1:length(yf{select})) = yf{select};
%         peak_freq
%         peak_power*1E18
%         semilogy(frequency,PSD_real(z,:)*1E18,frequency,PSD_fit(z,:),'r');
%         grid;
%         waitforbuttonpress;
    else
        % Stores zeros
        peak_power(z,:) = 0;
        int_power(z,:) = 0;
        linewidth(z,:) = 0;
        peak_freq(z,:) = 0;
        PSD_real(z,:) = 0;
        PSD_fit(z,:) = 0;
    end
    f_I_matrix(z,:) = clean_spectrum;
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

%% Peak grouping including regrouping routine
% Define variables
[loop,Max_peaks] = size(peak_freq);
group_no = 1;
output{loop+1,group_no+1} = {0};
first = 0;
flag = 0;
l = 1;
number = 1;
flag_break = 0;

if FM == 0
    % Build a cell matrix with the indices of the sorted peak frequencies
    for i = 1 : loop
        for j = 1 : Max_peaks
            % Check for valid peak frequencies
            if (peak_freq(i,j) ~= NaN) && (peak_freq(i,j) > frequency(1)) && (peak_freq(i,j) < frequency(points))
                % If it is the first, save without comparing
                if (first == 0)
                    first = 1;
                    output{i,group_no} = [i,j];
                % If the first has been already stored, check for similarities
                elseif i ~= 1
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
                            if abs(peak_freq(i,j) - peak_freq(output{i-l,k}(1),output{i-l,k}(2))) < 50e6
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
                        output{i,similar(min(select))} = [i,j];
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
    
    % Regroup
    if REGROUP == 1
        % Check if less than Max_peaks have been found
        if (j - 1) < Max_peaks
            Max_peaks = j - 1;
        end
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
    else
        regroup = j - 1;
        % Substitute zeros for NaN
        for i = 1 : loop
            for k = 1 : regroup
                if peak_freq_org2(i,k) == 0
                    peak_freq_org2(i,k) = NaN;
                    peak_power_org2(i,k) = NaN;
                    int_power_org2(i,k) = NaN;
                    linewidth_org2(i,k) = NaN;
                end
            end
        end
        peak_freq_org3 = peak_freq_org2;
        peak_power_org3 = peak_power_org2;
        int_power_org3 = int_power_org2;
        linewidth_org3 = linewidth_org2;
    end
else
    % FM dedicated variables
    select = find(int_power(1,:) > 0.01E-10);
    STOMod_freq = peak_freq(1,select);
    Mod_freq = RFF1_log(1);
    range = Mod_freq / 2;
    PeakMod_freq = zeros(loop+1,1);
    
    % Determine how many single sidebands should look for
    max_sideband = floor((STOMod_freq(2) - STOMod_freq(1))/2/Mod_freq);
    
    % Grouping routine
    for j = 1 : length(STOMod_freq)
        % Carry on the grouping for main peak and corresponding sidebands
        for jj = 1 : 2 * max_sideband + 1
            % Select upper and lower sideband each time
            l = (-1)^jj * floor(jj/2);
            
            % Assignment of variable carrier frequency
            PeakMod_freq(1) = STOMod_freq(j);
            
            for i = 1 : loop
                flag = 0;
                % Check for peak inside the range and in an existent frequency
                for k = 1 : Max_peaks
                    if (abs(peak_freq(i,k) - PeakMod_freq(i) - l*Mod_freq) < range) && (PeakMod_freq(i) + l*Mod_freq > frequency(1)) && (PeakMod_freq(i) + l*Mod_freq < frequency(length(frequency)))
                        peak_freq_org(i,group_no) = peak_freq(i,k);
                        peak_power_org(i,group_no) = peak_power(i,k);
                        int_power_org(i,group_no) = int_power(i,k);
                        linewidth_org(i,group_no) = linewidth(i,k);
                        flag = 1;
                        
                        % Update carrier frequency
                        if  jj == 1
                            PeakMod_freq(i+1) = peak_freq(i,k);
                        end
                        break;
                    end
                end
                
                % Check if any peak was found
                if flag == 0;
                    peak_freq_org(i,group_no) = 0;
                    peak_power_org(i,group_no) = 0;
                    int_power_org(i,group_no) = 0;
                    linewidth_org(i,group_no) = 0;
                    
                    % Keep the previous carrier frequency
                    PeakMod_freq(i+1) = PeakMod_freq(i);
                end
            end
            
            % Store next peak
            group_no = group_no + 1;
        end
    end
    
    % Substitute zeros for NaN
    regroup = group_no - 1;
    for i = 1 : loop
        for k = 1 : regroup
            if peak_freq_org(i,k) == 0
                peak_freq_org(i,k) = NaN;
                peak_power_org(i,k) = NaN;
                int_power_org(i,k) = NaN;
                linewidth_org(i,k) = NaN;
            end
        end
    end
    peak_freq_org3 = peak_freq_org;
    peak_power_org3 = peak_power_org;
    int_power_org3 = int_power_org;
    linewidth_org3 = linewidth_org;
end
    

%% Export variables
export_vector = [rfp1' rff1' current' resistance' field'];
for ii=1:regroup
    export_vector(:,ii*4) = peak_freq_org3(:,ii);
    export_vector(:,ii*4 + 1) = peak_power_org3(:,ii);
    export_vector(:,ii*4 + 2) = int_power_org3(:,ii);
    export_vector(:,ii*4 + 3) = linewidth_org3(:,ii);
end
end
end
%% Save results to a file (don't forget to modify the filename)
framecount=sprintf('7test_Iac_linear_DBL_FreqMod_EBL_Iac_0-1_5mA_Cu28_29_6mA_%s.dat',strrep(strrep(strcat(seqfig1_2name, num2str(seqfig1_2_unique(mmm)), seqfig1_2units, '_', seqfig1_1name, num2str(seqfig1_1_unique(jjj)), seqfig1_1units), ' ', ''), '=', ''));
% saveppt(framecount,'EBL_Cu28','-r000')
save(framecount,'export_vector','-ascii','-append');
% save('5NC_Figure4c.dat','export_vector','-ascii','-append');
end
end

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
plot(peak_freq/1e9, abs(rfp1'), '.');
figure; clf;
plot(abs(rfp1'), int_power/1e-12, '.');


figure; clf;
plot(peak_freq_org/1e9, abs(rfp1'), '.');

figure; clf;
plot(abs(rfp1'), int_power_org3(:,1:10)/1e-12, '.');

figure; clf;
plot(abs(rfp1'), linewidth_org3/1e6, '.');

figure; clf;
image(15*log10(PSD_fit)+50);

figure; clf;
image(30*log10(PSD_real*2e19+5))
