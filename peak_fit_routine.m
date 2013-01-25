function [peak_power_out,int_power_out,linewidth_out,peak_freq_out,PSDreal,error,y_fit_out]=...
    peak_fit_routine(frequency,spectrum,clean_spectrum_dBw,...
    gain,rbw,refl_corr,s,correction,Max_peak,threshold,noise_level,min_power,ml)

% Performs fitting of the spectra by first finding peaks and then fitting
% lorentzians. The low-power peak correction is found by forcing extra
% peaks and find the answer that reduces the absolute error of the fitting.
% The information of the location of the peak is randomly assigned from the
% mean of the previous peaks plus a Gaussian random index. The high-power
% correction forces less peaks and limits the results to a certain power
% level in order to reduce the influence of noisy data.
% 
% Inputs:
% frequency         = frequency axis in Hz
% spectrum          = spectrum in dB
% clean_spectrum    = clean spectrum in dB
% gain              = measurement gain
% rbw               = resolution bandwidth
% s                 = flag
% correction        = flag for the power correction type
% Max_peak          = maximum number of peaks to fit
% threshold         = input to the peak detection routine
% noise_level       = calculated noise level of the measurement
% min_power         = minimum power expected from the measurement.
%                     Estimated from the standard deviation of the data
% ml                = maximum linewidth to fit
% 
% Outputs:
% peak_power        = peak power of each fitted peak
% int_power         = integrated power of each fitted peak
% linewidth         = linewidth of each fitted peak
% PSDreal           = power spectrum in mW
% peak_freq         = freqeucny of each fitted peak
% error             = absolute fitting error


    %% Variables
    df = frequency(2)-frequency(1);
    delta_error = -1;
    count = 0;
    flag = 0;
    
    %% Initialize empty output variables
    peak_power_out = [];
    int_power_out = [];
    linewidth_out = [];
    peak_freq_out = [];
    PSDreal = [];
    error = Inf;
    y_fit_out = [];
    
    %% Peak detection routine
    if s~=0     
        % Calls the function that semiautomatically identifies multiple peaks
        % on each spectrum above a threshold. The output can be a vector.
        [f_peak,SA_peak,ind_peak_pre]=peak_finder_auto(frequency,smooth(clean_spectrum_dBw),threshold);
        
        % Sets the number of peaks found and selects the ones above the
        % minimum power
        dummy = 1;
        for i = 1 : length(ind_peak_pre)
            if (1/rbw)*refl_corr*10.^(clean_spectrum_dBw(ind_peak_pre(i))/10) >= min_power
                ind_peak(dummy) = ind_peak_pre(i);
                dummy = dummy + 1;
            end
        end
        if dummy == 1
            ind_peak = ind_peak_pre;
        end
        clear dummy
        peaks_no=length(ind_peak);
        
        % Raises a flag if no peak is found
        if isempty(ind_peak)==0
            j=1;
        else
            j=0;
        end
    end
    
    %% Peak fitting routine
    if j == 1
        % Perform peak fitting with the found peaks
        [peak_power,int_power,linewidth,peak_freq,PSDreal,error_old,y_fit]=...
            pow_and_linew_calc_lorentzian_3peaks(frequency,spectrum,...
            clean_spectrum_dBw,ind_peak,gain,rbw,refl_corr,noise_level,ml,min_power);
    else
            peak_freq = 0;
            peak_power = 0;
            int_power = 0;
            linewidth = 0;
            y_fit = 0;
            error_old = 10000000000;
    end
    
    % Power corrections
    switch correction
        case 0
            % Store results
            peak_freq_out = peak_freq;
            peak_power_out = peak_power;
            int_power_out = int_power;
            linewidth_out = linewidth;
            y_fit_out = y_fit;
            error = error_old;
            
        case 1
            % Low power peak correction. Here usually the peaks found are
            % fewer than the expected Max_peak
            
            % Loop to find solutions
            while ((delta_error <= 0) && (peaks_no <= Max_peak))
                % Save current fit
                if peaks_no ~= 0
                    if count == 0 
                        peak_freq_out(1:peaks_no) = peak_freq;
                        peak_power_out(1:peaks_no) = peak_power;
                        int_power_out(1:peaks_no) = int_power;
                        linewidth_out(1:peaks_no) = linewidth;
                        y_fit_out = y_fit;
                    end

                    % Update ind_peak with fitted values
                    ind_peak = ceil((peak_freq - min(frequency))/df)';
                    if ind_peak(peaks_no) > length(frequency)
                       ind_peak = length(frequency);
                    end
                end

                % Increase number of peaks
                peaks_no = peaks_no + 1;
                
                % Give a guess based on the mean of the previous peaks and
                % a random variation. The random variation depends on the
                % number of peaks so that with more detected peaks, the
                % frequency range is bigger
                if peaks_no > 1
                    ind_peak(peaks_no) = ceil(mean(ind_peak) +...
                        abs((peaks_no - 1 ) * unidrnd(ceil(2 * linewidth_out(peaks_no-1)/df))));
                else
                    ind_peak(peaks_no) = abs(unidrnd(floor(max(frequency)/df)));
                end
                % Correct if the guess is out of the available indices
                for jj = 1 : length(ind_peak)
                    if (ind_peak(jj) < 1)
                        ind_peak(jj) = 1;
                    elseif (ind_peak(jj) > length(frequency))
                            ind_peak(jj) = length(frequency);
                    end
                end

                % Calculate new spectra
                [peak_power,int_power,linewidth,peak_freq,PSDreal,error_new,y_fit]=...
                    pow_and_linew_calc_lorentzian_3peaks(frequency,spectrum,...
                    clean_spectrum_dBw,ind_peak,gain,rbw,refl_corr,noise_level,ml,min_power);

                % Check if it is a valid peak by comparing its power with
                % the noise floor.
                % If it cannot find a suitable peak in two
                % consecutive runs, the correction is cancelled.
                
                if peak_power(peaks_no) < min_power - noise_level
                    peaks_no = peaks_no - 1;
                    flag = 1;
                    min_power = min_power / 10;
                else
                    for ii=1:peaks_no - 1
                        if (peak_freq(peaks_no) <= peak_freq(ii) + 1*linewidth(ii)) && (peak_freq(peaks_no) >= peak_freq(ii) - 1*linewidth(ii)) && (flag == 0)
                            peaks_no = peaks_no - 1;
                            flag = 1;
                        end
                    end
                end
                if flag == 1
                    count = count + 1;
                else
                    count = 0;
                end
                flag = 0;
                if count == 5
                    break;
                end
                    
                % Compare error
                error = error_old;
                delta_error = abs(error_new) - abs(error_old);
                error_old = error_new;
            end
            
        case 2
            % High-power threshold. Here the peaks found are usually higher
            % than Max_peak.

            % Allocate empy matrix
            ind_new = [];
            peak_freq_out = 0;
            peak_power_out = 0;
            int_power_out = 0;
            linewidth_out = 0;
            y_fit_out = 0;

            % Check if the detected peaks have enough power compared with
            % the minimum power variable
            count = 1;
            for kk = 1 : peaks_no
                % Store only peaks with sufficient power
                if peak_power(kk) > min_power - noise_level
                    ind_new(count) = kk;
                    count = count + 1;
                end
            end
            
            % If there are more than the maximum number of peaks expected,
            % select the most powerful
            if length(ind_new) > Max_peak
                dummy = sort(peak_power(ind_new));
            
                % Choose the maximum power data
                for kk = 1 : Max_peak
                    ind_new2(kk) = find(peak_power == dummy(end - kk + 1));
                end
            else
                ind_new2 = ind_new;
            end
            
            % Correct fit if needed
            if isempty(ind_new2) == 0
                
                % Create new peak positions)
                for kk = 1 : length(ind_new2)
                    dummy = find(frequency >= peak_freq(ind_new2(kk)));
                    ind_peak_new(kk) = dummy(1);
                end
                
                % Repeat fitting with new conditions
                [peak_power,int_power,linewidth,peak_freq,PSDreal,error_old,y_fit]=...
                    pow_and_linew_calc_lorentzian_3peaks(frequency,spectrum,...
                    clean_spectrum_dBw,ind_peak_new,gain,rbw,refl_corr,noise_level,ml,min_power);
        
                % Store corrected fits
                peak_freq_out = peak_freq;
                peak_power_out = peak_power;
                int_power_out = int_power;
                linewidth_out = linewidth;
                y_fit_out = y_fit;
            end
            error = error_old;
    end
    