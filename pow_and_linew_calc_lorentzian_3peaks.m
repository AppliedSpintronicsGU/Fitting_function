function [peak_power,int_power,linewidth,peak_freq,PSDreal,error,y_fit]=...
    pow_and_linew_calc_lorentzian_3peaks(frequency,spectrum,clean_spectrum_dBw,...
    ind_peak,gain,rbw,refl_corr,noise_level,ml,min_power)

% Calculates the real PSD, accounting for amplifier gain and resolution
% bandwidth. It does not take into account reflections due to impedance
% mismatch between the sample and the transmission line.
% 
% Inputs:
% frequency         = frequency axis in Hz
% spectrum          = spectrum in dB
% clean_spectrum_dBw= clean spectrum in dBw
% ind_peak          = indices of the found peaks
% gain              = measurement gain
% rbw               = resolution bandwidth
% refl_corr         = reflection correction given in the main
% noise_level       = noise_level calculated in the main
% ml                = maximum linewidth to fit
% 
% Outputs:
% peak_power        = peak power of each fitted peak
% int_power         = integrated power of each fitted peak
% linewidth         = linewidth of each fitted peak
% PSDreal           = power spectrum in mW
% peak_freq         = freqeucny of each fitted peak
% error             = absolute fitting error

%% Initialization

% Convert the spectrum to linear scale
PSDreal = (1/rbw)*refl_corr*10.^(clean_spectrum_dBw/10);
PSDreal = PSDreal - noise_level;

% Create output vectors
peak_power = NaN(1,length(ind_peak));
int_power = peak_power;
linewidth = int_power;
peak_freq = int_power;
y_fit = int_power;
error=0;

%% Fitting procedure

if ~isempty(ind_peak)
    % Define variables
    fitting_factor = 1E15;
    xf = frequency/1E9;
    yf = PSDreal*fitting_factor;

    % Build string of fits and fitting options
    clear equation; clear coeff;clear Lower;clear Upper;clear Startpoint;
    equation = '(A1/pi)*(Gamma_1/2)/((x-x1)^2+(Gamma_1/2)^2)';
    coeff = {'A1','Gamma_1','x1'};
    Lower = [0, 0, min(xf)];
    Upper = [1E3, ml, max(xf)];
    Startpoint = [0 0.01 xf(ind_peak(1))+rand(1,1)/1000];
    for kk=2:length(ind_peak)
        equation = strcat(equation,sprintf('+(A%g/pi)*(Gamma_%g/2)/((x-x%g)^2+(Gamma_%g/2)^2)',kk,kk,kk,kk));
        coeff((kk-1)*3+1:(kk-1)*3+3) = {sprintf('A%g',kk),sprintf('Gamma_%g',kk),sprintf('x%g',kk)};
        Lower((kk-1)*3+1:(kk-1)*3+3) = [0, 0, min(xf)];
        Upper((kk-1)*3+1:(kk-1)*3+3) = [1E3, ml, max(xf)];
        Startpoint((kk-1)*3+1:(kk-1)*3+3) = [0 0.01 xf(ind_peak(kk))];
    end
    s = fitoptions('Method','NonlinearLeastSquares','Lower',Lower,'Upper',Upper,'Startpoint',Startpoint,...
        'MaxIter',1E8,'TolX',1e-10,'TolFun',1e-10,'DiffMinChange',10e-11);
    f=fittype(equation,'coefficients',coeff,'options',s);
    
    % Perform fitting
    L = fit(xf,yf,f);
    cL = coeffvalues(L);
    
    % Calculate ouput
    for kk=1:length(ind_peak)
        int_power(kk) = cL((kk-1)*3+1) * 1E9 / fitting_factor;
        linewidth(kk) = cL((kk-1)*3+2) * 1E9;
        peak_power(kk) = (2/pi) * cL((kk-1)*3+1) / fitting_factor / cL((kk-1)*3+2);
        peak_freq(kk) = cL((kk-1)*3+3) * 1E9;
    end
    
    % Calculate spectra
    y_fit = feval(L, xf);
    
    % Calculate error
%     error = sum((y_fit-abs(yf)) .* ((sign(y_fit-abs(yf))+1)/2));
    error = sum(abs(y_fit-yf) .* (sign(yf-min_power*fitting_factor)+1)/2);
    
    % Update fit amplitude
    y_fit = y_fit / fitting_factor;
end
