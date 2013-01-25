function [f_max,PSD_max,ind_max]=peak_finder_auto(f,PSD,threshold)

f=f/1E9;
ind_max=[];
f_max=[];
PSD_max=0;

[maxima,minima]=peakdet(PSD,threshold);

if isempty(maxima)==0
    f_max=f(maxima(:,1));
    PSD_max=PSD(maxima(:,1));
    ind_max=maxima(:,1);
end