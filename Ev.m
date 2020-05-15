 function [FMarco,Fmicro ,pp, r] = Ev( truedata,labelt,ub )
pp=0;
r=0;
f_measure=0;
precision=0;
alpha=length(truedata);
[truedata] = sort(truedata);
[labelt] = sort(labelt);
for i=1:ub-1 
idx = (truedata()==i);

p = length(truedata (:,idx));
nn = labelt(:,idx)==i;
if nn==0;
    n=0;
else
n=sum(nn~=0);
end
N = p+n;
%[truelabels, S] = sort(truelabels);
tp = sum(truedata(idx)==labelt(idx));
tn = sum(truedata(~idx)==labelt(~idx));
fp = n-tn;
fn = p-tp;

tp_rate = tp/p;
tn_rate = tn/n;

precision =  (tp/(n));
 if isnan(precision)
    precision=0;
end
pp=pp+precision;

recall = tp/(p);
 if isnan(recall)
    recall=0;
end
r=r+recall;

f1=(2*precision*recall)/(precision + recall);
 if isnan(f1)
    f1=0;
end 
f_measure =(f_measure+ f1); 
end 
FMarco=(f_measure/(ub-1));
pp=pp/(ub-1);
rr= r/(ub-1);
Fmicro= (2*pp*rr)/(pp + rr);
% fprintf('P_Macro: %4.4f %% \n', pp);
% fprintf('R_Macro: %4.4f %% \n', rr);
% fprintf('FMacro: %4.4f %% \n', FMarco);
% fprintf('F_Micro: %4.4f %% \n', Fmicro);
% 
 end 
