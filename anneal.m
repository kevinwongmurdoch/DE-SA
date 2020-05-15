function [minimum,fval] = anneal(UB,LB,col,row,image,parent)


newsol =@(parent)(parent + (rand(size(parent))-rand(size(parent)))*0.5);      % neighborhood space function
minT = 1e-8;         % stopping temp
cool = @(T) (.8*T);        % annealing schedule
minF = -Inf;
max_consec_rejections = 1000;
 max_try = 300;
 max_success = 20;
 MaxIter=Inf;
 report = 1;
k = 1;                           % boltzmann constant

% counters etc
itry = 0;
success = 0;
finished = 0;
consec = 0;
T = 1;
initenergy = objective_functionFS( UB,col,row,image, parent',1);
oldenergy = initenergy;
total = 0;
if report==2, fprintf(1,'\n  T = %7.5f, loss = %10.5f\n',T,oldenergy); end

while ~finished && total<MaxIter;
    fprintf('T: %g, %i, %i, %i\n',T,success,consec,itry)
    itry = itry+1; % just an iteration counter
    current = parent;    
    
    
    % % Stop / decrement T criteria
    if itry >= max_try || success >= max_success;
        
        if T < minT || consec >= max_consec_rejections;
            finished = 1;
            total = total + itry;
            break;
        else
            T = cool(T);  % decrease T according to cooling schedule
            if report==2, % output
                fprintf(1,'  T = %7.5f, loss = %10.5f\n',T,oldenergy);
            end
            total = total + itry;
            itry = 1;
            success = 1;
        end
    end
    
    newparam = fix(newsol(current));
 for iter= 1:row
            if (newparam(iter) < LB ||newparam(iter)>=UB)
                newparam(iter)=fix(LB +(UB-LB).*rand);
            end;
 end 
   newenergy = objective_functionFS( UB,col,row,image, newparam',1);
    
    if (newenergy < minF),
        parent = newparam; 
        oldenergy = newenergy;
        break
    end
    
    if (newenergy-oldenergy >0)
        parent = newparam;
        oldenergy = newenergy;
        success = success+1;
        consec = 0;
    else
        if (rand < exp( (newenergy-oldenergy)/(k*T) ));
            parent = newparam;
            oldenergy = newenergy;
            success = success+1;
        else
            consec = consec+1;
        end
    end
   % T=0.93*T;
end

minimum = parent;
fval = oldenergy;

if report;
    fprintf(1, '  Final temperature:       \t%g\n', T);
    fprintf(1, '  Consecutive rejections:  \t%i\n', consec);
    fprintf(1, '  Number of function calls:\t%i\n', total);
    fprintf(1, '  Total final loss:        \t%g\n', fval);
end

