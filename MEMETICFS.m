function  MEMETICFS(image,truelabels)

% % -------------------------------------------------------------------------
% % original parametets 
K=3;
ub=4;
lb=1;
cr = 0.5;   
F = 0.2;
Fmacro=zeros(10,1);
Fmicro=zeros(10,1);
%MaxF=zeros(20,1);
BFarray = zeros(10,1);
ADDC=zeros(10,1);
max_iter = 10;
% % -------------------------------------------------------------------------
%SET PARAMETERS
%nd = 8; %dimension of genes----------------
UB=2;
LB=0;
[row col]=size(image);
old_f=col;
a_hi = UB; %boundary values of variables (gene)
a_lo = LB;
%paramters for Differential Evolution
% -------------------------------------------------------------------------
 pop_size=col; % set population size


%GENERATE INITIAL POPULATION
a = fix((a_hi-a_lo).*rand(pop_size,20)+a_lo); 
par = [a];

[ro co]=size(par);
for iter= 1:co
costassignment(iter,:)=par(:,iter);
end;

for iter= 1:co
fitness(iter,1)=objective_functionFS( UB,col,row,image, costassignment(iter,:),iter);
end;

%FIND INITIAL FITNESS VALUES FOR ENTIRE POPULATION
for ii=1:max_iter
 for j=1:co
     
       % DE STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Pick 3 parents
       r(1) = floor(rand()* co) + 1;
       while (r(1)==j)
	     r(1) = floor(rand()* co) + 1;
       end
       r(2) = floor(rand()* co) + 1;
       while ((r(2)==r(1))||(r(2)==j))
         r(2) = floor(rand()* co) + 1;
       end
       r(3) = floor(rand()* co) + 1;
       while ((r(3)==r(2))||(r(3)==r(1))||(r(3)==j))
         r(3) = floor(rand()* co) + 1;
       end
       
%mutant vector
     %  v=(par(:,r(1))+F*(par(:,r(2))-par(:,r(3))));
       x1=par(:,r(1));
       x2=par(:,r(2));
       x3=xor(x1,x2);
       [a, ~]=(size(x1));
      
%       ra=randi([0, 1], 1,a);
%       ra;
      ra =x1 (randperm(length(x1)));
      x5=and(x3,ra);
      x6=not(x3);
      x7=and(x6,x1);
      x8=or (x5,x7);
      v=x8;
        for iter= 1:a
            if (v(iter,1)<=0)
              v(iter,1)=0; 
            else v(iter,1)=1;
%             if (v(iter,1) < a_lo ||v(iter,1)>a_hi)
%                 v(iter,1)=fix(a_lo +(a_hi-a_lo).*rand);
            end;
        end;

%recombination to generate offspring
      u=par(:,j);
       rv=rand(1,row);
       for i=1:row
           if rv(1,i)<cr        
            u(i,1)=v(i,:);
           end
       end
  
 %recompute fitness value of offsprings 
 costassignment(j,:)=u;
       new_fitness=objective_functionFS( UB,col,row,image, costassignment(j,:),j);
       %new_fitness =- tr_op(global_mean,B,C,im_size,image,u(1,1),u(1,2),u(1,3),u(1,4));
     if(new_fitness >= fitness(j))
         par(:,j)=u;
          fitness(j)=new_fitness;
        end
  end
     
    [fmax,index]=max(fitness);
                                 parent=par(:,index);
                                [new_parent,fval]=anneal(UB,LB,col,row,image, parent);
                                 if(fval > fmax)
                                     par(:,index)=new_parent;
                                      fitness(index)=fval;
                                 end 
                                   [fmax,index]=max(fitness);
   BFitness= fmax;
 
fprintf('Best MAD: %4.2f %% \n',BFitness);   
BFarray(ii) = BFitness;
      
            hold on
            pause(0.01)
            figure(1)
            title(' MAD Feature selection')
            xlabel('Iterations')
            ylabel('function value of MAD')
            plot(ii,BFitness,'--r.','LineWidth',2,...
                'MarkerEdgeColor','k')
              
end

T = table(BFarray, 'VariableNames',{'MAD'});
writetable(T,'C:\Users\32599642\Desktop\MAD.xls');

cs=sort(fitness,'ascend')';
    for iter= 1:co
    up_par(iter,:)= par(:,iter);
%     [C, label] = ind2cluster(up_par);
  %   labelt(1,:)= transpose(label);  
    end;
      [Fmin,index]=max(fitness);
       Xmin=up_par(index,:);
       BFitness = Fmin;
       BestGen=Xmin; 


%.............................
       newimage=zeros(row,col);
       ai(1,1)=sum(Xmin(1,:));
       for fs=1:col
            if Xmin(:,fs)==1
                newimage(:,fs) =image(:,fs);
            end
       end
       
fprintf('Best MAD: %4.2f % ',BFitness);  
%fprintf('new datasets: %4.2f % ',newimage);
%fprintf('feature selection: %4.2f % ',Xmin);
   
 
image= newimage; 
ximage=[];
pxx=1;
 for px=1:col
     
    if sum(newimage(:,px))~=0
        
    ximage(:,pxx)=newimage(:,px);
    pxx=pxx+1;
    end
 end
image=ximage;
[ro co]=size(image);

fprintf('number of NEW feature: %4.2f % ',co);   
fprintf('number of old feature: %4.2f % ',old_f);   
r=1-(co/old_f);
fprintf('number of old feature: %4.2f % ',r);   

%MEMATIC CLUSTERING 
% % -------------------------------------------------------------------------

for ii=1:10
BestGenc=fix(lb+((ub-lb).*rand(1,row)));
centroid = centroids( K+1,col,row,image,BestGenc);
        interleaved = functionKmeans(image, centroid, K);
        srrcost= objective_function( K+1,col,row,image,interleaved,1);
        %fprintf('BFitnessHS+k-means: %d\n', srrcost);
        %fprintf('HS+k-means: %d\n', interleaved);
   labelid=interleaved;
    [C, label] = ind2cluster(labelid);
    valid_errorate(label, truelabels);
    truedata(1,:)= transpose(truelabels);       
labelt(1,:)= transpose(label);  
[C, label] = ind2cluster(labelt);
ADDC(ii)=srrcost;
[EVAL1, EVAL2, ~, ~] = Ev( truedata,labelt,ub );     
 Fmacro(ii)=EVAL1;
 Fmicro(ii)= EVAL2;
 %MaxF(ii)=EVAL;
end
T = table(Fmacro,Fmicro, ADDC, 'VariableNames',{'Fmacro', 'Fmicro', 'ADDC'});
writetable(T,'C:\Users\32599642\Desktop\DE_SA.xls');
end

  
