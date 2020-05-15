function [cost] = objective_functionFS( UB,col,nvars,A,x,i)
    costu=0;
    [row col]=size(A);
    
    ai=0;
    xi=0;
    xii=0;
    mad=0;
    madi=0;
    
    
    for sumA=1:1                                      
        ai(sumA,1)=sum(x(sumA,:));
    end
    
    for ii=1:col                                       %for centroids
       w=0;
       if x(1,ii)==1
           xi=xi+A(i,ii);
       end
       
       
       
    end
    xii=xi/ai;
     for ii=1:col                                       %for centroids
       w=0;
       if x(1,ii)==1
           if A(i,ii)~=0
           mad=mad+abs(A(i,ii)-xii);
           end
        end
     end
     madi=mad/ai;
     
     
     
     
                                        %cost all clusters for first solution
     cost=madi;
end

