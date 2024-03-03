function [K]=fKxx1( A, B,  M,bandw, flag )
    
% Here, we use A refer to pre-data and
% B refer to post-data

switch flag
    
    
    case 1  
        
        
        D = []; 
        for i = 1:M
            temp = bsxfun (@minus, A(:,(i+1):M), B(:, i)); 
            D(i,(i+1):M) = sum(temp.^2,1);
           
        end

          D = D+D'; % D is the pairwise distance matrix
          K = exp(-D/2/bandw);
          K(logical(eye(size(K)))) = 0; % set the diagonal entries 0
                     
        
    case 2
        
         
        D = [];
        for i = 1:M
            temp = bsxfun (@minus, A(:,1:M), B(:, i));
            D(1:M,i) = sum(temp.^2,1);      
        end
         
         K = exp(-1/2/bandw * D);
         K(logical(eye(size(K)))) = 0; % set the diagonal entries 0
   
                
         
    case 3  
        
        K = [];
        temp = bsxfun (@minus, A, B);
        D = sum(temp.^2,1);
        K = exp(-1/2/bandw * D);
        K(M) = 0;
          
        
end

