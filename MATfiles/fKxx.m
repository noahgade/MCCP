function [K] = fKxx ( A, B, N, M, bandw, flag )
    
% Here, we use A refer to pre-change data and
% B refer to post-change data
% N is the number of blocks
% M is the block size




switch flag
    
    
    case 1    % Kxx_pre = fKxx(data_post, data_post,.... , flag = 1);
        
        
        D = [];
        
        for i = 1:M
          
            temp = bsxfun (@minus, A(:,(i+1):M), A(:, i)); 
            D(i, (i+1):M) = sum(temp.^2,1);
            
        end

        D = D+D'; % D is the pairwise distance matrix
        K = exp(-D/2/bandw);
        K(logical(eye(size(K)))) = 0; % set the diagonal entries 0
        
        
    case 2   % Kxx_pre=fKxx(data_pre, data_pre,.., flag = 2); 
        
          K = [];
        
        
        for j = 1:N
            
            D = [];
 
            for i = 1:M
                temp = bsxfun (@minus, A(:,(j-1)*M+i+1 : j*M ), A(:, (j-1)*M+i));
                D(i,(i+1):M) = sum(temp.^2,1);
            end

            D=D+D';
    
          temp_1 = exp(-1/2/bandw * D);
          temp_1(logical(eye(size(temp_1)))) = 0; % set the diagonal entries 0
   
          K = [K;temp_1];
          
        end
            
              
        
    case 3
        
         K = [];
         
         for j = 1:N
           
             D = [];
             
              for i = 1:M
                  temp = bsxfun (@minus, A(:,(j-1)*M+1 : j*M), B(:, i));                
                  D(i,1:M) = sum(temp.^2,1);
                
              end
         
            temp_1 = exp(-1/2/bandw * D);
            temp_1(logical(eye(size(temp_1)))) = 0; % set the diagonal entries 0
      
            K=[K;temp_1];
        
         end
        
        
end

