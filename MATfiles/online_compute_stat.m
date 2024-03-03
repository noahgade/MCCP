function M_stat = online_compute_stat (reference, data_test, b, bandw, S_var, M, N)

       [~, n] = size (data_test); % n is the length of the testin data
 
 %% initialization
       
       Pool = reference;

       index=M;  
       
       X = reference(:, 1 : N*M); % X is reference data
       X_sample = datasample( Pool, N*M, 2, 'Replace', false ); % X is sampled reference data      
       
       Y = data_test; % Y is testing data
 
       Kxx_post = fKxx1( Y(:, index-M+1:index), Y(:, index-M+1:index), M, bandw, 1); % M by M
       
       Kxx_pre=[];
       Kxx_cross=[];
       
       for j = 1:N     
         
        Kxx_pre = [Kxx_pre; fKxx1(X(:, (j-1)*M+1: j*M),X(:,(j-1)*M+1: j*M),M,bandw,1)]; %  N*M by M 
        Kxx_cross = [Kxx_cross; fKxx1(X(:,(j-1)*M+1: (j-1)*M+M), Y(:, index-M+1: index),M,bandw,2)]; % N*M by M
       
       end
       
       M_stat = zeros(1,n); 
       
       
 %% introducing sliding window to online monitor change-point   
       
  for index = M:size(data_test, 2)
         
          MMD = [];        
          temp1 = 1/M/(M-1) * sum( Kxx_post(:) );
          
         for j = 1:N  
             A = Kxx_pre( j*M-M+1:j*M,  1:M );
             C = Kxx_cross( j*M-M+1:j*M,  1:M);
             
             MMD(j) = 1/M/(M-1)*sum(A(:))+ temp1 - 2/M/(M-1)*sum(C(:));
         end 
       
         M_stat(index)=( mean(MMD)./sqrt(S_var)) ;  % compute statistic
         
     
        % online update given new data
        
         Pool = [ Pool, Y(:, index-M+1)];  % update "Pool" if no change 
   
         % index is the end of sliding window 
         
         % given new data, update Kxx_post
         
         Kxx_post(1:M-1, 1:M-1) = Kxx_post(2:M, 2:M);
         temp = fKxx1( Y(:, index-M+1:index), Y(:,  index), M, bandw, 3); % M by 1
         Kxx_post(:,M) = temp;
         Kxx_post(M,:) = temp';
            
        % given new data stream, we sometimes need to random sample data from Pool as
        % reference data
        
         r = mod(index, M); % if r==0, we need to sample the reference data blocks
                  
           if r == 0
               
               for j = 1:N
                 
            % update Kxx_pre
                  
               Kxx_pre( (j-1)*M+1:(j-1)*M+M-1, 1:M-1 ) = Kxx_pre((j-1)*M+2:(j-1)*M+M, 2:M);
               temp = fKxx1( X_sample(:,(j-1)*M+1:(j-1)*M+M), X_sample(:, (j-1)*M+M),M,bandw,3);
               Kxx_pre((j-1)*M+1:(j-1)*M+M, M) = temp;
               Kxx_pre((j-1)*M+M,:) = temp';
                                
           %  update Kxx_cross
                     
               Kxx_cross((j-1)*M+1:(j-1)*M+M-1, 1:M-1) = Kxx_cross((j-1)*M+2:(j-1)*M+M, 2:M);
               temp1 = fKxx1(X_sample(:,(j-1)*M+1:(j-1)*M+M), Y(:,index),M,bandw,3);
               temp2 = fKxx1( Y( :,index-M+1:index), X_sample(:,(j-1)*M+M),M,bandw,3);
               Kxx_cross((j-1)*M+1:(j-1)*M+M, M) = temp1;
               Kxx_cross((j-1)*M+M,:) = temp2';
              
               end
                
                X = X_sample;  
                
                X_sample = datasample( Pool, N*M, 2,  'Replace', false );   
             
           elseif r~=0
           
               for j=1:N
                 
            % update Kxx_pre
            
               Kxx_pre((j-1)*M+1:(j-1)*M+M-1, 1:M-1) = Kxx_pre((j-1)*M+2:(j-1)*M+M, 2:M);
               temp = fKxx1( [X(:, j*M-(M-r)+1:j*M), X_sample(:,(j-1)*M+1:(j-1)*M+r)], X_sample(:, (j-1)*M+r),M,bandw,3);
               Kxx_pre((j-1)*M+1:(j-1)*M+M, M) = temp;
               Kxx_pre((j-1)*M+M,:) = temp';
                                
           % update Kxx_cross
                     
               Kxx_cross((j-1)*M+1:(j-1)*M+M-1, 1:M-1) = Kxx_cross((j-1)*M+2:(j-1)*M+M, 2:M);
               temp1 = fKxx1([X(:, j*M-(M-r)+1:j*M), X_sample(:,(j-1)*M+1:(j-1)*M+r)] , Y(:,index),M,bandw,3);
               temp2 = fKxx1( Y( :,index-M+1:index), X_sample(:,(j-1)*M+r),M,bandw,3);
               Kxx_cross((j-1)*M+1:(j-1)*M+M, M) = temp1;
               Kxx_cross((j-1)*M+M,:) = temp2';
                           
              end   
               
           end
  
 
  end
  
end
