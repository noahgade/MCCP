
clear
clc

close all



%% load your reference (background) data & testing data here:

%  Your objective is to detect whether there is potential change-point in
%  your testing data.
%  Your reference data should contain NO change-point.

%% One dimensional example: testing data change from N(0,1) to N(0,2);
%                        change-point location = 250



 
     load reference_variance_1D.mat
     load data_test_variance_1D.mat
     figure; subplot(2,1,1); plot(data_test); xlabel('Location'); ylabel('Signal');


     
%% High dimensional example (Just uncomment the following lines) 
%                        change-point location = 100



%     load reference_randomgraph.mat
%     load data_test_randomgraph.mat
%     figure; subplot(2,1,1); [x, y] = find(data_test); scatter(y, x); xlabel('Time'); ylabel('Random Graph');

    
    
%% Set your parameters here:

% M: Length of your testing data; or you can set any length you want
%    Remeber M will influence your detection power. For weak signals, longer
%    M is needed.

% N: Number of blocks (5 to 10 would be a reasonable number)
%    Make sure your input reference data with length >M*N !!


[~,M] = size( data_test );
N = 5;


%% Estimate bandwidth from reference data

%  r = 1 by default (median trick);
%  You can also try other bandwidth by tuning r in comparision with median trick;

%  bandw1(.) is a function to estimate bandwidth from reference data using median trick

reference = reference(:, 1:N*M);

r = 1;
bandw = r * bandw1(reference);

%% Estimate Variance of the statistic (Lemma 1 in our paper)

%  est_var(.) is a function to estimate variacen of statistic under null (using lemma 1)

S_var = est_var (reference, bandw, M, N);


%% Compute kernel matrix  

%  fKxx(data1, data2, N, M, bandw, flag) is a function to oompute kernel matrix


Kxx_post = fKxx(data_test, data_test, N, M, bandw, 1);  
Kxx_pre = fKxx(reference, reference, N, M, bandw, 2);
Kxx_cross = fKxx(reference, data_test ,N, M, bandw, 3) ; 


%% Compute statistic

% M_stat is your computed statistic
 
M_stat = zeros(1,M);
 
 for B = 2:M  % scan your testing data by tuning the block size
   
         MMD = [];
         
         T = Kxx_post(M-B+1:M, M-B+1:M);
         
         
         for j = 1:N % compute MMD for all blocks and take an average to get the statistic 
           
             A = Kxx_pre( j*M-B+1:j*M,  M-B+1:M );
             C = Kxx_cross( j*M-B+1:j*M,  M-B+1:M);
             
             MMD(j) = 1/B/(B-1)*sum(A(:))+ 1/B/(B-1)*sum(T(:)) - 2/B/(B-1)*sum(C(:));
             
         end 
   
         M_stat(M-B+1) = mean(MMD) ./ sqrt(S_var(B-1)) ;
        
 end

%% Estimate the threshold given significance level $alpha$ (using Theorem 1)

%  You can change alpha = 0.1, 0.05, 0.01 ... for different significance level
%  tail_est(.) is a function to get threshold $b$ (theretical estimation)

alpha = 0.05;  
b = tail_est (alpha, M);

%% Plot computed statistic

subplot(2,1,2);
plot(M_stat); hold on

% localize the change-point location

loc = find(M_stat==max(M_stat));
plot(loc, M_stat(loc), 'ro'); hold on 
xlabel('Index');
ylabel('M-Statistic');

 
fprintf('-- The estimated threshold is %f given significance level %f.\n', b, alpha);

if M_stat(loc) > b
  
  fprintf('-- There exists a change-point in testing data at location %d. \n', loc);
  
else
  
  fprintf('There is no change-point in testing data. \n');
  
end
  











