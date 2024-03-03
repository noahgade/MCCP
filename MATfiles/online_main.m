
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

% M: Length of your fixed block size
%    Remeber M will be fixed in the online setting
%    By tuning M, we can control false alarm and detection delay
%    Small M would yield higher false alarm but lower detection delay
%    Big M would yield lower false alarm but higher detection delay
%    Therefore, M would play a key role in the detecting results.
%    M to be 50~200 would be a reasonable number


% N: Number of blocks (5 to 10 would be a reasonable number)

M = 100; % you can also try M = 20, 50, 100, 150, 200 to get a sense 
N = 5;


%% estimate bandwidth

%  r = 1 by default (median trick);
%  You can also try other bandwidth by tuning r in comparision with median trick;

%  bandw1(.) is a function to estimate bandwidth from reference data using median trick

[~, L] = size(reference); % Note, you can choose a segment of reference data to estimate bandwidth 

reference = reference(:, 1:L);
r = 1;
bandw = r * bandw1(reference);


%% estimate variance
          
reference = reference(:, 1:L); % Note, you can choose a segment of reference data to estimate variance
S_var = est_var_online(reference, bandw, M, N);


%% Estimate the threshold given average run length (ARL) (using Theorem 2)
%  Set ARL = 5000~10000 would be reasonable


ARL = 5000; 
b = find_thre(M, ARL); 


%% Online monitor the change; 
%  Note: The algorithm would stop once the statistic hits the pre-defined threshold

M_stat = online_compute_stat (reference, data_test, b, bandw, S_var, M, N);

    
%% Plot computed statistic

subplot(2,1,2);
plot(M_stat); % note: the algorithm stops whenever the change-point is detected


xlabel('Index');
ylabel('M-Statistic');


hold on
% localize the change-point location

loc = find(M_stat==max(M_stat));
plot(loc, M_stat(loc), 'ro'); hold on 


 
fprintf('-- The estimated threshold is %f given ARL %f.\n', b, ARL);

if M_stat(loc) > b
  
  fprintf('-- There exists a change-point in testing data at location %d. \n', loc);
  
else
  
  fprintf('There is no change-point in testing data. \n');
  
end
  

 

  



         
 