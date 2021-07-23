%%Set size and output matrices
n = 104;
ntrial = 2; % Number of independent PMF runs on each hidden increment (gives error bars)
ninc = 5; % Number of increments for fraction hidden data

%%
epsilon = 2.170330e+00;
lambda = 7.061362e-03;
momentum = 2.284520e-01;

%Read in all 60 files and append them to create matix T
%T is the collection of all 60 104x104 drug-drug interaction matrices,
%T(1-104, 1-104) is cell line 1, T(105-208, 1-104) is cell line 2, etc.
%NOTE: Make sure filepath leads to directory of CSV files

%%
cellnum = 60; %Number of celllines

%Output matrices for graphs
err_known = zeros(cellnum, ninc);
errbar_known = zeros(cellnum, ninc);
err_hidden = zeros(cellnum, ninc);
errbar_hidden = zeros(cellnum, ninc);
err_all = zeros(cellnum, ninc);
errbar_all = zeros(cellnum, ninc);
%%
files = dir('C:\Users\ronna\OneDrive\Documents\MATLAB\ComboDrugMat\*.csv'); %Replace with your directory
num_files = length(files);
results = cell(length(files), 1);
%%
for cellline = 1:cellnum
    cellline
    M = readmatrix(['C:\Users\ronna\OneDrive\Documents\MATLAB\ComboDrugMat\' files(cellline).name]);
    M = triu(M) + triu(M)'; %% Make data symmetric
    M(isnan(M))=0;

    A = M - mean(mean(M(M ~= 0)));
    A = A ./ std(M(M ~= 0),1,'all');
    normM = A; %%Normalize Values
    %Increment start and end indices for next iteration
    startindex = startindex + 104;
    endindex = endindex + 104;
    
    %Output arrays for storing error of predictions
    output_known = zeros(ninc,ntrial);
    output_hidden = zeros(ninc,ntrial);
    output_all = zeros(ninc,ntrial);
    
    %Start iterations of PMF
    for trial = 1:ntrial %%Middle Loop: performs repeated trials on selected subset
        for inc = 0:ninc-1 %%Inner Loop: performs the PMF algorithm on the repeated subset at a specific hidden percentage

            hidden = (inc/ninc); %Randomly select values to hide
            Hidden_Ind = rand(n);
            Hidden_Ind = triu(Hidden_Ind,1) + triu(Hidden_Ind)'; %Make Hidden_Ind Symetric
            Hidden_Ind(Hidden_Ind < hidden)=0;
            Hidden_Ind(Hidden_Ind >= hidden)=1;

            legend_str = ['Cellline ' num2str(cellline,'%02d')];
            input = normM .* Hidden_Ind; %If the index is hidden, set coresponding value in input to zero
            PMFMat = zeros(n^2,3);

            count = 1;
            %Create N^2 x 3 matrix for input
            for i = 1:n
                for j = 1:n
                    if M(i,j) ~= 0 && Hidden_Ind(i,j) == 1
                        PMFMat(count, 1) = i;
                        PMFMat(count, 2) = j;
                        PMFMat(count, 3) = input(i,j);
                        count = count + 1;
                    end
                end
            end
            PMFMat = PMFMat(1:count-1,:);
            [D,T] = pmfnest(PMFMat,n,n,n,epsilon,lambda,momentum,5);
            B = D*T';
            V = var(A(M ~= 0),1,'all');
            Err_all = immse(A(M ~= 0),B(M ~= 0))/V; %%Error ignoring hidden values
            Err_known = 0;
            Err_hidden = 0;
            for i=1:n
                for j=1:n
                    if M(i,j) ~= 0
                        Err_known = Err_known + (((input(i,j)-B(i,j)).^2)*Hidden_Ind(i,j));%Known Value only error
                        Err_hidden = Err_hidden + (((normM(i,j)-B(i,j)).^2)*(1-Hidden_Ind(i,j)));%Hidden Value only error
                    end
                end
            end
            %%Save data from trial
            Err_known = (Err_known/(nnz(Hidden_Ind(M ~= 0)==1)*V));
            Err_hidden = (Err_hidden/(nnz(Hidden_Ind(M ~= 0)==0)*V));
            output_known(inc+1,trial) = Err_known;
            output_hidden(inc+1,trial) = Err_hidden;
            output_all(inc+1,trial) = Err_all;
        end
    end 
    err_known(cellline, :) = mean(output_known,2);
    errbar_known(cellline, :) = std(output_known, [], 2);
    
    err_hidden(cellline, :) = mean(output_hidden,2);
    errbar_hidden(cellline, :) = std(output_hidden, [], 2);
    
    err_all(cellline, :) = mean(output_all,2);
    errbar_all(cellline, :) = std(output_all, [], 2);

%     if cellline == 1
%         str = legend_str;
%     else
%         str = [str, legend_str];
%     end
end
%%
subplot(2,2,1), hold on
xx = 0:1/ninc:1-1/ninc;
shadedErrorBar(xx,mean(err_known),mean(errbar_known)); %Shaded error bar function from Matlab forums
xlabel('Fraction Hidden');
ylabel('MSE');
title('Error Across Known Values');
lgd.FontSize = 18;
hold on

subplot(2,2,2)
shadedErrorBar(xx,mean(err_hidden),mean(errbar_hidden));
xlabel('Fraction Hidden');
ylabel('MSE');
title('Error Across Hidden Values');
lgd.FontSize = 18;
hold on

subplot(2,2,[3,4])
shadedErrorBar(xx,mean(err_all),mean(errbar_all));
hold on
xlabel('Fraction Hidden');
title('Error Across All Values');
ylabel('MSE');
lgd.FontSize = 18;
% legend(str{:});