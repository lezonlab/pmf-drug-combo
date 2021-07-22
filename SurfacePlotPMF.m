s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

csv = readmatrix('C:\Users\ronna\OneDrive\Documents\MATLAB\ComboDrugMat\786-0Mat'); %%Read in from csv file
csv = csv(:,all(~isnan(csv))); 
M = triu(csv) + triu(csv)'; %% Make data symmetric
M(isnan(M))=0;
%%
meanM = mean(mean(M(M ~= 0)));
stdM = std(M(M ~= 0),1,'all');
A = M - mean(mean(M(M ~= 0)));
A = A ./ std(M(M ~= 0),1,'all');
normM = A; %%Normalize Values
%%

%%Set size and output matrices
n = 104;
nsub = 1; % 
ntrial = 5; % Number of independent PMF runs on each hidden increment (gives error bars)
ninc = 20; % Number of increments for fraction hidden data
  
epsilon = 2.170330e+00;
lambda = 7.061362e-03;
momentum = 2.284520e-01;

%%
%%Create grids for 3d plot
hidden_inc = [0:1/ninc:1-1/ninc];
eta_inc = [0:2:100];
[hidden_grid,eta_grid] = meshgrid(hidden_inc,eta_inc);
output_all = zeros(size(hidden_grid));
%%
figure();
legend_str = '';
eta = 0;

for I = 1:size(hidden_grid,1) %%Loop over all combinations on the grid of the surface plot
    I
    for J = 1:size(hidden_grid,2)
        AUC_mean = 0;
        for trial = 1:ntrial %%Middle Loop: performs repeated trials on selected subset
            hidden = hidden_grid(I,J); %Randomly select values to hide
            Hidden_Ind = rand(n);
            Hidden_Ind = triu(Hidden_Ind,1) + triu(Hidden_Ind)'; %Make Hidden_Ind Symetric
            Hidden_Ind(Hidden_Ind < hidden)=0;
            Hidden_Ind(Hidden_Ind >= hidden)=1;
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
            %Run PMF on hidden subset
            [D,T] = pmfnest(PMFMat,n,n,n,epsilon,lambda,momentum,5);
            B = D*T';
            %Calculate AUROC
            orig = A(triu(Hidden_Ind,1)==0 & M ~= 0);
            pred = B(triu(Hidden_Ind,1)==0 & M ~= 0);
            %Transform cutoff
            eta_trans = (eta_grid(I,J)-meanM)/stdM;
            decision_values = zeros(length(orig),1);
            decision_values(orig(:)>eta_trans) = 1;
            [X,Y,TT,AUC] = perfcurve(decision_values, pred, 1);
            AUC_mean = AUC_mean + AUC;
        end
        AUC_mean = AUC_mean/ntrial;
        output_all(I,J) = AUC;
    end
end 
%%
%figure();
surf(hidden_grid, 100 - eta_grid, output_all); %eta is percent growth, minus 100
xlabel('Percetage Hidden');
ylabel('Efficacy cutoff');
zlabel('AUC');
colorbar
