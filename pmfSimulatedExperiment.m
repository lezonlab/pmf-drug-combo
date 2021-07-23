s = RandStream('mt19937ar','Seed',41); %Try setting an actual seed, not shuffle
RandStream.setGlobalStream(s);

csv = readmatrix('C:\Users\ronna\OneDrive\Documents\MATLAB\ComboDrugMat\COLO 205Mat'); %%Read in from csv file
csv = csv(:,all(~isnan(csv))); 
M = triu(csv) + triu(csv)'; %% Make data symmetric
M(isnan(M))=0;
%%
%%Standardize
meanM = mean(mean(M(M ~= 0)));
stdM = std(M(M ~= 0),1,'all');
A = M - mean(mean(M(M ~= 0)));
A = A ./ std(M(M ~= 0),1,'all');
normM = A; %%Normalize Values

%%
epsilon = 2.170330e+00;
lambda = 7.061362e-03;
momentum = 2.284520e-01;
%%
hidden = .9; %Randomly select values to hide for starting subset
n = 104;
ninc = 20;
ntrial = 18;

res = zeros(1,ntrial + 1);
known_synergy = zeros(1,ntrial + 1);
percentage_hidden = zeros(1,ntrial + 1);
%%
%Pick initial values for experiment
Hidden_Ind = rand(n);
Hidden_Ind = triu(Hidden_Ind,1) + triu(Hidden_Ind)'; %Make Hidden_Ind Symetric
Hidden_Ind(Hidden_Ind < hidden)=0;
Hidden_Ind(Hidden_Ind >= hidden)=1;
%%
%Test initial values for inital predictions, then Identify synergystic
%candidates
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
nnz(Hidden_Ind)
PMFMat = PMFMat(1:count-1,:);
[D,T] = pmfnest(PMFMat,n,n,n,epsilon,lambda,momentum,5);
B = D*T';
V = var(A(M ~= 0),1,'all');

orig = A(triu(Hidden_Ind,1)==0);
pred = B(triu(Hidden_Ind,1)==0);

zero_trans = (70-meanM)/stdM;
decision_values = zeros(length(orig),1);
decision_values(orig(:)>zero_trans) = 1;
[X,Y,TT,AUC] = perfcurve(decision_values, pred, 1);
res(1,1) = AUC;
%Perecentage of synergystic drugs in original matrix currently known
known_synergy(1,1) = nnz(A(Hidden_Ind==1)>zero_trans)/nnz(A>zero_trans); 
percentage_hidden(1,1) = count/n^2;
%%


%%
for trial = 2:ntrial+1
    x = B .* ~Hidden_Ind;
    x = triu(x,1);
    x(x == 0) = -Inf;
    [~,idx] = sort(x(:), 'descend');
    total = floor(((n^2)/ninc)/2); %5% of matrix
    if  total > nnz(triu(Hidden_Ind==0))
        total = nnz(triu(Hidden_Ind==0));
    end
    to_test = zeros(n);
    to_test(idx(1:total)) = x(idx(1:total));
    to_test = triu(to_test,1) + triu(to_test)';
    to_test = to_test~=0;
    
%     [row,col] = find(triu(Hidden_Ind==0));
%     ix = randperm(length(row));
%     total = 270;
%     to_test = zeros(n);
%     for i = 1:total
%         Ensure we are not adding along the diagonal
%         if row(ix(i)) ~= col(ix(i))
%             to_test(row(ix(i)), col(ix(i))) = 1;
%             to_test(col(ix(i)), row(ix(i))) = 1;
%         end
%     end
    
    %Values have now been tested, and are now known. Add them to
    %Hidden_Ind, and begin procedure once more
    Hidden_Ind = Hidden_Ind + to_test;
    
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
    nnz(Hidden_Ind)
    PMFMat = PMFMat(1:count-1,:);
    AUC_mean = 0;
    %%PMF is stochastic and depends on the initial conditions
    %%We can run it several times and pick the best AUC
    for epoch = 1:10
        AUC_best = Inf;
        [D,T] = pmfnest(PMFMat,n,n,n,epsilon,lambda,momentum,5);
        B = D*T';
        V = var(A(M ~= 0),1,'all');

        orig = A(triu(Hidden_Ind,1)==0);
        pred = B(triu(Hidden_Ind,1)==0);

        zero_trans = (70-meanM)/stdM;
        %zero_trans = .6;
        % zero_trans2 = (-5-minM)/maxM;
        decision_values = zeros(length(orig),1);
        decision_values(orig(:)>zero_trans) = 1;
        [X,Y,TT,AUC] = perfcurve(decision_values, pred, 1);
        if AUC < AUC_best
            AUC_best = AUC;
        end
    end
    res(1,trial) = AUC_best;
    known_synergy(1,trial) = nnz(A(Hidden_Ind==1)>zero_trans)/nnz(A>zero_trans); 
    percentage_hidden(1,trial) = count/n^2;
    %%
    %Add tested values to known indices
    
end

%%
xx = 1-hidden:.05:1;
figure
yyaxis left
plot(xx,res,'--o'); hold on;
ylabel('AUROC');
yyaxis right
plot(xx,xx,'--ok');
plot(xx,known_synergy,'--o');
xlabel('Percetage Hidden');
ylabel('Known Combinations with Efficacy > 70');
legend('PMF AUROC','Random Choice','Erdos-Reyni PMF Guided');
title('Experimental Designs vs Known Synergism');
xlim([.1 1]);

% xx = 1-hidden:.05:1;
% figure
% plot(xx,known_synergy,'--o'); hold on;
% plot(xx,xx,'--ok');
% xlabel('Percetage Hidden');
% ylabel('Known Combinations with Efficacy > 70');
% legend('PMF Guided','Random Choice');
% title('Experimental Designs vs Known Synergism');




%     [row,col] = find(triu(Hidden_Ind==0));
%     ix = randperm(length(row));
%     total = floor(((n^2)/20)/2)
%     if  total > nnz(triu(Hidden_Ind==0))
%         total = nnz(triu(Hidden_Ind==0));
%     end
%     for i = 1:total
%         %Ensure we are not adding along the diagonal
%         if row(ix(i)) ~= col(ix(i))
%             to_test(row(ix(i)), col(ix(i))) = 1;
%             to_test(col(ix(i)), row(ix(i))) = 1;
%         end
%     end
%%
% x = B .* ~Hidden_Ind;
% x = triu(x,1);
% x(x == 0) = -Inf;
% [~,idx] = sort(x(:), 'descend');
% total = floor(((n^2)/ninc)/2); %5% of matrix
% if  total > nnz(triu(Hidden_Ind==0))
%     total = nnz(triu(Hidden_Ind==0));
% end
% to_test = zeros(n);
% to_test(idx(1:total)) = x(idx(1:total));
% to_test = triu(to_test,1) + triu(to_test)';
% to_test = to_test~=0;
%%
