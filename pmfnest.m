function [D,T]=pmfnest(train_set,num_d,num_t,rank,epsilon,lambda,momentum,bsize)
%input training dataset and num_d total number of the drug in whole dataset and
%num_t total number of targets in whole dataset
%batch is the number of groups we wanna divide our training set

batch_item = bsize;%choose the number of batch we want to use, default was 5.
data = train_set;
[train_num,~] = size(data);%num of item in this dataset
remain = rem(train_num,batch_item);
train_num = train_num-remain;

D = 0.1*randn(num_d,rank); % initial drug matrix
T = 0.1*randn(num_t,rank); % initial target matrix

%step=0.005;%step size 0.005

Drug_incr = zeros(num_d,rank);% increase matrix to store the drug derivative
Target_incr = zeros(num_t,rank);%increase matrix to store the target derivative

% construct a batch dictionary
batch_dict = zeros(batch_item,2);
for uuu=1:batch_item-1
batch_dict(uuu,:)=[(uuu-1)*round((train_num/batch_item))+1 , uuu*round((train_num/batch_item))];
end
batch_dict(batch_item,:)=[(batch_item-1)*round((train_num/batch_item))+1,train_num+remain];

bat_num_use=batch_dict(1,2); % average number of items used in each batch
train_num=train_num+remain;% set the actual number back
for epoch=1:200 %number of iterations, default is 150.
rr = randperm(train_num); % randomly permute the train dataset
data_t = data(rr,:); % randomly updated order training dataset after each epoch
clear rr
for uu=1:batch_item
    Drug_d = zeros(num_d,rank);%drug derivative matrix
    Target_d = zeros(num_t,rank);%target derivative matrix


    Nesterov_D = D - momentum*Drug_incr; %Nesterov gradients
    Nesterov_T = T - momentum*Target_incr;
    for i=batch_dict(uu,1):batch_dict(uu,2) %items in each batch
        D_id = data_t(i,1);
        T_id = data_t(i,2);

        R = data_t(i,3);%real weight value drugbank default value is always 1
        R_hat = sum(Nesterov_D(D_id,:).*(Nesterov_T(T_id,:))); %this the dot plot of Drug vector and target vector to estimate weight value

        %======compute the gradient======%
        ii = repmat(2*(R_hat-R),1,rank);
        E_d = ii.*Nesterov_T(T_id,:)+lambda*Nesterov_D(D_id,:);%partial differentiation of Drug
        E_t = ii.*Nesterov_D(D_id,:)+lambda*Nesterov_T(T_id,:);%partial differentiation of Target

        %compute the updata
        Drug_d(D_id,:) = Drug_d(D_id,:)+E_d;
        Target_d(T_id,:) = Target_d(T_id,:)+E_t;
    end

    %update
    Drug_incr = momentum*Drug_incr+epsilon*Drug_d/bat_num_use;
    Target_incr = momentum*Target_incr+epsilon*Target_d/bat_num_use;
    D = Nesterov_D-Drug_incr;%update Drug matrix
    T = Nesterov_T-Target_incr;%update Target matrix
end
end
end