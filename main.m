clc
clear all

%% Input raw data
data=xlsread('D:\Desktop\D1\D1_datasets\NNGC1_dataset_D1_V1_003');%%%原始数值

%% I. Time series decomposition
[MObje_max MX Alpha_opt MIG_opt y]=l_1_filter(data); %% 1.1 Obtain an optimal decomposition

PPL=[];%% 1.2  Trend Screening
for ii=1:size(data,2)
z1=MX(:,ii);
test1=size(data,1);
i=1;l1=[];
while i+1<test1&&i+j<=test1
    key=1;
    k=z1(i+1)-z1(i);
    b=z1(i)-k;
    j=2;
    while key==1&&i+j<=test1
        if abs(z1(i+j)-(k*(j+1)+b))<1e-2
            j=j+1;
        else
            l1=[l1;i+j-1];
            key=2;
            i=i+j;
        end
    end
end
%Calculate data - regain trend data
l1=[1;l1;test1];
for i=1:size(l1,1)-1
    T=l1(i+1)-l1(i)+1;
    k=(z1(l1(i+1))-z1(l1(i)))/(T-1);
    b=z1(l1(i))-k;
    delta=norm(z1(l1(i):l1(i+1),1)-y(l1(i):l1(i+1),1),2)/sqrt(T);
    PL1(i,:)=[k b  delta T];
end
PPL=[PPL mat2cell(PL1,size(PL1,1),size(PL1,2))];
end

%% II. Scale equalization                                                                                                                                                                二、周期粒度均衡
T0=10; %%Length of fixed time window
T=T0:T0:length(data); 
[GN_stand MIG_Max_cellj0]= equalization(data, PPL, MX, T, T0);%%Determine granularity criteria

%% III. Four BPNNs for long-term forecasting

%% 3.1 Construct 4 corresponding datasets
Data_j17=[];
for j17=1:4
    Data_j16=[];
    for j16=1:size(MIG_Max_cellj0,1) 
        Data_j15=[];
        for j15=1:size(MIG_Max_cellj0,2) 
            data_j15=MIG_Max_cellj0{j16,j15}(:,j17);
            Data_j15=[Data_j15 data_j15'];
        end
        Data_j16=[Data_j16 Data_j15'];
    end
    Data_j16=Data_j16'; 
    Data_j17=[Data_j17  mat2cell(Data_j16,size(Data_j16,1),size(Data_j16,2))];%% 4 corresponding datasets
end

%% 3.2 BPNN

Data_j18=[];
for j18=1:size(Data_j17,2)

Toal_Data=Data_j17{1,j18};
input=Toal_Data(1:end-1,:);
output=Toal_Data(2:end,GN_stand*(size(data,2)-1)+1:end);

input=input';     
output=output';   

%Training and prediction data
input_train=input(:,1:end-1);
input_test=input(:,end)
output_train=output(:,1:end-1);
output_test=output(:,end);

%% Normalization
[inputn,inputps]=mapminmax(input_train,0,1);
[outputn,outputps]=mapminmax(output_train);
inputn_test=mapminmax('apply',input_test,inputps);

%% Get the number of nodes
inputnum=size(input,1);
outputnum=size(output,1);

MSE=1e+5; 
transform_func={'tansig','purelin'}; 
train_func='trainlm';  
for hiddennum=fix(sqrt(inputnum+outputnum))+1:fix(sqrt(inputnum+outputnum))+10
    
    %Build the network
    net=newff(inputn,outputn,hiddennum,transform_func,train_func);
    % Network parameters
    net.trainParam.epochs=10000;         % Number of trainings
    net.trainParam.lr=0.01;                   % Learning rate
    net.trainParam.goal=1e-6;        % Minimum error of training objectives
    % Network training
    net=train(net,inputn,outputn);
    an0=sim(net,inputn);  %Simulation results
    mse0=mse(outputn,an0);  %The mean square error of the simulation
     
    %Update the best implicit layer node
    if mse0<MSE  %%It's guaranteed not to "overfit."
        MSE=mse0;
        hiddennum_best=hiddennum;
    end
end

%% Building BP Neural Networks with Optimal Hidden Layer Nodes
net=newff(inputn,outputn,hiddennum_best,transform_func,train_func);

% 
net.trainParam.epochs=10000;         
net.trainParam.lr=0.01;                 
net.trainParam.goal=1e-6;         

net=train(net,inputn,outputn);
  
an=sim(net,inputn_test); 
test_simu=mapminmax('reverse',an,outputps);  

Data_j18=[Data_j18 test_simu'];
end
Data_pre=reshape(Data_j18, [GN_stand,4])

for j19=1:GN_stand
if Data_pre(j19,4)<=0
    Data_pre(j19,4)=0;
else 
    Data_pre(j19,4)=fix(Data_pre(j19,4));
end
end

X_pre=[];
for j20=1:GN_stand 
x_pre=Data_pre(j20,1).*[1:Data_pre(j20,4)]+Data_pre(j20,2);
X_pre=[X_pre x_pre];
end

PP_network=(X_pre(1,[1:T0]))';
A_value=data((size(T,2)-1)*T0+1:end,end);

K=10;
P_network=PP_network(1:K);
Actual_value=A_value(1:K);

Error=P_network-Actual_value
[P_network Actual_value]

%Evaluation indicators
RMSE_end=sqrt(mse(P_network-Actual_value)); 
MAE_end=mean(abs(Actual_value-P_network)); 
MAPE_end=mean(abs((Actual_value-P_network)./Actual_value))*100; 

att_end=[RMSE_end MAE_end MAPE_end]
