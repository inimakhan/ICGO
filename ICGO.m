
%__________________________________________________________________     %
%        MOCGO: Multi-objective Chaos Game Optimization (IOCGO)         %
%                                                                       %
%                                                                       %
%                  Developed in MATLAB R2023a (MacOs)                   %
%                                                                       %
%                      Author and programmer                            %
%                ---------------------------------                      %
%                Nima Khodadadi (ʘ‿ʘ)   University of Miami             %
%                         SeyedAli Mirjalili                            %
%                             e-Mail                                    %
%                ---------------------------------                      %
%                      Nima.khodadadi@miami.edu                         %
%                                                                       %
%                                                                       %
%                            Homepage                                   %
%                ---------------------------------                      %
%                    https://nimakhodadadi.com                          %
%                                                                       %
%  % I acknowledge that this version of ICGO has been written using     %
% a large portion of the following code:Chaos Game Optimization (CGO)by %   
% Siamak Talatahari                                                     %
%                                                                       %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------------------------------------------------------- %
 

clc;clear;close all
%% Get Required Problem Information
ObjFuncName = @(x) Sphere(x); % @YourCostFunction ;
Var_Number = 100 ; % Number of your variables ;
LB = -100 *ones(1,Var_Number) ;  % Lower bound of variable ;
UB = 100 *ones(1,Var_Number) ;  % Upper bound of variable ;
%% Get Required Algorithm Parameters
MaxIter = 100 ; % Maximum number of generations ;
Seed_Number = 25 ; % Maximum number of initial eligible points ;

%% Initialization
for i=1:Seed_Number
    Seed(i,:)=unifrnd(LB,UB);
    Fun_eval(i,:)=feval(ObjFuncName,Seed(i,:));
end
tic
%% Main Loop of The CGO
for Iter=1:MaxIter
    for i=1:Seed_Number
        % Update the best Seed
        [~,idbest]=min(Fun_eval);
        BestSeed=Seed(idbest,:);
        %% Generate New Solutions
        % Random Numbers
        I=randi([1,2],1,12); % Beta and Gamma
        Ir=randi([0,1],1,5);
        % Random Groups
        RandGroupNumber=randperm(Seed_Number,1);
        RandGroup=randperm(Seed_Number,RandGroupNumber);
        % Mean of Random Group
        MeanGroup=mean(Seed(RandGroup,:)).*(length(RandGroup)~=1)...
            +Seed(RandGroup(1,1),:)*(length(RandGroup)==1);
        % New Seeds
        Alfa(1,:)=rand(1,Var_Number);
        Alfa(2,:)= 2*rand(1,Var_Number)-1;
        Alfa(3,:)= (Ir(1)*rand(1,Var_Number)+1);
        Alfa(4,:)= (Ir(2)*rand(1,Var_Number)+(~Ir(2)));
        ii=randi([1,4],1,3);
        SelectedAlfa=Alfa(ii,:);

        NewSeed(1,:)=Seed(i,:)+SelectedAlfa(1,:).*(I(1)*BestSeed-I(2)*MeanGroup);
        NewSeed(2,:)=BestSeed+SelectedAlfa(2,:).*(I(3)*MeanGroup-I(4)*Seed(i,:));
        NewSeed(3,:)=MeanGroup+SelectedAlfa(3,:).*(I(5)*BestSeed-I(6)*Seed(i,:));
        NewSeed(4,:)=WAVELETMUTAION(BestSeed,LB,UB,Iter,MaxIter);
        for j=1:4
            % Checking/Updating the boundary limits for Seeds
            NewSeed(j,:)=min(max(NewSeed(j,:),LB),UB);
            % Evaluating New Solutions
            Fun_evalNew(j,:)=feval(ObjFuncName, NewSeed(j,:));
        end
        Seed=[Seed; NewSeed];
        Fun_eval=[Fun_eval; Fun_evalNew];

    end
    [Fun_eval, SortOrder]=sort(Fun_eval);
    Seed=Seed(SortOrder,:);
    [BestFitness,idbest]=min(Fun_eval);
    BestSeed=Seed(idbest,:);
    Seed=Seed(1:Seed_Number,:);
    Fun_eval=Fun_eval(1:Seed_Number,:);
    % Store Best Cost Ever Found
    Conv_History(Iter)=BestFitness;
    % Show Iteration Information
    disp(['Iteration ' num2str(Iter) ': Best Cost = ' num2str(Conv_History(Iter))]);
end
toc

%% Objective Function
function z=Sphere(x)
z=sum(x.^2);
end

function XWaveletMutaion=WAVELETMUTAION(Position,LB,UB,It,MaxIt)
s=randi(MaxIt);
alpha=s.*(1/s).^(1-It/MaxIt);
Indx=rand;
if Indx<0.5
    Sigma=exp(-(UB-Position)/alpha).^2/2.*cos(5*(UB-Position)/alpha);
    XWaveletMutaion=Position+Sigma;
elseif Indx>0.5
    Sigma=exp(-(Position-LB)/alpha).^2/2.*cos(5*Position-LB/alpha);
    XWaveletMutaion=Position+Sigma;
end
end