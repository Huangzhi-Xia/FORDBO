%% FORDBO
%_________________________________________________________________________%

function [Best_pos,Best_score,curve] = FORDBO(pop,MaxIter,lb,ub,dim,fobj)
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end
PballRolling = 6/30; % Proportion of ball-rolling dung beetles
PbroodBall = 6/30; % Proportion of brood balls
PSmall = 7/30; % Proportion of small dung beetles
Pthief = 11/30; % Proportion of thief dung beetles
BallRollingNum = pop*PballRolling; % Number of ball-rolling dung beetles
BroodBallNum = pop*PbroodBall; % Number of brood balls
SmallNum = pop*PSmall; % Number of small dung beetles
ThiefNum = pop*Pthief; % Number of thief dung beetles

X = initializationNewJ(pop,dim,ub,lb); %% Population initialization using good nodes set method
fitness = zeros([1,pop]);
ObjectiveFunction = fobj;
for i =1:pop
    fitness(i)=fobj(X(i,:));
end
%The population optimal value is recorded
[GbestScore,minIndex] = min(fitness);
GbestPos = X(minIndex,:);
Xl = X; %For recording X(t-1)
%Current generation population
cX = X;
cFit=fitness;
%fmin=min(fitness);
%fmax=max(fitness);

%S=zeros(1,pop);
for t=1:MaxIter
    %% Ball-rolling dung beetle
    %Population worst value
    [fmax,maxIndex]=max(fitness);
    Worst=X(maxIndex,:); 
    r2=rand;
    
    %% Reduction factor
        a=1.6; 
        p=rand;
        if p>0 %
            w=exp(-(0.3)/((0.1)^a)*(t)/MaxIter);
        else
            w=1;
        end
    for i=1:BallRollingNum
        if r2<0.9
            if rand>0.5
                alpha=1;
            else
                alpha=-1;
            end
            b=0.3;
            k=0.1;
            
            X(i,:)=w.*X(i,:)+b.*abs(cX(i,:)-Worst)+alpha.*k*(Xl(i,:));% No update occurs when theta is equal to these three degrees
            
        else
            theta=randi(180);
            if(theta == 0 || theta == 90 || theta==180)
                X(i,:)=cX(i,:);
            else
                theta = theta*pi/180;
                X(i,:)=w.*cX(i,:)+tan(theta).*abs(cX(i,:)-Xl(i,:));
            end
        end
        for j=1:dim
            if X(i,j)>ub(j)
                X(i,j)=ub(j);
            end
            if X(i,j)<lb(j)
                X(i,j)=lb(j);
            end
        end
        fitness(i)=fobj(X(i,:));
        if fitness(i)<GbestScore
            GbestScore=fitness(i);
            GbestPos=X(i,:);
        end
    end
    %The current iteration is the best
    [~,minIndex]=min(fitness);
    GbestB = X(minIndex,:);
    %% Brood ball
    R=1-t/MaxIter;
    X1=GbestB*(1-R);
    X2=GbestB*(1+R);
    Lb = max(X1,lb);
    Ub = min(X2,ub);

    %% Fractional order calculus optimization %%%%%%%%%%%%%%%%%%%%%%
       Xold{t} = X;%It is used to store historical data%
       %Calculate fractional-order correlation parameters
       dg = sum(sum(sqrt((GbestB-X).^2)));%The average distance from GbestB to other locations
       for i = 1:pop%
       d(i) =  sum(sum(sqrt((X(i,:) - X).^2)));%
       end%
       dmin = min(d);%
       dmax = max(d);%
       f=(dg-dmin)/(dmax-dmin);% Evolutionary factor
       v =1/(2*exp(-0.47*f));%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = BallRollingNum+1:BallRollingNum+BroodBallNum
        b1=rand;
        b2=rand;
         %% Fractional order calculus is introduced to improve the dynamic upper and lower boundaries using the location of previous generations of populations
       if t>3
              
              X1 = cell2mat(Xold(t-1));%Individuals from the first 1 generations
              X2 = cell2mat(Xold(t-2));%Individuals from the first 2 generations
              X3 = cell2mat(Xold(t-3));%Individuals from the first 3 generations
              r1=rand;
              if r1<0.9
                  Lb = v.*Lb-0.5*v.*(v-1).*X1(i,:) + v.*(v-1)*(v-2).*X2(i,:)/6-...
                  v.*(v-1).*(v-2).*(v-3).*X3(i,:)/24 + GbestB/MaxIter; %randn().*ones(1,dim)
                  Ub = v.*Ub-0.5*v.*(v-1).*X1(i,:) + v.*(v-1)*(v-2).*X2(i,:)/6-...
                  v.*(v-1).*(v-2).*(v-3).*X3(i,:)/24 - GbestB/MaxIter ; %
              else
                  Lb = v.*Lb-0.5*v.*(v-1).*X1(i,:) + v.*(v-1)*(v-2).*X2(i,:)/6-...
                  v.*(v-1).*(v-2).*(v-3).*X3(i,:)/24 + randn().*ones(1,dim); %randn().*ones(1,dim)
                  Ub = v.*Ub-0.5*v.*(v-1).*X1(i,:) + v.*(v-1)*(v-2).*X2(i,:)/6-...
                  v.*(v-1).*(v-2).*(v-3).*X3(i,:)/24 - randn().*ones(1,dim) ; %
              end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         
              
       end
        X(i,:) = GbestB+b1*(cX(i,:)-Lb)+b2*(cX(i,:)-Ub);
     
        for j=1:dim
            if X(i,j)>Ub(j)
                X(i,j)=Ub(j);
            end
            if X(i,j)<Lb(j)
                X(i,j)=Lb(j);
            end
        end
        fitness(i)=fobj(X(i,:));
        if fitness(i)<GbestScore
            GbestScore=fitness(i);
            GbestPos=X(i,:);
        end
    end
    %The current iteration is the best
    [~,minIndex]=min(fitness);
    GbestB = X(minIndex,:);
    %% Small dung beetle update
    R=1-t/MaxIter;
    X1=GbestPos*(1-R);
    X2=GbestPos*(1+R);
    Lb = max(X1,lb);
    Ub = min(X2,ub);
    
    %% Fractional order calculus optimization %%%%%%%%%%%%%%%%%%%%%%
       Xold{t} = X; %It is used to store historical data%
       %Calculate fractional-order correlation parameters
       dg = sum(sum(sqrt((GbestPos - X).^2)));%The average distance from GbestPos to other locations
       for i = 1:pop%
       d(i) =  sum(sum(sqrt((X(i,:) - X).^2)));%
       end%
       dmin = min(d);%
       dmax = max(d);%
       f=(dg-dmin)/(dmax-dmin);% Evolutionary factor
       v =1/(2*exp(-0.47*f));%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=BallRollingNum+BroodBallNum+1:BallRollingNum+BroodBallNum+SmallNum
        C1=rand;
        C2=rand;
        
        %% Fractional order calculus is introduced to improve the dynamic upper and lower boundaries using the location of previous generations of populations
       if t>3
              X1 = cell2mat(Xold(t-1));%Individuals from the first 1 generations
              X2 = cell2mat(Xold(t-2));%Individuals from the first 2 generations
              X3 = cell2mat(Xold(t-3));%Individuals from the first 3 generations
              r2=rand;
              if r2<0.9
                  Lb = v.*Lb-0.5*v.*(v-1).*X1(i,:) + v.*(v-1)*(v-2).*X2(i,:)/6-...
                  v.*(v-1).*(v-2).*(v-3).*X3(i,:)/24+GbestPos/MaxIter; %randn().*ones(1,dim)
                  Ub = v.*Ub-0.5*v.*(v-1).*X1(i,:) + v.*(v-1)*(v-2).*X2(i,:)/6-...
                  v.*(v-1).*(v-2).*(v-3).*X3(i,:)/24-GbestPos/MaxIter ; %
              else
                  Lb = v.*Lb-0.5*v.*(v-1).*X1(i,:) + v.*(v-1)*(v-2).*X2(i,:)/6-...
                  v.*(v-1).*(v-2).*(v-3).*X3(i,:)/24 + randn().*ones(1,dim); %randn().*ones(1,dim)
                  Ub = v.*Ub-0.5*v.*(v-1).*X1(i,:) + v.*(v-1)*(v-2).*X2(i,:)/6-...
                  v.*(v-1).*(v-2).*(v-3).*X3(i,:)/24 - randn().*ones(1,dim) ; %
              end
       end
        X(i,:) = GbestPos+C1*(cX(i,:)-Lb)+C2*(cX(i,:)-Ub);

        for j=1:dim
            if X(i,j)>Ub(j)
                X(i,j)=Ub(j);
            end
            if X(i,j)<Lb(j)
                X(i,j)=Lb(j);
            end
        end
        fitness(i)=fobj(X(i,:));
        if fitness(i)<GbestScore
            GbestScore=fitness(i);
            GbestPos=X(i,:);
        end
    end
    %% Thief dung beetle phase
    %% Define pathfinder dung beetle, repetitive renewal mechanism of the pathfinder dung beetle
        global_pop=GbestPos; % GbestPos remains unchanged in the current iteration.
        path_p=global_pop;
        path_old=path_p;
        u=-1+(1-(-1)).*rand(1,dim);
        r3=rand(1,dim);
        A=u.*exp(-2.*t/MaxIter);
        % Pathfinder dung beetle position update
        path_p(1,:) = path_old(1,:) + ...
                ((2).*r3(1,:)).*(path_old(1,:) - GbestB(1,:)) + 1.*A(1,:);
        % Boundary handling to prevent crossing the boundary
        Flag4ub=path_p(1,:)>ub(1,:); % check lower bound for pathfinder
        Flag4lb=path_p(1,:)<lb(1,:); % check upper bound for pathfinder
        path_p(1,:)=(path_p(1,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        if ObjectiveFunction(path_p(1,:))<GbestScore
            GbestScore=ObjectiveFunction(path_p(1,:));
            GbestPos=path_p(1,:);
        end
        
    
    % Calculate the average dung beetle position
    Mean = mean(X);
    
    for i=pop-ThiefNum:pop
        g=randn();
        S=0.5;
        r4=rand;
        if r4<1/2
            %Thief dung beetle update
             X(i,:) = GbestPos + (g.*S.*(abs(cX(i,:)-GbestB)+abs(cX(i,:)-GbestPos)));
        else
             X(i,:) = GbestPos + g.*S.*(Mean - cX(i,:));
        end
        global_pop=GbestPos; 
        path_p=global_pop;
        path_old=path_p(1,:);
        u2=-1+(1-(-1)).*rand(1,dim);
        r3=rand(1,dim);
        A=u2.*exp(-2.*t/MaxIter);
        % Pathfinder dung beetle position update
        path_p(1,:) = path_old(1,:) + ...
                ((2).*r3(1,:)).*(path_old(1,:) - GbestB(1,:)) + 1.*A(1,:);
            
        % Boundary handling to prevent crossing the boundary
        Flag4ub=path_p(1,:)>ub(1,:);                            % check lower bound for pathfinder
        Flag4lb=path_p(1,:)<lb(1,:);                               % check upper bound for pathfinder
        path_p(1,:)=(path_p(1,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        if ObjectiveFunction(path_p(1,:))<GbestScore
            GbestScore=ObjectiveFunction(path_p(1,:));
            GbestPos=path_p(1,:);
        end
        
        for j=1:dim
            if X(i,j)>ub(j)
                X(i,j)=ub(j);
            end
            if X(i,j)<lb(j)
                X(i,j)=lb(j);
            end
        end
        
        
        fitness(i)=fobj(X(i,:));
        if fitness(i)<GbestScore
            GbestScore=fitness(i);
            GbestPos=X(i,:);
        end
    end
    %The population of t generations was recorded
    Xl=cX;
    %Update the current generation population
    for i=1:pop
        if fitness(i)<cFit(i)
            cFit(i)=fitness(i);
            cX(i,:)=X(i,:);
        end
    end
    Best_score = GbestScore;
    Best_pos = GbestPos;
    curve(t)=GbestScore;
end
end

