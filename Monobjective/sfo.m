%   SUNFLOWER OPTIMIZATION (SFO) ALGORITHM FOR NONLINEAR UNCONSTRAINED
%                             OPTIMIZATION 
%
%
% Copyright (c) 2018, Guilherme Ferreira Gomes
% All rights reserved.
%
%
% Please cite this algorithm as:
%
% Gomes, G. F., da Cunha, S. S., & Ancelotti, A. C. 
% A sunflower optimization (SFO) algorithm applied to damage identification
% on laminated composite plates. Engineering with Computers, p. 1-8, 2018.
% DOI: https://doi.org/10.1007/s00366-018-0620-8


%% Test Function optimization examples:
%
% clear all
% close all
% clc
% format long
% 
% set(0,'DefaultAxesFontName', 'Times New Roman')
% set(0,'DefaultAxesFontSize', 14)
% set(0,'DefaultTextFontname', 'Times New Roman')
% set(0,'DefaultTextFontSize', 26)
% 
% n = 20;           %number of sunflowers
% p = 0.05;         %pollination rate best values 0.01 < p < 0.10
% m = 0.1;          %mortality rate, best values 0.01 > m < 0.10
% s = 1-(p+m);      %survival rate, best values 0.80 > s < 0.90
% d = 2;            %problem dimension
% LB = [-5 -5];     %lower bounds
% UB = [5 5];       %upper bounds
% n_iter = 100;     %max number os iterations/gerations
% 
% %Objective test functions
% Fun1 = @(x)  100*(x(2)-x(1)^2)^2+(1-x(1))^2;                 % rosenbrok (1,1)
% Fun2 = @(x) (x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;             % booth (1,3)
% Fun3 = @(x) (1.5-x(1)*(1-x(2)))^2+(2.25-x(1)*(1-x(2)^2))^2;  % beale's (3,0.5)
% Fun4 = @(x) 0.26*(x(1)^2+x(2)^2)-0.48*x(1)*x(2)              % matyas (0,0)
% 
% Fun = Fun4; %choose the objtective test function 
% 
% [x,fval,iter,state,population]=sfo(Fun,n,p,m,d,LB,UB,n_iter);


%% SunFLower Optimization Function
function [best,fmin,n_iter,state,population]=sfo(Fun,n,p,m,d,LB,UB,n_iter)

% n = number of plants
% p = pollination rate
% m = mortality rate
% d = problem dimension
% LB = lower bounds
% UB = upper bounds
% n_iter = max number os gerations

for i=1:n
    Plants(i,:)=LB+(UB-LB).*rand(1,d);
    Fitness(i)=Fun(Plants(i,:));
end

[fmin,I]=min(Fitness);
best=Plants(I,:);
S=Plants;

for t=1:n_iter
    for i=1:n

        % pollination 
        for j=1:(round(p*n))
        S(j,:) = (Plants(j,:)-Plants(j+1,:))*rand(1) + Plants(j+1,:);
        end
        
        % steps
        for j=(round(p*n)+1):(round(n*(1-m)))
        S(j,:)=Plants(j,:)+rand*((best-Plants(j,:))/(norm((best-Plants(j,:)))));
        end
        
        % mortality of m% plants
        for j = ((round(n*(1-m)))+1):n
        S(j,:)= (UB-LB)*rand+LB;
        end
        
        for j = ((round(n*(1-m)))+1):n
            for k=1:length(LB)  
            S(j,k)= (UB(k)-LB(k))*rand+LB(k);
            end
        end
        S(i,:)=bound_check(S(i,:),LB,UB);
         
        Fnew = Fun(S(i,:));
        if (Fnew <= Fitness(i))
            Plants(i,:) = S(i,:);
            Fitness(i) = Fnew;
        end
        if  Fnew <= fmin
            best = S(i,:);
            fmin = Fnew;
        end
    end
        state(t,:) = [n_iter best fmin];
        population{t} = S;
        
        figure(1)
        set(gcf,'color','w');
        plot(best(1,1),best(1,2),'ro',S(2:end,1),S(2:end,2),'b*','markers',12)
        axis([LB(1) UB(1) LB(2) UB(2)]);
        xlabel('x_1')
        ylabel('x_2')
        axis square
        legend('best','plants')
        
        
        %Use this pause fuction to view the sfo graphical interations 
        %pause(0.01)
        
        X = sprintf('Iteration %0.f  x1 = %.4f x2 = %.4f function=%.4f\n',t,best(1),best(2),fmin);
        disp(X);

    if round(t/100) == t/100
        best
        fmin
    end
end
disp(['Total number of function evaluation:',num2str(n_iter*n)]);
disp(['Best solution (Sun) =',num2str(best), 'Best Fitness',num2str(fmin)]);

function s = bound_check(s,LB,UB)
ns_tmp=s;
I=ns_tmp<LB;
ns_tmp(I)=LB(I);

J=ns_tmp>UB;
ns_tmp(J)=UB(J);
s=ns_tmp;
