%% FORDBO for CEC2005 benchmark functions
clear all
clc
SearchAgents_no = 30; % Number of search agents

Function_name='F1'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
Max_iteration=1000; % Maximum numbef of iterations
% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);  %Set bounds and optimize functions
%% The experiment was performed 30 times and the results were counted
for i = 1:30
    disp(['Experiment ',num2str(i),' ------']);
%DBOFOR
[Best_pos0(i,:),Best_score0(i),DBO3_curve0(i,:)]=FORDBO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); 
end

%% Comparison of results
figure

semilogy(mean(DBO3_curve0),'Color','r','MarkerFaceColor','r','linewidth',2,'MarkerSize',10,'MarkerIndices',1:100:1000)
hold on

title('\fontname{Times New Roman}Objective Space')
xlabel('\fontname{Times New Roman}Iteration');
ylabel('\fontname{Times New Roman}Average of Fitness Value');

axis tight
%The grid is set only in the X-axis direction
%y = rand(100,1);
%ay = gca;
%ay.YGrid = 'on';
%ay.XGrid = 'off';
grid on
box on
legend('\fontname{Times New Roman}FORDBO')

display(['FORDBO---The best fitness value of 30 experiments (Best) : ', num2str(min(Best_score0))]);
display(['FORDBO---The average fitness value of 30 experiments (mean) : ', num2str(mean(Best_score0))]);
display(['FORDBO---The worst fitness value of 30 experiments (worst) : ', num2str(max(Best_score0))]);
display(['FORDBO---The standard deviations of 30 experiments (std) : ', num2str(std(Best_score0))]);
