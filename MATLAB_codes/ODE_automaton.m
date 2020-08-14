%Hybrid ODE & discrete automaton model

%% ODE Model
%ODE- from Novak-Tyson (1993) model of Xenopus Cell Cycle
%Described in Sible & Tyson (2007)

%%% Constants 
global k1 k3
k1 = 1 ;
k3 = 0.005 ;

global ka Ka kb Kb kc Kc kd Kd
global ke Ke kf Kf kg Kg kh Kh
ka = 0.02 ;
Ka = 0.1 ;
kb = 0.1 ;
Kb = 1 ;
kc = 0.13 ;
Kc = 0.01 ;
kd = 0.13 ;
Kd = 1 ;
ke = 0.02 ;
Ke = 1 ;
kf = 0.1 ;
Kf = 1 ;
kg = 0.02 ;
Kg = 0.01 ;
kh = 0.15 ;
Kh = 0.01 ;

global v2_1 v2_2 v25_1 v25_2 vwee_1 vwee_2
v2_1 = 0.005 ;
v2_2 = 0.25 ;
% % To increase bistability range
% v25_1 = 0.017 ;
% v25_2 = 0.17 ;
v25_1 = 0.5*0.017 ;
v25_2 = 0.5*0.17 ;
vwee_1 = 0.01 ;
vwee_2 = 1 ;

global CDK_total cdc25_total wee1_total IE_total APC_total PPase
CDK_total = 100 ; %100nm is found in typical frog oocyte
% % To increase bistability range
% cdc25_total = 1 ;
cdc25_total = 5 ;
wee1_total = 1 ;
IE_total = 1 ;
APC_total = 1 ;
PPase = 1 ;

%%% Simulation time - minutes
tlast = 1840; % 1 day= 1440min + 400min 
%(first 400 minutes excluded for automaton-
%to avoid affect of initial conditions
%%% Initial conditions 
cyclin = 0 ;
MPF = 0 ;
preMPF = 0 ;
cdc25P = 0 ;
wee1P = 0 ;
IEP = 1 ;
APC = 1 ;

statevar_i = [cyclin,MPF,preMPF,cdc25P,wee1P,IEP,APC] ;

%%% Run ODE

[time,statevars] = ode15s(@dydt_NovakTyson,0:1:tlast,statevar_i) ;

cyclin = statevars(:,1) ;
MPF = statevars(:,2) ;
preMPF = statevars(:,3) ;
cdc25P = statevars(:,4) ;
wee1P = statevars(:,5) ;
IEP = statevars(:,6) ;
APC = statevars(:,7) ;

cyclin_tot = cyclin + MPF + preMPF ;

%%% Results -Implementing only stable limit cycles
% This ignores initial condition-dependent effects
dices = find(time > 400) ;
time = time(dices) - time(dices(1)) ;
cyclin_tot = cyclin_tot(dices) ;
MPF = MPF(dices) ;
APC = APC(dices);

%%% Figures 
%Plotted for 1440 minutes, time can be changed for clarity in the figure
figure(1)
hold on
plot(time,cyclin_tot,'k') 
plot(time,MPF,'r') 
%plot(time,preMPF,'r')
%plot(time,cdc25P,'g')
%plot(time,wee1P,'m')
%plot(time,IEP,'c')
plot(time,APC,'c')
set(gca,'TickDir','Out')

figure(2)
hold on
plot(MPF,cyclin_tot,'g','LineWidth',2.25)
xlabel('{MPF}')
ylabel('[Cyclin_t_o_t_a_l]')
set(gca,'TickDir','Out')

%% Sort concentrations into cell phases

%G1 phase= cyclin_tot <10 & MPF <2.5
%S phase= cyclin_tot >10 & MPF <7
%G2 phase= cyclin_tot >16 & 7<= MPF <10
%M phase = MPF > 2.5

G1MPF = find(MPF <2.5);
G1cyclin_tot = find(cyclin_tot<10);
G1= intersect(G1MPF, G1cyclin_tot);

SMPF = find(MPF< 7);
Scyclin_tot = find(cyclin_tot>=10);
S= intersect(SMPF, Scyclin_tot);
remove=  find(cyclin_tot./MPF <2.5);
remove= intersect(S, remove);
for i= 1: length(remove)
    idx = find(S==remove(i));
    S(idx)=[];
end
    
G2MPF = find(MPF >= 7 & MPF < 10);
G2cyclin_tot = find(cyclin_tot>=16);
G2= intersect(G2MPF, G2cyclin_tot);
remove=  find(cyclin_tot./MPF <2);
remove= intersect(G2, remove);
for i= 1: length(remove)
    idx = find(G2==remove(i));
    G2(idx)=[];
end
    
MMPF = find(MPF >=2.5);
M= MMPF;
for i= 1: length(G1)
    idx = find(M==G1(i));
    M(idx)=[];
end
for i= 1: length(G2)
    idx = find(M==G2(i));
    M(idx)=[];
end
for i= 1: length(S)
    idx = find(M==S(i));
    M(idx)=[];
end

extra=  find((cyclin_tot./MPF >2.5));
inS= intersect(M, extra);
for i= 1: length(inS)
    idx = find(M==inS(i));
    M(idx)=[];
end
S= [S; inS];
S= sort(S);

figure(3)
hold on
sz = 25;
scatter(MPF(M),cyclin_tot(M),sz,'b','filled')
scatter(MPF(G1),cyclin_tot(G1),sz,'k','filled')
scatter(MPF(S),cyclin_tot(S),sz,'r','filled')
scatter(MPF(G2),cyclin_tot(G2),sz,'g','filled')
legend('M phase','G1 phase', 'S phase', 'G2 phase')
xlabel('{MPF}')
ylabel('[Cyclin_t_o_t_a_l]')

%%
%Discrete cell automaton model
max_CellNo= 60000;       % total number of cells tracked
initial_CellNo= 10000;      % initial number of cells

logistic=1;         % logistic regulation of death (0=no, 1=yes)
ks=0.01;           % parameter used in logistic regulation of pdeath
no_cellss=30000;    % number of cells at steady state (used in logistic regulation of death)

initial_phase = 'equal'; % initial condition ('G1' or 'equal') 
%equal--> cells are randomly distributed with equal probability in the 5 phases

prob_death=0.018;

total_time = length(time);

CellDataMatrix =zeros(max_CellNo, 2);
%%%CellDataMatrix is the working matrix 
% every row= one cell
% col1= state of cell (0=no cell, 1=G1, 2=S, 3=G2, 4=M )
% col2 = cell marked for death (0= no, 1=yes)

Output_matrix =zeros(total_time,7);
%%% Output_matrix contains the output data
% each columns = one state ( col1= time, col2=G1, col3=S, col4=G2, col5=M, col6= No.dead, col7= No.marked_for_death)
% each row = number of cells in the corresponding state for each time step

fprintf('Initialisation: ')

if strcmp(initial_phase,'G1')        % cells are initially in G1
  CellDataMatrix(1:initial_CellNo,1)=1*ones(1,initial_CellNo); 
elseif strcmp(initial_phase,'equal') % cells are randomly distributed in the 4 phases
  CellDataMatrix(1:initial_CellNo, 1)=ceil(4*rand(1,initial_CellNo));
end

%%% Probability to die (/min)
probab_death= prob_death+(logistic==1)*ks*(initial_CellNo/no_cellss-1);   

%%% Initialization of time and counting variables
total_cells= initial_CellNo; % number of cells (alive + dead)
imax= initial_CellNo;     %initially set to living cells
inew=0;                   %row for new cell

No_dead=0;              % number of dead cells
No_marked=0;            % number of cells marked for death

fprintf('Complete \n\n')
fprintf('Run simulation: ')

for t= 1: total_time  % loop on time
    
    for i=1:imax   % loop on cells

        if(CellDataMatrix(i,1)>0)   % if cell exists and is alive
            if(rand(1)<probab_death)
                CellDataMatrix(i,2)=1;  % cell marked for death (random)
            end
        
            if (CellDataMatrix(i,1)==1) %cell in G1 
                if ismember(cyclin_tot(t), cyclin_tot(S)) && ismember(MPF(t), MPF(S)) % cell at end of G1
                    if (CellDataMatrix(i,2)==1) % if cell marked for death
                        No_dead=No_dead+1;   % update number of dead cells
                        CellDataMatrix(i,1)=0;  % cell dead
                    else
                        CellDataMatrix(i,1)=2; %Proceed to S phase G1->S 
                    end
                end
            
            elseif (CellDataMatrix(i,1)==2) %cell in S phase
                if ismember(cyclin_tot(t), cyclin_tot(G2)) && ismember(MPF(t), MPF(G2)) % cell at end of S
                    CellDataMatrix(i,1)=3; %Proceed to G2 phase S->G2 
                end
            
            elseif (CellDataMatrix(i,1)==3) %cell in G2 phase
                if ismember(cyclin_tot(t), cyclin_tot(M)) && ismember(MPF(t), MPF(M))% cell at end of G2
                    if (CellDataMatrix(i,2)==1) % if cell marked for death
                        No_dead=No_dead+1;   % update number of dead cells
                        CellDataMatrix(i,1)=0;  % cell dead
                    else
                       CellDataMatrix(i,1)=4; %Proceed to M phase G2->M
                    end
                end
                
            elseif (CellDataMatrix(i,1)==4)  % cell in M phase(division)
                CellDataMatrix(i,1)=1;             % back to G1
                CellDataMatrix(i,2)=CellDataMatrix(i,2); % marked for death if parent marked for death
                
                k=find(CellDataMatrix(:,1)==0);  % look for a free place in matrix data_matrix
                if k~=0
                    inew=k(1);
                else
                    break;
                end
                                
                CellDataMatrix(inew,1)=1;          % new cell starts in G1
                CellDataMatrix(inew,2)=CellDataMatrix(i,2);% marked for death if parent marked for death

            end 
            
        end    % end "if cell alive"

        imax=max(inew,imax);

    end   % end "loop on cells"

    L=CellDataMatrix(:,1);      % extract 1st col of data_matrix (phases of each cell)
    NG1=length(find(L==1));     % count number of cells in G1
    NS=length(find(L==2));      % count number of cells in S
    NG2=length(find(L==3));     % count number of cells in G2
    NM=length(find(L==4));      % count number of cells in M
    NLiving=NG1+NS+NG2+NM;            % total number of living cells
    LD=CellDataMatrix(:,2);              % extract lines of data_matrix (mark for death)
    No_marked=length(find(L>0 & LD==1));     % count number of cells marked for death
        
    Output_matrix(t,:)= [t NG1 NS NG2 NM No_dead No_marked];
       
    if (imax>= max_CellNo)    % Number of cells exceeds ncellmax
       fprintf('\n\nSTOP! Maximum number of cells has exceeded! \n')
       break;
    elseif (NG1+NS+NG2+NM == 0) % No more living cells
        fprintf('\n\nSTOP! There are no more living cells! \n')
        fprintf('tdeath = %g min\n',t)
        break;
    end
    
end     % end loop on time

fprintf('\n\n');

NTOT=NG1+NS+NG2+NM;  %After whole simulation
fprintf('At t=%g (min): \n\n',t);

fprintf('  Number of living cells: %g \n',NTOT);
fprintf('  Number of dead cells: %g \n',No_dead);
fprintf('  Number of cells condemned to death: %g \n\n',No_marked);

fprintf('  Number of cells in phase G1: %g (%g %%) \n',NG1,NG1*100/NTOT);
fprintf('  Number of cells in phase S: %g (%g %%) \n',NS,NS*100/NTOT);
fprintf('  Number of cells in phase G2: %g (%g %%) \n',NG2,NG2*100/NTOT);
fprintf('  Number of cells in phase M: %g (%g %%) \n\n',NM,NM*100/NTOT);

figure(4)  %No of cells in individual phases across time
plot(Output_matrix(:,1), Output_matrix(:,2:5),'LineWidth',1.25);   
xlabel('Time (min)','fontsize',18);
ylabel('Number of cells','fontsize',18);
legend('G1','S','G2','M');
title('No of living cells in individual phases across time')

figure(5)  %Total no. of living cells across time
No_alive =Output_matrix(:,2)+Output_matrix(:,3)+Output_matrix(:,4)+Output_matrix(:,5);
plot(Output_matrix(:,1), No_alive(:,1),'LineWidth',1.25);   
xlabel('Time (min)','fontsize',18);
ylabel('Total number of living cells','fontsize',18);
title('Total no. of living cells across time')

