% An ODE Model for yeast cell-cycle, including asymmetric division
% The parameter g indicates where 'mother' cells have a size advantage over
% new daughter cells
% Includes a calculation of Bud Index
% Includes realistic division and harvesting
% Includes logisitic growth with defined Carrying Capacity

% Figure exporting settings can be optionally toggled in the plotting section

close all

%% Model parameters, Signalling and Responsive Regions
% For a "perfectly clustered 1:2 solution" set:
% r1=0.25, r2=0.55, s1=0.6, s2=0.9, beta=-0.6, gen=0.47, bud=0.65
% Not Perfect: same except gen=0.457
%
% For a 2:3 PCS: r1=.35, r2=.55, s1=.6,s2=.8, gen=.33, bud=.75
% For a 2:2 near PCS: r1=.21,r2 =.55, s1=0.6,s2=.9, gen=0.01, bud=.65

r1 = .25;     % beginning of responsive region
r2 = .55;     % end of responsive region
s1 =  0.60;     % beginning of S phase
s2 = 0.90;     % end of S Phase
gen = 0.457; % generational difference cell designated mother
bud = .65;  % position at which budding takes place
beta = -.6;  % recruitment(+)/retardation(-) rate for linear feedback.
mu= .9; % harvest rate
Cap=10020; % Carrying capacity
gamma=.51; % coefficient of growth nonlinearity

% Simulation parameters
ncells = Cap;         % initial number of cells
bins = 98;          % number of bins used for histogram
Cycles = 200;       % number of cell cycles
h = .004;          % step size
steps = ceil(Cycles/h);   % number of time steps
sigma = .0005;  % noise level
sqrh = sqrt(h);   % square root h, used for Euler with noise
hbeta = h*beta;   % used repeatedly, calculate once
harvest_step=(1-exp(-mu*h)); %fraction to be harvested each step
%
% Used for the histograms
    edges = 0:0.01:1;
    H = zeros(100,length(edges)-1);
%
t = 0;   % start time at 0
Ncells = [];
BI = []; % record of bud index
T = [];  % record of time
u = unifrnd(0,1,ncells,1); % initial position of cells in cell cycle
% u = linspace(0,1-1/ncells,ncells)';  % an alternate initial distribution
D = []; M1 =[]; M2 =[]; M3 =[]; M4 =[]; M5=[]; M6=[]; % Space for record
%
% Assign initial cells to generations
    ug = zeros(size(u));
    ug(2:2:end) = 1;
    ug(4:4:end) = 2;
    ug(8:8:end) = 3;
    ug(16:16:end) = 4;
for i = 1:steps
%     if C < ncells;
   % u= u+h;
%     else
    u = u + h*(1 + gamma*(Cap-ncells)/Cap);
    %u = u + h;  % normal growth w/noise
    t = t + h;  % step forward time, t is not used until the end
    S = find(s1<=u & u<=s2); % indices of cells in S phase
    I = max(size(S))/Cap;  % fraction of cells in S phase
    J = find( r1 <=u & u <= r2);  %indices of cells in R phase
    %J = (u>= r1 & u <= r2);
    u(J) = u(J) + hbeta*I*ones(size(J)); % apply feedback to cells in R 
    u=u+sqrh*sigma*(rand(1,ncells)-.5)'; % add noise
    Div = find(u>1);  % find cells that have reached division at 1
    if isempty(Div) == 0  % Did any divide? 
        u_d=u(Div)-1; % locations of the new daughters
        u(Div) = u(Div)+gen-1; % location of new mothers
        ug(Div) = ug(Div)+1;
        u=[u;u_d]; % adding the new daughters to u  
        ug_new = zeros(size(Div));
        ug = [ ug; ug_new];
    end
    keep = rand(size(u)) > harvest_step;
    u = u(keep); % keep the cells that were not harvested.
    ug=ug(keep);
   ncells=length(u);
% Record the last two cycles   
    if t >= Cycles-2.0  % last 2 cycles
        B = find(u>bud);  % find the budded cells
        bi = length(B)/ncells;  % calculate the current Bud Index
        BI = [BI; bi];    % record Bud Index
        Ncells=[Ncells, ncells];
        T = [T; t];       % record Time
     % Find and count number in each generation
        d  = sum(ug==0)/ ncells;
        m1 = sum(ug==1)/ ncells;
        m2 = sum(ug==2)/ ncells;
        m3 = sum(ug==3)/ ncells;
        m4 = sum(ug==4)/ ncells;
        m5 = sum(ug==5)/ncells;
        m6 = sum(ug>=6)/ncells;
            D = [ D d];
            M1 = [M1 m1];
            M2 = [M2 m2];
            M3 = [M3 m3];
            M4 = [M4 m4];
            M5 = [M5 m5];
            M6 = [M6 m6];
            M = M1+M2+M3+M4+M5+M6;
    end 
end
%%
% Plot the final distribution
    imom = find(ug>0);
    ichild = find(ug==0);
   %
    figure  
    histogram(u(ichild),edges,FaceColor='g',FaceAlpha = 1) % Histogram of density in final state
    xlim([0 1])
    xlabel('x')
    ylabel('count')
    ax = gca;
    ax.FontSize = 20;
    hold on
    histogram(u,edges,FaceColor='y',FaceAlpha = .4)
    N = histcounts(u,edges);
    plot([r1 r2],[0 0],'r','LineWidth',6)
    plot([s1 s2],[0 0],'b','LineWidth',6)
    plot(gen,0,'kx','markersize',18,'linewidth',2)
    legend('daughters','mothers','R','S','g','location','northwest','FontSize',18)
    exportgraphics(gcf,'nPCS_12_snap.png')
    % exportgraphics(gcf,'4_5_PCS.eps')
    % savefig('4_5_PCS.fig')
    hold off
%    
% Plot the generational distribution vs. time for the last two cycles.
%
figure
plot(T,D,'k--',T,M,'b-.',T,M1,T,M2,T,M3,T,M4,T,M5,'linewidth',3)
xlabel('time');
ylabel('fraction of cells in each generation')
legend('D','M','M1','M2','M3','M4','M5','fontsize',17)
ax = gca;
ax.FontSize = 20;
exportgraphics(gcf,'nPCS_12_gen.png')
%
% Plot the time series of Density and Bud Index over the last two cycles.
%
figure
yyaxis left
plot(T,BI,'bd-')
xlabel('time')
ylabel('fraction budded')
ylim([0 1])
ax = gca;
ax.FontSize = 20;
hold on
% Plot inverse golden ratio
iphi = 0.618034; %  phi-1 or 1/phi where phi is the golden ratio
plot([T(1) T(end)], [iphi  iphi],'b')
yyaxis right
plot(T,Ncells,'r*-')
ylim([6000 14000])
ylabel('cells')
legend('Bud Idx','1/phi','Density','location','northwest','FontSize',18)
exportgraphics(gcf,'nPCS_12_time.png')
hold off


