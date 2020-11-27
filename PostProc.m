%*************************************************************************%
%  TU BS - Technical University of Braunschweig
%  script PostProc
%*************************************************************************%
%  Post-Processing for fully developed channel flow simulation using
%  openFOAM and a DNS test case
%  Includes comparisons of distinct turbulence models with DNS data
%  Models: kOmega (kOmg), kEpsilon (kEps), SST, Spalart-Allmaras (SA)
%*************************************************************************%
%  Autoren: Henrik Wüstenberg, ...
%  Erstellt am:        Jun.  2019
%  Letzte Änderung am: 06.06.2019 durch H. Wüstenberg                     
%*************************************************************************%
%% ------------------------------------------------------------------------
%	DATA INPUT
%--------------------------------------------------------------------------
% kOmega
kOmg = read_data('kOmega','k-w');

% kEpsilon
kEps = read_data('kEpsilon','k-e');

% SST-kOmega
sst = read_data('SSTkOmega','SST');

% Spalart-Allmaras
sa = read_data('SA','SA');


% DNS
DNS = struct();
fileID = fopen('DNS\DNS_Database.dat','r');
DNS.all = fscanf(fileID,'%f %f %f %f',[4 inf]);
fclose(fileID);

DNS.yp = DNS.all(1,:);
DNS.up = DNS.all(2,:);
DNS.kp = DNS.all(3,:);
DNS.Rp = DNS.all(4,:);
DNS.abbr = 'DNS';


% Constants
rho     = 1.223;            % Density               [kg/m^3]
nu      = 2e-5;             % Kinematic viscosity   [m^2/s]
u_b     = 0.1335;           % Bulk velocity         [m/s]
Re      = 13350;            % Reynoldsnumber        [-]
h       = 1;                % Half channel height   [m]

yp      = 1;                % Normalized mesh height[-]
ny      = length(kOmg.U);   % Amount of grid points [-]
ny_half = ny/2-1;           % Grid points for half channel width    [-]
ny_plot = 200;              % Guessed value for y+ > 400 [-]



%% ------------------------------------------------------------------------
%	NON-DIMENSIONALISE DATA
%--------------------------------------------------------------------------
% kOmega
kOmg = normalise(kOmg,nu,rho);

% kEpsilon
kEps = normalise(kEps,nu,rho);

% SST-kOmega
sst = normalise(sst,nu,rho);

% Spalart-Allmaras
sa = normalise(sa,nu,rho);



%% ------------------------------------------------------------------------
%	PLOT RESULTS
%--------------------------------------------------------------------------
% Graph u+/y+
comp_Up(DNS,kOmg,sst,sa,kEps,'ALL',ny_half)


% Graph k+/y+
comp_kp(DNS,kOmg,sst,sa,kEps,'ALL',ny_plot)


% Graph -R+/y+
comp_Rp(DNS,kOmg,sst,sa,kEps,'ALL',ny_plot)






%% ------------------------------------------------------------------------
%	ADDITIONAL FUNCTIONS
%--------------------------------------------------------------------------
function strct = read_data(filename,abbr)
    strct = struct();
    % Read U
    fileID = fopen(strcat(filename,'\Uf.xy'),'r');
    strct.U = fscanf(fileID,'%f %f',[2 inf]);
    fclose(fileID);
    
    % Read k
    fileID = fopen(strcat(filename,'\k.xy'),'r');
    strct.k = fscanf(fileID,'%f %f',[2 inf]);
    fclose(fileID);
    
    % Read R
    fileID = fopen(strcat(filename,'\Rfw.xy'),'r');
    strct.R = fscanf(fileID,'%f %f',[2 inf]);
    fclose(fileID);
    
    % Read y
    strct.y = strct.U(1,:);
    
    % Set abbreviation for legend
    strct.abbr = abbr;
    
    % Add zero coordinates
    strct.y = [0          strct.y 0         ];
    strct.U = [zeros(2,1) strct.U zeros(2,1)];
    strct.k = [zeros(2,1) strct.k zeros(2,1)];
    strct.R = [zeros(2,1) strct.R zeros(2,1)];
end



function strct = normalise(strct,nu,rho)
    % Determine u_tau
    dudy_w  = strct.U(2,2)/strct.y(2);
    tau_w   = nu*(dudy_w);
    u_tau   = sqrt(tau_w/rho);

    % Normalise all quantities
    strct.yp     = strct.y .* u_tau ./ nu;    % y+
    strct.up     = strct.U ./ u_tau;          % U+
    strct.kp     = strct.k ./ u_tau^2;        % k+
    strct.Rp     = strct.R ./ u_tau^2;        % R+ = mean(u'v')+
    
    % Save quantities for comparison
    strct.u_tau  = u_tau;
    strct.Re_tau = u_tau/nu;
end



function comp_Up(Data1,Data2,Data3,Data4,Data5,name,nPlot)
    figure('name',strcat(name,' u+/y+'));
    semilogx(Data1.yp,Data1.up,'-k')
    hold on
    semilogx(Data2.yp(1:nPlot),Data2.up(2,1:nPlot),'-r')
    semilogx(Data3.yp(1:nPlot),Data3.up(2,1:nPlot),'--r')
    semilogx(Data4.yp(1:nPlot),Data4.up(2,1:nPlot),'-.r')
    semilogx(Data5.yp(1:nPlot),Data5.up(2,1:nPlot),'.r')
    hold off
    xlabel('y+')
    ylabel('u+')
    xlim([0 200])
    ylim([0 30])
    legend(Data1.abbr,Data2.abbr,Data3.abbr,Data4.abbr,Data5.abbr)
    axUp = gca;
    set(axUp,'FontSize',12)
    print(strcat(name,'_u+_y+'),'-dpng')
end



function comp_kp(Data1,Data2,Data3,Data4,Data5,name,nPlot)
    figure('name',strcat(name,' k+/y+'));
    plot(Data1.yp,Data1.kp,'-k')
    hold on
    plot(Data2.yp(1:nPlot),Data2.kp(2,1:nPlot),'-r')
    plot(Data3.yp(1:nPlot),Data3.kp(2,1:nPlot),'--r')
    plot(Data4.yp(1:nPlot),Data4.kp(2,1:nPlot),'-.r')
    plot(Data5.yp(1:nPlot),Data5.kp(2,1:nPlot),'.r')
    hold off
    xlabel('y+')
    ylabel('k+')
    xlim([0 400])
    ylim([0 4])
    legend(Data1.abbr,Data2.abbr,Data3.abbr,Data4.abbr,Data5.abbr)
    axkp = gca;
    set(axkp,'FontSize',12)
    print(strcat(name,'_k+_y+'),'-dpng')
end



function comp_Rp(Data1,Data2,Data3,Data4,Data5,name,nPlot)
    figure('name',strcat(name,' -R+/y+'));
    plot(Data1.yp,-Data1.Rp,'-k')
    hold on
    plot(Data2.yp(1:nPlot),-Data2.Rp(2,1:nPlot),'-r')
    plot(Data3.yp(1:nPlot),-Data3.Rp(2,1:nPlot),'--r')
    plot(Data4.yp(1:nPlot),-Data4.Rp(2,1:nPlot),'-.r')
    plot(Data5.yp(1:nPlot),-Data5.Rp(2,1:nPlot),'.r')
    hold off
    xlabel('y+')
    ylabel('-R+')
    xlim([0 400])
    ylim([0 1.2])
    legend(Data1.abbr,Data2.abbr,Data3.abbr,Data4.abbr,Data5.abbr)
    axRp = gca;
    set(axRp,'FontSize',12)
    print(strcat(name,'_R+_y+'),'-dpng')
end





