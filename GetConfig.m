function config = GetConfig()
%returns a structure of the configuration
config.secFactor = 1e9; % time scaling factor 1e9 means that time is in 
%nano seconds and frequency in GHz
config.c = 3e8/config.secFactor; %light speed m/ns
config.risPhaseMethod = "random"; % in the set {"random","directional"}
%% ofdm
config.f = 30e9/config.secFactor; %frequency Ghz
config.Df= 120e3/config.secFactor;%bandwidth Ghz
config.Nsc = 3e3;% Subcarrier number
config.T = 256;% Subcarrier number
config.lambda = config.c/config.f;
%% geometry
config.bsPos = [1,1,0]; % Bs position
config.risPos = [0,0,0]; % RIS position
config.Mc = 64; % number of RIS elements in a row
config.risElementDist = .5*config.lambda; % RIS element distance
config.TxPower = 100e-3; % transmit power
config.NPSD = 1e-3*10^(-174/10)*config.secFactor; % noise spectral density
config.Noise_Factor = 10^(8/10); % noise factor
% ------------------------
config.sigmaPriorInf =1; %the initial position uncertenty
config.PointsNum=30; % number of points that is  considered in FIG 1
config.xvec = [-10:-1,-.5,.5,1:10];
config.yvec = [0.5,1:10];
[X,Y] = meshgrid(config.xvec,config.yvec);
config.xyz = [X(:).';Y(:).';-Y(:).'];%contians xyz corrdinates of the points
end

