clc
clear
format long g
addpath('C:\Users\nurul\Desktop\GPS Receiver Psition\GPS Receiver Position\EKF')
addpath('C:\Users\nurul\Desktop\GPS Receiver Psition\GPS Receiver Position\data')
addpath('C:\Users\nurul\Desktop\GPS Receiver Psition\GPS Receiver Position\data\record')

% NURULLAH - determines how many epochs will be got in observation
% and navigation data.
nepoch = 1000;   % max 2880

% read observation file to get orbit elements
[obs_raw,rec_xyz_actual] = read_rinex_obs('madr2000.06o', nepoch);
% [obs,rec_xyz] = read_rinex_obs('nor10040.23o', nmeasurement);
rec_xyz = zeros(3,1);

% NURULLAH - moved before ephemeris datas
epochs = unique(obs_raw.data(:, obs_raw.col.TOW));

% NURULLAH - find ephemeris data if only epochs of observation data is
% found
% read navigation file to get orbit elements
[ephemeris_raw, matched_epochs] = read_rinex_nav('brdc2000.06n', epochs);
% ephemeris = read_rinex_nav('nor10040.23d', nmeasurement);

% NURULLAH - filter observation data by selecting only epochs found in
% navigation data
j = 1;
obs_epoch_filtered.col = obs_raw.col;
for i = 1 : size(obs_raw.data,1)
    if ~isempty(find(obs_raw.data(i,2) == matched_epochs, 1)) 
        obs_epoch_filtered.data(j,:) = obs_raw.data(i,:);
        j = j+1;
    end
end

% NURULLAH - filter observation data regarding a satellite is
% available in both observation and navigation data.
k = 1;
obs.col = obs_epoch_filtered.col;
for i = 1 : size(obs_epoch_filtered.data,1)
    for j = 1 : size (ephemeris_raw,1)
        if ~isempty(find(obs_epoch_filtered.data(i,3) == ephemeris_raw(j,1), 1)) ...
            && obs_epoch_filtered.data(i,2) == ephemeris_raw(j,20) 
            
            obs.data(k,:) = obs_epoch_filtered.data(i,:);
            k = k+1;
        end
    end
end

% NURULLAH - filter navigation data regarding a satellite is
% available in both observation and navigation data.
k = 1;
for i = 1 : size(ephemeris_raw,1)
    for j = 1 : size (obs.data,1)
        if ~isempty(find(ephemeris_raw(i,1) == obs.data(j,3), 1)) ...
            && ephemeris_raw(i,20) == obs.data(j,2)
            
            ephemeris(k,:) = ephemeris_raw(i,:);
            k = k+1;
        end
    end
end

% NURULLAH - fixed to number of matched epochs
% TimeSpan=epochs(1:2880);
TimeSpan=matched_epochs(1:size(matched_epochs,1));

% Broadcast Orbit
satOrbits.XS=zeros(1,length(TimeSpan));
satOrbits.YS=zeros(1,length(TimeSpan));
satOrbits.ZS=zeros(1,length(TimeSpan));
satOrbits.VXS=zeros(1,length(TimeSpan));
satOrbits.VYS=zeros(1,length(TimeSpan));
satOrbits.VZS=zeros(1,length(TimeSpan));
satOrbits.clk=zeros(1,length(TimeSpan));
satOrbits.Rel=zeros(1,length(TimeSpan));

% GPS Satellite Measurements
c = 2.99792458e8 ; % speed of light (m/s)
fL1 = 1575.42e6;   % L1 frequency (Hz)
fL2 = 1227.6e6;    % L2 frequency (Hz)
B=fL2^2/(fL2^2-fL1^2);
A=-B+1;
satOrbits.C1=zeros(1,length(TimeSpan));
satOrbits.L1=zeros(1,length(TimeSpan));
satOrbits.P2=zeros(1,length(TimeSpan));
satOrbits.L2=zeros(1,length(TimeSpan));
satOrbits.P3=zeros(1,length(TimeSpan)); % Iono free pseudorange
satOrbits.CorrP1=zeros(1,length(TimeSpan)); % Corrected Pseudorange from broadcast orbit
satOrbits.CorrP2=zeros(1,length(TimeSpan)); % Corrected Pseudorange from precise orbit
satOrbits.TOW=TimeSpan';
satOrbits.PRN=0;

satOrbits = repmat(satOrbits,1,32);
for ii=1:32
    satOrbits(ii).PRN=ii;
end

% Initialize User Position
userPos=zeros(length(TimeSpan),4);

% NURULLAH - JUST COMMENT - calculates pseudorange and satpos for each
% epoch
for ii=1:length(TimeSpan)    
    this_TOW = TimeSpan(ii);
    index = find(obs.data(:,obs.col.TOW) == this_TOW);
    curr_obs.data = obs.data(index, :);
    curr_obs.col = obs.col;

    % NURULLAH - JUST COMMENT - calculates pseudoranges
    for jj=1:size(curr_obs.data,1)        
        PRN_obs.data = curr_obs.data(jj,:);
        PRN_obs.col = curr_obs.col;
        
        % Record Measurements
        satOrbits(PRN_obs.data(PRN_obs.col.PRN)).C1(ii)=PRN_obs.data(PRN_obs.col.C1);
        satOrbits(PRN_obs.data(PRN_obs.col.PRN)).L1(ii)=PRN_obs.data(PRN_obs.col.L1);
        satOrbits(PRN_obs.data(PRN_obs.col.PRN)).P2(ii)=PRN_obs.data(PRN_obs.col.P2);
        satOrbits(PRN_obs.data(PRN_obs.col.PRN)).L2(ii)=PRN_obs.data(PRN_obs.col.L2);
        
        % Calculate Iono Free Measurement
        P1 = satOrbits(PRN_obs.data(PRN_obs.col.PRN)).C1(ii);
        P2 = satOrbits(PRN_obs.data(PRN_obs.col.PRN)).P2(ii);
        P3=A*P1+B*P2;
        satOrbits(PRN_obs.data(PRN_obs.col.PRN)).P3(ii)=P3;
    end
    
    stop = 10;
    while stop ~= 1
        for jj=1:size(curr_obs.data,1)
            PRN_obs.data = curr_obs.data(jj,:);
            PRN_obs.col = curr_obs.col;
            
            % Obtain the broadcast orbits
            % NURULLAH - JUST COMMENT - sents ephemeris of only one
            % satellite for only one epoch
            PRN_obs = get_broadcast_orbits(PRN_obs,ephemeris,rec_xyz');
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).XS(ii)=PRN_obs.data(PRN_obs.col.XS);
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).YS(ii)=PRN_obs.data(PRN_obs.col.YS);
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).ZS(ii)=PRN_obs.data(PRN_obs.col.ZS);
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).VXS(ii)=PRN_obs.data(PRN_obs.col.VXS);
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).VYS(ii)=PRN_obs.data(PRN_obs.col.VYS);
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).VZS(ii)=PRN_obs.data(PRN_obs.col.VZS);
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).clk(ii)=PRN_obs.data(PRN_obs.col.satClkCorr);
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).Rel(ii)=PRN_obs.data(PRN_obs.col.Rel);        
            % NURULLAH
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).Weight(ii)=mod((PRN_obs.data(PRN_obs.col.C1))*1000,10);        

            % Calculate corrected pseudorange based on broadcast orbit
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).CorrP1(ii)=...
                satOrbits(PRN_obs.data(PRN_obs.col.PRN)).P3(ii)+...
                satOrbits(PRN_obs.data(PRN_obs.col.PRN)).clk(ii)+satOrbits(PRN_obs.data(PRN_obs.col.PRN)).Rel(ii);            
        end
        
        % Calculate User Position
        [broadcast_obs,~]=createObs(this_TOW,satOrbits);    %satpos, pseudorange, PRN
        % NURULLAH - edit parameters as to fit posCalc functions
        satpos = zeros(3, size(broadcast_obs.data, 1));
        for mm = 1: size(broadcast_obs.data, 1)
            for nn = 1:3
                satpos(nn, mm) = broadcast_obs.data(mm, nn);
            end
            n = 1;
            obs_pseudo(mm) = broadcast_obs.data(mm, 4);
            weight(mm) = broadcast_obs.data(mm, 7);
            if weight(mm) == 0
                weight(mm) = 1;
            end
        end

        delta_xyz = comp_pos(broadcast_obs,rec_xyz');
        rec_xyz = rec_xyz + delta_xyz(1:3);
        
        stop=stop-1;
    end

    userPos(ii,1:4) = [rec_xyz; delta_xyz(4)]';
    [lon1(ii),lat1(ii),alt1(ii)] = Geodetic(rec_xyz);

    % NURULLAH - calculate position error for each axis
    err_x(ii) = userPos(ii,1) - rec_xyz_actual(1);
    err_y(ii) = userPos(ii,2) - rec_xyz_actual(2);
    err_z(ii) = userPos(ii,3) - rec_xyz_actual(3);

    % NURULLAH
    
    % NURULLAH - calculates least squares 
    [pos_least_squares(ii,1:4)] = leastSquarePos(satpos, obs_pseudo);  
    
    % % NURULLAH - calculates weighted least squares
    % [pos_weighted_ls(ii,1:4), el, az, dop_weighted_ls] = weightedLeastSquarePos(satpos, obs_pseudo, weight);
    
    % NURULLAH - calculates multi lateration
    % [pos_multilateration(ii,1:4)] = algebraicGPSequations(satpos, obs_pseudo);
    
    % NURULLAH - calculates extended kalman filter
    % pos_extended_kalman_filter(:,ii) = extended_kalman(satpos, obs_pseudo);

end


fprintf('mean of user position:');
fprintf('%16f%16f%16f\n',mean(userPos(:,1)),mean(userPos(:,2)),mean(userPos(:,3)));
fprintf('delta_xyz:');
fprintf('%16f\n',mean(userPos(:,4)));

% NURULLAH
% plot_err(rec_xyz_actual, pos_least_squares, pos_weighted_ls, ...
%     pos_multilateration, pos_extended_kalman_filter);
plot_err(rec_xyz_actual, matched_epochs, pos_least_squares, pos_multilateration);


