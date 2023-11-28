function [userPos] = originalCodePos(satOrbits, this_TOW, rec_xyz)

stop = 10;
    while stop ~= 1
        % Calculate User Position
        [broadcast_obs,~]=createObs(this_TOW,satOrbits);    %satpos, pseudorange, PRN
        delta_xyz = comp_pos(broadcast_obs,rec_xyz');
        rec_xyz = rec_xyz + delta_xyz(1:3);

        stop=stop-1;
    end

    % userPos(ii,1:4) = [rec_xyz; delta_xyz(4)]';
    % [lon1(ii),lat1(ii),alt1(ii)] = Geodetic(rec_xyz);

    userPos(1:4) = [rec_xyz; delta_xyz(4)]';
    [lon1,lat1,alt1] = Geodetic(rec_xyz);
