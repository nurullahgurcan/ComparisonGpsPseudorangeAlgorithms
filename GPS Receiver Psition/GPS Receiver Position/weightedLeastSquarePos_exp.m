function [pos] = weightedLeastSquarePos_exp(satpos, obs, weight)
%Function calculates the Least Square Solution.
%
%[pos, el, az, dop] = leastSquarePos(satpos, obs, settings);
%
%   Inputs:
%       satpos      - Satellites positions (in ECEF system: [X; Y; Z;] -
%                   one column per satellite)
%       obs         - Observations - the pseudorange measurements to each
%                   satellite:
%                   (e.g. [20000000 21000000 .... .... .... .... ....])
%       settings    - receiver settings -> deleted
%
%   Outputs:
%       pos         - receiver position and receiver clock error 
%                   (in ECEF system: [X, Y, Z, dt]) 
%       el          - Satellites elevation angles (degrees)
%       az          - Satellites azimuth angles (degrees)
%       dop         - Dilutions Of Precision ([GDOP PDOP HDOP VDOP TDOP])

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%--------------------------------------------------------------------------
%Based on Kai Borre
%Copyright (c) by Kai Borre
%Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
%
% CVS record:
% $Id: leastSquarePos.m,v 1.1.2.12 2006/08/22 13:45:59 dpl Exp $
%==========================================================================

    %=== Initialization =======================================================
    nmbOfIterations = 5;
    useTropCorr = 0;    % NURULLAH
    
    dtr     = pi/180;
    pos     = zeros(4, 1);
    X       = satpos;
    nmbOfSatellites = size(satpos, 2);
    
    A       = zeros(nmbOfSatellites, 4);
    omc     = zeros(nmbOfSatellites, 1);
    az      = zeros(1, nmbOfSatellites);
    el      = az;
    R       = [zeros(nmbOfSatellites)];
    c = 299792458;
    
    % weight = reweight(weight);
    s = zeros(size(weight,2), size(weight,2));



    % exponential model weighting
    a = 0;
    b = 1;
    k = 0.3;
    for i = 1 : size(weight,2)
        s(i,i) = a + b * exp(k * (min(weight) - weight(1,i)));
    end

    carpim_s = 1;
    for i = 1 : size(weight,2)
        carpim_s = carpim_s * s(i,i); 
    end
    
    for i = 1 : size(weight,2)
        sigma(i,i) = s(i, i) / carpim_s^(1/size(weight,2));
    end

    R = sigma;  
    % R = inv(sigma);  

    % sigma model weighting
    % a = 0;
    % b = 1.e5;
    % for i = 1 : size(weight,2)
    %     s(i,i) = sqrt(a + b * 10 ^ (weight(i)/10) );
    % end



    % for i = 1 : size(weight)
    %     weight(1,i) = 1 / weight(1,i);
    % end
    % R = diag((weight));
    
    
    %=== Iteratively find receiver position ===================================
    for iter = 1:nmbOfIterations
    
        for i = 1:nmbOfSatellites
            if iter == 1
                %--- Initialize variables at the first iteration --------------
                Rot_X = X(:, i);
                % trop = 0;   % NURULLAH
                trop = 2;
            else
                %--- Update equations -----------------------------------------
                rho2 = (X(1, i) - pos(1))^2 + (X(2, i) - pos(2))^2 + ...
                       (X(3, i) - pos(3))^2;
                traveltime = sqrt(rho2) / c ;
    
                % NURULLAH - corrections are still made
                %--- Correct satellite position (do to earth rotation) --------
                % Rot_X = e_r_corr(traveltime, X(:, i));
                Rot_X = X(:, i);
    
                %--- Find the elevation angel of the satellite ----------------
                [az(i), el(i), dist] = topocent(pos(1:3, :), Rot_X - pos(1:3, :));
    
                % NURULLAH - pseudrange corrections are still made.
                if (useTropCorr == 1)
                    %--- Calculate tropospheric correction --------------------
                    trop = tropo(sin(el(i) * dtr), ...
                                 0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
                else
                    % Do not calculate or apply the tropospheric corrections
                    trop = 0;
                end
            end % if iter == 1 ... ... else 
    
            %--- Apply the corrections ----------------------------------------
            omc(i) = (obs(i) - norm(Rot_X - pos(1:3), 'fro') - pos(4) - trop);
    
            %--- Construct the A matrix ---------------------------------------
            A(i, :) =  [ (-(Rot_X(1) - pos(1))) / obs(i) ...
                         (-(Rot_X(2) - pos(2))) / obs(i) ...
                         (-(Rot_X(3) - pos(3))) / obs(i) ...
                         1 ];
            E(i) = (omc(i))- A(i, :) * pos;
            
            % NURULLAH
            % PRN = settings.activePRN(i);
            % R(i,i) = settings.acqPeakMetrics(PRN);
            % R(i,i) = weight(i);
            
        end % for i = 1:nmbOfSatellites
    
        % These lines allow the code to exit gracefully in case of any errors
        if rank(A) ~= 4
            pos     = zeros(1, 4);
            return
        end
    

        
        %--- Find position update ---------------------------------------------
        % x   = A \ omc;
        % x = inv(A' * A) * A' * omc;
        x = inv(A' * inv(R) * A) * A' * inv(R) * omc;
        % x = (A' * inv(R) * A) * A' * inv(R) * omc;        % BLUE
        %--- Apply position update --------------------------------------------
        pos = pos + x;
        
    end % for iter = 1:nmbOfIterations
    
    pos = pos';
    
    %=== Calculate Dilution Of Precision ======================================
    if nargout  == 4
        %--- Initialize output ------------------------------------------------
        dop     = zeros(1, 5);
        
        %--- Calculate DOP ----------------------------------------------------
        Q       = inv(A'*A);
        
        dop(1)  = sqrt(trace(Q));                       % GDOP    
        dop(2)  = sqrt(Q(1,1) + Q(2,2) + Q(3,3));       % PDOP
        dop(3)  = sqrt(Q(1,1) + Q(2,2));                % HDOP
        dop(4)  = sqrt(Q(3,3));                         % VDOP
        dop(5)  = sqrt(Q(4,4));                         % TDOP
    end
end

function re_weight = reweight(weight)
    % madr2000.06o
    for i = 1 : size(weight,2)
        if weight(1,i) > 316
            re_weight(1,i) = 9;
        elseif weight(1,i) > 100
            re_weight(1,i) = 8;
        elseif weight(1,i) > 31.6
            re_weight(1,i) = 7;
        elseif weight(1,i) > 10
            re_weight(1,i) = 6;
        elseif weight(1,i) > 3.2
            re_weight(1,i) = 5;
        elseif weight(1,i) > 0
            re_weight(1,i) = 4;
        elseif weight(1,i) == 0
            re_weight(1,i) = 1;
        end
    end

    % elevation angle
    % for i = 1 : size(weight,2)
    %     if weight(1,i) > 60
    %         re_weight(1,i) = 4;
    %     elseif weight(1,i) > 50
    %         re_weight(1,i) = 5;
    %     elseif weight(1,i) > 40
    %         re_weight(1,i) = 6;
    %     elseif weight(1,i) > 30
    %         re_weight(1,i) = 7;
    %     elseif weight(1,i) > 20
    %         re_weight(1,i) = 8;
    %     elseif weight(1,i) > 10
    %         re_weight(1,i) = 9;
    %     elseif weight(1,i) > 0
    %         re_weight(1,i) = 10;
    %     end
    % end
    
end
