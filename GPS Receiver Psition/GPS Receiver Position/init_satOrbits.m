function satOrbits = init_satOrbits(TimeSpan)
    % Broadcast Orbit
    satOrbits.XS=zeros(1,length(TimeSpan));
    satOrbits.YS=zeros(1,length(TimeSpan));
    satOrbits.ZS=zeros(1,length(TimeSpan));
    satOrbits.VXS=zeros(1,length(TimeSpan));
    satOrbits.VYS=zeros(1,length(TimeSpan));
    satOrbits.VZS=zeros(1,length(TimeSpan));
    satOrbits.clk=zeros(1,length(TimeSpan));
    satOrbits.Rel=zeros(1,length(TimeSpan));
    
    satOrbits.C1=zeros(1,length(TimeSpan));
    satOrbits.L1=zeros(1,length(TimeSpan));
    satOrbits.P2=zeros(1,length(TimeSpan));
    satOrbits.L2=zeros(1,length(TimeSpan));
    satOrbits.P3=zeros(1,length(TimeSpan)); % Iono free pseudorange
    satOrbits.CorrP1=zeros(1,length(TimeSpan)); % Corrected Pseudorange from broadcast orbit
    satOrbits.CorrP2=zeros(1,length(TimeSpan)); % Corrected Pseudorange from precise orbit
    satOrbits.TOW=TimeSpan';
    satOrbits.PRN=0;
    satOrbits.S1 = zeros(1,length(TimeSpan)); % NURULLAH
    satOrbits.S2 = zeros(1,length(TimeSpan)); % NURULLAH

    satOrbits = repmat(satOrbits,1,32);

    for ii=1:32
        satOrbits(ii).PRN=ii;
    end
end
