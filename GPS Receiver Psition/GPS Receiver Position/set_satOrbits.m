function satOrbits = set_satOrbits(curr_obs, ephemeris, rec_xyz, satOrbits, A, B, ii)
    for jj=1:size(curr_obs.data,1)        
        PRN_obs.data = curr_obs.data(jj,:);
        PRN_obs.col = curr_obs.col;
        
        % NURULLAH - JUST COMMENT - calculates pseudoranges
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
        
        % NURULLAH - calculates sat positions
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
        % satOrbits(PRN_obs.data(PRN_obs.col.PRN)).S1(ii)=mod((PRN_obs.data(PRN_obs.col.C1))*1000,10);        
        satOrbits(PRN_obs.data(PRN_obs.col.PRN)).S1(ii)=(PRN_obs.data(PRN_obs.col.S1));   
        satOrbits(PRN_obs.data(PRN_obs.col.PRN)).S2(ii)=(PRN_obs.data(PRN_obs.col.S2));  
        % Calculate corrected pseudorange based on broadcast orbit
        satOrbits(PRN_obs.data(PRN_obs.col.PRN)).CorrP1(ii)=...
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).P3(ii)+...
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).clk(ii)+...
            satOrbits(PRN_obs.data(PRN_obs.col.PRN)).Rel(ii);  
    end
end