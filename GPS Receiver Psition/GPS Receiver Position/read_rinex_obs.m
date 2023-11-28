function [rinex,rec_xyz] = read_rinex_obs(fname, nepoch, nlines )
%*******************************************************
% function [ rinex ] = read_rinex_obs(fname)
%
% DESCRIPTION:
%     This function reads a RINEX format GPS data
%     file and returns the data in an array.
%  
% ARGUMENTS:
%     fname (str) - RINEX file
%     NURULLAH nmeasurement - to restrict number of measurements by user
%  
% OUTPUT:
%     rinex - rinex data (v3 struct)
%
% CALLED BY:
%     process_rinex.m
%
% FUNCTIONS CALLED:
%     read_rinex_header.m
%
% MODIFICATIONS:    
%     XX-XX-03  :  Jan Weiss - Original
%     03-14-06  :  Jan Weiss - Cleanup
%               :  See SVN log for further modifications
%     10-26-06  :  simplified for class positioning project & allows for
%     limited number of lines read (PA)
% 
% Colorado Center for Astrodynamics Research
% Copyright 2006 University of Colorado, Boulder
%*******************************************************

if (nargin < 3)
    nlines = 1e6;
end

% Initialize variables
rinex_data = [];
% line_count = 1; % NURULLAH - Initialized at read_rinex_header function
epoch_count = 0; % NURULLAH - Counts how many epochs are done

% Read header
[ fid, rec_xyz, observables, line_count ] = read_rinex_header(fname);
num_obs = length(observables);

% Get the first line of the observations.
current_line = fgetl(fid);
line_count = line_count+1;
    
% If not at the end of the file, search for the desired information.
while current_line ~= -1 & line_count < nlines & epoch_count < nepoch

    % NURULLAH - skip if line is comment or wrong
    while(isspace(current_line(2:3)) | ~isempty(strfind(current_line,'COMMENT')))
        current_line = fgetl(fid);
        line_count = line_count + 1;
    end

    % Get the time for this data epoch.
    current_time = [ str2num(current_line(2:3)) ; str2num(current_line(5:6)) ; ...
            str2num(current_line(8:9)) ; str2num(current_line(11:12)) ; ...
            str2num(current_line(14:15)) ; str2num(current_line(17:27)) ]';
    
    % How many SV's are there?
    current_num_sv = str2num(current_line(30:32));    
    
    % NURULLAH - Edited for 13 or more num_sv. Rinex 2.11 adds one more row
    % for the rest.
    % Figure out which PRN's there are.
    iii = 1;
    iiii = 1;
    for ii=1:current_num_sv  
        if ii == 13 || ii == 25
            current_line = fgetl(fid);
            line_count = line_count + 1;
        end

        if ii >= 13 && ii <= 24
            current_prn_satType(ii) = current_line( 30+3*iii );
            current_prn(ii) = str2num(current_line( 31+3*iii : 32+3*iii ));
            iii = iii + 1;
        elseif ii > 24
            current_prn_satType(ii) = current_line( 30+3*iiii );
            current_prn(ii) = str2num(current_line( 31+3*iiii : 32+3*iiii ));
            iiii = iiii + 1;
        else
            current_prn_satType(ii) = current_line( 30+3*ii );
            current_prn(ii) = str2num(current_line( 31+3*ii : 32+3*ii ));
        end
    end


    % Get the data for all SV's in this epoch.
    for ii=1:current_num_sv
                
        % Get the next line.
        current_line = fgetl(fid);
        line_count = line_count + 1;
        
        % Check the length of the line and pad it with zeros to
        % make sure it is 80 characters long.
        current_line = check_rinex_line_length(current_line);
        
        % Get the observables on this line.
        current_obs = [ str2num(current_line(1:14)) ; str2num(current_line(17:30)) ; ...
                str2num(current_line(33:46)) ; str2num(current_line(49:62)) ; str2num(current_line(65:78)) ];
        
        % If there are > 5 observables, read another line to get the rest of the observables for this SV.
        if num_obs > 5            
             % Get the next line.
             current_line = fgetl(fid);
             line_count = line_count + 1;

             % Check the length of the line and pad it with zeros to
             % make sure it is 80 characters long.
             current_line = check_rinex_line_length(current_line);
           
            % Append the data in this line to the data from previous line.
            current_obs = [ current_obs ; str2num(current_line(1:14)) ; ...
                            str2num(current_line(17:30)) ; str2num(current_line(33:46)) ; ...
                            str2num(current_line(49:62)) ; str2num(current_line(65:78)) ];
               
        end  % if num_obs > 5        
       
        % Format the data for this PRN as Date/Time, PRN, Observations.
        current_data = [ current_time , current_prn(ii) , current_obs' ];       
       
        % Keep only data for the specified PRNs
        if nargin == 3 & PRN_list & isempty(find(PRN_list == current_prn(ii)))
            continue
        end           
       
        if current_prn_satType(ii) ~= 'G'

        else
            %Append to the master rinex data file.
            rinex_data = [ rinex_data ; current_data ];   
        end     % if not G type
 
    end  % for ii=1:current_num_sv
   
    % Get the next line.
    current_line = fgetl(fid);
    line_count = line_count + 1;  

    % NURULLAH
    epoch_count = epoch_count+1;
   
end  % while current_line ~= -1

% Convert time format
[ gpswk, gpssec ] = cal2gpstime(rinex_data(:,1:6));
rinex.data = [ gpswk gpssec rinex_data(:, 7:end) ];

% Define columns
rinex = define_cols(rinex, observables);

% Convert CP to meters
rinex = convert_rinex_CP(rinex);


% =========================================================================
function rinex = define_cols(rinex, observables)

% Defaults
rinex.col.WEEK = 1;
rinex.col.TOW = 2;
rinex.col.PRN = 3;

col_offset = 3;

for ii=1:length(observables)
    switch observables{ii}
        case {'L1'}
            rinex.col.L1 = ii + col_offset;
        case {'L2'}
            rinex.col.L2 = ii + col_offset;
        case {'P1'}
            rinex.col.P1 = ii + col_offset;
        case {'P2'}
            rinex.col.P2 = ii + col_offset;
        case {'C1'}
            rinex.col.C1 = ii + col_offset;    
        case {'D1'}
            rinex.col.D1 = ii + col_offset;
        case {'S1'}
            rinex.col.S1 = ii + col_offset;
        case {'S2'}
            rinex.col.S2 = ii + col_offset;
    end  % switch    
end  % for ii=1:length(observables)


% =========================================================================
function [ rinex ] = convert_rinex_CP(rinex)

set_constants;

if rinex.col.L1 ~= 0
    rinex.data(:, rinex.col.L1) = rinex.data(:, rinex.col.L1) * LAMBDA_L1;
end


% =========================================================================
function [ fid, rec_xyz, observables, line_count ] = read_rinex_header( file_name )

% NURULLAH - Initialize line_count
line_count = 0;

% Initialize the observables variable.
observables={};

% Assign a file ID and open the given header file.
fid=fopen(file_name);

% If the file does not exist, scream bloody murder!
if fid == -1
    display('Error!  Header file does not exist.');
else    
    % Set up a flag for when the header file is done.
    end_of_header=0;
    
    % Get the first line of the file.
    current_line = fgetl(fid);
    line_count = line_count+1;
    
    % If not at the end of the file, search for the desired information.
    while end_of_header ~= 1        
        % Search for the approximate receiver location line.
        if strfind(current_line,'APPROX POSITION XYZ')
                        % Read xyz coordinates into a matrix.
            [rec_xyz] = sscanf(current_line,'%f');
        end
        
        % NURULLAH - Edited for 10 or more number of types. Rinex 2.11 adds
        % one more row for the rest.
        % Search for the number/types of observables line.
        if strfind(current_line,'# / TYPES OF OBSERV')            
            % Read the non-white space characters into a temp variable.
            [observables_temp] = sscanf(current_line,'%s');            
            
            % Read the number of observables space and then create
            % a matrix containing the actual observables.
            [~, tf] = str2num(observables_temp(2));   %Total  number is 2 or 1 digit?
            if tf == false      %Total number is 1-digit.
                for ii = 1:str2num(observables_temp(1))                
                    observables{ii} = observables_temp( 2*ii : 2*ii+1 );
                end
            else                %Total number is 2-digit.
                for ii = 1:(str2num(observables_temp(1)) * 10 + str2num(observables_temp(2)))
                    % Max 9 types are involved in 1-row, so jump to other row to
                    % read the rest.
                    observables{ii} = observables_temp( mod(2*ii+1,20) : mod(2*ii+2,20) );
                    if ii == 9
                        observables{ii} = observables_temp( 19 : 20 );
                        current_line = fgetl(fid);
                        line_count = line_count+1;
                        [observables_temp] = sscanf(current_line,'%s');
                    end
                end   
            end
        end
                  
        % Get the next line of the header file.
        current_line = fgetl(fid);
        line_count = line_count+1;
        
        %Check if this line is at the end of the header file.
        if strfind(current_line,'END OF HEADER')
            end_of_header=1;
        end        
    end
end

