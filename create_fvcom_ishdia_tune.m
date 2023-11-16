%How to use?
%First;install TMD(TMD_Matlab_toolbox) from internet and install model data by visiting TPXO website
%Second; set path
%for tuning salinity and temprature
%
% Which system am I using?
% CHANGE THESE TO YOUR OWN DIRECTOR
    basedir = '/Users/ishid/';
      % Insert your mapped drive to \\stoe\login\your_name here
    addpath([basedir,'']);
    addpath([basedir,'Github/fvcom-toolbox-main/fvcom_prepro/']);
    addpath([basedir,'Github/fvcom-toolbox-main/utilities/']);
    addpath([basedir,'Github/fvcom-toolbox-main/']);
    addpath([basedir,'Github/OceanMesh2D/work/']);
    addpath([basedir,'Github/fvcom-toolbox-main/TMD_Matlab_Toolbox_v2.5/']);
    addpath([basedir,'Github/fvcom-toolbox-main/TMD_Matlab_Toolbox_v2.5/TMD/']);
    addpath([basedir,'Github/fvcom-toolbox-main/TMD_Matlab_Toolbox_v2.5/TMD/DATA/']);
    addpath([basedir,'Github/OceanMesh2D/datasets/RiverDischarge'])
    addpath([basedir,'Github/OceanMesh2D/work2/output/'])
    
% Base folder location - where do you want your forcing input files to end
% up?
inputConf.base = [basedir,'fvcominputs/input_inner_JST'];
inputConf.outbase = [inputConf.base];

%inputConf.AMM_folder = '\\store\projectsa\ISO_Modelling\POLCOMS\OPERATIONAL\OUTPUT\NETCDF\S12\AMM.hourly.';
activate_Tide =0;
activate_River = 0;
activate_HYCOM=0;
activate_Meteo=1; %現在２倍
% Which version of FVCOM are we using (for the forcing file formats)?
inputConf.FVCOM_version = '4.4.2';

% Case name for the model inputs and outputs
% Change this to whatever you want
inputConf.casename ='TokyoBay';
Mobj = read_sms_mesh('2dm','TokyoBay_inner_wide_slope01.2dm');
%Mobj = read_sms_mesh('2dm','highrs.2dm');
write_admesh_mesh(Mobj,'filename','out');
inputConf.grid = ['out.14'];
%Mobj = read_sms_mesh('2dm','Tokyo_bay_fix.2dm')
%write_admesh_mesh(Mobj,'filename','Tokyo_bay2')
%inputConf.grid = Mobj;
% output coordinates (FVCOM only likes cartesian at the moment)
inputConf.coordType = 'cartesian'; % 'spherical' or 'cartesian'
% input coordinates (what's my input bathy in?)
inputConf.coordInput = 'spherical'; % 'spherical' or 'cartesian'
% Input grid UTM Zone (if applicable)
inputConf.utmZone = {'54 N'};

% vertical coordinates type: sigma or hybrid
inputConf.verticalCoordType = 'sigma';

% Sponge layer parameters
inputConf.spongeRadius = -1; % in metres, or -1 for variable
inputConf.spongeCoeff = 0.001;

% z0 value in metres
inputConf.bedRoughness = 0.025; % or 0.015, 0.025 or 0.03 - Davies and Furnes (1980) shelf model

% Estimated velocity (m/s) and tidal range (m) for time step estimate
inputConf.estVel = 1.5;
inputConf.estRange = 2.0;

% Uniform temperature and salinity values
%inputConf.temperature = 10;
%inputConf.salinity = 35;

% Model time ([Y,M,D,h,m,s])
inputConf.modelYear = 2020;
%inputConf.startDate = [inputConf.modelYear,1,1,0,00,00];
inputConf.startDate = [2020,1,1,0,00,00];
inputConf.endDate = [2021,1,1,0,00,00];


% How many tidal constituents do we actually want to use at the model
% boundaries? Case sensitive (M2 != m2).
% Only relevant if using TPXO.
inputConf.tidalComponents = {'M2','S2','N2','K2','K1','O1','P1','Q1'};

% Give some names to the boundaries. This must match the number of node
% strings defined in SMS. Ideally, the order of the names should match the
% order in which the boundaries were made in SMS.
inputConf.boundaryNames = {'obc'};

% Smooth the bathymetry if desired.
%inputConf.smoothBathy = 'yes'; % 'yes' or 'no'.
%if strcmpi(inputConf.smoothBathy, 'yes')
    % Set the smoothing factor and number of iterations (see smoothmesh).
%    inputConf.smoothFactors = [0.5, 4]; % [factor, iterations]
%end

%%%------------------------------------------------------------------------
%%%                      END OF INPUT CONFIGURATION
%%%------------------------------------------------------------------------

%% Prepare the data

% Convert times to Modified Julian Date
inputConf.startDateMJD = greg2mjulian(inputConf.startDate(1),inputConf.startDate(2),inputConf.startDate(3),inputConf.startDate(4),inputConf.startDate(5),inputConf.startDate(6));
inputConf.endDateMJD = greg2mjulian(inputConf.endDate(1),inputConf.endDate(2),inputConf.endDate(3),inputConf.endDate(4),inputConf.endDate(5),inputConf.endDate(6));
%inputConf.inputTimeTS = inputConf.startDateMJD:inputConf.dtTS:inputConf.endDateMJD;

% Read the input mesh and bathymetry. Also creates the data necessary for
% the Coriolis correction in FVCOM.
Mobj = read_grid_mesh('grid',inputConf.grid,...
    'coordinate',inputConf.coordType,'in_coord',inputConf.coordInput,...
    'project',true,'zone',inputConf.utmZone,'addCoriolis',true);

% Parse the open boundary nodes and add accordingly
% Add the sponge nodes
for i=1:size(Mobj.read_obc_nodes,2)
    nodeList = double(cell2mat(Mobj.read_obc_nodes(i)));
    Mobj = add_obc_nodes_list(Mobj,nodeList,inputConf.boundaryNames{i},1);
    if inputConf.spongeRadius < 0    % if we want a variable sponge radius
        if i==1
            % Create an array to store the radii
            Mobj.sponge_rad = zeros(size(Mobj.sponge_nodes));
        end
        % calculate the sponge radius
        spongeRadius = calc_sponge_radius(Mobj,nodeList);
        % Add the sponge nodes to the list
        Mobj = add_sponge_nodes_list(Mobj,nodeList,...
            [inputConf.boundaryNames{i},' sponge'],spongeRadius,...
            inputConf.spongeCoeff);
    else
        Mobj = add_sponge_nodes_list(Mobj,nodeList,...
            [inputConf.boundaryNames{i},' sponge'],inputConf.spongeRadius,...
            inputConf.spongeCoeff);
    end
    clear nodeList
end

clear i

%シグマ座標の取得
% Get the sigma depths in order to interpolate from the POLCOMS depths
if exist(fullfile(inputConf.outbase, 'sigma.dat'),'file')
    % If the sigma.dat file exists, read it
    Mobj = read_sigma(Mobj, fullfile(inputConf.base, 'sigma.dat'));
else
    % If we can't find the sigma.dat file, print an error message and
    % finish
    error(['sigma.dat not found. Please put your sigma.dat file into ',...
        fullfile(inputConf.outbase),' and try again.'])
end

% Do the bed roughness
Mobj.z0 = ones(1,Mobj.nElems)*inputConf.bedRoughness;


% Grid
write_FVCOM_grid(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_grd.dat']));

% Bathymetry
write_FVCOM_bath(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_dep.dat']));

% Coriolis
write_FVCOM_cor(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_cor.dat']));

% Open boundaries
write_FVCOM_obc(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_obc.dat']));

% Sponge file
write_FVCOM_sponge(Mobj,fullfile(inputConf.outbase,[inputConf.casename,'_spg.dat']));

% Bed roughness (constant or variable (see above))
write_FVCOM_z0(Mobj.z0,fullfile(inputConf.outbase,[inputConf.casename,'_z0.nc']),'bottom roughness');

%%


%%---------------------------------------------------
%                        tides
%---------------------------------------------------


%if strcmpi(inputConf.obcForcing, 'phase-amp')
    % Use elevation harmonics, not currents from TPXO
%    inputConf.extractType = 'z'; 
    % Need to cd to TPXO directory or it doesn't work
    % (Yes, this is inelegant but it's the easiest way for now)
%%    here = pwd; % store the current working directory to return later
 %   tpxo_dir = which('TMD');    % find the TPXO directory
 %   tpxo_dir = tpxo_dir(1:end-5);   % remove TPXO.m
 %   cd(tpxo_dir)    % go to TPXO directory
    % Location of the TMD model description file
    %inputConf.Model = [tpxo_dir,'DATA/Model_tpxo7.2'];
%elseif strcmpi(inputConf.obcForcing, 'z')
    % Need to cd to TPXO directory or it doesn't work
    % (Yes, this is inelegant but it's the easiest way for now)
 %%  tpxo_dir = which('TMD');    % find the TPXO directory
  %  tpxo_dir = tpxo_dir(1:end-5);   % remove TPXO.m
  %  cd(tpxo_dir)    % go to TPXO directory
     %Location of the TMD model description file
 %   inputConf.Model = ['Model_PO'];
%elseif strcmpi(inputConf.obcForcing, 'model-output')
    % Use NOC Operational Tide Surge Model output
%    inputConf.extractType = 'm';
%end


if activate_Tide == 1
    inputConf.obcForcing = 'z'; 
    inputConf.datetide = 1/24;
    inputConf.tidesMJD = inputConf.startDateMJD:inputConf.datetide:inputConf.endDateMJD;
    %inputConf.tidesMJD(:10) %ishid
    if strcmpi(inputConf.obcForcing, 'z')
        % Need to cd to TPXO directory or it doesn't work
        % (Yes, this is inelegant but it's the easiest way for now)
        here = pwd; % store the current working directory to return later
        tpxo_dir = which('TMD');    % find the TPXO directory
        tpxo_dir = tpxo_dir(1:end-5);   % remove TPXO.m
        cd(tpxo_dir)    % go to TPXO directory
        %Location of the TMD model description file
        inputConf.Model = ['Model_Tokyo'];

        % Use tmd_tide_pred to predict surface elevations for a given time
        % range.
        
        % Add the tidal components to the Mobj.
        Mobj.Components = inputConf.tidalComponents;
        
        % Create a time series in MATLAB datenum format with ten minute
        % inputs
        inputConf.inputTimeTideZ = datenum(inputConf.startDate):1/144:datenum(inputConf.endDate);
        % Also do Modified Julian Day for the output to NetCDF
        inputConf.inputTimeTideZMJD = datenum(inputConf.startDateMJD):1/144:datenum(inputConf.endDateMJD);
        % Get the indices to use the tidal constituents defined in
        % inputConf.tidalComponents for TPXO (which requires a numerical
        % array of the constituents to be used). The order of the TPXO
        % constituents is M2, S2, N2, K2, K1, O1, P1, Q1, MF, MM, M4, MS4,
        % MN4.
        tpxoConsts = {'2n2','M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', ...
            'MF', 'MM', 'M4', 'MS4', 'MN4'};
        tIndUse = nan(length(Mobj.Components), 1);
        tInd = 1:length(tpxoConsts);
        for i=1:length(Mobj.Components)
            tPos = tInd(strcmp(Mobj.Components{i}, tpxoConsts));
            if ~isempty(tPos)
                tIndUse(i) = tPos;
            else
                warning('Supplied constituent (%s) is not present in the TPXO data', Mobj.Components{i}) %#ok<WNTAG>
            end
        end
        % Tidy up a bit
        clear c tpxoConsts tPos tInd
        % We can't just use tmd_tide_pred to do all the surface elevations
        % at once. Instead, the useful approaches are:
        %
        %   1. Time series at a single location
        %   2. Map of a given time step at all locations
        %
        % Since I'm likely to have many more time steps than locations,
        % it's probably best to do the time series at each location than
        % all the locations and a single time step.
        %
        % The order of the surface elevations in Mobj.surfaceElevation
        % should reflect the order of the open boundary node IDs as FVCOM
        % assumes they just map directly. So, rather than iterate through
        % each position, we need to get the position based on the list of
        % node IDs (Mobj.obc_nodes, without the zeros and in order
        % of each boundary).
        tmpObcNodes = Mobj.obc_nodes';
        ObcNodes = tmpObcNodes(tmpObcNodes~=0)';
        clear tmpObcNodes
    %         obc_lat = Mobj.lat(Mobj.obc_nodes(Mobj.obc_nodes~=0));
    %         obc_lon = Mobj.lon(Mobj.obc_nodes(Mobj.obc_nodes~=0));
        %Mobj.surfaceElevation = nan(size(ObcNodes,2), size(inputConf.inputTimeTideZ,2));
        %%for i=1:size(ObcNodes,2)
        %   currLon = Mobj.lon(ObcNodes(i));
        %   currLat = Mobj.lat(ObcNodes(i));
        %   fprintf('%,%',currrat,currLon);
        %end

    % for i=1:size(ObcNodes,2)
    %     currLon = Mobj.lon(ObcNodes(i));
    %     currLat = Mobj.lat(ObcNodes(i));
    %     fprintf('Position %i of %i (%.3f %.3f)... \n', i, size(ObcNodes,2), currLon, currLat)
    % end
    %tidedata = readmatrix('C:/Users/ishid/data/tide/output/tide_2020_MR.csv');
    tidedata = readmatrix('C:/Users/ishid/data/tide/output/tide_2020_MR_obs.csv');

    %tidedata = readmatrix('C:/Users/ishid/data/tide/output/tide_2020_abura.csv');
    surfaceElevation = nan(Mobj.nObcNodes, size(inputConf.tidesMJD, 2), length(inputConf.boundaryNames));
    for i=1:length(inputConf.boundaryNames)
        for j=1:Mobj.nObcNodes
            % Get the current location (from the node ID)
            currLon = Mobj.lon(Mobj.obc_nodes(i,j));
            currLat = Mobj.lat(Mobj.obc_nodes(i,j));
            %if ftbverbose
                fprintf('Position %i of %i (%.3f %.3f)... \n', j, Mobj.nObcNodes, currLon, currLat);
            %end
        %if UTC
           % [surfaceElevation(j,:,i), ~] = tmd_tide_pred(inputConf.Model, ...
           %     inputConf.tidesMJD+678942.000000, currLat, currLon, 'z', tIndUse);
        %else if JST
        %    [surfaceElevation(j,:,i), ~] = tmd_tide_pred(inputConf.Model, ...
        %  inputConf.tidesMJD+678942.000000-0.375, currLat, currLon, 'z', tIndUse); %0.375=UTC→JST
                %ここを変えます！！
        % else if USE observation data
            surfaceElevation(j,:,i) = tidedata((2:size(inputConf.tidesMJD, 2)+1),3)-0.18;
            %スタートが1/1/2020,すべての開境界ノードに同一のデータを入れるという前提
            if isnan(surfaceElevation(j,:))
                % Try the global model instead.
                [surfaceElevation(j,:,i), ~] = tmd_tide_pred(inputConf.Model, ...
                inputConf.tidesMJD, currLat, currLon, 'z', tIndUse);
            end
        end
    end
    
    Mobj.surfaceElevation = surfaceElevation;
    % Tidy up some more
    clear i j tIndUse obc_lat obc_lon currLon currLat surfaceElevation
    cd(here);
    %save('varb/Mobj_01.mat','Mobj','-v7.3','-nocompression');

        %for i=1:size(ObcNodes,2)
            % Get the current location (from the node ID)
        %    currLon = Mobj.lon(ObcNodes(i));
        %   currLat = Mobj.lat(ObcNodes(i));
        %    fprintf('Position %i of %i (%.3f %.3f)... \n', i, size(ObcNodes,2), currLon, currLat)
        %    [Mobj.surfaceElevation(i,:), ~] = tmd_tide_pred(inputConf.Model, inputConf.inputTimeTideZ, currLat, currLon, 'z', tIndUse);
        %end
        % Tidy up a bit
        %clear tIndUse obc_lat obc_lon ObcNodes currLon currLat
        write_FVCOM_elevtide(Mobj, ...
        inputConf.tidesMJD,...
        fullfile(inputConf.outbase, [inputConf.casename, 'm018_julian_obc.nc']),...
        'Model surface elevation boundary input',...
        'floattime', true,...
        'julian', true);
        fprintf('done.\n')
    else
        % No idea what has been supplied
        error('Unrecognised open boundary node forcing type. Choose phase')
    end
end

cd('C:/Users/ishid/Github/fvcom-toolbox-main/')



%%
%%%------------------------------------------------------------------------
%%%                     River discharge and output
%%%------------------------------------------------------------------------


Mobj.xc = nodes2elems(Mobj.x, Mobj);
Mobj.yc = nodes2elems(Mobj.y, Mobj);
Mobj.lonc = nodes2elems(Mobj.lon, Mobj);
Mobj.latc = nodes2elems(Mobj.lat, Mobj);



inputConf.riverForcing = 'FLUX';

inputConf.river.infos = {...
     'Tamagawa',...
     'Arakawa',...
     'Sumidagawa',...
     'Edogawa',...
     'Tsurumigawa',...
     'Ebigawa',...
     'Mamagawa',...
     'Yorogawa',...
     'Obitsugawa'...
     'koitogawa',...
     'Muratagawa',...
     'Hanamigawa'};

     inputConf.river.infos = {...
     'Tamagawa','Sumidagawa','Edogawa','Tsurumigawa','WestArakawa','EastArakawa','Mamagawa','Ebigawa',...
     'Yorogawa','Obitsugawa','koitogawa','Muratagawa','Hanamigawa'...
     };


inputConf.river.infos={...
'EastArakawa', 'CenterArakawa', 'WestArakawa', 'SouthArakawa',...
       'FirstSumidagawa', 'SecondSumidagawa', 'ThirdSumidagawa', 'OneEdogawa',...
       'TwoEdogawa', 'ThreeEdogawa', 'IchiTamagawa', 'NiTamagawa',...
       'SanTamagawa', 'ATsurumigawa', 'BTsurumigawa', 'Mamagawa', 'Ebigawa',...
       'Yorogawa', 'Obitsugawa', 'koitogawa', 'Muratagawa', 'Hanamigawa'...
};

%}
% Location of river file
%only arakawa ver.
 %inputConf.river.flux = [basedir,'data/river/scr/riverflux_arakawa.csv'];
 %inputConf.river.temp = [basedir,'data/river/scr/rivertemp_arakawa.csv'];
 %inputConf.river.salt = [basedir,'data/river/scr/riversalt_arakawa.csv'];
 %inputConf.river.location = [basedir,'data/river/scr/riverloc_arakawa.csv'];
%{
 inputConf.river.flux = [basedir,'data/river/scr/river_flux2_2019.csv'];
 inputConf.river.temp = [basedir,'data/river/scr/river_temp2_2019.csv'];
 inputConf.river.salt = [basedir,'data/river/scr/river_salt2_2019.csv'];
 inputConf.river.location = [basedir,'data/river/scr/river_location2_2019.csv'];

inputConf.river.flux = [basedir,'data/river/scr/river_flux2020all.csv'];
inputConf.river.salt = [basedir,'data/river/scr/river_salt2020.csv'];
inputConf.river.temp = [basedir,'data/river/scr/river_temp2020.csv'];
inputConf.river.location = [basedir,'data/river/scr/river_loc2020.csv'];
%}
%inputConf.river.flux = [basedir,'data/river/scr/river_flux2020DB_octfix.csv'];
% inputConf.river.flux = [basedir,'data/river/scr/river_flux2020DB_nodev_sumifix.csv'];
 inputConf.river.flux = [basedir,'data/river/scr/river_flux2020_final.csv'];
 inputConf.river.temp = [basedir,'data/river/scr/river_temp2020DB_nodev_2d.csv'];
 inputConf.river.salt = [basedir,'data/river/scr/river_salt2020DB_nodev.csv'];
 inputConf.river.location = [basedir,'data/river/scr/river_loc2020DB_nodev.csv'];

%{
inputConf.river.flux = [basedir,'data/river/scr/river_flux2020DB_Aradiv.csv'];
inputConf.river.temp = [basedir,'data/river/scr/river_temp2020_Aradiv.csv'];
inputConf.river.salt = [basedir,'data/river/scr/river_salt2020_Aradiv.csv'];

inputConf.river.location = [basedir,'data/river/scr/river_loc2020_Aradiv.csv'];
%}
develop_mode = 1 ;
if activate_River == 0
    disp('Skip creating river discharge files.')
else
if strcmpi(inputConf.riverForcing, 'FLUX')
    if develop_mode == 3
    fprintf('Loading Model objet file...')
        load('varb/Mobj_04.mat');
        fprintf('Done!\n');
    else
        fprintf('Reading river forcing file...');
        Mobj = get_COSTUM_river_location(inputConf, Mobj);
        Mobj = get_COSTUM_river_variable(Mobj, ...
            {inputConf.river.flux,inputConf.river.temp,inputConf.river.salt}, ...
            {'flux','temp','salt'},...
            inputConf.river.infos,...
            'time', true);
        Mobj.nRivers = Mobj.river.number;
        inc = 1;
        Mobj.rivermouth = cell(1); % don't preallocate as we don't know how many we'll have
        for s = 1:Mobj.nRivers
            [node, ~] = find_nearest_pt(Mobj.river.location(s, 3), Mobj.river.location(s, 4), Mobj);
            [~, elem] = min(abs(sqrt((Mobj.xc - Mobj.river.location(s, 3)).^2 + Mobj.yc - Mobj.river.location(s, 4)).^2));
            Mobj.rivermouth{inc} = {inc, Mobj.river.location(s, 3), Mobj.river.location(s, 4), node, Mobj.h(node), Mobj.river.name(s), elem};
            riverList(s,1) = Mobj.rivermouth{inc}(1,4);
            Mobj.riverList = cell2mat(riverList);
            inc = inc + 1;
        end
        clear node elem inc s
        Mobj = add_river_nodes_list(Mobj,Mobj.riverList,Mobj.river.name,1);
        save('../00_data/Mobj_04.mat','Mobj','-v7.3','-nocompression');
        fprintf('Done!\n');
        clear riverList;
    end
    % river file
    
    fprintf('Writing river forcing file...\n');
    %{
    write_FVCOM_river(fullfile(inputConf.outbase,...
        [inputConf.casename,'small_river.nc']),...
         inputConf.river.infos,...
         Mobj.river.timeMJD,...
         Mobj.river.flux/3600,...
         Mobj.river.temp,...
         Mobj.river.salt,...
        'Tokyo Bay rivers',...
        'Model river boundary input');
  %}
    write_FVCOM_river(fullfile(inputConf.outbase,...
        [inputConf.casename,'final_river.nc']),...
         inputConf.river.infos,...
         Mobj.river.timeMJD,...
         Mobj.river.flux*1.0,...
         Mobj.river.temp,...
         Mobj.river.salt,...
        'Tokyo Bay rivers',...
        'Model river boundary input');
      
        write_FVCOM_river(fullfile(inputConf.outbase,...
        [inputConf.casename,'final12_river.nc']),...
         inputConf.river.infos,...
         Mobj.river.timeMJD,...
         Mobj.river.flux*1.2,...
         Mobj.river.temp,...
         Mobj.river.salt,...
        'Tokyo Bay rivers',...
        'Model river boundary input');
 
    write_FVCOM_river(fullfile(inputConf.outbase,...
        [inputConf.casename,'final14_river.nc']),...
         inputConf.river.infos,...
         Mobj.river.timeMJD,...
         Mobj.river.flux*1.4,...
         Mobj.river.temp,...
         Mobj.river.salt,...
        'Tokyo Bay rivers',...
        'Model river boundary input');
        
       
    write_FVCOM_river(fullfile(inputConf.outbase,...
        [inputConf.casename,'final16_river.nc']),...
         inputConf.river.infos,...
         Mobj.river.timeMJD,...
         Mobj.river.flux*1.6,...
         Mobj.river.temp,...
         Mobj.river.salt,...
        'Tokyo Bay rivers',...
        'Model river boundary input');

    write_FVCOM_river(fullfile(inputConf.outbase,...
    [inputConf.casename,'final18_river.nc']),...
        inputConf.river.infos,...
        Mobj.river.timeMJD,...
        Mobj.river.flux*1.8,...
        Mobj.river.temp,...
        Mobj.river.salt,...
    'Tokyo Bay rivers',...
    'Model river boundary input');

    write_FVCOM_river(fullfile(inputConf.outbase,...
    [inputConf.casename,'final20_river.nc']),...
        inputConf.river.infos,...
        Mobj.river.timeMJD,...
        Mobj.river.flux*2.0,...
        Mobj.river.temp,...
        Mobj.river.salt,...
    'Tokyo Bay rivers',...
    'Model river boundary input');
 
    write_FVCOM_river_nml(Mobj, ...
        fullfile(inputConf.outbase,'RIVERS_NAMELIST.nml'), ...
        [inputConf.casename,'_river.nc'],...
        '''uniform''');
      
    fprintf('Done!\n');
end
end
%%
%clear ans tpxo_dir


%%
%%%------------------------------------------------------------------------
%%%                    HYCOM S&T forcing and staff
%%%------------------------------------------------------------------------
tic

% Open boundary temperatures and salinities (string for source or number for constant).
% The data is avaliable from 1992-10-02 00:00:00
inputConf.obc_temp = 'HYCOM';
inputConf.obc_salt = 'HYCOM';
inputConf.obc_u = 'None';
inputConf.obc_v = 'None';
inputConf.obctsMJD = [inputConf.startDateMJD, round(inputConf.endDateMJD)+2]; %+1toha?
%inputConf.obctsMJD = [inputConf.startDateMJD-1, inputConf.endDateMJD-1];%こうしないと一日後のデータスタートになってしまう



%about HYCOM
%GLBy0.08 grid is 0.08 deg lon x 0.04 deg lat that covers 80S to 90N.
%** GLBv0.08 Discontinued on 2020-Feb-18 ** grid is 0.08 deg lon x 0.08 deg lat between 40S-40N. 
%Poleward of 40S/40N, the grid is 0.08 deg lon x 0.04 deg lat. It spans 80S to 90N.
%3hours(standard) 1hours(sur)

% Increment used for open boundary ST (days)
inputConf.dateobs = 1/24;

develop_mode = 2;
%clear hycom_s hycom_t
% Now we need some boundary temperature and salinity conditions.
if activate_HYCOM==1
    if strcmpi('HYCOM', {inputConf.obc_temp, inputConf.obc_salt})
        % Use HYCOM data for the boundary forcing.
        % Offset the times to give us a bit of wiggle room.
        if develop_mode == 3
            fprintf('Loading Model objet file...\n')
            load('varb/Mobj_02.mat');
            fprintf('Done!\n');
        else
            if develop_mode == 1
                             % modelTime = inputConf.obctsMJD;varlist={'temperature', 'salinity'}
                fprintf('Downloading daliy open boundary S&T forcing from HYCOM...\n');
                % T
                for i = 1:20
                    try
                        hycom = get_HYCOM_forcing(Mobj, inputConf.obctsMJD, {'temperature','salinity'}); 
                        %hycom_t.startDateMJD = inputConf.startDateMJD;
                        %hycom_t.endDateMJD = inputConf.endDateMJD;
                        break;  % Break out of the i-loop on success
                    catch ME
                        disp(ME);
                        fprintf('Retrying...\n');
                    end
                end
                save(['../00_data/hycom_2020_test.mat'],'hycom','-v7.3','-nocompression')
                Mobj = get_HYCOM_tsobc(Mobj, hycom,inputConf.obctsMJD ,{'temperature'});

                Mobj.backup_temp = Mobj.temperature; 
                Mobj = get_HYCOM_series(Mobj, inputConf.dateobs, 'temperature',true);
               
                %clear hycom_*
               % save(['../00_data/hycom_2017_test.mat'],'hycom','-v7.3','-nocompression')
                %save(['../00_data/hycom_t','_',fname,num2str(inputConf.startDate(2)),'_',num2str(inputConf.startDate(3)),'to',num2str(inputConf.endDate(2)),'_',num2str(inputConf.endDate(3)),'.mat'], 'hycom_t','-v7.3','-nocompression');
              %save(['../00_data/hycom_t','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(2)),'.mat'],...
               %     'hycom_t','-v7.3','-nocompression');
              % '../00_data/hycom_t','_',num2str(inputConf.startDate(2)),'_',num2str(inputConf.startDate(3)),'to',...
              %  num2str(inputConf.endDate(2),'_',num2str(inputConf.endDate(3)),'.mat'
                % S
                %for i = 1:10
                %    try
                %        hycom_s = get_HYCOM_forcing(Mobj, inputConf.obctsMJD, {'salinity'});
                %        hycom_s.startDateMJD = inputConf.startDateMJD;
                %        hycom_s.endDateMJD = inputConf.endDateMJD;
                %        break;  % Break out of the i-loop on success
                %    catch ME
                %        disp(ME);
                %        fprintf('Retrying...\n');
                %    end
                %end
                %save(['../00_data/hycom_2018_2022.mat'],'hycom_s','-v7.3','-nocompression')
                %save(['../00_data/hycom_s','_',fname,num2str(inputConf.startDate(2)),'_',num2str(inputConf.startDate(3)),'to',...
                %num2str(inputConf.endDate(2)),'_',num2str(inputConf.endDate(3)),'.mat'],...
                %    'hycom_s','-v7.3','-nocompression');
                %}
                fprintf('Downloading daliy open boundary meanflow from HYCOM...\n');
                if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
                    fprintf('Writing daliy open boundary meanflow file.\n')
                    % u
                    for i = 1:10
                        try
                            hycom_u = get_HYCOM_forcing(Mobj, inputConf.obctsMJD, {'u'}); 
                            break;  % Break out of the i-loop on success
                        catch ME
                            disp(ME);
                            fprintf('Retrying...\n');
                        end
                    end
                    save(['../00_data/hycom_u','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat'],...
                        'hycom_u','-v7.3','-nocompression');
                    % v
                    for i = 1:10
                        try
                            hycom_v = get_HYCOM_forcing(Mobj, inputConf.obctsMJD, {'v'}); 
                            break;  % Break out of the i-loop on success
                        catch ME
                            disp(ME);
                            fprintf('Retrying...\n');
                        end
                    end
                    save(['../00_data/hycom_v','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat'],...
                        'hycom_v','-v7.3','-nocompression');
                end
                fprintf('Downloading daliy open boundary S&T forcing from HYCOM...Done!\n');
            elseif develop_mode == 2
                fprintf('Loading daliy open boundary S&T forcing from local HYCOM database...\n')
                %load(['../00_data/hycom_t','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(2)),'.mat']);
                %load(['../00_data/hycom_s','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(2)),'.mat']);
                load('../00_data/hycom2020_rev.mat')
                fprintf('Downloading daliy open boundary S&T forcing from HYCOM...Done!\n');
                if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
                    load(['../00_data/hycom_u','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat']);
                    load(['../00_data/hycom_v','_',num2str(inputConf.startDate(1)),'_',num2str(inputConf.endDate(1)),'.mat']);
                end
            end
            %2020_revのときだけ（保存する変数名を間違えた）
            hycom = hycom_t;
            % Interpolate the 4D HYCOM data on the FVCOM vertical grid at the open boundaries.
            Mobj = get_HYCOM_tsobc(Mobj, hycom,inputConf.obctsMJD, {'temperature'});
            Mobj = get_HYCOM_tsobc(Mobj, hycom, inputConf.obctsMJD,{'salinity'});
            if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
                % finding nesting region shows a nand of elements,
                % that is not we wanted, we need obc faces, instead.
                % Nested = find_nesting_region(inputConf, Mobj);
                % Mobj.read_obc_elems = Nested.read_obc_elems;
                Mobj = find_boundary_elements(Mobj);
                Mobj = get_HYCOM_tsobc(Mobj, hycom_u,inputConf.obctsMJD,{'u'});
                Mobj = get_HYCOM_tsobc(Mobj, hycom_v,inputConf.obctsMJD, {'v'});
            end
            %hycom.temperature.data(5,7,:,1)
            clear hycom_*
            % backup daliy data
            Mobj.backup_temp = Mobj.temperature; 
            Mobj.backup_salt = Mobj.salt;
            Mobj.backup_tstm = Mobj.ts_times;

            if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
                Mobj.backup_mflu = Mobj.u;
                Mobj.backup_mflv = Mobj.v;
                Mobj.backup_mftm = Mobj.mf_times;
            end
            % recover daliy data
            
            Mobj.temperature = Mobj.backup_temp; Mobj.salt = Mobj.backup_salt;
           % Mobj.u = Mobj.backup_mflu;           Mobj.v = Mobj.backup_mflv;
            Mobj.ts_times = Mobj.backup_tstm;   % Mobj.mf_times = Mobj.backup_mftm;
            
            % Interpolate the 4D HYCOM data on the hourly time series
            Mobj = get_HYCOM_series(Mobj, inputConf.dateobs,'temperature',true);
            Mobj = get_HYCOM_series(Mobj, inputConf.dateobs,'salinity',true);
            if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
                Mobj = get_HYCOM_series(Mobj, inputConf.dateobs, 'u',true);
                Mobj = get_HYCOM_series(Mobj, inputConf.dateobs, 'v',true);
                % Generate boundary mean flow values at the centroids of the open boundary elements. 
                % For the various output files, we need both u and v as well as depth averaged velocity 
                % at each boundary node position. 
                % An example: get_POLCOMS_meanflow uses PML POLCOMS-ERSEM NetCDF files 
                % to interpolate mean flow to the FVCOM open boundary elements and vertical grid. 
                % The velocity data should be saved in Mobj.velocity of size [nElements, nTime] 
                % and the u and v components in Mobj.meanflow_u and Mobj.meanflow_v as arrays 
                % of size [nElements, nSiglay, nTime].
                % Mobj.meanflow = zeros(numel(Mobj.read_obc_elements{1}), numel(Mobj.mf_times));
                % Stick the values in the mesh structure.
                Mobj.meanflow_u = Mobj.u;
                Mobj.meanflow_v = Mobj.v;
                % Now we have the 3D arrays, create depth averaged velocities too
                Mobj.meanflow_ubar = squeeze(mean(Mobj.meanflow_u, 2));
                Mobj.meanflow_vbar = squeeze(mean(Mobj.meanflow_v, 2));
                % Depth averaged velocity
                Mobj.velocity = squeeze(mean(sqrt(Mobj.meanflow_u.^2 + Mobj.meanflow_v.^2), 2));
            end
            save('varb/Mobj_02.mat','Mobj','-v7.3','-nocompression');
        end
    
    elseif strcmpi('FRA-JCOPE', {inputConf.obc_temp, inputConf.obc_salt})
        % [data,header]=read_grads(file_name,var_name,varargin)
        % 
        % file_name = ['C:\Users\Yulong WANG\Documents\GitHub\jcope-convert\fra_jcope\el.ctl'];
        % var_name = ['all'];
        % [data,header]=read_grads('C:\Users\Yulong WANG\Documents\GitHub\jcope-convert\fra_jcope\t.ctl','all'); 
    end
    fprintf('Open boundary ST forcing making time: %.2f minutes\n', toc / 60)

    tic
    % Write the temperature and salinity.
    %size(Mobj.temperature, 2)
    
    if strcmpi('HYCOM', {inputConf.obc_temp, inputConf.obc_salt})
        fprintf('Writing daliy open boundary S&T forcing file.\n')
        write_FVCOM_tsobc(fullfile(inputConf.outbase, inputConf.casename), ...
            Mobj.ts_times, ...
            size(Mobj.temperature, 2), ...
            Mobj.temperature, ...
            Mobj.salt, ...
            Mobj, ...
            'floattime', true,...
            'julian', true);
    end
    clear hycom_*
    

    % Write the meanflow.
    if strcmpi('HYCOM', {inputConf.obc_u, inputConf.obc_v})
        fprintf('Writing daliy open boundary meanflow file.\n')
        write_FVCOM_meanflow(Mobj, ...
            fullfile(inputConf.outbase,[inputConf.casename,'_mfobc.nc']), ...
            Mobj.velocity);
    end

    % plot temp and salt
    %plot(datetime(Mobj.ts_times+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),squeeze(Mobj.temperature(1,1,:)));
    %plot(datetime(Mobj.ts_times+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),squeeze(Mobj.salt(1,1,:)));

end
    %%
%{
%%%------------------------------------------------------------------------
%%%                     Meteorology and output
%%%------------------------------------------------------------------------
% Get the surface heating data.
%inputConf.obctsMJD = [inputConf.startDateMJD, inputConf.endDateMJD +1];
inputConf.dateforcing = 1/24;
inputConf.doForcing = 'NCEP';
inputConf.forceMJD = inputConf.startDateMJD:inputConf.dateforcing:inputConf.endDateMJD;
if strcmpi(inputConf.doForcing, 'GWO')
    inputConf.forceMJD = inputConf.startDateMJD:inputConf.dateforcing:inputConf.endDateMJD
elseif strcmpi(inputConf.doForcing, 'NCEP')
    inputConf.forceMJD = [inputConf.startDateMJD, inputConf.endDateMJD+1];
end
activate_Meteo=0;
develop_mode = 1;
if activate_Meteo==1
    if strcmpi(inputConf.doForcing, 'NCEP')
        
        % Use the OPeNDAP NCEP script to get the following parameters:
        %     - Potential evaporation rate (pevpr)[W/m^2]                   : for extracting land mask
        %陸上でしかpevprは与えられない、lhtflから導ける
        %lftflはreanalysis2では使えない
        %     - u wind component (uwnd)[m/s]                                : for wind
        %     - v wind component (vwnd)[m/s]                                : for wind
        %     - Precipitation rate (prate)[Kg/m^2/s]                        : for precipitation
        %     - Downward Longwave Radiation Flux at Surface (dlwrf) [W/m^2] : for heat flux
        %     - Downward Solar Radiation Flux at Surface (dswrf) [W/m^2]    : for heat flux
        %     - Upward Longwave Radiation Flux at Surface (ulwrf) [W/m^2]   : for heat flux
        %水面からの長波放射、水温から計算できるのでいらない。暴走は防げるかも
        %     - Upward Solar Radiation Flux at Surface (uswrf) [W/m^2]      : for heat flux
        %     - Sea level pressure (pres) [Pa]                              : for air pressure
        %     - Air temperature at 2 m (air) [Kelvins]
        %     - Relative Humidity on Pressure Levels (rhum) [%]
        % The script also calculate the following parameters:
        %     - Momentum flux (tau)
        %     - Net solar radiation surface (nswrs = uswrf - dswrf)
        %     - Net longwave radiation surface (nlwrs = ulwrf - dlwrf)
        if develop_mode == 3
            fprintf('Loading Model objet file...\n')
            load('varb/Mobj_03.mat');
            load('forcing_ncep_interp.mat');
            fprintf('Done!\n');
        else
            if develop_mode == 1
                % The script converts the NCEP data from the OPeNDAP server from longitudes 
                % in the 0 to 360 range to the latitudes in the -180 to 180 range. 
                % It also subsets for the right region (defined by Mobj.lon and Mobj.lat).
                % Uncomment the variables in get_NCEP_forcing as the varlist shows.
                fprintf('Downloading NCEP forcing from OPeNDAP server database...\n')
                %modelTime=inputConf.forceMJD; varargin
                inputConf.forceMJD = [inputConf.startDateMJD, inputConf.endDateMJD+1 ];
                forcing_ncep = get_NCEP_forcing(Mobj, inputConf.forceMJD, ...
                    'varlist', {...
                     'dswrf', 'dlwrf',...
                    'nlwrs','nswrs','lhtfl','shtfl','uswrf','ulwrf'...
                    },...
                    'source', 'reanalysis2');
                   % 'dswrf', 'dlwrf',...
                    %'nlwrs','nswrs','lhtfl','shtfl','uswrf','ulwrf'...
                    %   'uwnd', 'vwnd',...
                    %   'prate', 'pres', 'air', 'rhum','pevpr',...
                %forcing_ncep = get_NCEP_forcing(Mobj, inputConf.forceMJD, ...
                %    'varlist', {'dswrf'},...
                %    'source', 'reanalysis2');   
                %forcing_ncep.nlwrs.data(1,1)
                forcing_ncep.domain_cols = length(forcing_ncep.lon);
                forcing_ncep.domain_rows = length(forcing_ncep.lat);
                if isfield(forcing_ncep, 'rhum')||isfield(forcing_ncep, 'pres')
                    forcing_ncep.domain_cols_alt = length(forcing_ncep.rhum.lon);
                    forcing_ncep.domain_rows_alt = length(forcing_ncep.rhum.lat);
                end
                % Convert the small subdomain into cartesian coordinates. We need
                % to do this twice because some of the NCEP data are on different
                % grids (e.g. sea level pressure, relative humidity etc.).
                tmpZone = regexpi(inputConf.utmZone,'\ ','split');
                [tmpLon, tmpLat] = meshgrid(forcing_ncep.lon, forcing_ncep.lat);
                [forcing_ncep.x, forcing_ncep.y] = wgs2utm(tmpLat(:), tmpLon(:), str2double(char(tmpZone{1}(1))), char(tmpZone{1}(2)));
                if isfield(forcing_ncep, 'rhum')||isfield(forcing_ncep, 'pres')
                    [tmpLon2, tmpLat2] = meshgrid(forcing_ncep.rhum.lon, forcing_ncep.rhum.lat);
                    [forcing_ncep.xalt, forcing_ncep.yalt] = wgs2utm(tmpLat2(:), tmpLon2(:), str2double(char(tmpZone{1}(1))), char(tmpZone{1}(2)));
                end
                clear tmpLon tmpLat tmpLon2 tmpLat2 tmpZone
                % Create arrays of the x and y positions.
                forcing_ncep.x = reshape(forcing_ncep.x, forcing_ncep.domain_rows, forcing_ncep.domain_cols);
                forcing_ncep.y = reshape(forcing_ncep.y, forcing_ncep.domain_rows, forcing_ncep.domain_cols);
                if isfield(forcing_ncep, 'rhum')||isfield(forcing_ncep, 'pres')
                    forcing_ncep.xalt = reshape(forcing_ncep.xalt, forcing_ncep.domain_rows_alt, forcing_ncep.domain_cols_alt);
                    forcing_ncep.yalt = reshape(forcing_ncep.yalt, forcing_ncep.domain_rows_alt, forcing_ncep.domain_cols_alt);
                end
                [forcing_ncep.lon, forcing_ncep.lat] = meshgrid(forcing_ncep.lon, forcing_ncep.lat);
                forcing_ncep = rmfield(forcing_ncep, {'domain_rows', 'domain_cols'});
                if isfield(forcing_ncep, 'rhum')||isfield(forcing_ncep, 'pres')
                    forcing_ncep = rmfield(forcing_ncep, {'domain_rows_alt', 'domain_cols_alt'});
                end
                fprintf('Saving NCEP forcing from OPeNDAP server database...')
                save('forcing_ncep.mat','forcing_ncep','-v7.3','-nocompression');
                fprintf('Done\n')
                
                % % Have a look at some data.
                %{
                for i=1:240%size(forcing_ncep.air.data, 3)
                    figure1 = figure(1);
                    ax = axes('Parent',figure1);
                    clf
                    %           1      2     3      4     ???5      6     7     8
                    varb = {'uwnd','vwnd','air','rhum','prate','pres'};
                    %value = sqrt(forcing_ncep.(varb{1}).data(:, :, i).^2 + forcing_ncep.(varb{2}).data(:, :, i).^2);
                    value = forcing_ncep.(varb{3}).data(:,:,i);
                    [X, Y] = meshgrid(forcing_ncep.lon, forcing_ncep.lat);
                    s = pcolor(forcing_ncep.lon, forcing_ncep.lat, value);
                    shading flat
                    s.FaceColor = 'interp';
                    axis('equal','tight')
                    ylabel('Latitude (degree)','FontSize',12);
                    xlabel('Longtitude (degree)','FontSize',12);
                    pause(0.01);
                end
                clear ans ax figure1 i s value varb
                plot(squeeze(forcing_ncep.nswrs.data(2,4,1:36)));
                hold on;
                plot(squeeze(forcing_ncep.nlwrs.data(2,4,1:36)));
                hold on;
                plot(squeeze(forcing_ncep.shtfl.data(2,4,1:36)));
                hold on;
                plot(squeeze(forcing_ncep.lhtfl.data(2,4,1:36)));
                % The result of reanalysis2 shows the sensible and 
                % latent is totaly wrong.
                % reanalysis2 shtfl and lhtfl can not be used.
                %}
            elseif develop_mode == 2
                fprintf('Loading NCEP forcing from the local database...');
                load('forcing_ncep.mat');
                fprintf('Done!\n');
            end %developmode1,elseif 2 ,end
            % Interpolate the data onto the FVCOM unstructured grid.
            interpfields = {'time', 'lon', 'lat', 'x', 'y',...
          'nlwrs','nswrs','lhtfl','shtfl','dswrf', 'dlwrf'... 
                     };
                     % 'nlwrs','nswrs','lhtfl','shtfl','dswrf', 'dlwrf'...
               %    'uwnd', 'vwnd',...    %    'air', 'rhum','pres',...
             %       'prate', 'pres', 'air', 'rhum','pevpr'
             %       'uswrf', 'ulwrf', 'dswrf', 'dlwrf',...
            %interpfields = {'time', 'lon', 'lat', 'x', 'y','dswrf'};
            forcing_ncep_interp = grid2fvcom(Mobj,inputConf.doForcing,interpfields, forcing_ncep,...
                'add_elems', false);
            
           % forcing_ncep_interp.prate
            %size(forcing_ncep_interp.prate.node)
            %forcing_ncep_interp
        end %developmode3かそれ以外かの分岐

        %write out to nc
        
            write_FVCOM_forcing_calculated(Mobj, ...
            fullfile(inputConf.outbase,[inputConf.casename]),...
            forcing_ncep_interp, ...
            [inputConf.doForcing, ' atmospheric forcing data'],...
            inputConf.FVCOM_version, ...
            'floattime', true,...
            'julian', true);

        % (Mobj, fileprefix, data, infos, fver, varargin)

             save('Mobj_03.mat','Mobj','-v7.3','-nocompression');
             save('forcing_ncep_interp.mat','forcing_ncep_interp','-v7.3','-nocompression');
    

        % elseif strcmpi(inputConf.doForcing, 'NCEP-CALCULATED')
        % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.uwnd.node(1,:));
        % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.vwnd.node(1,:));
        % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.prate.node(1,:));
        % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.evap.node(1,:));
        % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.nswrs.node(1,:));
        % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.dlwrf.node(1,:));
        % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.air.node(1,:));
        % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.pres.node(1,:));
        % plot(datetime(forcing_interp_calculated.time+678942,'ConvertFrom','datenum','TimeZone','Asia/Tokyo'),forcing_interp_calculated.rhum.node(1,:));
        % write_FVCOM_forcing_calculated(Mobj, ...
        %     fullfile(inputConf.outbase,inputConf.casename),...
        %     forcing_interp_calculated,...
        %     [inputConf.doForcing, 'atmospheric forcing data'],...
        %     inputConf.FVCOM_version, ...
        %     'floattime', true,...
        %     'julian', true);
        % fileprefix=fullfile(inputConf.outbase,inputConf.casename);
        % data=forcing_interp_calculated;
        % infos=[inputConf.doForcing, 'atmospheric forcing data'];
        % fver=inputConf.FVCOM_version;
    %inputConf.doForcing = 'GWO';
    %inputConf.forceMJD = inputConf.startDateMJD:inputConf.dateforcing:inputConf.endDateMJD;
    %develop_mode = 1;
    elseif strcmpi(inputConf.doForcing, 'GWO')
        clear forcing_gwo forcing_gwo_interp;
    end
end   %meteo end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputConf.dateforcing = 1/24;
inputConf.doForcing = 'GWO';
inputConf.forceMJD = inputConf.startDateMJD:inputConf.dateforcing:inputConf.endDateMJD;
if strcmpi(inputConf.doForcing, 'GWO')
    inputConf.forceMJD = inputConf.startDateMJD:inputConf.dateforcing:inputConf.endDateMJD;
elseif strcmpi(inputConf.doForcing, 'NCEP')
    inputConf.forceMJD = [inputConf.startDateMJD, inputConf.endDateMJD];
end

develop_mode = 1;
if activate_Meteo==1
 %   if strcmpi(inputConf.doForcing, 'NCEP')

%    elseif strcmpi(inputConf.doForcing, 'GWO')
        %日本の気象庁のデータ
        % Get the surface heating data.
        % Use the GWO data to get the following parameters:
        % - u wind component (uwnd)[m/s]                       : for wind
        % - v wind component (vwnd)[m/s]                       : for wind
        % - Precipitation rate (prate)[m/s]                    : for precipitation
        % - Solar Radiation Flux (dsw) [W/m^2]                 : for heat flux
        % - Sea level pressure (pres) [Pa]                     : for air pressure
        % - Air temperature at ? m (air) [deg C]               : for temperature
        % - Relative Humidity on Pressure Levels (rhum) [%]    : for relative humidity
        % - Cloud coverage (cloud)

        % Generate the following variables for FVCOM net sea surface heat flux and sea
        % surface shortwave radiation. For this purpose, heat_calculated and custom
        % created dlw_calculated options are actived.
        %  &NML_SURFACE_FORCING
        %  WIND_ON                         = T,
        %  WIND_TYPE                       = 'speed',
        %  WIND_FILE                       = 'tokyobay_hfx.nc', 
        %  &NML_HEATING_CALCULATED
        %  HEATING_CALCULATE_ON            = T,
        %  HEATING_CALCULATE_TYPE          = 'flux',
        %  HEATING_CALCULATE_FILE          = 'tokyobay_hfx.nc',
        %  COARE_VERSION                   = 'BULKALGORITHM',
        %  The air pressure is not turnning on.
        %  This is used for non calculated heat flux and hurrican model.
        if develop_mode == 3
            fprintf('Loading Model objet file...\n')
            load('varb/Mobj_03.mat');
            % load('forcing_gwo.mat');
            load('forcing_gwo_interp.mat');

        else
            %inputConf.forceMJD = [inputConf.startDateMJD, inputConf.endDateMJD+1 ]
            Mobj.gwo.time = inputConf.forceMJD';
            %Mobj.gwo.time = Mobj.gwo.time(1:length(inputConf.forceMJD')-2);
    
            if develop_mode == 1
                clear forcing_gwo_interp
                fprintf('Start GWO data reading,develop_mode = 1');
                forcing_gwo = get_GWO_forcing(Mobj.gwo.time);
                UTMzone = regexpi(inputConf.utmZone,'\ ','split');
                for s = 1:length(forcing_gwo.lon)
                    [forcing_gwo.x(s),forcing_gwo.y(s),~,~] = wgs2utm(...
                        forcing_gwo.lat(s),forcing_gwo.lon(s),...
                        str2double(char(UTMzone{1}(1))),char(UTMzone{1}(2)));
                end
                [forcing_gwo.x,forcing_gwo.y] = meshgrid(forcing_gwo.x,forcing_gwo.y);
                clear s UTMzone
                save('forcing_gwo.mat','forcing_gwo','-v7.3','-nocompression');
                save('varb/Mobj_03.mat','Mobj','-v7.3','-nocompression');
            elseif develop_mode == 2
                fprintf('Loading GWO forcing from the local database...');
                load('forcing_gwo.mat');
                fprintf('Done!\n');
            end
            %size(forcing_gwo.air.data)
            % air

            interpfields = {'time', 'lon', 'lat', 'x', 'y',...
                'air','cld','dsw','dlwrf','hum','prs','rin'};
            % vars = interpfields; data = forcing;
            forcing_gwo_inter_node = grid2fvcom(Mobj,inputConf.doForcing,interpfields, forcing_gwo,...
                'add_elems', false);
            %dlwrf = readmatrix('C:/Users/ishid/Github/00_data/dlwrf_2020.csv');
            %for i = 1:3210
            %    forcing_gwo_inter_node.dlwrf.node(:,i) =  dlwrf(2:length(dlwrf),2)
            %end
            save('forcing_gwo_inter_node.mat','forcing_gwo_inter_node','-v7.3','-nocompression');
            clear forcing_gwo_inter_node;
            %{
            % cld
            interpfields = {'time', 'lon', 'lat', 'x', 'y',...
                'cld'};
            % vars = interpfields; data = forcing;
            forcing_gwo_inter_cld = grid2fvcom(Mobj,inputConf.doForcing,interpfields, forcing_gwo,...
                'add_elems', false);
            save('forcing_gwo_inter_cld.mat','forcing_gwo_inter_cld','-v7.3','-nocompression');
            clear forcing_gwo_inter_cld;
            
            % dsw
            interpfields = {'time', 'lon', 'lat', 'x', 'y',...
                'dsw'};
            % vars = interpfields; data = forcing;
            forcing_gwo_inter_dsw = grid2fvcom(Mobj,inputConf.doForcing,interpfields, forcing_gwo,...
                'add_elems', false);
            save('forcing_gwo_inter_dsw.mat','forcing_gwo_inter_dsw','-v7.3','-nocompression');
            clear forcing_gwo_inter_dsw;

                % dlwrf
            interpfields = {'time', 'lon', 'lat', 'x', 'y',...
                'dlwrf'};
            % vars = interpfields; data = forcing;
            forcing_gwo_inter_dlwrf = grid2fvcom(Mobj,inputConf.doForcing,interpfields, forcing_gwo,...
                'add_elems', false);
            save('forcing_gwo_inter_dlwrf.mat','forcing_gwo_inter_dlwrf','-v7.3','-nocompression');
            clear forcing_gwo_inter_dlwrf;

            % hum
            interpfields = {'time', 'lon', 'lat', 'x', 'y',...
                'hum'};
            % vars = interpfields; data = forcing;
            forcing_gwo_inter_hum = grid2fvcom(Mobj,inputConf.doForcing,interpfields, forcing_gwo,...
                'add_elems', false);
            save('forcing_gwo_inter_hum.mat','forcing_gwo_inter_hum','-v7.3','-nocompression');
            clear forcing_gwo_inter_hum;
            
            % prs
            interpfields = {'time', 'lon', 'lat', 'x', 'y',...
                'prs'};
            % vars = interpfields; data = forcing;
            forcing_gwo_inter_prs = grid2fvcom(Mobj,inputConf.doForcing,interpfields, forcing_gwo,...
                'add_elems', false);
            save('forcing_gwo_inter_prs.mat','forcing_gwo_inter_prs','-v7.3','-nocompression');
            clear forcing_gwo_inter_prs;
            %}
            % wnd
            interpfields = {'time', 'lon', 'lat', 'x', 'y',...
                'uwnd','vwnd'};
            % vars = interpfields; data = forcing_gwo;
            forcing_gwo_inter_wnd = grid2fvcom(Mobj,inputConf.doForcing,interpfields, forcing_gwo,...
                'add_elems', true);
            save('forcing_gwo_inter_wnd.mat','forcing_gwo_inter_wnd','-v7.3','-nocompression');
            clear forcing_gwo_inter_wnd;
            %{
            % rin
            interpfields = {'time', 'lon', 'lat', 'x', 'y',...
                'rin'};
            % vars = interpfields; data = forcing;
            forcing_gwo_inter_rin = grid2fvcom(Mobj,inputConf.doForcing,interpfields, forcing_gwo,...
                'add_elems', false);
            save('forcing_gwo_inter_rin.mat','forcing_gwo_inter_rin','-v7.3','-nocompression');
            clear forcing_gwo_inter_rin;
            
            % evpr
            interpfields = {'time', 'lon', 'lat', 'x', 'y',...
                'evpr'};
            % vars = interpfields; data = forcing;
            forcing_gwo_inter_rin = grid2fvcom(Mobj,inputConf.doForcing,interpfields, forcing_gwo,...
                'add_elems', false);
            save('forcing_gwo_inter_evpr.mat','forcing_gwo_inter_evpr','-v7.3','-nocompression');
            clear forcing_gwo_inter_evpr;
            %}

        end
    
 %   end
            clear forcing_gwo;
            clear interpfields;
            %{
            load('forcing_gwo_inter_air.mat'); forcing_gwo_interp.air = forcing_gwo_inter_air.air; clear forcing_gwo_inter_air;
            load('forcing_gwo_inter_cld.mat'); forcing_gwo_interp.cld = forcing_gwo_inter_cld.cld; clear forcing_gwo_inter_cld;
            load('forcing_gwo_inter_dsw.mat'); forcing_gwo_interp.dsw = forcing_gwo_inter_dsw.dsw; clear forcing_gwo_inter_dsw;
            load('forcing_gwo_inter_hum.mat'); forcing_gwo_interp.hum = forcing_gwo_inter_hum.hum; clear forcing_gwo_inter_hum;
            load('forcing_gwo_inter_prs.mat'); forcing_gwo_interp.prs = forcing_gwo_inter_prs.prs; clear forcing_gwo_inter_prs;
            load('forcing_gwo_inter_rin.mat'); forcing_gwo_interp.rin = forcing_gwo_inter_rin.rin; clear forcing_gwo_inter_rin;
            %load('forcing_gwo_inter_evpr.mat'); forcing_gwo_interp.evpr = forcing_gwo_inter_evpr.evpr; clear forcing_gwo_inter_evpr;
            load('forcing_gwo_inter_dlwrf.mat'); forcing_gwo_interp.dlwrf = forcing_gwo_inter_dlwrf.dlwrf; clear forcing_gwo_inter_dlwrf;
            %}
            load('forcing_gwo_inter_wnd.mat'); 
            load('forcing_gwo_inter_node.mat'); 
            forcing_gwo_interp.uwnd =forcing_gwo_inter_wnd.uwnd; 
            forcing_gwo_interp.vwnd =forcing_gwo_inter_wnd.vwnd; 
            forcing_gwo_interp.time = forcing_gwo_inter_wnd.time;
            forcing_gwo_interp.lon  = forcing_gwo_inter_wnd.lon; 
            forcing_gwo_interp.lat  = forcing_gwo_inter_wnd.lat; 
            forcing_gwo_interp.x    = forcing_gwo_inter_wnd.x; 
            forcing_gwo_interp.y    = forcing_gwo_inter_wnd.y;  
            forcing_gwo_interp.air =forcing_gwo_inter_node.air;
            forcing_gwo_interp.cld =forcing_gwo_inter_node.cld;
            forcing_gwo_interp.dsw =forcing_gwo_inter_node.dsw;
            forcing_gwo_interp.dlwrf =forcing_gwo_inter_node.dlwrf;
            forcing_gwo_interp.hum =forcing_gwo_inter_node.hum;
            forcing_gwo_interp.prs =forcing_gwo_inter_node.prs;
            forcing_gwo_interp.rin =forcing_gwo_inter_node.rin;

            clear forcing_gwo_inter_wnd;
            
            save('forcing_gwo_interp.mat','forcing_gwo_interp','-v7.3','-nocompression');
        %end
        
        write_FVCOM_gwo_forcing(Mobj, ...
            forcing_gwo_interp, ...
            fullfile(inputConf.outbase,[inputConf.casename]),...
            [inputConf.doForcing, ' atmospheric forcing data'],...
            inputConf.FVCOM_version, ...
            'floattime', true,...
            'julian', true);


        clear forcing_gwo forcing_gwo_interp;
        
    % Have a look at some data.
%{
    for i=1:240%size(forcing.air.data, 3)
        figure1 = figure(1);
        ax = axes('Parent',figure1);
        clf
        %           1      2     3     4     5     6     7     8
        varb = {'uwnd','vwnd','air','hum','rin','prs','cld','dsw'};
        %value = sqrt(forcing_gwo.(varb{1}).data(:, :, i).^2 + forcing_gwo.(varb{2}).data(:, :, i).^2);
        value = forcing_gwo.(varb{8}).data(:,:,i);
        s = pcolor(forcing_gwo.lon, forcing_gwo.lat, value');
        shading flat
        s.FaceColor = 'interp';
        axis('equal','tight')
        ylabel('Latitude (degree)','FontSize',12);
        xlabel('Longtitude (degree)','FontSize',12);
        pause(0.01);
    end
    clear ans ax figure1 i s value varb
%}

end   %meteo end




fprintf('All done!\n')

