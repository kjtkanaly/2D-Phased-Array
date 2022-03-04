
% Triangluar Phased Array
% [1] - Radar Handbook - Merrill Skolnik

%% Constants and Array Parameters

c					= physconst('lightspeed');	% Speed of light
frequency	= 5.7e9;								% Operating Frequency
lambda		= c/frequency;						% Operating wavelength

ArrayXDim	= 15;									% Number of elements in the X-Dim
ArrayYDim   = ArrayXDim*2;					% Number of elements in the Y-Dim, multipled by two to account for the triangular spacing

dx = 28.45e-3;											% Spacing between elements in the X-Dim
dy = 32.77e-3/2;										% Spacing between elements in the Y-Dim, divided by two to account for the triangular spacing

u =  -1:0.001:1;										% Observation locations in the X-Dim
v = -1:0.001:1;											% Observation locations in the Y-Dim
Tx = (2*pi / lambda) .* dx .* u;					% The phase difference between elements in the X-Dim due to position
Ty = (2*pi / lambda) .* dy .* v.';				% The phase difference between elements in the Y-Dim due to position

Us = 0;													% The beam steering X position
Vs = 0;														% The beam steering Y position
Txs = (2*pi / lambda) .* dx .* Us;				% The element-element phase shifting for the given beam steering X position
Tys = (2*pi / lambda) .* dy .* Vs;				% The element-element phase shifting for the given beam steering Y position

% Type of Element Tapering 
%     0 - Uniform,				1 - Vertical Fan Beam,		2 - Cosine Tapering, 
%     3 - HanningTaper,	4 - Symmetrical Uniform,	
taperType = 0;

% Desired Plots
%		0 - None,		1 - Az/El Principle Plane Cuts,	2- 3-D Surf Plot,		3 - All Plots,
desiredPlots = 1;

deleteProgress = '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b';

%% Tapering Matrix

% Uniform Illumination
if (taperType == 0)

	Amn = ones(ArrayYDim, ArrayXDim);
	Amn(1:2:ArrayYDim-1,1:2:ArrayXDim-1) = 0;
	Amn(2:2:ArrayYDim,2:2:ArrayXDim) = 0;

	taperTypeLabel = "Uniform";


% Vertical Fan Illumination
elseif (taperType == 1)

	Amn(2:2:ArrayYDim, 6:2:8) = 1;

	taperTypeLabel = "Vertical Fan Beam";


% Cosine Tapering
elseif (taperType == 2)

	padding = [0, 8];	% Set to 8 for a 32x32 array

	temp = cosineTaperTriangleSpacing(ArrayXDim - padding(1), ArrayYDim - padding(2));

	Amn = zeros(ArrayYDim, ArrayXDim);
	Amn(padding(2)/2 + 1:padding(2)/2 + size(temp,1),padding(1)/2 + 1: padding(1)/2 + size(temp,2)) = temp;

	taperTypeLabel = "Cosine";


% Hanning Taper
elseif (taperType == 3)

	padding = [0, 8];	% Set to [0, 8] for a 32x32 array to be symmetrical in Az and El

	temp = hannTriangleSpacing(ArrayXDim - padding(1), ArrayYDim - padding(2));

	Amn = zeros(ArrayYDim, ArrayXDim);
	Amn(padding(2)/2 + 1:padding(2)/2 + size(temp,1),padding(1)/2 + 1: padding(1)/2 + size(temp,2)) = temp;
	
	taperTypeLabel = "Hanning";


% Symmetrical Uniform,	
elseif (taperType == 4)

	Amn = ones(ArrayYDim, ArrayXDim);
	Amn(2:2:ArrayYDim,1:2:ArrayXDim-1) = 0;
	Amn(1:2:ArrayYDim-1,2:2:ArrayXDim) = 0;
	Amn(1:4, 1:ArrayXDim) = 0;
	Amn(ArrayYDim-3:ArrayYDim,1:ArrayXDim) = 0;

	taperTypeLabel = "Symmetrical Uniform";

end

%% 2-D Array Factor Calculation [1]

AF = zeros(length(v), length(u));

currentPercent = 0;
fprintf('Calculate Array Factor Progress:%6.2f%%',currentPercent);

% Iterate through the elements in the array
for m = 1:ArrayXDim
	for n = 1:ArrayYDim

		AF = AF + (  abs(Amn(n,m)) .* exp(1i.*( m.*(Tx - Txs) + n.*(Ty - Tys) ))  );

	end

	currentPercent = ( m / ArrayXDim) * 100;
	
	fprintf(deleteProgress);
	fprintf('Calculate Array Factor Progress:%6.2f%%',currentPercent);
end

currentPercent = 100;
fprintf(deleteProgress);
fprintf('Calculate Array Factor Progress:%6.2f%%\n',currentPercent);

% AF = flip(rot90(AF));

% Calculating Azimuth BeamWidth
index	= find(mag2db(abs(AF(round(length(Tx)/2),:) ./ max(AF,[],'all'))) >= -3);
if ~isempty(index)
	thetaAz = rad2deg(asin(u(index)));
	thetaAz = abs(thetaAz(end) - thetaAz(1));
else
	thetaAz = "N/A";
end

% Calculating Elevation BeamWidth
index	= find(mag2db(abs(AF(:,round(length(Tx)/2)) ./ max(AF,[],'all'))) >= -3);
if ~isempty(index)
	thetaEl = rad2deg(asin(v(index)));
	thetaEl = abs(thetaEl(end) - thetaEl(1));
else
	thetaEl = "N/A";
end

% Rounding the Az and El beamwidths for visual aesthetic
format bank;
thetaAz = round(thetaAz*100)/100;
thetaEl  = round(thetaEl*100)/100;
format default;

%% Plotting

if (desiredPlots == 2 || desiredPlots == 3)

	fig = figure;
	fig.Name = "3-D Array Factor";
	fig.NumberTitle = 'off';
	srf = surf(u,v,mag2db(abs(AF ./ max(AF,[],'all'))));
	srf.FaceColor = "Interp";
	srf.LineStyle = ":";
	srf.LineWidth = 0.25;
	xlabel('U','FontSize',18)
	ylabel('V','FontSize',18)
	zlabel('dB','FontSize',18)
	title([strcat(string(ArrayXDim),"x",string(ArrayYDim/2)," Triangular Array Factor") strcat("Tapering: ", taperTypeLabel, " | Az/El Beamwidth: ", string(thetaAz), "/", string(thetaEl),"\circ")],'FontSize',18)
	zlim([-80 0])
	caxis([-80 -20])
	colorbar

end

if (desiredPlots == 1 || desiredPlots == 3)

	fig = figure;
	fig.Name = "Principal Plane Cuts";
	fig.NumberTitle = 'off';
	subplot(1,2,1)
	plot(u,mag2db(abs(AF(round(length(Tx)/2),:) ./ max(AF,[],'all'))),'LineWidth',1.5)
	% plot(rad2deg(asin(u)),mag2db(abs(AF(501,:) ./ max(AF,[],'all'))))
	grid on;
	ylim([-60 0])
	xlim([-1 1])
	xlabel('U','FontSize',18)
	ylabel('dB','FontSize',18)
	title(strcat("Azimuth Cut: \theta_3 = ", string(thetaAz)),'FontSize',16)
	
% 	t = text(0,0,'\theta_{3} = ' + string(thetaAz) + '\circ');
% 	t.Position = [-0.959537572254336,-2.147619047619015,1.4e-14];
% 	t.FontSize = 16;
% 	t.EdgeColor = 'k';
% 	t.LineWidth = 1;
	
	subplot(1,2,2)
	plot(v,mag2db(abs(AF(:,round(length(Tx)/2)) ./ max(AF,[],'all'))),'LineWidth',1.5)
	% plot(rad2deg(asin(v)),mag2db(abs(AF(:,501) ./ max(AF,[],'all'))))
	grid on;
	ylim([-60 0])
	xlim([-1 1])
	xlabel('V','FontSize',18)
	ylabel('dB','FontSize',18)
	title(strcat("Elevation Cut: \theta_3 = ", string(thetaEl)),'FontSize',16)

% 	t = text(0,0,'\theta_{3} = ' + string(thetaEl) + '\circ');
% 	t.Position = [-0.959537572254336,-2.147619047619015,1.4e-14];
% 	t.FontSize = 16;
% 	t.EdgeColor = 'k';
% 	t.LineWidth = 1;

	sgtitle(strcat(string(ArrayXDim),"x",string(ArrayYDim/2)," Triangular Array Factor"),'FontSize',18)

end


%% Useful Functions


% [2-D Cosine Angles and Normalizing Values come from J]
function taper = cosineTaperTriangleSpacing(cols,rows)

	x = linspace(deg2rad(-80),deg2rad(80),cols);     % This should be the length of the number of horizontal elements
	y = linspace(deg2rad(-80),deg2rad(80),rows);   % This should be the length of twice the number of vertical elements
                                                	
	map = cos(y.') .* cos(x);					% Create the 2-D cosine
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Normalizing the Cosine Taper map to fit the old cosine taper
	P1 = min(map,[],'all');	% Min value in the 2-D cosine map
	P1Prime = 1175;			% Min taper value from previous version of code
	P2Prime = 32000;		% Max taper value from previous version of code
	maxValue = 32767;		% Max Value of an unsigned in

	map = map - P1;
	map = map ./ max(map,[],'all');
	map = map .* (P2Prime - P1Prime);
	map = map + P1Prime;
	map = round(map) ./ maxValue;			% Round before normalzing so that the percentages are whole bytes
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Sampling the taper map in a fashion that reflects the array's actual
	% checkerboard distribution of elements
	taper						= map;
	taper(1:2:rows-1,1:2:cols-1)	= 0;
	taper(2:2:rows,    2:2:cols)		= 0;

end


% Hanning 2-D Tapering Function with Triangular Sampling
% [Equation from Matlab's Hann Fx and Normalizing Values come from J]
function taper = hannTriangleSpacing(cols,rows)
	
	N = 1;
	
	x = linspace(0.1, 0.9,cols);
	y = linspace(0.1, 0.9,rows);
	
	map = 0.5.*(1 - cos(2.*pi.*(x./N))) .* (0.5.*(1 - cos(2.*pi.*(y./N)))).';

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Normalizing the Cosine Taper map to fit the old cosine taper
	P1 = min(map,[],'all');	% Min value in the 2-D cosine map
	P1Prime = 1175;			% Min taper value from previous version of code
	P2Prime = 32000;		% Max taper value from previous version of code
	maxValue = 32767;		% Max Value of an unsigned in

	map = map - P1;
	map = map ./ max(map,[],'all');
	map = map .* (P2Prime - P1Prime);
	map = map + P1Prime;
	map = round(map) ./ maxValue;			% Round before normalzing so that the percentages are whole bytes

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Sampling the taper map in a fashion that reflects the array's actual
	% checkerboard distribution of elements
	taper						= map;
	taper(1:2:rows-1,1:2:cols-1)	= 0;
	taper(2:2:rows,    2:2:cols)		= 0;

end












