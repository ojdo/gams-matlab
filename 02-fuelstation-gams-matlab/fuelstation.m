% fuelstation.m
% sets
t      = GAMS.set('t', 1:8760);
i      = GAMS.set('i', {'pv', 'windon', 'windoff'});

% parameters
cs     = GAMS.param('cs',100); % cost of storage (EUR/MWh)
c      = GAMS.param('c',[3000 1500 2500],i.uels); % cost of plant (EUR/MWh)
demand = GAMS.param('demand',rand(8760,1),t.uels);

% renewable timeseries
values = [ ... 
    min(max(0, sin((1:8760)'/24*3.14/2).^4+0.15*randn(8760,1)), 1), ...
    min(max(0, rand(8760,1)), 1), ...
    min(max(0, rand(8760,1).^0.25), 1) ];
onset = [ t.uels i.uels ];
cf = GAMS.param('cf', values, onset);
clear values onset;

% write to GDX file
GAMS.putGDX('input.gdx',t,i,c,cs,demand,cf);

% run GAMS model
g = GAMS(struct('model','fuelstation.gms'));
g.run; % executes "gams.exe fuelstation.gms -GDX=result.gdx"

% read result variable x if run successful
if g.status == 0
	x = GAMS.getGDX('result.gdx','x');
	x = GAMS.rectify(x, i.uels);

    % create bar chart of installed plant capacities
	bar(1000*x.val);
	set(gca,'XTickLabel',x.uels{1});
	ylabel('Installed capacity (kW)');
end
