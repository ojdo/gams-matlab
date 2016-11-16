classdef FS < GAMS
    properties
        % input data
        set_t      % timesteps
        set_i      % technologies
        db_cs      % cost of storage (€/MWh)
        db_c       % cost of plant (€/MWh)
        ts_demand  % demand timeseries (1)
        ts_cf      % renewable input timeseries (1)
        
        % result data
        Z       % total cost (k€)
        X       % plant sizes per technology (MW)
        S       % storage size (MWh)
    end
    
    methods
        function obj = FS()
            % Call GAMS constructor
            obj = obj@GAMS((struct('model','fuelstation.gms')));
            
            % Set values for input data
            obj.set_t = GAMS.set('t', 1:8760);
            obj.set_i = GAMS.set('i', {'pv', 'windon', 'windoff'});
            
            obj.db_cs     = GAMS.param('cs',100); 
            obj.db_c      = GAMS.param('c',[3000 1500 2500],obj.set_i.uels); 
            obj.ts_demand = GAMS.param('demand',rand(8760,1),obj.set_t.uels);
            
            values = [ ... 
                min(max(0, sin((1:8760)'/24*3.14/2).^4+0.15*randn(8760,1)), 1), ...
                min(max(0, rand(8760,1)), 1), ...
                min(max(0, rand(8760,1).^0.25), 1) ];
            onset = [ obj.set_t.uels obj.set_i.uels ];
            obj.ts_cf = GAMS.param('cf', values, onset);
        end
        
        function writeInputs(obj)
            GAMS.putGDX('input.gdx', obj.set_t, obj.set_i, ...
                obj.db_cs, obj.db_c, obj.ts_demand, obj.ts_cf);
        end
        
        function readResults(obj)
            obj.Z = GAMS.getGDX(obj.path.result, 'z');
            obj.X = GAMS.getGDX(obj.path.result, 'x');
            obj.S = GAMS.getGDX(obj.path.result, 's');
            
            obj.X = GAMS.rectify(obj.X, obj.set_i.uels);
        end
        
        function plot(obj)
            bar(1000*obj.X.val);
            set(gca,'XTickLabel',obj.X.uels{1});
            ylabel('Installed capacity (kW)');
            grid on;
        end
    end
end
