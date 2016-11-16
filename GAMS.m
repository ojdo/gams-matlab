classdef GAMS < handle
%GAMS Interface class to call GAMS models from within Matlab
%   This class can be used for quick scenario generation and result
%   analysis for GAMS models. Call 'help GAMS.method_name' for information
%   on individual methods. Get started by: set, param, getGDX, putGDX and
%   getXLS. 
%   
%   Methods
%       g = GAMS                  constructor, set paths
%       g.setPath                 set paths to GAMS, model, result
%       g.setMail                 set preferences
%       g.run                     execute GAMS
%
%   Static methods
%     
%       GAMS.set                  create set
%       GAMS.param                create parameter
%
%       GAMS.getGDX               read a symbol from GDX file
%       GAMS.putGDX               write set or parameter to GDX file
%       GAMS.globGDX              read symbol from multiple GDX files
%       GAMS.getXLS               read entity from XLS file
%       GAMS.putXLS               write sets, params, variables to XLS file
%
%       GAMS.rectify              make data structure compliant to given uels
%       GAMS.merge                merge contents of two data structures into one
%       GAMS.sum                  summarize multi-dimensional data structures
%
%       GAMS.to_uels              creates uels {'1' '2' '3'} from [1:3]
%       GAMS.uel_to_ids           create ids from uels
%       GAMS.full_to_sparse       convert value matrix form
%       GAMS.sparse_to_full       convert value matrix from
%
%       GAMS.GDX2MDB              convert GDX file to MS Access MDB
%
%   Example
%       [set_node att_node db_node] = GAMS.getXLS('input.xlsx', 'Node');
%       GAMS.putGDX('node.gdx',set_node,att_node,db_node);
%       g = GAMS(struct('model','fuelstation.gms'));
%       g.run();
%       total_cost = GAMS.getGDX('result.gdx','z');
%

    properties
        path   % paths for GAMS binary and model file
        status % return code of last GAMS run
        result % stdout of last GAMS run
    end

    methods
        function obj = GAMS(path) % constructor initialises  variables
            %GAMS.GAMS constructs a GAMS object
            %
            %   Usage
            %       obj = GAMS(path)
            %
            %   Parameters
            %       path        struct with the following fields
            %           gams        path to GAMS executable ('gams.exe')
            %           model       path to model file ('model.gms')
            %           result      path to result file ('result.gdx')
            %                   missing fields are replaced by the default
            %                   values (in parenthesis above)
            %
            %   Returns
            %       obj         GAMS object
            %
            
            % create empty struct if called without input parameter
            if nargin == 0, path = struct(); end
            
            obj.setPath(path);
            
        end
        
        function [status, result] = run(obj, verbose)
            %GAMS.run calls GAMS and collects return code
            %   Executes 'gams model.gms', where both path to gams
            %   and to model file are taken from the objects property path.
            %
            %   Usage
            %       [status, result] = G.run(verbose)
            %
            %   Parameters
            %       verbose     (optional) integer, specifies how many
            %                   lines of the GAMS stdout to be displayed in
            %                   the command window
            %                   default: 0 (=off)
            %
            %   Returns
            %       status      Return code. Non-zero values indicate errors
            %       result      STDOUT GAMS output, different from *.lst!
            %
            if nargin < 2, verbose = 0; end
            
            % create run command
            if isfield(obj.path,'result')
                gams_command = sprintf('"%s" "%s" -GDX="%s" -lo=3\n', obj.path.gams, obj.path.model, obj.path.result);
            else
                gams_command = sprintf('"%s" "%s" -lo=3\n', obj.path.gams, obj.path.model);
            end
            
            % execute command
            tstart = tic;            
            [status, result] = system(gams_command);
            telapsed = toc(tstart);
            
            % catch status and result
            dlmwrite(strrep(obj.path.model,'.gms','.log'), result, 'delimiter', '', 'newline', 'pc');
            obj.result = regexp(result, '[\f\n\r]', 'split');
            obj.status = status;
            result = char(obj.result);
            
            % optional: display last lines of output
            if verbose > 0
                verbose = min(size(result,1)-1, verbose);
                disp(result(end-verbose:end,:)); 
            end

        end
        
        function setPath(obj, path)
            %GAMS.setPath sets path property
            %
            %   Usage
            %       g.setPath(path)
            %
            %   Parameters
            %       path        struct with the following fields
            %           gams        path to GAMS executable ('gams.exe')
            %           model       path to model file ('model.gms')
            %           result      path to result file ('result.gdx')
            %                   missing fields are replaced by the default
            %                   values (in parenthesis above)
            %
            %   Returns
            %       nothing
            %            
            
            % set up struct with default values
            defaultpath = struct(...
                'gams','gams.exe',...
                'model','model.gms',...
                'result','result.gdx');
            
            % now fill up missing values in path structure with defaults
            fns = fieldnames(defaultpath);
            for fn=fns'
                fn=char(fn); %#ok<FXSET>
                % checks, whether fn exists as field in parameter path
                if sum(strcmp(fieldnames(path),fn)) == 0
                    % if not, copy default value
                    path.(fn) = defaultpath.(fn);
                    
                    % warn if no model filename was given
                    if strcmp(fn,'model_file'), warning(WrongNumberOfArguments,['No model file given, default ''' defaultpath.model_file ''' is used.']); end
                end
            end
            
            % assign path structure to created object that is returned
            obj.path = path;
        end       
    end
    
    methods (Static)

        
        function Gset   = set(name, vals, onsets, form)
            %GAMS.set creates GAMS set data structure
            %   Used to create one- and multi-dimensional sets from MATLAB
            %   data structure
            %   Example
            %       g_sites = GAMS.set('Site', {'AT', 'DE'});
            %
            %   Usage
            %       S = GAMS.set(name, vals, onsets, form)
            %
            %   Parameters
            %       name        name of set (string)
            %       vals        list of values (cell array of elements;
            %                   element may be a string [one-dimensional]
            %                   or a cell array of strings [multi-dimensional])
            %                   alternative: logical incidence matrix
            %       onsets      list of domain sets for multi-dimensional
            %                   (cell array of cell arrays like vals)
            %       form        (optional) 'full' or 'sparse' (default: 'full')
            %                   try sparse when using extremely huge
            %                   datasets with very few existing
            %                   (onset-)combinations
            %
            %   Returns
            %       S.name      name of variable, equation, parameter
            %       S.type      'set'
            %       S.val       incidence matrix (nd-array)
            %       S.form      'full' or 'sparse'
            %       S.uels      cell array of dimension labels
            %       S.ids       structures with dimension labels as fieldnames
            %
            %   Advanced examples
            %       % multi-dimensional sets
            %       node = {'AT' 'DE' 'FR' 'ES'};
            %       edge = {{'AT' 'DE'} {'DE' 'FR'} {'FR' 'ES'}};
            %       g_node = GAMS.set('Node',node);
            %       g_edge = GAMS.set('Edge',edge,{node node});
            %       GAMS.putGDX('grid.gdx',g_node,g_edge);
            %
            %       % sets with numeric entries are possible, but must
            %       % be converted to cell arrays when used as onsets
            %       t = 1:5;  % main set
            %       t0 = 1;   % subset
            %       g_time = GAMS.set('Time',t);
            %       g_t0   = GAMS.set('T0',t0,g_time.uels); % domain check
            %       GAMS.putGDX('time.gdx',g_time,g_t0);
            %
            
            % sparse option checking
            if nargin < 4, form = 'full'; end
            
            if sum(strcmp({'full' 'sparse'}, form)) ~= 1
                error('GAMS:set:WrongForm', 'Optional rgument form must be ''full'' (default) or ''sparse''');
            end
            
            % convert numeric inputs to cell array of strings
            if isnumeric(vals)
                vals = reshape(vals,numel(vals),1);
                vals = cellstr(num2str(vals,'%-g'))';
            end
            
            if nargin > 2 % multi-dimensional sets or one-dimensional subsets
                
                % handle dimensions defined by number of onsets
                dimensionality = cellfun(@(x) size(x,2), onsets);
                if numel(dimensionality) == 1
                    dimensionality = [1 dimensionality];
                end
                
                % create onset ids
                ids = cell(1,length(onsets));
                for d=1:length(onsets), ids{d} = GAMS.uel_to_ids(onsets{d}); end
                
                if prod(dimensionality) > 1e8 || strcmp(form,'sparse')
                    % create sparse matrix
                    
                    % incidence_matrix
                    incidence_matrix = zeros(numel(vals), length(dimensionality));
                    
                    
                    for k = 1:numel(vals)
                        el = vals(k);
                        if iscell(el{1}) % if element is multi-dimensional
                            el = el{:};  % unpack cell array wrapping around the element
                        end
                        
                        indices = zeros(1,length(el));

                        for dim = 1:length(el) % for each dimension of the element
                            % determine subscripts of el in corresponding onset
                            matching_index = find(strcmp(onsets{dim},el{dim}));

                            % basic domain checking
                            if isempty(matching_index)
                                error(['GAMS.set: No matching ''' el{dim} ''' found in onset dimension ' num2str(dim) '!']);
                            end

                            indices(dim) = matching_index;
                        end
                        
                        incidence_matrix(k,:) = indices;
                    end
                    
                    Gset = struct(...
                        'name',name,...
                        'type','set',...
                        'val',incidence_matrix,...
                        'form','sparse',...
                        'dim',length(onsets),...
                        'uels',{onsets},...
                        'ids',{ids}...
                    );
                else
                    % create full incidence matrix
                
                    % prepare incidence matrix
                    incidence_matrix = zeros(dimensionality);

                    if islogical(vals)
                        % special case: vals is already an incidence matrix
                        incidence_matrix = double(vals); 
                    else
                        % now cycle through values and fill incidence matrix
                        for el = vals
                            if iscell(el{1}) % if element is multi-dimensional
                                el = el{:};  %#ok<FXSET> % unpack cell array wrapping around the element
                            end

                            indices = zeros(1,length(el));

                            for dim = 1:length(el) % for each dimension of the element
                                % determine subscripts of el in corresponding onset
                                matching_index = find(strcmp(onsets{dim},el{dim}));

                                % basic domain checking
                                if isempty(matching_index)
                                    error(['GAMS.set: No matching ''' el{dim} ''' found in onset dimension ' num2str(dim) '!']);
                                end

                                indices(dim) = matching_index;
                            end

                            % and, because Matlab does not support indexing by vector
                            % of subscripts, convert the vector to a numeric cell array 
                            % and use that to convert indices to a linear index using 
                            % sub2ind.
                            indices_as_cell = num2cell(indices);
                            indices_linear  = sub2ind(size(incidence_matrix),indices_as_cell{:});

                            incidence_matrix(indices_linear) = 1;
                        end
                    end

                
                    Gset = struct(...
                        'name',name,...
                        'type','set',...
                        'val',incidence_matrix,...
                        'form','full',...
                        'dim',length(onsets),...
                        'uels',{onsets},...
                        'ids',{ids}...
                    );
                end
                
            else % one-dimensional sets 

                % make every input a row vector of elements (=flat list)
                vals = reshape(vals,1,numel(vals));
                
                ids = GAMS.uel_to_ids(vals);
                
                % create set data structure
                Gset = struct(...
                    'name',name,...
                    'type','set',...
                    'val',ones(size(vals)),...
                    'form','full',...
                    'dim',1,...
                    'uels',{{vals}},...
                    'ids',{{ids}}...
                );
            end
        end
        
        function Gparam = param(name, vals, onsets, form)
            %GAMS.set creates GAMS set data structure
            %   Used to create one- and multi-dimensional sets from MATLAB
            %   data structure
            %   Example
            %       g_param = GAMS.param('myParam', 5);
            %
            %   Usage
            %       P = GAMS.param(name, vals, onsets, form)
            %
            %   Parameters
            %       name        name of parameter (string)
            %       vals        matrix of values 
            %       onsets      list of domain set(s)
            %                   (cell array of cell arrays)
            %       form        (optional) 'full' or 'sparse' (default: 'full')
            %                   try sparse when using extremely huge
            %                   datasets with very few non-zero
            %                   (onset-)combinations of values
            %
            %   Returns
            %       P.name      name of variable, equation, parameter
            %       P.type      'parameter'
            %       P.val       value matrix (full: nd-array, sparse: matrix)
            %       P.uels      cell array of dimension labels
            %       P.ids       structures with dimension labels as fieldnames
            %       P.form      'full' or 'sparse'
            %       P.dim       number of dimensions
            %        
            %
            %   More examples
            %       % scalar parameter
            %       scalar = GAMS.param('scalar',5);
            %
            %       % 1D parameter
            %       domain = {'a' 'b' 'c' 'd'};
            %       vals   = [ 1   2   3   4];
            %       param = GAMS.param('param',vals,{domain});
            %
            %       % 2D parameter on (node,atts)
            %       node = {'AT' 'DE' 'FR' 'ES'};
            %       atts = { 'demand' 'price' };
            %       vals = [ 100 5; 200 4; 300 2; 400 1 ];
            %       db_node = GAMS.param('db_node',vals,{node atts});
            %
            %       % 3D parameter on (node,node,atts)
            %       edge = {{'AT' 'DE'} {'DE' 'FR'} {'FR' 'ES'}};
            %       atts = { 'length' 'capacity' };
            %       vals = [ 600 1; 800 3; 200 10 ];
            %       db_edge = GAMS.param('db_edge',vals,{edge atts});
            %
            %   Last changed
            %       2011-08-26 13:00 GMT+2      added this documentation
            
            if nargin < 4, form = 'full'; end
            
            if sum(strcmp({'full' 'sparse'}, form)) ~= 1
                error('GAMS:param:WrongForm', 'Optional rgument form must be ''full'' (default) or ''sparse''');
            end
            
            if nargin == 2 % scalar parameer
                Gparam = struct(...
                    'name',name,...
                    'type','parameter',...
                    'val', vals, ...
                    'form', 'full', ...
                    'dim', 0 ...
                    );
                return
            end

            if ~iscell(onsets{1}{1}) % elementary onsets
                % arbitrary number of dimensions (<20) allowed
                value_matrix = vals;
                
                ids = cell(1,length(onsets));
                for d=1:length(onsets)
                    ids{d} = GAMS.uel_to_ids(onsets{d});
                end
                
                Gparam = struct(...
                    'name',name,...
                    'type','parameter',...
                    'val',value_matrix,...
                    'form','full',...
                    'dim',length(onsets),...
                    'uels',{onsets},...
                    'ids',{ids}...
                );
            else % multi-dimensional onsets{1}
                % then onsets{2} is required and must be elementary
                
                % initialise uels
                Ndim = length(onsets{1}{1});
                uels = cell(1,Ndim+1); % +1 for attribute dimension (from onsets{2})
                
                % declare sparse 2D value matrix
                [Nrows,Ncols] = size(vals);
                value_matrix = zeros(numel(vals),Ndim+1);
                
                % extract uels from onsets{1} by dimension, extract unique
                % identifiers and assign positional references (using 
                % ismember) used in value_matrix
                for k = 1:Ndim
                    uels_for_vals = cellfun(@(x) x{k},onsets{1},'UniformOutput',false);
                    uels{k} = unique(uels_for_vals);
                    
                    % assign positional references to uels
                    %  uels_for_vals = {'DE' 'AT' 'AT' 'DE' 'FR' }
                    %  uels{k} = {'AT' 'DE' 'FR}
                    %  --> loc = [2 1 1 2 3]
                    [~,loc] = ismember(uels_for_vals, uels{k});
                    
                    % repeat each entry of loc Ncols times (one time for 
                    % each attribute)
                    loc = repmat(loc, Ncols, 1);
                    value_matrix(:,k) = loc(:);
                end
                
                % increment Ndim by one to include the additional dimension
                % spawned by attributes, add them to the uels cell array
                % and add references (just repeating 1:Ncols for each row
                % of vals) to the value_matrix
                Ndim = Ndim + 1;
                uels{Ndim} = onsets{2};
                value_matrix(:, Ndim) = repmat(1:Ncols,1,Nrows);
                
                % finally: copy attribute values to value_matrix
                % reshape required for keeping right order (line by line,
                % left to right), while a simple vals(:) would concatenate 
                % column by column
                value_matrix(:, Ndim+1) = reshape(vals', [], 1);
                                
                % convert uels to ids
                ids = cell(1,length(uels));
                for d=1:length(uels)
                    ids{d} = GAMS.uel_to_ids(uels{d});
                end
                
                % create output data structure
                Gparam = struct(...
                    'name',name,...
                    'type','parameter',...
                    'val',value_matrix,...
                    'form','sparse',...
                    'dim',Ndim,...
                    'uels',{uels},...
                    'ids',{ids}...
                );
                
                % and convert to full form if desired (default behaviour)
                % and feasible (dimensionality in range)
                if strcmp(form, 'full')
                    if prod(cellfun(@(x) size(x,2), Gparam.uels)) < 1e8
                        Gparam = GAMS.sparse_to_full(Gparam);
                    elseif nargin > 3
                        % warn if full data structure was desired
                        warning('GAMS:param:TooBigForFull',['Data structure ''' name ''' too big for full value matrix, fallback to sparse output.']);
                    end
                end
            end
        end
        
        function Gdata = rectify(Gdata, uels)
            %GAMS.rectify makes variable or parameter conform to given uels
            %   Data delivered from GAMS.getGDX often has missing elements
            %   in some dimensions due to the sparse data structure that is
            %   handed back by rgdx. This helper functions automatically
            %   adds or removes values from Gdata in order to make it the
            %   size specified by the parameter uels. Missing entries are
            %   filled up with zeros in the value matrix, while superfluous
            %   entries are removed, issuing a warning.
            %       This function also sorts entries in value matrices, so 
            %   that uels for all dimensions match after rectification.
            %   Example
            %       estocon = GAMS.getGDX('result.gdx','e_sto_con');
            %       estoin  = GAMS.getGDX('result.gdx','e_sto_in');
            %       % estoin.uels{1} (timesteps) has missing elements, so
            %       % dimensions of value matrices do not match
            %       estoin = GAMS.rectify(estoin,estcon.uels);
            %       % now estoin.uels is identical to estocon.uels
            %
            %   Usage
            %       Gdata = GAMS.rectify(Gdata, uels)
            %
            %   Parameters
            %       Gdata   a GAMS set, parameter or variable
            %       uels    cell array of cell arrays of the desired uels;
            %               usually taken from (an)other Gdata object that 
            %               is deemed complete in all dimensions
            %   
            %   Returns
            %       Gdata   the rectified input parameter
            %
            
            if nargin ~= 2, error('GAMS:rectify:WrongNumberOfArguments','Wrong number of arguments.'); end
            if length(uels) ~= length(Gdata.uels), error('GAMS:rectify:DifferentDimensions','Gdata and uels must have same number of dimensions.'); end
            
            % initialize variables
            Ndim = Gdata.dim;
            [uels_to_copy, uels_to_remove, new_idx, old_idx, ids]  = deal(cell(1,Ndim));
            [dims, delete_count] = deal(zeros(1,Ndim)); % dimensionality vector, delete count
            
            % Loop through each dimension
            for d=1:Ndim
                % calculate size of new value matrix
                dims(d) = length(uels{d});
                
                % convert uels to ids
                ids{d} = GAMS.uel_to_ids(uels{d});
                
                % determine which uels from Gdata will be kept and removed
                uels_to_copy{d}      = intersect(Gdata.uels{d}, uels{d});
                uels_to_remove{d}    = setdiff(Gdata.uels{d}, uels{d});
                                
                % from uels_to_copy, derive indexes for both old (Gdata)
                % and new value matrix
                new_idx{d} = zeros(1,length(uels_to_copy{d}));
                for i=1:length(new_idx{d}), new_idx{d}(i) = find(strcmp(uels_to_copy{d}(i),uels{d})); end
                
                old_idx{d} = zeros(1,length(uels_to_copy{d}));
                for i=1:length(old_idx{d}), old_idx{d}(i) = find(strcmp(uels_to_copy{d}(i),Gdata.uels{d})); end
                
                % from uels_to_remove, count
                delete_count(d) = sum(cellfun(@(x) ~isempty(x),uels_to_remove{d}));
            end
            
            if sum(delete_count) > 0
                if sum(delete_count) == 1
                    warning('GAMS:rectify:EntriesRemoved',['1 uel is removed from data structure ''' Gdata.name '''.']); 
                else
                    warning('GAMS:rectify:EntriesRemoved',[num2str(sum(delete_count)) ' uels are removed from data structure ''' Gdata.name '''.']); 
                end
            end
            
            % for vectors, value matrix becomes column vector
            if length(dims) == 1, dims = [dims 1]; end
            
            A = zeros(dims); % initialize new value matrix
            A(new_idx{:}) = Gdata.val(old_idx{:}); % copy values
            
            % create output data structure
            Gdata.val  = A;
            Gdata.uels = uels;
            Gdata.ids  = ids;      
        end
        
        function Gdata = merge(g1, g2)
            %GAMS.merge merge contents of two data structures into one
            %   Data from two sets, variables or parameters are copied into
            %   one data structure that contains the contents of both
            %   inputs.
            %       Both inputs must be full (as in "not sparse") and have
            %   the same number of dimensions and type (set, parameter,
            %   variable or equation). Values of g2 overwrite values of g1.
            %       For scalar (zero-dimensional) data structures, merge
            %   falls back to adding the values of g1 and g2.
            %   Usage
            %       Gdata = GAMS.merge(g1, g2)
            %
            %   Parameters
            %       g1      first GAMS data structure
            %       g2      second GAMS data structure
            %
            %   Returns
            %       Gdata   result GAMS data structure
            
            % skip merge if one argument is empty
            if isempty(g1) || isempty(g2)
                if isempty(g2), Gdata = g1; end
                if isempty(g1), Gdata = g2; end
                return;
            end
            % error checking
            if ~strcmp(g1.type, g2.type)
                error('GAMS:merge:TypeMismatch','Both inputs must be of same type.');
            end
            if g1.dim ~= g2.dim
                error('GAMS:merge:DimMismatch','Both inputs must have the same number of dimensions.');
            end
            if ~strcmp(g1.form,'full') || ~strcmp(g2.form,'full')
                error('GAMS:merge:OnlyFull','Only full data structures can be merged. Please contact GAMS class maintainer if you need merging for sparse structures.');
            end
            
            % initialise temporary variables
            [new_uels, new_ids, g1_idx, g2_idx] = deal(cell(1, g1.dim));
            new_dims = zeros(1, g1.dim);
            
            for d=1:g1.dim
                % determine result uels, ids and size of value matrix
                new_uels{d} = union(g1.uels{d}, g2.uels{d});                
                new_ids{d}  = GAMS.uel_to_ids(new_uels{d});
                new_dims(d) = length(new_uels{d});
                
                
                if isnumeric(new_ids{d})
                    % if ids are a matrix, uels are numeric. In that case, sort
                    % them by value, not by literal
                    [new_ids{d}, sort_order] = sort(new_ids{d});
                    new_uels{d} = new_uels{d}(sort_order);
                    
                    % determine indices in matrix A for copying values
                    [~, g1_idx{d}] = intersect(new_uels{d}, g1.uels{d});
                    [~, g2_idx{d}] = intersect(new_uels{d}, g2.uels{d});
                    
                    % determine correct sort order
                    [~, g1_sort] = sort(new_ids{d}(g1_idx{d}));
                    [~, g2_sort] = sort(new_ids{d}(g2_idx{d}));
                    
                    % and sort g1/g2_idx accordingly
                    g1_idx{d} = g1_idx{d}(g1_sort);
                    g2_idx{d} = g2_idx{d}(g2_sort);
                    
                    % and clean up those sort vectors
                    clear sort_order g1_sort g2_sort;
                else
                    % default behaviour for textual uels
                    % determine indices in matrix A for copying values
                    [~, g1_idx{d}] = intersect(new_uels{d}, g1.uels{d});
                    [~, g2_idx{d}] = intersect(new_uels{d}, g2.uels{d});
                end

            end
            
            % for vectors, value matrix becomes column vector
            if length(new_dims) == 1
                new_dims = [new_dims 1]; 
                g1.val = g1.val(:);
                g2.val = g2.val(:);
            end
            
            
            if isempty(new_dims) % i.e. new_dims has length 0
                % scalars are merged by adding their values
                A = g1.val + g2.val;
            else
                % create value matrix
                A = zeros(new_dims);

                % copy values from g1
                A(g1_idx{:}) = g1.val;

                % copy values from g2
                if strcmp(g1.type,'set')
                    % for sets, perform logical or "+" on set elements
                    A(g2_idx{:}) = A(g2_idx{:}) + g2.val;
                    A = min(1, A); % remove value 2 from incidence matrix
                else
                    % for parameters, variables and equations,
                    % g2.val replaces an identical g1.val
                    A(g2_idx{:}) = g2.val;
                end
            end
            
            % determine new name for data structure
            if strcmp(g1.name, g2.name)
                % keep original name if both are identical
                new_name = g1.name;
            else
                % append names if names don't match
                new_name = [g1.name '_' g2.name];
            end
            
            % construct result data structure
            Gdata.name = new_name;
            Gdata.type = g1.type;
            Gdata.dim  = g1.dim;
            Gdata.val  = A;
            Gdata.uels = new_uels;
            Gdata.ids  = new_ids;
            Gdata.form = g1.form;
        
            % append field 'field' if inputs are variable or equation
            if sum(strcmp(g1.type,{'variable' 'equation'})) > 0
                Gdata.field = g1.field;
            end
        end
        
        function Gdata = sum(Gdata, dims)
            %GAMS.sum calculates sum for given dimensions
            %   Calculates sum of values for a given GAMS data structure.
            %   uels and ids are automatically adapted to match the
            %   result.
            %   Usage
            %       Gdata = GAMS.sum(Gdata, dims)
            %
            %   Example
            %       % CO2Out(time, site, pro, coin, coout)
            %       co2_by_site = GAMS.sum(CO2Out, [1 3:5]);
            %
            %   Parameters
            %       Gdata       original data structure
            %       dims        vector of dimensions
            %
            %   Returns
            %       Gdata       summarized data structure
            %
            if ~strcmp(Gdata.form,'full')
                error('GAMS:sum:OnlyFull','Only full data structures can be summed. Please contact your GAMS class maintainer if you need summing for sparse structures.');
            end
            
            if max(dims) > Gdata.dim
                error('GAMS:sum:WrongDimensions', ...
                    'Parameter dims contains wrong dimensions.');
            end
            
            % sort dimensions in descending order
            dims = sort(dims(:),'descend')';
            
            % vector of remaining dimensions
            remaining_dims = setdiff(1:Gdata.dim, dims);
            
            % sum value matrix...
            for d=dims
                Gdata.val = sum(Gdata.val, d);
            end
            % ... squeeze out any singleton dimensions...
            Gdata.val  = squeeze(Gdata.val);
            
            % ... and change dim, uels and ids accordingly
            Gdata.dim  = length(remaining_dims);
            Gdata.uels = Gdata.uels(remaining_dims);
            Gdata.ids  = Gdata.ids(remaining_dims);            
            
        end
        
        function ids = uel_to_ids(uels)
            %GAMS.uel_to_ids creates lookup structure for non-numeric uels
            %   This internal helper function creates a struct from a one-
            %   dimensional cell array, using its contents as fieldnames
            %   and their position 1..N as values. In case of numeric uels,
            %   this function returns the numeric array.
            %   Usage
            %       ids = GAMS.uel_to_ids(uels)
            %
            %   Example
            %       ids = GAMS.uel_to_ids({'a' 'b' 'c'})
            %       % returns struct('a',1,'b',2,'c',3)
            %
            %   Parameters
            %       uels    one-dimensional cell array of uels
            %
            %   Returns
            %       ids     struct with uels as fieldnames or numeric array
            %
            if sum(isnan(str2double(uels))) == 0
                % convert numeric uels
                ids = str2double(uels);
            else % provide id structure for easier access to named entries
                uels = regexprep(uels,'[^A-Za-z0-9]','_');
                uels = regexprep(uels,'^(\d|_)','x\1'); % if uel starts with a digit (\d) or underscore, prepend it with a character (here: x)
                ids = cell2struct(num2cell(1:length(uels)),uels,2);
            end
        end
        
        function uels = to_uels(vector)
            %GAMS.to_uels converts a numeric vector to cell array of strings
            %   Usage
            %       uels = GAMS.to_uels(vector)
            %
            %   Example
            %       uels = GAMS.to_uels(1:4)
            %       % returns {'1', '2', '3', '4'}
            %
            vals = reshape(vector,numel(vector),1);
            uels = cellstr(num2str(vals,'%-g'))';
        end
            
        
        function Gparam = sparse_to_full(Gparam)
            %GAMS.sparse_to_full converts sparse value matrix to full form
            %   This internal helper function is needed to convert sparse
            %   data structures as delivered from rgdx to the compact
            %   internal representations with nd-arrays.
            %   Usage
            %       param_full = GAMS.sparse_to_full(param_sparse)
            %
            %   Parameters
            %       param_sparse    sparse data structure
            %
            %   Returns
            %       param_full      full data structure
            %
            if ~(strcmp(Gparam.form,'sparse') && strcmp(Gparam.type,'parameter'))
                error('GAMS:sparse_to_full:WrongType','Argument either no parameter or not sparse.');
            end
            
            % determine dimensionality of output 
            Ndim = length(Gparam.uels);
            dims = zeros(1,Ndim);

            % initialise value matrix
            for k=1:Ndim
                dims(k) = length(Gparam.uels{k});
            end
            value_matrix = zeros(dims);

            % loop through values
            for v=1:size(Gparam.val,1)
                % convert subscript vector to cell array
                idx = num2cell(Gparam.val(v,1:Ndim));
                % use subscript vector for position addressing in value matrix
                value_matrix(idx{:}) = Gparam.val(v,Ndim+1);
            end
            
            % create output
            Gparam.val = value_matrix;
            Gparam.form = 'full';
        end
        
        function Gparam = full_to_sparse(Gparam)
            %GAMS.full_to_sparse converts full value matrix to sparse form
            %   This function converts the value matrix of a (default) full 
            %   GAMS data structure to the 2-dimension sparse form.
            %   Supported data structures are sets, parameters and
            %   variables. Equations might work as well (untested).
            %   Usage
            %       Gdata_sparse = GAMS.full_to_sparse(Gdata_full)
            %
            %   Parameters
            %       Gdata_full      full data structure
            %
            %   Returns
            %       Gdata_sparse    sparse data structure
            %
            if ~(strcmp(Gparam.form,'full'))
                error('GAMS:full_to_sparse:WrongType','Argument not full.');
            end
            
            % fix for row vector value matrices that acutally should be 
            % column vectors according to its uels. Example:
            %     name: 'x'
            %     type: 'parameter'
            %      val: [1 2 3]
            %     form: 'full'
            %      dim: 2
            %     uels: {{'a' 'b' 'c'}  {'one'}}
            % In this case and this case only, the value matrix is brought
            % to column vector form.
            if Gparam.dim == 2 && size(Gparam.val,1) == 1 && length(Gparam.uels{2}) == 1
                Gparam.val = Gparam.val(:); 
            end
            
            % find positions (linear index) of non-zero entries
            non_zero_entries = find(Gparam.val);
            
            % derive size of sparse value matrix
            Ndim = Gparam.dim;
            Nval = numel(non_zero_entries);
            dims = size(Gparam.val);
            idx  = cell(1,Ndim);
            value_matrix = zeros(Nval, Ndim+1);
            
            % loop through all non-zero values
            for v=1:Nval
                % create subscript vector idx from linear index
                [idx{:}] = ind2sub(dims, non_zero_entries(v));
                
                % write idx to first Ndim columns of value matrix
                value_matrix(v,1:Ndim) = cell2mat(idx);
                
                % add value to last column
                value_matrix(v,Ndim+1) = Gparam.val(non_zero_entries(v));
            end
            
            % special set treatment
            if strcmp(Gparam.type, 'set')
                % remove incidence column
                Gparam.val = value_matrix(:,1:end-1);
            else % parameter, variable, equation
                Gparam.val = value_matrix;
            end
            
            % create output
            Gparam.form = 'sparse';
            Gparam.val  = sortrows(Gparam.val);
        end
        
        % GDX import/export
        function Gdata = getGDX(filename, var, form, field)
            %GAMS.getGDX Reads GDX file to compact value array
            %   GAMS.getGDX is a wrapper function for rgdx, returning a 
            %   structure read from a GDX file.
            %   Example
            %       EprOut = GAMS.getGDX('result.gdx', 'eprout')
            %
            %   Usage
            %       G = GAMS.getGDX(filename, name, form, field)
            %
            %   Parameters
            %       filename    file name of GDX file to be read
            %       name        name of variable, equation, set or parameter
            %       form        (optional) form of output data structure
            %                   possible: 'full' or 'sparse'
            %                   default: 'full' if possible, 'sparse' else
            %       field       (optional) field to be read. l is level or
            %                   actual value, m is marginal value, lo and
            %                   up are the boundaries of a value
            %                   possible: 'l', 'm', 'lo', 'up'
            %                   default: 'l'
            %
            %   Returns
            %       G.name      name of variable, equation, parameter
            %       G.type      type ('variable', 'equation', ...) of name
            %       G.dim       array of dimension size, same as SIZE(G.val)
            %       G.val       nd-array of values
            %       G.uels      cell array of dimension labels
            %       G.ids       structures with dimension labels as fieldnames
            %       G.form      form of value matrix:  'full' or 'sparse'
            %       G.field     field of value matrix: 'l', 'm', 'lo' or 'up'
            %
            if nargin < 2, error('Wrong number of arguments: getGDX(filename[, var][, form][, field])'); end
            if isempty(dir(filename)), error('File not found.'); end
            if length(strfind(filename,'.')) ~= 1
                % rgdx() cannot handle filenames with more than one '.'.
                % Those files are temporarily copied in order to generate
                % a save name ... (1)
                tempdire = tempname('C:');
                mkdir(tempdire);
                tempfile = [tempname(tempdire) '.gdx'];
                copyfile(filename, tempfile);
                rgdxfn = tempfile;
            else
                rgdxfn = filename;
            end
            
            if nargin < 3 || strcmp(form,''), form = 'full'; end
            if nargin < 4, field = 'l'; end
                        
            % read GDX file using rgdx()
            rgdxopt = struct('name', var, 'compress', true);
            if nargin > 3, rgdxopt.field = field; end
            g = rgdx(rgdxfn, rgdxopt);

            if length(strfind(filename,'.')) ~= 1
                % (1) ... and deleted afterwards
                delete(tempfile);
                rmdir(tempdire);
            end
                        
            % UELS AND IDENTIFIERS
            C = cell(1, g.dim); % uels
            I = cell(1, g.dim); % ids
            dims = zeros(1, g.dim); % vector of lengths per dimension
            for k=1:g.dim
                dims(k) = length(g.uels{k});
                C{k} = g.uels{k};
                if sum(isnan(str2double(g.uels{k}))) == 0 % convert numeric uels
                    I{k} = str2double(g.uels{k});
                else % provide id structure for easier access to named entries
                    I{k} = GAMS.uel_to_ids(g.uels{k});
                end
            end

            % VALUES
            if isempty(dims) % handle scalar stuff
                A = g.val;
            else
                if length(dims) == 1 % handle vectorial (one-dimensional) stuff
                    if strcmp(g.type,'set')
                        A = g.val; % one-dimensional sets
                    else
                        A = g.val(:,2); % one-dimensional equations, parameters and variables
                    end
                else
                    if prod(dims) > 1e8 || strcmp(form,'sparse')
                        % sparse
                        A = g.val;
                        
                        if strcmp(form,'full')
                            warning('getGDX:TooBigForFull',['Data structure ''' g.name ''' too big for full value matrix, fallback to sparse output.']);
                            form = 'sparse';
                        end
                    else
                        
                        A = zeros(dims);
                        % convert g.val to n-dimensional array A
                        for v=1:size(g.val,1) % surprisingly, a for loop works quite well (=fast) here
                            idx = num2cell(g.val(v,1:end-1)); % convert subscript vector to cell array
                            A(idx{:}) = g.val(v,end);
                        end
                    end
                end
            end
            
            % create output data structure
            Gdata.name = g.name;
            Gdata.type = g.type;
            Gdata.dim = g.dim;
            Gdata.val = A;
            Gdata.uels = C;
            Gdata.ids = I;
            Gdata.form = form;
            
            % append field 'field' if variable or equation
            if sum(strcmp(g.type,{'variable' 'equation'})) > 0
                Gdata.field = g.field;
            end
            
            % remove uels/ids for scalars
            if Gdata.dim == 0
                Gdata = rmfield(Gdata, {'uels' 'ids'});
            end
        end
        
        function putGDX(filename, varargin) 
            %GAMS.putGDX Writes data structures to GDX file
            %   Takes an arbitrary number of arguments. Each argument must
            %   be a data structure like the ones created by GAMS.set() or
            %   GAMS.param().
            %   Example
            %       GAMS.putGDX('input.gdx', g_cost, g_param);
            %
            %   Usage
            %       GAMS.putGDX(filename, var1, var2, ..., varN)
            %
            %   Parameters
            %       filename        output GDX filename
            %       var1 to varN    GAMS data structures (Gset, Gparam)
            %
            %   Returns
            %       Nothing
            %
            %   Last changed
            %       2011-08-25 12:42 GMT+2      full_to_sparse for singleton
            %       2011-08-24 16:43 GMT+2      added this documentation
            
            % wgdx filename work-around, part 1
            if length(strfind(filename,'.')) ~= 1
                % wgdx cannot handle filenames with more than one '.'.
                % Those files are temporarily written with a temporary name ...
                tempdire = tempname('C:');
                mkdir(tempdire);
                tempfile = [tempname(tempdire) '.gdx'];
                rgdxfn = tempfile;
            else
                rgdxfn = filename;
            end
            
            for k=1:length(varargin)                
                % varargin{k} is a parameter in full (nd-array) form 
                % and its last dimension is singleton. Convert to
                % sparse form for writing because of a bug in wgdx.                
                if strcmp(varargin{k}.form, 'full') && varargin{k}.dim > 0
                    if iscell(varargin{k}.uels{end}) && length(varargin{k}.uels{end}) == 1
                        varargin{k} = GAMS.full_to_sparse(varargin{k});
                    end
                end
                
                % remove unused fields from struct
                unused_fieldnames = {'ids' 'field'};
                found = isfield(varargin{k},unused_fieldnames);
                varargin{k} = rmfield(varargin{k},unused_fieldnames(found));
                
                % automatic variable-to-parameter conversion
                if strcmp(varargin{k}.type, 'variable')
                    varargin{k}.type = 'parameter';
                    warning('GAMS:putGDX:VariableToParameter', ...
                        ['Changed type of variable ''' varargin{k}.name ...
                         ''' to parameter.']);
                end
            end
            
            % call wgdx
            wgdx(rgdxfn, varargin{:});
            
            % wgdx filename work-around, part 2
            if length(strfind(filename,'.')) ~= 1
                % ... and moved afterwards to the desired filename
                copyfile(tempfile,filename);
                delete(tempfile);
                rmdir(tempdire);
            end
        end
        
        function symbols = globGDX(var, files)
            %GAMS.globGDX Reads a symbol from multiple GDX files
            %   Reads a symbol from multiple GDX files and returns a
            %   cell array of GAMS data structures. This function is 
            %   a convenient way of comparing multiple scenario results
            %   programmatically. 
            %       See advanced example for a real-world
            %   application. Most usage will depend on function cellfun.
            %
            %   Examples
            %       co2 = GAMS.globGDX('PrOutEnv', {'a.gdx' 'b.gdx'})
            %       eprout = GAMS.globGDX('EprOut', 'result')
            %
            %   Advanced example
            %       co2 = GAMS.globGDX('PrOutEnv', 'result')
            %       co2_sum = cellfun(@(x) GAMS.sum(x, [1:6]), co2, 'UniformOutput', false)
            %       co2_val = cellfun(@(x) x.val, co2_sum)
            %       bar(co2_val)
            %
            %   Usage
            %       symbols = GAMS.globGDX(var, files)
            %
            %   Parameters
            %       var         name of symbol; must exist in all files
            %       files       cell array of GDX filenames or
            %                   name of a directory; in this case, all
            %                   GDX files in that directory will be read
            %
            %   Returns
            %       symbols     cell array of GAMS data structures
            %
            if ischar(files)
                % check if it is a directory, glob GDX files
                if isdir(files)
                    directory_name = files;
                    % add trailing slash
                    if directory_name(end) == '/' || directory_name(end) == '\'
                        directory_name(end) = '/';
                    else
                        directory_name = [directory_name '/'];
                    end
                    % glob files and create relative path
                    files = dir([directory_name '*.gdx']);
                    files = {files.name};
                    files = cellfun(@(x) [directory_name x], files, 'UniformOutput', false);
                % else pack into a cell array of size one
                else
                    files = {files}; 
                end
            elseif iscell(files)
                % do nothing, I don't have to check for everything!
            else
                % neither string nor cell array given: glob working
                % directory
                files = dir('*.gdx');
                files = {files.name};
            end
            
            all_vars = cell(1, length(files));
            for k=1:length(files)
                all_vars{k} = GAMS.getGDX(files{k},var);
            end
            
            symbols = all_vars;
        end
        
        % XLS import/export
        function [set_entity, att_entity, db_entity, s_onsets] = getXLS(filename, entity, type, form)
            %GAMS.getXLS Reads data structures from XLS files
            %   Reads an entity from its Excel spreadsheet and creates
            %   several GAMS.set and one GAMS.param (db_entity) that
            %   can be directly written using GAMS.putGDX.
            %
            %   An entity is a standardised form of spreadsheet which
            %   follows the following structure. Set names must begin with
            %   a capital lettre, attribute names with a miniscule lettre.
            %   Set elements and attribute names may contain lettres,
            %   numbers and dashes (-) or underscores (_) only. Spaces are
            %   possible, but might cause headaches in later data handling.
            %
            %   | Set1 Set2 Set3 | att1 att2 att3 att4 |
            %   +----------------+---------------------+
            %   | IT   Sun  pv   |  10   0.3  0.1 2500 |
            %   | IT   Sun  csp  |  20   0.1  0.9 3800 |
            %   | IT   Wind wt   |  25   0.4  ...      |
            %   | IT   Coal pp   | ...                 |
            %   | AT   Wind wt   |                     |
            %   | CH   Wind wt   |                     |
            %
            %   This function converts the Excel table to the following
            %   GAMS data structures:
            %   
            %   set_entity = { {'IT' 'Sun' 'pv'} {'IT 'Sun' 'csp'} ... }
            %   att_entity = { 'att1' 'att2' 'att3' 'att4' }
            %   db_entity  = [ 10 0.3 0.1 2500; 20 0.1 0.9 3800; ... ]
            %   s_onsets   = { Set1: {'AT' 'CH' 'IT'} 
            %                  Set2: {'Coal' 'Sun' 'Wind'}
            %                  Set3: {'csp' 'pp' 'pv' 'wt'} }
            %
            %   If the optional third parameter type is set to 'timeseries', 
            %   the following table structure is assumed:
            %
            %   | t | a.1 a.2 a.3 b.1 b.2 b.3 c.1 c.4 |
            %   +---+---------------------------------+
            %   | 1 |   4   3   1   1   0   2   1 0.4 |
            %   | 2 |   7   2   1   0   5   2   1 0.3 |
            %   | 3 |   9   1   1   0   0   1   0 0.2 |
            %
            %   The first column must contain (not necessarily numeric)
            %   timestep identifiers, while the columns are one- or
            %   multidimensional labels (separated by dots, like GAMS
            %   tuples) that define the onsets for the parameter.
            %
            %   Usage
            %       [set_entity att_entity db_entity s_onsets] = ...
            %           GAMS.getXLS(filename, name, 'entity', form)
            %     or
            %       [ts t cols onsets] = ...
            %           GAMS.getXLS(filename, name, 'timeseries')
            %
            %   Parameters
            %       filename        filename of an Excel spreadsheet
            %       name            name of worksheet
            %       type            'entity' (default) or 'timeseries'
            %       form            'full' (default) or 'sparse'
            %
            %   Returns
            %       set_entity      set of multi- or one-dimensional elements
            %       att_entity      one-dimensional set of attribute names
            %       db_entity       parameter of attribute values
            %       s_onsets        cell array of indivudal 1-dim onsets
            %     or
            %       ts              parameter of timeseries values
            %       t               set of timesteps
            %       cols            set of columns (one- or multi-dimensional)
            %       onsets          cell array of individual 1-dim onsets
            %
            if nargin < 4, form = 'full'; end
            if nargin < 3, type = 'entity'; end
            if sum(strcmp(type, {'entity' 'timeseries'})) ~= 1
                error('GAMS:getXLS:WrongType', ...
                     'Wrong value for parameter ''type'' given.'); 
            end

            % Read Excel file
            [~, ~, xlsdata] = xlsread(filename, entity);
            first_row = xlsdata(1,:);
            first_col = xlsdata(:,1);
            
            % Find first empty cell in first row/col
            empty_cols = find(cellfun(@(x) any(isnan(x)),first_row));
            if isempty(empty_cols)
                Ncols = length(first_row);
            else
                Ncols = empty_cols(1)-1;
            end
            empty_rows = find(cellfun(@(x) any(isnan(x)),first_col));
            if isempty(empty_rows)
                Nrows = length(first_col);
            else
                Nrows = empty_rows(1)-1;
            end
            
            % And truncate table there
            first_row = first_row(1:Ncols);
            first_col = first_col(1:Nrows);
            
            if strcmp(type, 'entity')
                % finds those columns, whose caption begins with a 
                % (captial/miniscule) lettre
                attribute_columns = strcmp(cellfun(@(str) {lower(str(1))},first_row), cellfun(@(str) {str(1)},first_row));
                set_columns       = strcmp(cellfun(@(str) {upper(str(1))},first_row), cellfun(@(str) {str(1)},first_row));

                % and create sets of onsets and attributes
                onset_names = first_row(set_columns);
                att_names   = first_row(attribute_columns);

                % determine size of submatrices
                Nsets = length(onset_names);
                Natts = length(att_names);
                Nvals = Nrows-1;

                % initialise cell array for set elements
                elements = cell(1,Nvals);
                
                % convert numeric set elements to strings
                numeric_set_entries = cellfun(@isnumeric, xlsdata(:,1:Nsets));
                xlsdata(numeric_set_entries) = cellfun(...
                        @(x) num2str(x,'%-g'), ...
                        xlsdata(numeric_set_entries), ...
                        'UniformOutput', false);
                
                % fill set element values
                if Nsets == 1
                    for k=1:Nvals
                        elements{k} = xlsdata{k+1,1}; % note curly braces
                    end
                else
                    for k=1:Nvals
                        elements{k} = xlsdata(k+1,1:Nsets);
                    end
                end

                % Onsets
                [s_onsets, onsets] = deal(cell(1,Nsets));
                for k = 1:Nsets
                    onsets{k}   = unique(xlsdata(2:Nrows,k)');
                    s_onsets{k} = GAMS.set(lower(onset_names{k}), onsets{k});
                end

                % Values
                values_in_cells = xlsdata(2:Nrows,Nsets+1:Ncols);
                non_numeric = cellfun(@isstr, values_in_cells);
                values_in_cells(non_numeric) = cellfun(@str2num, values_in_cells(non_numeric),'UniformOutput',false);

                values = cell2mat(values_in_cells);

                % create output data structures

                if Nsets == 0
                    set_entity = {};
                elseif Nsets == 1
                    set_entity = GAMS.set(['set_' lower(entity)], elements);
                else
                    set_entity = GAMS.set(['set_' lower(entity)], elements, onsets, form);
                end
                
                if Natts > 0
                    att_entity = GAMS.set(['att_' lower(entity)],att_names);
                    db_entity  = GAMS.param(['db_' lower(entity)],values,{elements att_names}, form);
                else
                    att_entity = {};
                    db_entity  = {};
                end

                % special case: if only one attribute exists and its name is 
                % 'value', omit attributes completely from .uels and .val 
                % for easier addressing in GAMS: db_entity(a,b,c,'value') 
                % becomes db_entity(a,b,c) 
                if Natts == 1 && strcmp(att_names,'value')
                    att_entity = {};
                    db_entity.dim = db_entity.dim -1;
                    db_entity.uels = db_entity.uels(1:end-1);
                    db_entity.ids  = db_entity.ids(1:end-1);
                    
                    % remove uels/ids if empty
                    if isempty(db_entity.uels)
                        db_entity = rmfield(db_entity,{'uels' 'ids'});
                    end
                end
                
            elseif strcmp(type, 'timeseries');
                
                if ~strcmp(first_row(1), 't'), warning('GAMS:getXLS:UncommonName', ['First column in timeseries table ''' entity ''' is not labeled ''t''. Please inform GAMS maintainer if your timeseries set is called differently.']); end
                
                % find commodity captions
                % and split them by the dot (.) character
                captions = first_row(2:end);
                captions = regexp(captions, '\.', 'split');
                
                % error check: do all column title have same number of
                % dimensions?
                element_lengths = cellfun(@(x) length(x),captions);
                if ~all(element_lengths==element_lengths(1)), error('GAMS:getXLS:DimensionMismatch', ['One column in timeseries table ''' entity ''' has wrong number of dimensions.']); end
                
                % determine onsets for column title substrings
                Nsets = element_lengths(1);
                [onsets, s_onsets] = deal(cell(1,Nsets));
                for k = 1:Nsets
                    onsets{k}   = unique(cellfun(@(x) x(k),captions));
                    s_onsets{k} = GAMS.set(['on_' lower(entity) '_' num2str(k)], onsets{k});
                end
                
                % extract timesteps
                t = cell2mat(first_col(2:end));
                set_t = GAMS.set(['t_' lower(entity)], t);
                
                % and create set of column captions
                set_cols = GAMS.set(['col_' lower(entity)], captions, onsets);
                
                % now create uels/ids for db_entity
                uels = [set_t.uels(1) onsets];
                ids  = cellfun(@(x) GAMS.uel_to_ids(x),uels,'UniformOutput',false);
                
                % determine size of value matrix
                dimensionality = cellfun(@(x) length(x),uels);
                
                % create scaffold of db_entity
                ts.name = ['ts_' lower(entity)];
                ts.type = 'parameter';
                ts.form = 'full';
                ts.dim  = length(uels);
                ts.val  = zeros(dimensionality);
                ts.uels = uels;
                ts.ids  = ids;
                
                % pre-process timeseries in xlsdata
                values_in_cells = xlsdata(2:Nrows,2:Ncols);
                non_numeric = cellfun(@isstr, values_in_cells);
                values_in_cells(non_numeric) = cellfun(@str2num, values_in_cells(non_numeric),'UniformOutput',false);

                values = cell2mat(values_in_cells);
                
                % and fill values column by column
                for k = 1:length(captions)
                    el = captions{k};
                    
                    % find idx vector to address 'el' in db_entity.val
                    idx = cell(1,Nsets);
                    for d=1:Nsets
                        idx{d} = find(strcmp(el{d}, onsets{d}));
                    end
                    
                    % and copy data
                    ts.val(:,idx{:}) = values(:,k);
                end
                
                % output
                % [ts t cols onsets]
                set_entity = ts;
                att_entity = set_t;
                db_entity  = set_cols;
                % s_onsets is already set
                
            end
                
        end
        
        function putXLS(filename, varargin)
            %GAMS.putXLS Writes data structures to XLS file
            %   Takes an arbitrary number of arguments. Each argument must
            %   be a data structure like the ones created by GAMS.set,
            %   GAMS.param or GAMS.getXLS.
            %       By default, the last dimension of each data structure
            %   is used for columns. If that is not desired, give the
            %   keyword 'sparse' as the second argument.
            %   Example
            %       GAMS.putXLS('result.xls', g_set, g_cost);
            %
            %   Usage
            %       GAMS.putXLS(filename, var1, var2, ...)
            %       GAMS.putXLS(filename, 'sparse', var1, var2, ...)
            %
            %   Parameters
            %       filename        output XLS filename
            %       var1, var2      GAMS data structures
            %
            %   Returns
            %       Nothing
            %
            
            % determine whether new file is to be created
            if ~exist(filename, 'file'), new_file = 1; else new_file = 0; end
            
            % suppress warning "added specified worksheet"
            warning('off','MATLAB:xlswrite:AddSheet');
            
            % write_mode: sparse or full?
            if length(varargin) > 1 && ischar(varargin{1}) && strcmp(varargin{1},'sparse')
                varargin = varargin(2:end);
                write_mode = 'sparse';
            else
                write_mode = 'full';
            end
            
            for k=1:length(varargin)
                % var is the variable to be written
                var = varargin{k};
                
                % skip if variable is empty and issue a warning
                if isempty(var) 
                    warning('GAMS:putXLS:EmptyVariable',['Argument nr. ''' num2str(k) ''' is empty and not written to XLS file.']); 
                    continue;
                end
                if isempty(var.val)
                    warning('GAMS:putXLS:EmptyVariable',['Data structure ''' var.name ''' is empty and not written to XLS file.']); 
                    continue;
                end
                
                % initialize output data structure
                xls.name = var.name;
                
                % protection against legacy data structures:
                % if var uses "plain uels" (one-dimensional), pack them in
                if isfield(var,'uels') && var.dim > 0 && ~iscell(var.uels{1})
                    var.uels = {var.uels};
                end
                
                if strcmp(var.type, 'set')
                    % convert data structures from full to sparse
                    if ~strcmp(var.form, 'sparse')
                        var = GAMS.full_to_sparse(var);
                    end

                    % initialize output cell array
                    Nvals = size(var.val,1);
                    Ndims = size(var.val,2);
                    xls.val = cell(Nvals+1, Ndims);

                    % fill xls.val
                    for n=1:Ndims
                        % write column headings
                        xls.val{1,n} = ['Onset' num2str(n)];

                        % write set element names
                        uel_indices = var.val(:,n); % extract uel indices
                        xls.val(2:end,n) = var.uels{n}(uel_indices); % and select uels
                    end
                else
                    % parameters, variables, possibly even equations
                    if var.dim > 1
                        % two- and more-dimensional params/variables

                        % convert data structures from full to sparse
                        if ~strcmp(var.form, 'sparse')
                            var = GAMS.full_to_sparse(var);
                        end
                        

                        if strcmp(write_mode, 'full')
                            % write_mode full means that the
                            % last onset is used for column headings
                            % (like attribute names)
                            
                            % determine number of onsets as number of
                            % columns minus 2: last column is value, second
                            % last column is attribute, rest onsets
                            Ndims = size(var.val,2);
                            Nons  = Ndims - 2;

                            % onsets
                            % unique index combinations in first Nons
                            % columns become the rows of xls.val
                            [onsets, ~, ic] = unique(var.val(:,1:Nons),'rows');

                            % declaration of xls.val
                            Nrows = 1 + size(onsets,1);
                            Ncols = Nons + length(var.uels{end});
                            xls.val = cell(Nrows, Ncols);

                            % index handling
                            % determine index vector idx, so that all
                            % values of xls.val can be accessed in a single
                            % assignment like (in pseudo-code) this:
                            %     xls.val(idx) = var.val(:,end)
                            % For this, first determine row (set_idx) and
                            % column (att_idx) indices from the unique
                            % onset rows (while creating them!), while
                            % column numbers can be directly obtained from
                            % attribute numbers in column end-1, offset by
                            % Nons (number of offset columns in xls.val).
                            set_idx = ic+1;
                            att_idx = var.val(:,end-1) + Nons;
                            idx = sub2ind(size(xls.val), set_idx, att_idx);
                            
                            % fill xls.val

                            % onsets
                            for n = 1:Nons
                                xls.val(1,    n) = {['Onset' num2str(n)]};
                                xls.val(2:end,n) = var.uels{n}(onsets(:,n));
                            end

                            % attribute column titles
                            for n = Nons+1:Ncols
                                xls.val(1,n) = var.uels{end}(n-Nons);
                            end

                            % values
                            % first convert numeric values to cell array
                            % of scalar values. Then replace infinite
                            % values to the string 'inf'. Finally initialise 
                            % all values in xls.val to zero, then assign 
                            % non-zero values using the linear index vector
                            values = num2cell(var.val(:,end));
                            is_infinite = cellfun(@isinf, values);
                            values(is_infinite) = {'inf'};
                            xls.val(2:end, Nons+1:end) = num2cell(0);
                            xls.val(idx) = values;
                            
                        else % sparse
                            % write_mode sparse means that the
                            % all onsets are used for row indexing
                            % (like in GDXviewer)
                            
                            % initialize output cell array
                            Nvals = size(var.val,1); %#ok<NASGU>
                            Ndims = size(var.val,2);
                            Nons  = Ndims - 1;

                            onsets = unique(var.val(:,1:Nons),'rows');

                            Nrows = 1 + size(onsets,1);
                            Ncols = Nons + 1;
                            xls.val = cell(Nrows, Ncols);

                            % fill xls.val

                            % onsets
                            for n = 1:Nons
                                xls.val(1,    n) = {['Onset' num2str(n)]};
                                xls.val(2:end,n) = var.uels{n}(onsets(:,n));
                            end
                            
                            % values
                            xls.val{1,Nons+1} = 'value';
                            value_cell = num2cell(var.val(:,end)); % convert to cell array
                            value_cell(isinf(var.val(:,end))) = deal({'inf'}); % change infinite values to string 'inf'
                            if ~isempty(value_cell)
                                xls.val(2:end, Nons+1) = value_cell;
                            else
                                xls.val(2:end, Nons+1) = num2cell(0);
                            end
                        end
                            

                    elseif var.dim == 1
                        % one-dimensional 

                        % initialize xls.val
                        Nrows = 1 + length(var.val);
                        Ncols = var.dim + 1;
                        xls.val = cell(Nrows, Ncols);

                        % fill xls.val
                        xls.val(1,    :) = {'Onset' 'value'};
                        xls.val(2:end,1) = var.uels{1};
                        xls.val(2:end,2) = num2cell(var.val);

                    else
                        % zero-dimensional
                        xls.val = {'value'; var.val};                    
                    end
                    
                end % if strcmp(var.type, 'set')
                    
                xlswrite(filename,xls.val,xls.name);
            end
            
            if new_file % delete default three sheets
                excelFileName = filename;
                excelFilePath = pwd; % Current working directory.
                defaultSheetNames = {'Sheet' 'Tabelle'}; % EN: Sheet, DE: Tabelle, etc. (Lang. dependent)

                % Open Excel file.
                objExcel = actxserver('Excel.Application');
                objExcel.Workbooks.Open(fullfile(excelFilePath, excelFileName)); % Full path is necessary!

                % Delete default sheets.
                for sheetName=defaultSheetNames
                    for s=1:3
                        try
                            % Throws an error if the sheets do not exist.
                            objExcel.ActiveWorkbook.Worksheets.Item([sheetName{:} num2str(s)]).Delete;
                        catch  %#ok<CTCH>
                            ; %#ok<NOSEM> % do nothing
                        end
                    end
                end

                % Save, close and clean up.
                objExcel.ActiveWorkbook.Save;
                objExcel.ActiveWorkbook.Close;
                objExcel.Quit;
                objExcel.delete;
            end
            
            % turn warning back on
            warning('on','MATLAB:xlswrite:AddSheet');
        end
        
        function status = GDX2MDB(filename)
            %GAMS.GDX2MDB convert GDX file to Access DB
            %   Uses the executable 'gdx2access' in the GAMS installation
            %   folder to convert a GDX file to a Microsoft Access MDB
            %   file.
            %   Example
            %       GAMS.GDX2MDB('result\result.gdx');
            %
            %   Usage
            %       status = GAMS.GDX2MDB(filename)
            %
            %   Parameters
            %       filename    path to GDX file
            %
            %   Returns
            %       status      return value of gdx2access
            %
            command = sprintf('gdx2access "%s"\n', filename);
            status = system(command);
        end
    end
end

