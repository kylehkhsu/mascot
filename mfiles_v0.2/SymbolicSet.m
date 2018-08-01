classdef SymbolicSet < handle 
% Access the grid points stored in a BDD produced by 
%
% scots::write_to_file(const Cudd& manager, const SymbolicSet& set, const BDD& bdd, const std::string& filename, char mode='B')
% 
% USAGE:
% 
% set = SymbolicSet('filename')  reads the SymbolicSet from file
%  
% X = set.points;                X is a matrix containing all grid points 
%
% X = set.points([2 3]);         X is a matrix containing all grid points projected 
%                                on the dimensions 2 and 3
%
% Y = set.restriction(x);        Y is a matrix containing all grid points such
%                                that for each grid point y in Y, the grid point (x,y)
%                                is an element of the set
%
% Y = set.restriction(x,dim);    same as before, but the entries of x are placed
%                                at the indices specified by dim (see also the
%                                documentation in ./doc/html of SymbolicSet::restriction)
% 
%
  properties (SetAccess=private)
    filename  % name of file which contains the SymbolicSet
    h         % C++ object handles
    dim       % dimension 
  end
  methods 
    function obj=SymbolicSet(filename)
    % the constructor opens the file 
      if(ischar(filename))
        obj.filename=filename;
      else
        error('filname is not a string');
      end

      [h, dim]=mexSymbolicSet('init',filename);
      obj.filename=filename;
      obj.h = h;
      obj.dim = dim;
    end
    function disp(obj)
      disp(['Matlab object to access the points cointained in the SymbolicSet stored in ', obj.filename])
      disp(' ')
    end
    function delete(obj)
      if(~isempty(obj.h))
        mexSymbolicSet('delete', obj.h);
      end
    end

    function points=restriction(obj,x,varargin)
      x=x(:);
      n=length(x);
      if(n>=obj.dim)
        error('SymbolicSet: vector equals or exceeds dimension');
      end
      if(isempty(varargin))
        r_dim=0:1:length(x)-1;
      else 
        r_dim=varargin{1}(:);
        if(length(r_dim)~=n)
          error('SymbolicSet: vector and project dimensions must be of the same length');
        end
        if(max(r_dim)>obj.dim)
          error('SymbolicSet: project dimensions exceed dimension');
        end
        r_dim = r_dim-1;
      end
      r_dim=r_dim(:);
      points=mexSymbolicSet('restriction',obj.h,x,r_dim);
      if(isempty(points))
        error(['SymbolicSet: no grid points found for x=[',num2str(x'),']']);
      end
    end

    function points=points(obj,varargin)
    % read the grid points from file
    % optional projection indices
    % get.point('project',[1 3]) loads the grid points in the symbolic set
    % projected onto the dimensions [1 3]
      if(isempty(varargin))
        p_dim=0:1:obj.dim-1;
        p_dim=p_dim(:);
      else 
        p_dim=varargin{1}(:);
        if(length(p_dim)>=obj.dim)
          error('SymbolicSet: vector of project dimensions exceeds dimension');
        end
        if(max(p_dim)>obj.dim)
          error('SymbolicSet: project dimensions exceed dimension');
        end
        p_dim = p_dim-1;
      end
      points=mexSymbolicSet('gridpoints',obj.h,p_dim);
    end
  end
end
