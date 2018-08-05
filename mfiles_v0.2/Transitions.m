classdef Transitions < handle 
% Acces the controller stored by
%
% scots::write_to_file(const scots::TransitionFunction&, const std::string&) 
%
% USAGE:
% 
% con = Transitions('filename')  reads the TransitionFunction from file
%  
% U = con.control(x);                 U is a matrix containing all valid control
%                                     inputs at x
%
% X = con.domain;                     X is a matrix containing all centers of
%                                     cells that are in the winning domain
% 
%
  properties (SetAccess=private)
    filename  % name of file which contains the TransitionFunction
    layer     % the abstraction layer
    domain    % the gird points for which there exist valid inputs
    handle    % handle to C++ TransitionFunction (created in the mex file)
  end
  methods 
    function obj=Transitions(filename,varargin)
    % the constructor opens the file and creates an instance of 
    % a TransitionFunction with the data read from file 'filename'

      if(isstr(filename))
        obj.filename=filename;
      else
        error('filname is not a string');
      end
      obj.filename=filename;

      h=mexTransitions('init',filename);
      if(isempty(h))
        error(['Transitions: could not read TransitionFunction from ', filename])
      end
      obj.handle=h;
      obj.layer=varargin{1};
      
%        domain=mexTransitions('domain', obj.handle);
%        obj.domain=domain;
    end

%     function delete(obj)
%       if(~isempty(obj.handle))
%         mexTransitions('delete', obj.handle);
%       end
%     end

    function disp(obj)
      disp(['Matlab class to access the TransitionFunction computed with SCOTS stored in ', obj.filename])
      disp(' ')
    end
    
%     function u=control(obj,x)
%       % return control inputs associated with grid point x  
%       u=mexStaticController('control',obj.handle,x(:));
%       if(isempty(u))
%         error(['scots::StaticController: state ',...
%                 mat2str(x(:)), ' is out of winning domain: no progress possible.'])
%       end
%     end

    function domain=get.domain(obj)
      % return the grid points that are in the domain of the
      % TransitionFunction
      if ~isempty(obj.domain)
        domain=obj.domain;
      else
        domain=mexTransitions('domain',obj.handle,obj.layer);
        obj.domain=domain;
      end
    end
  end
end
