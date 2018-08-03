classdef StaticController < handle 
% Acces the controller stored by
%
% scots::write_to_file(const scots::StaticController&, const std::string&) 
%
% USAGE:
% 
% con = StaticController('filename')  reads the StaticController from file
%  
% U = con.control(x);                 U is a matrix containing all valid control
%                                     inputs at x
%
% X = con.domain;                     X is a matrix containing all centers of
%                                     cells that are in the winning domain
% 
%
  properties (SetAccess=private)
    filename  % name of file which contains the StaticController
    domain    % the gird points for which there exist valid inputs
    tau       % time sampling associated with the controller
    handle    % handle to C++ StaticController (created in the mex file)
  end
  methods 
    function obj=StaticController(filename,varargin)
    % the constructor opens the file and creates an instance of 
    % a StaticController with the data read from file 'filename'

      if(isstr(filename))
        obj.filename=filename;
      else
        error('filname is not a string');
      end
      obj.filename=filename;

      h=mexStaticController('init',filename);
      if(isempty(h))
        error(['StaticController: could not read StaticController from ', filename])
      end
      obj.handle=h;
      
      t = mexStaticController('parameters',obj.handle);
      obj.tau = t;
    end

%     function delete(obj)
%       if(~isempty(obj.handle))
%         mexStaticController('delete', obj.handle);
%       end
%     end

    function disp(obj)
      disp(['Matlab class to access the StaticController computed with SCOTS stored in ', obj.filename])
      disp(' ')
    end

    function u=control(obj,x)
      % return control inputs associated with grid point x  
      u=mexStaticController('control',obj.handle,x(:));
      if(isempty(u))
        error(['scots::StaticController: state ',...
                mat2str(x(:)), ' is out of winning domain: no progress possible.'])
      end
    end

    function domain=get.domain(obj)
      % return the grid points that are in the domain of the controller
      if ~isempty(obj.domain)
        domain=obj.domain;
      else
        domain=mexStaticController('domain',obj.handle);
        obj.domain=domain;
      end
    end
  end
end
