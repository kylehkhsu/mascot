classdef Goal < handle
% ReadGoal
% Read the vector of goal states from files in Z/ and G/
% generated during controller synthesis
%
%   Kaushik Mallik, July 2018

properties (SetAccess=private)
    filename
    dim
    eta
    points
    handle
end
methods
    function obj=Goal(filename,varargin)
    % the constructor opens the file and creates an instance of 
    % a StaticController with the data read from file 'filename'

      if(ischar(filename))
        obj.filename=filename;
      else
        error('filname is not a string');
      end
      obj.filename=filename;

      h=mexGoal('init',filename);
      if(isempty(h))
        error(['Goal: could not read Goal from ', filename])
      end
      obj.handle=h;
      
      [dim, eta] = mexGoal('parameters',obj.handle);
      obj.dim = dim;
      obj.eta = eta;
%       function points=get.points(obj)
      % return the grid points that are in the domain of the controller
%       if ~isempty(obj.points)
%         points=obj.points;
%       else
%         points=mexGoal('points',obj.handle);
%         obj.points=points;
%       end
%     end
    end
    
    function disp(obj)
      disp(['Matlab class to access the (intermediate) goals computed with MASCOT stored in ', obj.filename])
      disp(' ')
    end
    
    function points=get.points(obj)
      % return the grid points that are in the domain of the controller
      if ~isempty(obj.points)
        points=obj.points;
      else
        points=mexGoal('points',obj.handle);
        obj.points=points;
      end
    end
    
    function out = isElement(obj,x)
      x=x(:);
      if (mexGoal('checkstate',obj.handle,x))
        out=1;
      else
        out=0;
      end
    end

%     function obj = Goal(filename, varargin)
%         if ~ischar(filename)
%             error('invalid filename');
%         end
% 
%         fid = fopen(filename);
%         if fid==-1
%             error('file cannot be opened');
%         end
% 
%         obj.filename = filename;
%         tline = fgetl(fid);
%         while ischar(tline)
%             if contains(tline, "ETA")
% %                 obj.eta = str2double(extractAfter(tline, ":"));
%                 tline = fgetl(fid);
%                 str = extractAfter(tline, ":");
%                 dim = str2double(str);
%                 obj.eta = zeros(1, dim);
%                 break;
%             end
% %             if contains(tline, "BEGIN")
% %                 break;
% %             end
%             tline = fgetl(fid);
%         end
%         tline = fgetl(fid);
%         index = 1;
%         while ~contains(tline, "END")
%             obj.eta(index) = str2double(tline);
%             index = index+1;
%             tline = fgetl(fid);
%         end
%         
%         while ischar(tline)
%             if contains(tline, "TARGET_STATES")
%                 tline = fgetl(fid);
%                 newStr = extractAfter(tline, ":");
%                 size = str2double(newStr);
%                 obj.points = zeros(size, 1);
%                 break;
%             end
%             tline = fgetl(fid);
%         end
%         
%         tline = fgetl(fid);
%         index = 1;
%         while ~contains(tline, "END")
%             obj.points(index) = str2double(tline);
%             index = index +1;
%             tline = fgetl(fid);
%         end
%     end
end
end


