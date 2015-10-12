classdef toastRegul < handle
    % The toastRegul class represents a regularisation object in the toast
    % toolbox. Regularisation objects can be used to add an image-based
    % penalty term to optimisation-based image reconstruction algorithms.
    %
    % See also:
    %    toast2
    methods
        function obj = toastRegul(varargin)
            % Create a new regularisation object.
            %
            % Syntax: regul = toastRegul(prm,x0)
            %         regul = toastRegul(method,basis,x0,tau,...)
            %
            % Parameters:
            %         prm [struct]:
            %             structure containing the regularisation
            %             parameters (see notes)
            %         method [string]:
            %             Regularisation method. Available options are:
            %             'TV', 'TK0', 'TK1', 'Huber', 'PM', 'QPM', 'Tukey'
            %         basis [object]:
            %             basis mapping object
            %         x0 [real array]:
            %             initial solution estimate
            %         tau [real]:
            %             regularisation hyperparameter
            %         Depending on the regularisation method selected,
            %         additional parameters may be required (see notes to
            %         individual regularisation methods)
            %
            % Return values:
            %         regul [object]:
            %             new regularisation object
            %
            % Notes:  When using toastRegul with a parameter structure, the
            %         following fields must be present in the structure:
            %
            %         field        value
            %         -----------------------------------------------------
            %         prm.method   a string defining the regularisation
            %                      method. Supported are: 'TV' (total
            %                      variation), 'TK0' (0-th order Tikhonov),
            %                      'TK1' (1st order Tikhonov), 'Huber',
            %                      'PM' (Perona-Malik), 'QPM' (quadratic
            %                      Perona-Malik), 'Tukey'
            %         prm.tau      a scalar corresponding to the 'tau'
            %                      parameter
            %         prm.x0       a vector corresponding to the 'x0'
            %                      parameter
            %
            %         If the regularisation object is to support a
            %         structural prior, an additional 'prior' field
            %         containing a sub-structure must be added:
            %
            %         field             value
            %         -----------------------------------------------------
            %         prm.prior.refimg  reference image from which the edge
            %                           prior is generated. Must be of the
            %                           same dimension as the grid basis
            %                           defined by prm.basis.
            %         prm.prior.refname file name of a reference image to
            %                           be read. Cannot be used together
            %                           with refimg.
            %       
            %         For specific regularisation methods, additional
            %         parameters may be required. These are to be added as
            %         a sub-structure using the regularisation name
            %
            %         field         value
            %         -----------------------------------------------------
            %         prm.tv.beta   scalar TV regularisation parameter
            %         prm.tk0.xs    scaling vector for TK0 regulariser
            %         prm.huber.eps scalar regularisation parameter
            %         prm.pm.t      scalar regularisation parameter
            %         prm.qpm.t     scalar regularisation parameter
            %         prm.tukey.t   scalar regularisation parameter
            %         
            %         When using toastRegul without a parameter structure,
            %         the compulsory parameters are supplied directly.
            %         Additional parameters for structural priors or
            %         specific regularisation methods are supplied with
            %         key/value pairs e.g. ... 'beta', 0.1, ...
            %
            % See also:
            %         toastRegul
            if isstruct(varargin{1}) % uses a structure for parameters
                prm = varargin{1};
                x0 = varargin{2};
                if isfield(prm,'basis')
                    % need to unpack the basis handle here because we can't
                    % get at it from C++
                    hbasis = prm.basis.handle;
                else
                    error('toastRegul: Required field not found: basis');
                end
                obj.handle = toastmex(uint32(44),prm,hbasis,x0,0);
            else % uses parameters directly
                typestr = varargin{1};
                basis = varargin{2};
                x0 = varargin{3};
                tau = varargin{4};
                if nargin >= 5
                    extras = varargin(5:end);
                else
                    extras = [];
                end
    	        obj.handle = toastmex(uint32(44),typestr,basis.handle,x0,tau,extras);
            end
        end
        
        function delete(obj)
            if obj.handle > 0
                toastmex(uint32(45),obj.handle);
            obj.handle = 0;
            end
        end

        function val = Value (obj, x)
            % Return the regularisation value of an image.
            %
            % Syntax: v = regul.Value(x)
            %
            % Parameters:
            %         x [real array n]:
            %             parameter distribution
            %
            % Return values:
            %         v [real]:
            %             regularisation value
            %
            % Notes:  Given a parameter distribution x, this method returns
            %         the corresponding regularisation value.
            %
            % See also:
            %         toastRegul
            val = toastmex(uint32(46),obj.handle,x);
        end

        function grad = Gradient (obj, x)
            % Return the regularisation gradient for an image
            %
            % Syntax: g = regul.Gradient(x)
            %
            % Parameters:
            %         x [real array n]:
            %             coefficient vector
            %
            % Return values:
            %         g [real array n]:
            %             regularisation gradient vector
            %
            % Notes: Given a parameter distribution x, this method returns
            %        the vector representing the corresponding
            %        regularisation gradient.
            %
            % See also:
            %         toastRegul
            grad = toastmex(uint32(47),obj.handle,x);
        end

        function hdiag = HDiag (obj, x)
            % Return the diagonal of the Hessian of the regularisation
            % operator.
            %
            % Syntax: hdiag = regul.HDiag(x)
            %
            % Parameters:
            %         x [real array n]:
            %             coefficient vector
            %
            % Return values:
            %         hdiag [real array n]:
            %             diagonal of Hessian
            %
            % Notes:  Given a parameter distribution x, this method returns
            %         the vector representing the diagonal of the Hessian
            %         of the regularisation operator.
            %
            % See also:
            %         toastRegul
            hdiag = toastmex(uint32(48),obj.handle,x);
        end
        
        function hess = Hess (obj, x)
            % Return the Hessian of the regularisation operator as a sparse
            % matrix.
            %
            % Syntax: hess = regul.Hess(x)
            %
            % Parameters:
            %         x [real array n]:
            %             coefficient vector
            %
            % Return values:
            %         hess [real csr matrix n x n]:
            %             Hessian of regularisation operator
            %
            % See also:
            %         toastRegul
            hess = toastmex(uint32(49),obj.handle,x);
        end
        
        function kappa = Kappa(obj, x)
            kappa = toastmex(uint32(50),obj.handle,x);
        end

        function SetLocalScaling(obj, scale)
            toastmex(uint32(96),obj.handle,scale);
        end
    end

    properties
        handle = uint64(0);
    end
end
