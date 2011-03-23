function [PREC, options] = AMGfactor(A, options)
% [PREC, options] = AMGfactor(A, options)
% [PREC, options] = AMGfactor(A)
%
%
% Computes ILUPACK preconditioner PREC according to the given options.
% For details concerning `options' see `AMGinit'
%
% input
% -----
% A           nxn nonsingular matrix
% options     parameters. If `options' is not passed then the default options
%             from `AMGinit' will be used
%
% output
% ------
% PREC        ILUPACK multilevel preconditioner
%             PREC is a structure of `nlev=length(PREC)' elements indicating
%             the number of levels.
%             For every level l we have
%             PREC(l).n           size of level l
%             PREC(l).nB          size of the leading block of level l
%
%             PREC(l).L           (block) lower triangular matrix
%             PREC(l).D           (block) diagonal matrix
%             PREC(l).U           if present: (block) upper triangular matrix
%                                 for symmetrically structured matrices only
%                                 L is present.
%                                 L*D^{-1}*U 
%                                 is the approximate LU decomposition of the 
%                                 leading block of A after rescaling and 
%                                 reordering
%
%             PREC(l).E           except for l<nlev this refers to the lower
%                                 left block of A after rescaling and 
%                                 reordering and after the approximate LU
%                                 decomposition has been computed for the
%                                 leading block
%             PREC(l).F           except for l<nlev this refers to the upper
%                                 right block of A after rescaling and 
%                                 reordering and after the approximate LU
%                                 decomposition has been computed for the
%                                 leading block. Not used in the case
%                                 of symmetrically structured matrices
% 
%             PREC(l).rowscal     row scaling
%             PREC(l).colscal     column scaling
%
%             PREC(l).p           permutation for the rows of A
%             PREC(l).invq        inverse permutation for the columns of A
%
%             PREC(l).param       internal data to communicate with 'AMGsolver',
%             PREC(l).ptr         'AMGsol' and `AMGdelete'. DO NOT TOUCH!
%
%             PREC(l).isreal      flags to indicate which preconditioner has 
%             PREC(l).isdefinite  been computed, internal data to communicate
%             PREC(l).issymmetric with 'AMGsolver', `AMGsol' and `AMGdelete'. 
%             PREC(l).ishermitian DO NOT TOUCH!
%             PREC(l).A_H         coarse grid system
%             PREC(l).errorL      error estimate for L
%             PREC(l).errorU      error estimate for U
%             PREC(l).errorS      error estimate for coarse grid system A_H
%             PREC(l).isblock     block structured ILU
%
% options     updated parameters

% make sure that shifted system and A have same type
if isfield(options,'shiftmatrix')
   shiftreal=1;
   shifthermitian=1;
   shiftsymmetric=1;

   if ~isreal(options.shiftmatrix)
      shiftreal=0;
   end  
   if isfield(options,'shift0')
      if ~isreal(options.shift0)
	 shiftreal=0;
      end
   end
   if isfield(options,'shiftmax')
      if ~isreal(options.shiftmax)
	 shiftreal=0;
      end
   end
   if isfield(options,'shifts')
      if ~isreal(options.shifts)
	 shiftreal=0;
      end
   end

   if (norm(options.shiftmatrix-options.shiftmatrix',1)==0)
      shifthermitian=1;
      % make sure that the shifts are real
      if isfield(options,'shift0')
	 if ~isreal(options.shift0)
	    shifthermitian=0;
	 end
      end
      if isfield(options,'shiftmax')
	 if ~isreal(options.shiftmax)
	    shifthermitian=0;
	 end
      end
      if isfield(options,'shifts')
	 if ~isreal(options.shifts)
	    shifthermitian=0;
	 end
      end
   else
      shifthermitian=0;
   end
   if (norm(options.shiftmatrix-options.shiftmatrix.',1)==0)
      shiftsymmetric=1;
   else
      shiftsymmetric=0;
   end
   
   % A and shift do not match
   if isreal(A) & ~shiftreal
      % A becomes complex hermitian
      if norm(A-A',1)==0 & shifthermitian
	 A(1,2)=A(1,2)+norm(A,1)*eps^2*sqrt(-1);
	 A(1,2)=A(2,1)';
      % A becomes complex symmetric
      elseif norm(A-A',1)==0 & shiftsymmetric
	 A(1,1)=A(1,1)+norm(A,1)*eps^2*sqrt(-1);
      % A becomes complex unsymmetric
      else
	 A(1,2)=A(1,2)+norm(A,1)*eps^2*sqrt(-1);
      end
   elseif isreal(A) & shiftreal
      % A becomes unsymmetric
      if norm(A-A',1)==0 & ~shiftsymmetric
	 A(1,2)=A(1,2)+norm(A,1)*eps^2;
      end
   elseif ~isreal(A)
      % A becomes complex hermitian
      if norm(A-A',1)==0 & shifthermitian
	 % perfect
      % A becomes complex symmetric
      elseif norm(A-A.',1)==0 & shiftsymmetric
	 % perfect
      % A becomes complex unsymmetric
      else
	 A(1,2)=A(1,2)+norm(A,1)*eps^2*sqrt(-1);
      end
   end
   
end

if nargin==1
   options = AMGinit(A);
elseif nargin==2
   if isfield(options, 'isdefinite')
      myoptions.isdefinite=options.isdefinite;
      myoptions = AMGinit(A,myoptions);
   else
      myoptions = AMGinit(A);
   end

   % complete missing data
   if ~isfield(options, 'matching')
      options.matching=myoptions.matching;
   end
   if ~isfield(options, 'ordering')
      options.ordering=myoptions.ordering;
   end
   if ~isfield(options, 'droptol')
      options.droptol=myoptions.droptol;
   end
   if ~isfield(options, 'droptolS')
      options.droptolS=myoptions.droptolS;
   end
   if ~isfield(options, 'condest')
      options.condest=myoptions.condest;
   end
   if ~isfield(options, 'elbow')
      options.elbow=myoptions.elbow;
   end
   if ~isfield(options, 'lfil')
      options.lfil=myoptions.lfil;
   end
   if ~isfield(options, 'lfilS')
      options.lfilS=myoptions.lfilS;
   end
   
   if ~isfield(myoptions,'nrestart')
      options.nrestart=myoptions.nrestart;
   end

   if ~isfield(options,'typetv')
      options.typetv=myoptions.typetv;
   end
   if ~isfield(options,'tv')
      options.tv=myoptions.tv;
   end
   if ~isfield(options,'ind')
      options.ind=myoptions.ind;
   end
   if ~isfield(options,'mixedprecision')
      options.mixedprecision=myoptions.mixedprecision;
   end
   if ~isfield(options,'amg')
      options.amg=myoptions.amg;
   end
   if ~isfield(options,'npresmoothing')
      options.npresmoothing=myoptions.npresmoothing;
   end
   if ~isfield(options,'npostsmoothing')
      options.npostsmoothing=myoptions.npostsmoothing;
   end
   if ~isfield(options,'ncoarse')
      options.ncoarse=myoptions.ncoarse;
   end
   if ~isfield(options,'presmoother')
      options.presmoother=myoptions.presmoother;
   end
   if ~isfield(options,'postsmoother')
      options.postsmoother=myoptions.postsmoother;
   end
   if ~isfield(options,'FCpart')
      options.FCpart=myoptions.FCpart;
   end
   if ~isfield(options,'typecoarse')
      options.typecoarse=myoptions.typecoarse;
   end
   if ~isfield(options, 'solver')
      options.solver=myoptions.solver;
   end
   if ~isfield(options, 'damping')
      options.damping=myoptions.damping;
   end

else
   error('wrong number of input arguments');
end % if 

if nargout~=2
   error('wrong number of output arguments');
end % if 

myoptions=options;

if isreal(A)
   myoptions.isreal=1;
   if norm(A-A',1)==0
      myoptions.issymmetric=1;
      myoptions.ishermitian=1;
   else
      myoptions.issymmetric=0;
      myoptions.ishermitian=0;
   end
else
   myoptions.isreal=0;
   if norm(A-A',1)==0
      myoptions.ishermitian=1;
   else
      myoptions.ishermitian=0;
   end
   if ~myoptions.ishermitian
      if norm(A-A.',1)==0
	 myoptions.issymmetric=1;
      else
	 myoptions.issymmetric=0;
      end
   else
      myoptions.issymmetric=0;
   end
end


if ~isfield(options,'isdefinite')
   myoptions.isdefinite=0;
end


if myoptions.isreal
   if myoptions.issymmetric
      if myoptions.isdefinite
	 if ~strcmp(myoptions.typetv,'static') & ~strcmp(myoptions.typetv,'none')
	    rcomflag=-1;
	    tv=myoptions.tv;
	    PREC=0;
	    while rcomflag~=0
	          [PREC, myoptions,rcomflag,S,tv]=DSPDilupackfactor(A,myoptions,PREC,tv);
		  if rcomflag~=0
		     D=spdiags(diag(S),0,size(S,1),size(S,1));
		     S=(S-D)+S';
		     tv=feval(myoptions.typetv,S,tv);
		  end % if
	    end % while
	 else
	    [PREC, myoptions]=DSPDilupackfactor(A,myoptions);
	 end % if
	 if ~strcmp(myoptions.amg,'ilu')
	    for i=1:length(PREC)-2
	        n=size(PREC(i).A_H,1);
		D=spdiags(diag(PREC(i).A_H),0,n,n);
		PREC(i).A_H=(PREC(i).A_H-D)+PREC(i).A_H';
		Dr=spdiags(1./PREC(i+1).rowscal.',0,n,n);
		Dc=spdiags(1./PREC(i+1).colscal.',0,n,n);
		PREC(i).A_H=Dr*PREC(i).A_H*Dc;
	    end
	    i=length(PREC)-1;
	    n=0;
	    if i>0
	       n=size(PREC(i).A_H,1);
	    end
	    if n>0
	       D=spdiags(diag(PREC(i).A_H),0,n,n);
	       PREC(i).A_H=(PREC(i).A_H-D)+PREC(i).A_H';
	       Dr=spdiags(1./PREC(i+1).rowscal.',0,n,n);
	       Dc=spdiags(1./PREC(i+1).colscal.',0,n,n);
	       PREC(i).A_H=Dr*PREC(i).A_H*Dc;
	    elseif i>0
	       PREC(i).A_H=PREC(i).L*inv(PREC(i).D)*PREC(i).L';
	    end % if-else
	 end % if
      else
	 if ~strcmp(myoptions.typetv,'static') & ~strcmp(myoptions.typetv,'none')
	    rcomflag=-1;
	    tv=myoptions.tv;
	    PREC=0;
	    while rcomflag~=0
	          [PREC, myoptions,rcomflag,S,tv]=DSYMilupackfactor(A,myoptions,PREC,tv);
		  if rcomflag~=0
		     D=spdiags(diag(S),0,size(S,1),size(S,1));
		     S=(S-D)+S';
		     tv=feval(myoptions.typetv,S,tv);
		  end % if
	    end % while
	 else
	    [PREC, myoptions]=DSYMilupackfactor(A,myoptions);
	 end % if
	 % adapt final decomposition
	 nlev=length(PREC);
	 if ~issparse(PREC(nlev).L)
	    PREC(nlev).L=PREC(nlev).L*PREC(nlev).D;
	 end
	 if ~strcmp(myoptions.amg,'ilu')
	    for i=1:length(PREC)-2
	        n=size(PREC(i).A_H,1);
		D=spdiags(diag(PREC(i).A_H),0,n,n);
		PREC(i).A_H=(PREC(i).A_H-D)+PREC(i).A_H';
		Dr=spdiags(1./PREC(i+1).rowscal.',0,n,n);
		Dc=spdiags(1./PREC(i+1).colscal.',0,n,n);
		PREC(i).A_H=Dr*PREC(i).A_H*Dc;
	    end
	    i=length(PREC)-1;
	    n=0;
	    if i>0
	       n=size(PREC(i).A_H,1);
	    end
	    if n>0
	       D=spdiags(diag(PREC(i).A_H),0,n,n);
	       PREC(i).A_H=(PREC(i).A_H-D)+PREC(i).A_H';
	       Dr=spdiags(1./PREC(i+1).rowscal.',0,n,n);
	       Dc=spdiags(1./PREC(i+1).colscal.',0,n,n);
	       PREC(i).A_H=Dr*PREC(i).A_H*Dc;
	    elseif i>0
	       PREC(i).A_H=PREC(i).L*inv(PREC(i).D)*PREC(i).L';
	    end % if-else
	 end % if
      end % if
   else
      [PREC, myoptions]=DGNLilupackfactor(A,myoptions);
      if ~strcmp(myoptions.amg,'ilu')
	 for i=1:length(PREC)-1
	    n=size(PREC(i).A_H,1);
	    if n>0
	       Dr=spdiags(1./PREC(i+1).rowscal.',0,n,n);
	       Dc=spdiags(1./PREC(i+1).colscal.',0,n,n);
	       PREC(i).A_H=Dr*PREC(i).A_H*Dc;
	    end
	 end
	 i=length(PREC)-1;
	 n=0;
	 if i>0
	    n=size(PREC(i).A_H,1);
	 end
	 if n==0 & i>0
	    PREC(i).A_H=PREC(i).L*inv(PREC(i).D)*PREC(i).U;
	 end % if
      end
   end
else
   if myoptions.ishermitian
      if myoptions.isdefinite
	 if ~strcmp(myoptions.typetv,'static') & ~strcmp(myoptions.typetv,'none')
	    rcomflag=-1;
	    tv=myoptions.tv;
	    PREC=0;
	    while rcomflag~=0
	          [PREC, myoptions,rcomflag,S,tv]=ZHPDilupackfactor(A,myoptions,PREC,tv);
		  if rcomflag~=0
		     D=spdiags(diag(S),0,size(S,1),size(S,1));
		     S=(S-D)+S';
		     tv=feval(myoptions.typetv,S,tv);
		  end % if
	    end % while
	 else
	    [PREC, myoptions]=ZHPDilupackfactor(A,myoptions);
	 end % if
	 if ~strcmp(myoptions.amg,'ilu')
	    for i=1:length(PREC)-2
	        n=size(PREC(i).A_H,1);
		D=spdiags(diag(PREC(i).A_H),0,n,n);
		PREC(i).A_H=(PREC(i).A_H-D)+PREC(i).A_H';
		Dr=spdiags(1./PREC(i+1).rowscal',0,n,n);
		Dc=spdiags(1./PREC(i+1).colscal',0,n,n);
		PREC(i).A_H=Dr*PREC(i).A_H*Dc;
	    end
	    i=length(PREC)-1;
	    n=0;
	    if i>0
	       n=size(PREC(i).A_H,1);
	    end
	    if n>0
	       D=spdiags(diag(PREC(i).A_H),0,n,n);
	       PREC(i).A_H=(PREC(i).A_H-D)+PREC(i).A_H';
	       Dr=spdiags(1./PREC(i+1).rowscal.',0,n,n);
	       Dc=spdiags(1./PREC(i+1).colscal.',0,n,n);
	       PREC(i).A_H=Dr*PREC(i).A_H*Dc;
	    elseif i>0
	       PREC(i).A_H=PREC(i).L*inv(PREC(i).D)*PREC(i).L';
	    end % if-else
	 end % if
      else
	 if ~strcmp(myoptions.typetv,'static') & ~strcmp(myoptions.typetv,'none')
	    rcomflag=-1;
	    tv=myoptions.tv;
	    PREC=0;
	    while rcomflag~=0
	          [PREC, myoptions,rcomflag,S,tv]=ZHERilupackfactor(A,myoptions,PREC,tv);
		  if rcomflag~=0
		     D=spdiags(diag(S),0,size(S,1),size(S,1));
		     S=(S-D)+S';
		     tv=feval(myoptions.typetv,S,tv);
		  end % if
	    end % while
	 else
	    [PREC, myoptions]=ZHERilupackfactor(A,myoptions);
	 end % if
	 if ~strcmp(myoptions.amg,'ilu')
	    for i=1:length(PREC)-2
	        n=size(PREC(i).A_H,1);
		D=spdiags(diag(PREC(i).A_H),0,n,n);
		PREC(i).A_H=(PREC(i).A_H-D)+PREC(i).A_H';
		Dr=spdiags(1./PREC(i+1).rowscal',0,n,n);
		Dc=spdiags(1./PREC(i+1).colscal',0,n,n);
		PREC(i).A_H=Dr*PREC(i).A_H*Dc;
	    end
	    i=length(PREC)-1;
	    n=0;
	    if i>0
	       n=size(PREC(i).A_H,1);
	    end
	    if n>0
	       D=spdiags(diag(PREC(i).A_H),0,n,n);
	       PREC(i).A_H=(PREC(i).A_H-D)+PREC(i).A_H';
	       Dr=spdiags(1./PREC(i+1).rowscal',0,n,n);
	       Dc=spdiags(1./PREC(i+1).colscal',0,n,n);
	       PREC(i).A_H=Dr*PREC(i).A_H*Dc;
	    elseif i>0
	       PREC(i).A_H=PREC(i).L*inv(PREC(i).D)*PREC(i).L';
	    end % if-else
	 end % if
      end % if
   elseif myoptions.issymmetric
      if ~strcmp(myoptions.typetv,'static') & ~strcmp(myoptions.typetv,'none')
	 rcomflag=-1;
	 tv=myoptions.tv;
	 PREC=0;
	 while rcomflag~=0
	       [PREC, myoptions,rcomflag,S,tv]=ZSYMilupackfactor(A,myoptions,PREC,tv);
	       if rcomflag~=0
		  D=spdiags(diag(S),0,size(S,1),size(S,1));
		  S=(S-D)+S.';
		  tv=feval(myoptions.typetv,S,tv);
	       end % if
         end % while
      else
	 [PREC, myoptions]=ZSYMilupackfactor(A,myoptions);
      end % if
      if ~strcmp(myoptions.amg,'ilu')
	 for i=1:length(PREC)-2
	    n=size(PREC(i).A_H,1);
	    D=spdiags(diag(PREC(i).A_H),0,n,n);
	    PREC(i).A_H=(PREC(i).A_H-D)+PREC(i).A_H.';
	    Dr=spdiags(1./PREC(i+1).rowscal.',0,n,n);
	    Dc=spdiags(1./PREC(i+1).colscal.',0,n,n);
	    PREC(i).A_H=Dr*PREC(i).A_H*Dc;
	 end
	 i=length(PREC)-1;
	 n=0;
	 if i>0
	    n=size(PREC(i).A_H,1);
	 end
	 if n>0
	    D=spdiags(diag(PREC(i).A_H),0,n,n);
	    PREC(i).A_H=(PREC(i).A_H-D)+PREC(i).A_H.';
	    Dr=spdiags(1./PREC(i+1).rowscal.',0,n,n);
	    Dc=spdiags(1./PREC(i+1).colscal.',0,n,n);
	    PREC(i).A_H=Dr*PREC(i).A_H*Dc;
	 elseif i>0
	    PREC(i).A_H=PREC(i).L*inv(PREC(i).D)*PREC(i).L.';
	 end % if-else
      end % if
   else
      [PREC, myoptions]=ZGNLilupackfactor(A,myoptions);
      if ~strcmp(myoptions.amg,'ilu')
	 for i=1:length(PREC)-1
	    n=size(PREC(i).A_H,1);
	    if n>0
	       Dr=spdiags(1./PREC(i+1).rowscal.',0,n,n);
	       Dc=spdiags(1./PREC(i+1).colscal.',0,n,n);
	       PREC(i).A_H=Dr*PREC(i).A_H*Dc;
	    end
	 end
	 i=length(PREC)-1;
	 n=0;
	 if i>0
	    n=size(PREC(i).A_H,1);
	 end
	 if n==0 & i>0
	    PREC(i).A_H=PREC(i).L*inv(PREC(i).D)*PREC(i).U;
	 end % if
      end
   end
end % if


% update options
options.matching      =myoptions.matching;
options.ordering      =myoptions.ordering;
options.droptol       =myoptions.droptol;
options.droptolS      =myoptions.droptolS;
options.condest       =myoptions.condest;
options.elbow         =myoptions.elbow;
options.lfil          =myoptions.lfil;
options.lfilS         =myoptions.lfilS;
options.typetv        =myoptions.typetv;
options.amg           =myoptions.amg;
options.FCpart        =myoptions.FCpart;
options.mixedprecision=myoptions.mixedprecision;
options.typecoarse    =myoptions.typecoarse;
