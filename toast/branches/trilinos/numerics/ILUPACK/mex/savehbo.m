function savehbo(filename, A,b,x0,x)
% savehbo(filename, A)
% savehbo(filename, A,b)
% savehbo(filename, A,b,x0)
% savehbo(filename, A,b,x0,x)
% 
% save matrix A and optionally right hand side b, initial guess x0 and
% exact solution x


if (nargin<2)
   error('at least two input parameters required!');
end
[m,n]=size(A);

if nargin>=3
   [p,q]=size(b);
   if p~=m
      error('b must have as many rows as A');
   end
end
if nargin>=4
   [r,s]=size(x0);
   if r~=m
      error('x0 must have as many rows as A');
   end
   if s~=q
      error('x0 must have as many columns as b');
   end
end
if nargin>=5
   [r,s]=size(x);
   if r~=m
      error('x must have as many rows as A');
   end
   if s~=q
      error('x must have as many columns as b');
   end
end

if isreal(A)
   if norm(A-A',1)==0
      if nargin==2
	 DSYMsavehbo(filename, A);      
      elseif nargin==3
	 DSYMsavehbo(filename, A,b);      
      elseif nargin==4
	 DSYMsavehbo(filename, A,b,x0);      
      elseif nargin==5
	 DSYMsavehbo(filename, A,b,x0,x);      
      end % if
   else
      if nargin==2
	 DGNLsavehbo(filename, A);      
      elseif nargin==3
	 DGNLsavehbo(filename, A,b);      
      elseif nargin==4
	 DGNLsavehbo(filename, A,b,x0);      
      elseif nargin==5
	 DGNLsavehbo(filename, A,b,x0,x);      
      end % if
   end
else % A is complex
   if norm(A-A',1)==0
      if nargin==2
	 ZHERsavehbo(filename, A);      
      elseif nargin==3
	 ZHERsavehbo(filename, A,b);      
      elseif nargin==4
	 ZHERsavehbo(filename, A,b,x0);      
      elseif nargin==5
	 ZHERsavehbo(filename, A,b,x0,x);      
      end % if
   elseif norm(A-A.',1)==0
      if nargin==2
	 ZSYMsavehbo(filename, A);      
      elseif nargin==3
	 ZSYMsavehbo(filename, A,b);      
      elseif nargin==4
	 ZSYMsavehbo(filename, A,b,x0);      
      elseif nargin==5
	 ZSYMsavehbo(filename, A,b,x0,x);      
      end % if
   else
      if nargin==2
	 ZGNLsavehbo(filename, A);      
      elseif nargin==3
	 ZGNLsavehbo(filename, A,b);      
      elseif nargin==4
	 ZGNLsavehbo(filename, A,b,x0);      
      elseif nargin==5
	 ZGNLsavehbo(filename, A,b,x0,x);      
      end % if
   end
end % if
