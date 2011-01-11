function [mua, mus, ref] = assignMatValues(p, t, mua3, mus3, ref3)
	mua = zeros(size(t, 1), 1);
	mus = zeros(size(t, 1), 1);
	ref = zeros(size(t, 1), 1);
	tol = 1e-6;
	for i = 1 : size(t, 1)
		a = sqrt(p(t(i ,1), 1)^2 + p(t(i, 1), 2)^2 + p(t(i, 1), 3)^2);
		b = sqrt(p(t(i ,2), 1)^2 + p(t(i, 2), 2)^2 + p(t(i, 2), 3)^2);
		c = sqrt(p(t(i ,3), 1)^2 + p(t(i, 3), 2)^2 + p(t(i, 3), 3)^2);
		d = sqrt(p(t(i ,4), 1)^2 + p(t(i, 4), 2)^2 + p(t(i, 4), 3)^2);
		
		%if a - 10 <= tol & b - 10 <= tol & c - 10 <= tol & d - 10 <= tol
		%	mua(i, 1) = mua3(3);
		%	mus(i, 1) = mus3(3);
		%	ref(i, 1) = ref3(3);
		%	continue;
		%elseif a-15 <= tol & b-15 <= tol & c-15 <= tol & d-15 <= tol 
		%	mua(i, 1) = mua3(2);
		%	mus(i, 1) = mus3(2);
		%	ref(i, 1) = ref3(2); 
		%	continue;
		%elseif a-20 <= tol & b-20 <= tol & c-20 <= tol & d-20 <= tol 
		%	mua(i, 1) = mua3(1);
		%	mus(i, 1) = mus3(1);
		%	ref(i, 1) = ref3(1);
		%	continue;
		%else
		%   	i
		%	disp('dunno what to do ...')
		%	pause	
		%end;
		%if fix(a) <= 10 & fix(b) <= 10 & fix(c) <= 10 & fix(d) <= 10
		%	mua(i, 1) = mua3(3);
		%	mus(i, 1) = mus3(3);
		%	ref(i, 1) = ref3(3);
		%	continue;
		%elseif fix(a) <= 15  & fix(b) <= 15  & fix(c) <= 15 & fix(d) <= 15 
		%	mua(i, 1) = mua3(2);
		%	mus(i, 1) = mus3(2);
		%	ref(i, 1) = ref3(2); 
		%	continue;
		%elseif fix(a) <= 20  & fix(b) <= 20 & fix(c) <= 20 & fix(d) <= 20 
		%	mua(i, 1) = mua3(1);
		%	mus(i, 1) = mus3(1);
		%	ref(i, 1) = ref3(1);
		%	continue;
		%else
		%   	i
		%	disp('dunno what to do ...')
		%	pause	
		%end;

		if a <= 10 | (b) <= 10 | (c) <= 10 | (d) <= 10
			mua(i, 1) = mua3(3);
			mus(i, 1) = mus3(3);
			ref(i, 1) = ref3(3);
			continue;
		elseif (a) <= 15  | (b) <= 15  | (c) <= 15 | (d) <= 15 
			mua(i, 1) = mua3(2);
			mus(i, 1) = mus3(2);
			ref(i, 1) = ref3(2);
			continue;
		elseif (a) <= 20  | (b) <= 20 | (c) <= 20 | (d) <= 20 
			mua(i, 1) = mua3(1);
			mus(i, 1) = mus3(1);
			ref(i, 1) = ref3(1);
			continue;
		else
		   	[i a b c d]
			
			disp('dunno what to do ...')
			pause	
		end;

		
	  
	end;
