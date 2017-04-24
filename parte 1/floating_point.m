% MAC0420 - EP1 - Parte 1: Aritmética e Ponto Flutuante
%
% Arthur Coser Marinho               NUSP: 7210629 
% Ludmila Ferreira Vicente e Silva   NUSP: 7557136
%
%! /bin/octave -qf

% transform a binarry array into an integer
function retval = binary_array_to_int (x)
	s = size(x)(2);
	retval = 0;
	for i = 1:s
		if (x(i) == 1)
			retval += pow2(s-i);
		endif
	endfor
endfunction

% ulp = 2^-23
function retval = x_plus_ulp (x)
	retval = x(1:32);
	carry = 1;
	i = 32;
	while (carry == 1 && i > 9)
		sum = retval(i) + carry;
		carry = (sum > 1);
		retval(i) = mod(sum,2);
		i--;
	endwhile
	if (i == 9)
		carry = 1;
		while (carry == 1 && i > 1)
			sum = retval(i) + carry;
			carry = (sum > 1);
			retval(i) = mod(sum,2);
			i--;
		endwhile
	endif
endfunction

% round the number according to the chosen method
function retval = round_number (x)
	global method;
	sign = x(1);
	if (sign == 0)
		xminus = x;
		xplus = x_plus_ulp(x);
	else
		xminus = x_plus_ulp(x);
		xplus = x;
	endif
	if (strcmp(method, 'down'))
		retval = xminus;
	elseif (strcmp(method, 'up'))
		retval = xplus;
	elseif (strcmp(method, 'towards zero'))
		retval = x;
	else
		if (x(32+1) == 1) % guard bit = 1
			if (x(32+2) == 0 && x(32+3) == 0) %tie
				if (xminus(32) == 0)
					retval = x;
				else
					retval = x_plus_ulp(x);
				endif
			else 
				retval = x_plus_ulp(x);
			endif
		else
			retval = x;
		endif
	endif
	retval = retval(1:32);
endfunction

% simulates the ieee float representation with an array 
% [1: sign, 2:9 expoent, 10:32 mantissa, 33:34 guard bits, 35: sticky bit]
function retval = ieee_float (x)
	if (x == 0)
		retval = zeros(1,32+3);
		return;
	endif

	sign = 0;
	if (x < 0)
		x = -x;
		sign = 1;
	endif

	aux = x;
	e = 0;
	if (x >= 1)
		while (aux >= 2)
			aux /= 2;
			e++;
		endwhile
	else
		while (aux < 1)
			aux *= 2;
			e--;
		endwhile
	endif

	if (e > 127)
		expoent(1:8) = 1;
		mantissa(1:23+3) = 0;
		retval = [sign, expoent, mantissa];
	else
		aux = e + 127;
		expoent(1:8) = 0;
		if (aux > 0)
			for i = 1:8
				p = pow2(8 - i);
				if (aux - p >= 0)
					aux = aux - p;
					expoent(i) = 1;
				endif
			endfor
		endif
		x = x - pow2(e); % accounting for the hidden bit = 1
		mantissa(1:23+3) = 0; % extra 3 bits are the 2 guard bits and the sticky bit
		floating_point = false; % it will be used to check if the number is a floating point
		for i = 1:23
			p = pow2(e - i);
			if (x - p >= 0)
				x = x - p;
				mantissa(i) = 1;
			endif
			if (x == 0)
				floating_point = true;
			endif
		endfor
		p = pow2(e - (23+1));
		if (x - p >= 0)
			x = x - p;
			mantissa(32+1) = 1;
		endif
		if (x == 0)
			mantissa(32+2) = 0;
		endif

		% subnormal (expoent = [0,0,0,0,0,0,0,0])
		if (e < -126)
			mantissa = shift_right(x, -(e+126), 1);
		endif

		retval = [sign, expoent, mantissa];

		% if floating_point == true, the number is a floating point (it's perfectly represented and round(x) = x)
		% if floating_point == false, we need to round it
		if (!floating_point)
			retval = round_number(retval);
		endif
	endif
endfunction

% shift only the mantissa, r spaces to the right, index 33 and 34 are 'guard bits', index 35 is the 'sticky bit'
function retval = shift_right (x, r, hidden_bit)
	if (r > 23+3)
		x(32+3) = any([hidden_bit, x(10:32+3)]);
		x(10:32+2) = 0;
	elseif (r > 0)
		% sticky bit
		x(32+3) = any(x(10+r:32+3));
		% shift
		x(10+r:32+2) = x(10:32+2-r);
		x(10+r-1) = hidden_bit;
		x(10:10+r-2) = 0;
	endif
	retval = x;
endfunction

% perform the addition or the subtraction
function retval = round_operation (option, x, y)
	% guard bits: index 33 and 34 | sticky bit: index 35
	x(33:35) = 0;
	y(33:35) = 0;

	ex = binary_array_to_int(x(2:9)) - 127;
	ey = binary_array_to_int(y(2:9)) - 127;

	hidden_bit_x = 1;
	hidden_bit_y = 1;
	if (ex < -126)
		hidden_bit_x = 0;
	endif
	if (ey < -126)
		hidden_bit_y = 0;
	endif
	
	r = ex - ey;

	if (r > 0)
		y = shift_right(y, r, hidden_bit_y);
		y(2:9) = x(2:9);
		hidden_bit_y = 0;
	elseif (r < 0)
		x = shift_right(x, -r, hidden_bit_x);
		x(2:9) = y(2:9);
		hidden_bit_x = 0;
	endif

	sign_x = x(1);
	sign_y = y(1);

	if (option == 1)
		if (sign_x == 1 && sign_y == 0) % (-)x + (+)y = (+)y - (+)x
			retval = subtraction(y, hidden_bit_y, [0, x(2:35)], hidden_bit_x);
		elseif (sign_x == 0 && sign_y == 1) % (+)x + (-)y = (+)x - (+)y
			retval = subtraction(x, hidden_bit_x, [0, y(2:35)], hidden_bit_y);
		else
			retval = addition(x, hidden_bit_x, y, hidden_bit_y);
		endif
	else
		if (sign_x == 1 && sign_y == 0) % (-)x - (+)y = (-)x + (-)y
			retval = addition(x, hidden_bit_x, [1, y(2:35)], hidden_bit_y);
		elseif (sign_x == 0 && sign_y == 1) % (+)x - (-)y = +(x) + (+)y
			retval = addition(x, hidden_bit_x, [0, y(2:35)], hidden_bit_y);
		elseif (sign_x == 1 && sign_y == 1) % (-)x - (-)y = (+)y - (+)x
			retval = subtraction([0, y(2:35)], hidden_bit_y, [0, y(2:35)], hidden_bit_x);
		else
			retval = subtraction(x, hidden_bit_x, y, hidden_bit_y);
		endif
	endif

	% operation is done, now round it if it's not floating point
	if (retval(33) == 1 || retval(34) == 1 || retval(35) == 1) % result is not a floating point
		retval = round_number(retval);
	endif
endfunction

% addition of two ieee arrays (assumes x and y have the same sign)
function retval = addition (x, hidden_bit_x, y, hidden_bit_y)
	retval = x;
	carry = 0;

	i = 35;
	while (i > 9)
		sum = x(i) + y(i) + carry;
		carry = (sum > 1);
		retval(i) = mod(sum, 2);
		i--;
	endwhile

	if (hidden_bit_x == 1 || hidden_bit_y == 1) % else hidden_bit of the result will be 'carry', must add 'carry' to expoent
		sum = hidden_bit_x + hidden_bit_y + carry;
		carry = (sum > 1);
		if (sum > 1)
			retval = shift_right(retval, 1, sum - 2);
		endif
	endif

	% plus 1 to the expoent if carry == 1
	while (carry == 1 && i > 1)
		sum = retval(i) + carry;
		carry = (sum > 1);
		retval(i) = mod(sum,2);
		i--;
	endwhile
	if (i == 1) % expoent equals 11111111
		retval(10:35) = 0;
	endif
endfunction

% subtraction of two ieee arrays (assumes x and y are positive and expoent of both are the same)
function retval = subtraction (x, hidden_bit_x, y, hidden_bit_y)
	retval = x;
	borrow = 0;

	% if x > y, we do x - y and mantain the sign
	% if y > x, we do y - x and change the sign
	switch_x_y = false;
	if (hidden_bit_x == hidden_bit_y)
		for i = 10:32
			if (x(i) > y(i)) % x > y
				break;
			elseif (y(i) > x(i)) % y > x
				switch_x_y = true;
				break;
			endif
			i++;
		endfor
	elseif (hidden_bit_x == 0 && hidden_bit_y == 1)
		switch_x_y = true;
	endif
	if (switch_x_y)
		retval = y;
		aux = hidden_bit_y;
		y = x;
		hidden_bit_y = hidden_bit_x;
		x = retval;
		hidden_bit_x = aux;
		retval(1) = 1; % sign change to negative
	endif

	% sticky
	retval(35) = any([x(35), y(35)]);
	borrow = (x(35) == 0 && y(35) == 1);

	i = 34;
	while (i > 9)
		sub = (x(i) - borrow) - y(i);
		borrow = (sub < 0); 
		retval(i) = mod(sub, 2);
		i--;
	endwhile

	% if the numbers are equal, subtraction must be zero
	if (hidden_bit_x == hidden_bit_y && !any(retval(10:34)))
		retval(2:35) = 0;
		return;
	endif

	% hidden_bit_x == 0 and hidden_bit_y == 1 is not a possibility
	% hidden_bit_x == 1 and hidden_bit_y == 1 and borrow == 1 is not a possibility
	hidden_bit = (hidden_bit_x - borrow) - hidden_bit_y;

	subnormal = !any(retval(2:9));
	if (hidden_bit == 0 && !subnormal) % must shift left
		r = 1; % r is the number of shifts to the left
		while (retval(r+9) == 0 && r+9 <= 35)
			r++;
		endwhile
		for i = 1:r
			% shift one to the left
			retval(10:34) = retval(11:35); 
			retval(35) = 0;
			% subtracts 1 from expoent
			borrow = 1;
			j = 9;
			while(j > 1)
				sub = retval(j) - borrow;
				borrow = (sub < 0);
				retval(j) = mod(sub, 2);
				j--;
			endwhile
			if (borrow) % expoent equals to [0,0,0,0,0,0,0,0], number is subnormal, stop shifting to the left
				break;
			endif
		endfor
	endif
endfunction

% print the array of the ieee floating point as decimal
function print_ieee_number (x)
	if (!any(x(2:32)))
		result = 0;
	else 
		e = binary_array_to_int(x(2:9)) - 127;
		result = pow2(e);
		for i = 10:32
			if (x(i) == 1)
				result += pow2(e-(i-9));
			endif
		endfor
		if (x(1) == 1)
			result *= -1;
		endif
	endif
	printf('%e\n', result);
endfunction

% print the array of the ieee floating point as binary
function print_ieee_binary (x)
	for i = 1:32
		if (i == 2 || i == 10)
			printf('|');
		endif;
		printf('%d', x(i));
	endfor
	printf('\n');
endfunction

% print result
function print_result (option, x, y, result)
	if (option == 1)
		operator = '+';
	else
		operator = '-';
	endif
	printf('----------------------------------------\n')
	printf('x   = '); print_ieee_number(x);
	printf('y   = '); print_ieee_number(y);
	printf('x%cy = ', operator); print_ieee_number(result);
	printf('----------------------------------------\n')
	printf('x   = '); print_ieee_binary(x);
	printf('y   = '); print_ieee_binary(y);
	printf('x%cy = ', operator);print_ieee_binary(result);
	printf('----------------------------------------\n')
endfunction

% global variable for the round method chosen by the user
global method;

% start of the program, asks user for the round method
printf('MAC0210 EP1 - Part 1: Floating Point Arithmetic\n\nSelect the round method:\n    [0] Exit\n    [1] Round down\n    [2] Round up\n    [3] Round towards zero\n    [4] Round to nearest\n');
option = input('> ');
while (option != 0:4)
	printf('Invalid option.\n');
	option = input('> ');
endwhile
if (option == 0)
	return;
elseif (option == 1)
	method = 'down';
elseif (option == 2)
	method = 'up';
elseif (option == 3)
	method = 'towards zero';
else
	method = 'to nearest';
endif

% examples for report
printf('\nExamples for the report:\n\n');

printf('a) 2 + 3\n');
x = ieee_float(2);
y = ieee_float(3);
result = round_operation(1, x, y);
print_result(1, x, y, result);
printf('\n\n');

printf('b) 1 + 2^-24\n');
x = [0, [0,1,1,1,1,1,1,1], zeros(1,23)]; % 1
y = [0, [0,1,1,0,0,1,1,1], zeros(1,23)]; % (1.0) * 2^-24
result = round_operation(1, x, y);
print_result(1, x, y, result);
printf('\n\n');

printf('c) 1.0 - (1.1...1)*2^-1\n');
x = [0, [0,1,1,1,1,1,1,1], zeros(1,23)]; % 1
y = [0, [0,1,1,1,1,1,1,0], ones(1,23)]; % (1.1...1) * 2^-1
result = round_operation(2, x, y);
print_result(2, x, y, result);
printf('\n\n');

printf('d) 1.0 - (1.0...01)*2^-25\n');
x = [0, [0,1,1,1,1,1,1,1], zeros(1,23)]; % 1
y = [0, [0,1,1,0,0,1,1,0], zeros(1,22), 1]; % (1.0...01) * 2^-25
result = round_operation(2, x, y);
print_result(2, x, y, result);
printf('\n\n');

% asks user for the operation (addition or subtraction)
printf('Using the "round %s" method.\nSelect the operation:\n    [0] Exit\n    [1] Addition (x + y)\n    [2] Subtraction (x - y)\n', method);
option = input('> ');
while (option != 0:2)
	printf('Invalid option.\n');
	option = input('> ');
endwhile
while (option != 0)
	x = input('x = ');
	while (!isscalar(x))
		printf('Invalid value.\n');
		x = input('x = ');
	endwhile
	y = input('y = ');
	while (!isscalar(y))
		printf('Invalid value.\n');
		y = input('y = ');
	endwhile

	x = ieee_float(x);
	y = ieee_float(y);
	
	result = round_operation(option, x, y);

	print_result(option, x, y, result);

	printf('\nUsing the "round %s" method.\nSelect the operation:\n    [0] Exit\n    [1] Addition (x + y)\n    [2] Subtraction (x - y)\n', method);
	option = input('> ');
	while (option != 0:2)
		printf('Invalid option.\n');
		option = input('> ');
	endwhile
endwhile
return;
