function [coeff, weights, IndexF] = tri(arg1, arg2)
% This function creates integration points and weights for triangles 
% Syntax:
% [coeff, weights, IndexF] = tri(3) -   will do barycentric subdivision with
%                                       3*3 = 9 integration points
% [coeff, weights, IndexF] = tri(7,5) - will use Gaussian quadrature of fifth 
%                                       order with seven integration points
% Description:
% Uses Gaussian quadratures (with author's permission) from 
% http://www.cs.kuleuven.ac.be/~nines/research/ecf/ecf.html
% Formulas of third (4 integration points), fifth (7 integration points)
% seventh (13 integration points), and tenth (25 integration points) order
% of accuracy may be created.
%
% Uses the "edge" method for barycentric subdivision of arbitrary order,
% where the edges of smaller triangles (similar to the original one) are
% equally subdivided. This gives the desired barycentric points.

%   Copyright SNM 2002-2020

    nargs = nargin; 

    if  nargs < 1
        error('Requires at least 1 input');
    end

    if nargs == 1
    %   Barycentric triangle subdivision - coefficients for vertexes only
    %   M - subdivision order (number of subtriangles is M*M)
        M = arg1;	
        if M < 2 
            coeff(:,1) = [1/3 1/3 1/3]';
            weights(1) =   1;
        end

        coeff   = zeros(3, M*M);
        k       = 1;
        scale   = 1;
        eps     = 2 + 1e-9;   % scaling

        if mod(M,3) == 0
            N   =   M/3*2;
        elseif mod(M,3) == 2
            N   =   M/3*2 - 1/3;
        else
            N   =   (M-1)/3*2;
        end    

        %   Border loop - starts with the outer border of integration points
        %   and then goes inside - "triangle" by "triangle" 
        for m   = 1 : N
            div         = M - m - floor(m/eps);          %   integer - edge is divided into div (XX-jump-XX-jump)
            scale       = div/M;                         %   real - relative
            alpha       = (1 + 2*scale)/3;
            beta        = (1 - scale)/3;    
            coeff1      = [alpha beta beta]';    %   p1 new 
            coeff2      = [beta alpha beta]';    %   p2 new
            coeff3      = [beta beta alpha]';    %   p3 new
            %   first edge
            for n = 1 : div
                vector      = coeff1*(div - n + 1)/div + coeff2*(n - 1)/div;
                coeff(:,k)  = vector;
                k = k + 1;
            end
            %   second edge
            for n = 1 : div
                vector      = coeff2*(div - n + 1)/div + coeff3*(n - 1)/div;
                coeff(:,k)  = vector;
                k = k + 1;
            end
            %   third edge
            for n = 1 : div
                vector      = coeff3*(div - n + 1)/div + coeff1*(n - 1)/div;
                coeff(:,k)  = vector;
                k = k + 1;
            end  
        end

        %   Center point 
        if (3*floor(M/3) ~= M)
            coeff(:,k) = [1/3 1/3 1/3]';
        end	
        weights = 1/size(coeff,2)*ones(1,size(coeff,2)); 
    end

    if nargs == 2
    %   Gaussian quadrature formulae
    %   arg - number of integration points 
        if arg1 == 1    % first order (center)               
            coeff(:,1)  = [1/3 1/3 1/3]';
            weights(1)  =   1;
        end

        if (arg1 == 3) & (arg2 == 2)    %   second order (sides)                
            coeff(:,1)  = [1/2 1/2 0]';
            coeff(:,2)  = [0 1/2 1/2]';
            coeff(:,3)  = [1/2 0 1/2]';
            weights(1)  =   1/3;
            weights(2)  =   1/3;
            weights(3)  =   1/3;
        end
        if (arg1 == 4) & (arg2 == 3)    % third order                
            a1          =   0.6;
            b1          =   0.2;       
            coeff(:,1)  = [1/3 1/3 1/3]';
            coeff(:,2)  = [a1  b1  b1]';
            coeff(:,3)  = [b1  a1  b1]';
            coeff(:,4)  = [b1  b1  a1]';       
            weights(1)  =   -27/48;
            weights(2)  =    25/48;
            weights(3)  =    25/48;
            weights(4)  =    25/48;    	
        end
        if (arg1 == 6) & (arg2 == 3)   % third order (sides)     
            coeff(:,1)  = [1/2 1/2 0]';
            coeff(:,2)  = [0 1/2 1/2]';
            coeff(:,3)  = [1/2 0 1/2]';
            weights(1)  =   0.016666666;
            weights(2)  =   0.016666666;
            weights(3)  =   0.016666666;        
            a1          =   0.666666666;
            b1          =   0.166666666;       
            coeff(:,4)  = [a1  b1  b1]';
            coeff(:,5)  = [b1  a1  b1]';
            coeff(:,6)  = [b1  b1  a1]';       
            weights(4)  =    0.15;
            weights(5)  =    0.15;
            weights(6)  =    0.15;    	
            weights     = weights*2;
        end    

        if (arg1 == 7) & (arg2 == 5)    % fifth order                
            a1          =   0.797426985353087;
            b1          =   0.101286507323456;       
            a2          =   0.059715871789770;
            b2          =   0.470142064105115;	
            coeff(:,1)  = [1/3 1/3 1/3]';
            coeff(:,2)  = [a1  b1  b1]';
            coeff(:,3)  = [b1  a1  b1]';
            coeff(:,4)  = [b1  b1  a1]';       
            coeff(:,5)  = [a2  b2  b2]';
            coeff(:,6)  = [b2  a2  b2]';
            coeff(:,7)  = [b2  b2  a2]';    
            weights(1)    =   0.2250000;
            weights(2)    =   0.1259392;
            weights(3)    =   0.1259392;
            weights(4)    =   0.1259392;    
            weights(5)    =   0.1323942;
            weights(6)    =   0.1323942;
            weights(7)    =   0.1323942;        
       end
       if (arg1 == 9) & (arg2 == 5)    % fifth order (sides)
            coeff(:,1)  = [1 0 0]';
            coeff(:,2)  = [0 1 0]';
            coeff(:,3)  = [0 0 1]';
            coeff(:,4)  = [1/2 1/2 0]';
            coeff(:,5)  = [0 1/2 1/2]';
            coeff(:,6)  = [1/2 0 1/2]';
            a1          =   0.62283903060711;
            b1          =   0.18858048469644;       
            coeff(:,7)  = [a1  b1  b1]';
            coeff(:,8)  = [b1  a1  b1]';
            coeff(:,9)  = [b1  b1  a1]';
            weights(1)  =   0.01027006767296;
            weights(2)  =   0.01027006767296;
            weights(3)  =   0.01027006767296;        
            weights(4)  =   0.03098774943413;
            weights(5)  =   0.03098774943413;
            weights(6)  =   0.03098774943413;                          
            weights(7)  =   0.12540884955956;  
            weights(8)  =   0.12540884955956;
            weights(9)  =   0.12540884955956; 
            weights = weights*2;
       end
       if (arg1 == 13) & (arg2 == 7)    % seventh order
           a1 = 0.4793080678;
           b1 = 0.2603459660;
           a2 = 0.8697397941;
           b2 = 0.0651301029;
           a3 = 0.6384441885;			
           b3 = 0.3128654960;
           c3 = 0.0486903154;
           coeff(:,1)  = [1/3 1/3 1/3]';
           coeff(:,2)  = [a1  b1  b1]';
           coeff(:,3)  = [b1  a1  b1]';
           coeff(:,4)  = [b1  b1  a1]';       
           coeff(:,5)  = [a2  b2  b2]';
           coeff(:,6)  = [b2  a2  b2]';
           coeff(:,7)  = [b2  b2  a2]';    			       
           coeff(:,8) = [a3  b3  c3]'; 
           coeff(:,9) = [a3  c3  b3]'; 
           coeff(:,10) = [b3  a3  c3]'; 
           coeff(:,11) = [b3  c3  a3]'; 
           coeff(:,12) = [c3  a3  b3]'; 
           coeff(:,13) = [c3  b3  a3]'; 
           weights(1)  =    -0.14957004; 
           weights(2)  =	 0.1756152574; 
           weights(3)  =	 0.1756152574; 
           weights(4)  =	 0.1756152574; 
           weights(5)  =	 0.0533472356; 
           weights(6)  =	 0.0533472356; 
           weights(7)  =	 0.0533472356; 
           weights(8)  =	 0.0771137608; 
           weights(9)  =	 0.0771137608; 
           weights(10) =	 0.0771137608; 
           weights(11) =	 0.0771137608; 
           weights(12) =	 0.0771137608; 
           weights(13) =	 0.0771137608;       
       end
       if (arg1 == 25) & (arg2 == 10)    % tenth order       
            a1   =   0.1498275788;
            b1   =   0.4250862106;		
            a2   =   0.9533822650;
            b2   =   0.0233088675;		
            a3   =   0.6283074002; 
            b3   =   0.2237669736; 
            c3   =   0.1479256262; 		
            a4   =   0.6113138262; 
            b4   =   0.3587401419; 
            c4   =   0.0299460319;		
            a5   =   0.8210720699; 
            b5   =   0.1432953704; 
            c5   =   0.0356325597;    
            coeff(:,1)  = [1/3 1/3 1/3]';
            coeff(:,2)  = [a1  b1  b1]';
            coeff(:,3)  = [b1  a1  b1]';
            coeff(:,4)  = [b1  b1  a1]';       
            coeff(:,5)  = [a2  b2  b2]';
            coeff(:,6)  = [b2  a2  b2]';
            coeff(:,7)  = [b2  b2  a2]';    			       		
            coeff(:,8)  = [a3  b3  c3]'; 
            coeff(:,9)  = [a3  c3  b3]'; 
            coeff(:,10) = [b3  a3  c3]'; 
            coeff(:,11) = [b3  c3  a3]'; 
            coeff(:,12) = [c3  a3  b3]'; 
            coeff(:,13) = [c3  b3  a3]';        		
            coeff(:,14) = [a4  b4  c4]'; 
            coeff(:,15) = [a4  c4  b4]'; 
            coeff(:,16) = [b4  a4  c4]'; 
            coeff(:,17) = [b4  c4  a4]'; 
            coeff(:,18) = [c4  a4  b4]'; 
            coeff(:,19) = [c4  b4  a4]';        		
            coeff(:,20) = [a5  b5  c5]'; 
            coeff(:,21) = [a5  c5  b5]'; 
            coeff(:,22) = [b5  a5  c5]'; 
            coeff(:,23) = [b5  c5  a5]'; 
            coeff(:,24) = [c5  a5  b5]'; 
            coeff(:,25) = [c5  b5  a5]';        		
            weights(1) =     0.03994725237;
            weights(2) =     0.03556190112;
            weights(3) =     0.03556190112;
            weights(4) =     0.03556190112;        
            weights(5) =     0.00411190935;
            weights(6) =     0.00411190935;            
            weights(7) =     0.00411190935;                 
            weights(8) =     0.02271529614;
            weights(9) =     0.02271529614;
            weights(10) =    0.02271529614;                
            weights(11) =    0.02271529614;
            weights(12) =    0.02271529614;
            weights(13) =    0.02271529614;        
            weights(14) =    0.01867992812;
            weights(15) =    0.01867992812;
            weights(16) =    0.01867992812;                
            weights(17) =    0.01867992812;
            weights(18) =    0.01867992812;
            weights(19) =    0.01867992812;        
            weights(20) =    0.01544332844;
            weights(21) =    0.01544332844;
            weights(22) =    0.01544332844;
            weights(23) =    0.01544332844;
            weights(24) =    0.01544332844;
            weights(25) =    0.01544332844;
            weights     = weights*2;
        end
    end

    IndexF = size(coeff,2);
end