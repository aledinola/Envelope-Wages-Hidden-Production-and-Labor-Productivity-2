function [xopt,fopt,nacc,nfcnev,nobds,ier,t,vm]= SIMULANS_may(fcn,x,option,lb,ub,c,vm,param)
% [xopt, fopt, nacc, nfcnev, nobds, ier, t, vm] = simulan(fcn,x,option,lb,ub,c,vm,param);
%
% Input Parameters:
%   Note: The suggested values generally come from Corana et al. To
%         drastically reduce runtime, see Goffe et al., pp. 90-1 for
%         suggestions on choosing the appropriate RT and NT.
%   fcn - Name of a function which returns the value of function
%         to be optimized.
%   x - The starting values for the variables of the function to be
%       optimized. (Nx1)
%
%   option:
%
%	 option(1):
%   max - Denotes whether the function should be maximized or
%         minimized. A value =1 denotes maximization while a
%         value =0 denotes minimization. Intermediate output (see IPRINT)
%         takes this into account.
%
%   option(2):
%   rt - The temperature reduction factor. The value suggested by
%        Corana et al. is .85. See Goffe et al. for more advice.
%
%   option(3):
%   eps - Error tolerance for termination. If the final function
%         values from the last neps temperatures differ from the
%         corresponding value at the current temperature by less than
%         EPS and the final function value at the current temperature
%         differs from the current optimal function value by less than
%         EPS, execution terminates and IER = 0 is returned.
%
%   option(4):
%   ns - Number of cycles. After NS*N function evaluations, each
%        element of VM is adjusted so that approximately half of
%        all function evaluations are accepted. The suggested value
%        is 20.
%
%   option(5):
%   nt - Number of iterations before temperature reduction. After
%        NT*NS*N function evaluations, temperature (T) is changed
%        by the factor RT. Value suggested by Corana et al. is
%        max(100, 5*n). See Goffe et al. for further advice.
%
%   option(6):
%   neps - Number of final function values used to decide upon termi-
%          nation. See EPS. Suggested value is 4.
%
%   option(7):
%   maxevl - The maximum number of function evaluations. If it is
%            exceeded, IER = 1.
%
%   option(8):
%   iprint - controls printing inside SA.
%            Values: 0 - Nothing printed.
%                    1 - Function value for the starting value and
%                        summary results before each temperature
%                        reduction. This includes the optimal
%                        function value found so far, the total
%                        number of moves (broken up into uphill,
%                        downhill, accepted and rejected), the
%                        number of out of bounds trials, the
%                        number of new optima found at this
%                        temperature, the current optimal X and
%                        the step length VM. Note that there are
%                        N*NS*NT function evalutations before each
%                        temperature reduction. Finally, notice is
%                        is also given upon achieveing the termination
%                        criteria.
%                    2 - Each new step length (VM), the current optimal
%                        X (XOPT) and the current trial X (X). This
%                        gives the user some idea about how far X
%                        strays from XOPT as well as how VM is adapting
%                       to the function.
%                    3 - Each function evaluation, its acceptance or
%                        rejection and new optima. For many problems,
%                        this option will likely require a small tree
%                        if hard copy is used. This option is best
%                        used to learn about the algorithm. A small
%                        value for MAXEVL is thus recommended when
%                        using IPRINT = 3.
%            Suggested value: 1
%            Note: For a given value of IPRINT, the lower valued
%                  options (other than 0) are utilized.
%
%   option(9):
%   t - the initial temperature. See Goffe et al. for advice.
%
%   lb - The lower bound for the allowable solution variables. (N)
%   ub - The upper bound for the allowable solution variables. (N)
%        If the algorithm chooses X(I) < LB(I) or X(I) > UB(I),
%        I = 1, N, a point is from inside is randomly selected. This
%        This focuses the algorithm on the region inside UB and LB.
%        Unless the user wishes to concentrate the search to a par-
%        ticular region, UB and LB should be set to very large positive
%        and negative values, respectively. Note that the starting
%        vector X should be inside this region. Also note that LB and
%        UB are fixed in position, while VM is centered on the last
%        accepted trial set of variables that optimizes the function.
%   c - Vector that controls the step length adjustment. The suggested
%       value for all elements is 2.0. (N)
% Input/Output Parameters:
%   t - On input, the initial temperature. See Goffe et al. for advice.
%       On output, the final temperature.
%   vm - The step length vector. On input it should encompass the
%        region of interest given the starting value X. For point
%        X(I), the next trial point is selected is from X(I) - VM(I)
%        to  X(I) + VM(I). Since VM is adjusted so that about half
%        of all points are accepted, the input value is not very
%        important (i.e. is the value is off, SA adjusts VM to the
%        correct value). (N)
%
% Output Parameters:
%   xopt - The variables that optimize the function. (N)
%   fopt - The optimal value of the function.
%   nacc - The number of accepted function evaluations.
%   nfcnev - The total number of function evaluations. In a minor
%            point, note that the first evaluation is not used in the
%            core of the algorithm; it simply initializes the
%            algorithm.
%   nobds - The total number of trial function evaluations that
%           would have been out of bounds of LB and UB. Note that
%           a trial point is randomly selected between LB and UB.
%
%   ier - The error return number.
%         Values: 0 - Normal return; termination criteria achieved.
%                 1 - Number of function evaluations (NFCNEV) is
%                     greater than the maximum number (MAXEVL).
%                 2 - The starting value (X) is not inside the
%                     bounds (LB and UB).
%                 3 - The initial temperature is not positive.
%                 99 - Should not be seen; only used internally.
%
% © E.G. Tsionas
% Department of Economics
% University of Toronto
% 150 St George Str.
% Toronto Ontario M5S 1A1
% Canada
% (416) 978 4188 (office)
% (416) 978 6713 (fax)
% etsionas@epas.utoronto.ca (Internet)
%
% Adapted for Matlab © F. Collard (1998)
%

[n,dum]=size(x);
xopt=zeros(n,1);
xp=zeros(n,1);
%
% from options
%
if isempty(option);
  maxim=0;         % 1 if maximize, 0 otherwise
  rt=0.85;         % temperature decreasing ratio
  eps=1e-6;        % error tolerance for termination
  ns=20;           % number of cycles
  nt=max(100,5*n); % number of iterations before temperature reduction
  neps=4;          % number of final functions values used to decide upon termination
  maxevl=1000000;  % maximum number of function evaluations
  iprint=1;        % control inside printing
  t=1;             % initial temperature
else;
  if (option(1)~=1);
    maxim=0;
  else
    maxim=option(1);
  end;

  if (option(2)<=0);
    rt=0.85;
  else
    rt=option(2);
  end;
  if (option(3)<=0);
    eps=1e-6;
  else;
    eps=option(3);
  end;
  if (option(4)<=0);
    ns=20;
  else;
    ns=option(4);
  end;
  if (option(5)<=0);
    nt=max(100,5*n);
  else
 	 nt=option(5);
  end;
  if (option(6)<=0);
    neps=4;
  else;
    neps=option(6);
  end;
  if (option(7)<=0);
    maxevl=100000;
  else;
 	 maxevl=option(7);
  end;
  if (option(8)<=0);
    iprint=1;
  else;
    iprint=option(8);
  end;
  if (option(9)<=0);
   disp('The initial temperature is not positive. Reset the variable t');
   return;
  else;
    t=option(9);
  end;
end;
%
% Initial values
%
nacc=0;
nobds=0;
nfcnev=0;
ier=99;
xopt=x;
nacp=zeros(n,1);
fstar=1e20*ones(neps,1);
quit1=0;
%
%If the initial temperature is not positive, notify the user and abort.
%
%
% If the initial value is out of bounds, notify the user and abort.
%
if (sum(x>ub)+ sum(x<lb))>0;
  disp('initial condition out of bounds');
  ier = 2;
  quit1=1;
  return;
end;
%
% Evaluate the function with input x and return value as f.
%
if nargin>7
  f=feval(fcn, x, param); 
  %eval(['f=',fcn,'(x,param);']);
else
  f=feval(fcn,x); 
  %eval(['f=',fcn,'(x);']);
end;
%
% If the function is to be minimized, switch the sign of the function.
% Note that all intermediate and final output switches the sign back
% to eliminate any possible confusion for the user.
%
if (maxim == 0);
  f = -f;
end;
nfcnev = nfcnev + 1;
fopt = f;
fstar(1) = f;
if (iprint > 1);
  disp(' ');
  disp('initial x');disp(x(:)')
  if maxim;
    disp(sprintf('initial f : %g',f));
  else;
    disp(sprintf('initial f : %g',-f));
  end;
end;
%
% Start the main loop. Note that it terminates if:
%       (i) the algorithm succesfully optimizes the function or
%		 (ii) there are too many function evaluations (more than maxevl).
%
while ~quit1;
  nup = 0;
  nrej = 0;
  nnew = 0;
  ndown = 0;
  lnobds = 0;

  for m=1:nt;
    for j=1:ns
      for h=1:n;
%
% Generate xp, the trial value of x. Note use of vm to choose xp.
%
        for i=1:n;
          if i==h;
            xp(i) = x(i) + (rand(1,1)*2.- 1.) * vm(i);
          else;
            xp(i) = x(i);
          end;
%
% If xp is out of bounds, select a point in bounds for the trial.
%
          if ((xp(i)<lb(i))|(xp(i)>ub(i)));
            xp(i)=lb(i)+(ub(i)-lb(i))*rand(1,1);
            lnobds=lnobds+1;
            nobds=nobds+1;
            if (iprint>=3);
              disp(' ');
              disp('current x');disp(x(:)');
              if (maxim);
                disp(sprintf('current f : %g',f));
              else;
                disp(sprintf('current f : %g',-f));
              end;
              disp('trial x');disp(xp(:)');
              disp('point rejected since out of bounds');
            end;
          end;
        end; % end of loop on i
%
% Evaluate the function with the trial point xp and return as fp.
%
        if nargin>7
          %eval(['fp=',fcn,'(xp,param);']);
          fp=feval(fcn, xp, param); 
        else
          %eval(['fp=',fcn,'(xp);']);
          fp=feval(fcn, xp); 
        end;
        if(maxim==0);
          fp=-fp;
        end;
        nfcnev=nfcnev+1;
        if (iprint >= 3);
          disp(' ');
          disp('current x');disp(x(:)');
          if maxim;
            disp(sprintf('current f : %g',f));
            disp('trial x :');disp(xp(:)');
            disp(sprintf('resulting f : %g',fp));
          else;
            disp(sprintf('current f : %g',-f));
            disp('trial x :');disp(xp(:)');
            disp(sprintf('resulting f : %g',-fp));
          end;
        end;
%
% If too many function evaluations occur, terminate the algorithm.
%
        if(nfcnev >= maxevl);
          ier = 1;
          quit1=1;
          disp('Too many function evaluations; consider');
          disp('increasing maxevl or eps, or decreasing');
          disp('nt or rt. These results are likely to be poor');
          if (maxim == 0); fopt = -fopt; end;
          return;
        end;
%
% Accept the new point if the function value increases.
%
        if (fp>=f);
          if (iprint>=3);
            disp('point accepted');
          end;
          x = xp;
          f = fp;
          nacc = nacc + 1;
          nacp(h) = nacp(h) + 1;
          nup = nup + 1;
%
% If greater than any other point, record as new optimum.
%
          if (fp > fopt);
            if(iprint >= 3);
              disp('new optimum');
            end;
            xopt = xp;
            fopt = fp;
            nnew = nnew + 1;
          end;
%
% If the point is lower, use the Metropolis criteria to decide on
% acceptance or rejection.
%
        else;
          if (((fp - f)/t) > 709); p = 8.2184e+307;
          elseif (((fp - f)/t) < -708); p = 0;
          else; p = exp(((fp - f)/t));
          end;
          pp = rand(1,1);

          if (pp < p);
            if (iprint >= 3);
              if (maxim);
                disp('Though lower, point accepted');
              else;
                disp('though higher, point accepted');
              end;
            end;
            x = xp;
            f = fp;
            nacc = nacc + 1;
            nacp(h) = nacp(h) + 1;
            ndown = ndown + 1;
          else;
            nrej = nrej + 1;
            if (iprint >= 3);
              if (maxim);
                disp('lower point rejected');
              else;
                disp('higher point rejected');
              end;
            end;
          end;
        end;
      end; % end of loop on h
    end;   % end of loop on j
%
% Adjust vm so that approximately half of all evaluations are accepted.
%
        for i=1:n;
          ratio = nacp(i)/ns;
          if (ratio > .6);
            vm(i)=vm(i)*(1+c(i)*(ratio - .6)/.4);
          elseif (ratio < .4);
            vm(i)=vm(i)/(1+c(i)*((.4 - ratio)/.4));
          end;
          if (vm(i) > (ub(i)-lb(i)));
            vm(i)=ub(i)-lb(i);
          end;
        end;

        if(iprint >= 2);
          disp('intermediate results after step length adjustment');
          disp('new step length (vm)');disp(vm(:)');
          disp('current optimal x');disp(xopt(:)');
          disp('current x');disp(x(:)');
          disp(' ');
        end;

      nacp = zeros(n,1);
  end; % end of loop on m

  if(iprint >= 1);
    diary simdiary.txt;
    disp('Intermediate results before next temperature reduction');
    disp(sprintf('current temperature         : % g',t));
    if (maxim);
      disp(sprintf('max function value so far   : %g',fopt));
      disp(sprintf('total moves                 : %g',nup+ndown+nrej));
      disp(sprintf('uphill                      : %g',nup));
      disp(sprintf('accepted downhill           : %g',ndown));
      disp(sprintf('rejected downhill           : %g',nrej));
      disp(sprintf('out of bounds trials        : %g',lnobds));
      disp(sprintf('new maximima this temperature : %g',nnew));
    else;
      disp(sprintf('min function value so far   : %g',-fopt));
      disp(sprintf('total moves                 : %g',nup+ndown+nrej));
      disp(sprintf('downhill                    : %g',nup));
      disp(sprintf('accepted uphill             : %g',ndown));
      disp(sprintf('rejected uphill             : %g',nrej));
      disp(sprintf('out of bounds trials        : %g',lnobds));
      disp(sprintf('new minima this temperature : %g',nnew));
    end;
    disp('current optimal x');disp(xopt(:)');
    disp('strength (vm)');disp(vm(:)');
    diary off;
  end;
%
% Check termination criteria.
%
  quit1=0;
  fstar(1)=f;
  if ((fopt-fstar(1))<= eps);quit1=1;end;
  if ( sum(abs(f-fstar)>eps)>0);quit1=0;end;
  if ier~=99;
    quit1=1;
  end;
%
% If termination criteria are not met, prepare for another loop.
%
  t = rt*t;
  for i=neps:-1:2;
    fstar(i) = fstar(i-1);
  end
  f = fopt;
  x = xopt;
end;
if ier==1;
  disp('The initial temperature is not positive. Reset the variable t');
end;
x = xopt;
ier = 0;
if (maxim == 0); fopt = -fopt; end;
if(iprint >= 1);
  disp('SA achieved termination criteria. ier = 0');
end;
