% MEAN2OSC converts a set of mean Keplerian orbital elements to osculating
% Keplerian orbital elements perturbed by Earth's oblateness.
% 
%   Inputs:
%     mean_elem - vector of mean Keplerian orbital elements:
%                 a - semi-major axis [m]
%                 e - eccentricity [-]
%                 i - inclination [rad]
%                 O - right ascension of the ascending node [rad]
%                 w - argument of periapsis [rad]
%                 M - mean anomaly [rad]
%     J2_flag   - flag indicating whether or not to consider J2
%                 1 -> J2 is enabled, calculate according to algorithm
%                 0 -> J2 is disabled, osculating elements = mean elements
%                 DEFAULT J2_flag = 1
% 
%   OUTPUTS:
%     osc_elem - vector of osculating Keplerian orbital elements
%                a - semi-major axis [m]
%                e - eccentricity [-]
%                i - inclination [rad]
%                O - right ascension of the ascending node [rad]
%                w - argument of periapsis [rad]
%                M - mean anomaly [rad]

function osc_elem = mean2osc(mean_elem, J2_flag)

    % Check inputs
    if (nargin < 2) || isempty(J2_flag)
        J2_flag = 1;
    end
    if (nargin < 1) || isempty(mean_elem)
        error('Must input mean elements set');
    end
    
    % Format input to column vector and set tolerance
    mean_elem = mean_elem(:);

    % With J2, run method
    if J2_flag == 1
        
        % Convert to mean equinoctial elements
        mean_equi_elem = koe2equioe(mean_elem);

        % Convert to osculating equinoctial elements
        [~, osc_equi_elem] = mean_osc_closed_equi(mean_equi_elem,J2_flag);

        % Convert to keplerian elements
        osc_elem = equioe2koe(osc_equi_elem)';
        
        % Format output
        osc_elem = osc_elem(:);
        
    % Without J2, elements are equal
    else
        osc_elem = mean_elem;
    end
    
end
