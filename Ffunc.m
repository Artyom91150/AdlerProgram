function F = Ffunc(Phi, D)

if (D(1) >= 0)
    if (mod(Phi, 2*pi) >= D(1)) && (mod(Phi, 2*pi) <= D(2))
        F = 0;
    else
        F = 1;
    end
else if (D(1) >= -pi/2)
        if (mod(Phi, 2*pi) >= 0) && (mod(Phi, 2*pi) <= D(2)) || (mod(Phi, 2*pi) >= 2*pi + D(1)) && (mod(Phi, 2*pi) <= 2*pi)
            F = 0;
        else
            F = 1;
        end
    else
       F = 0; 
    end
end

% %Отклонение
% P = 0.01;
% 
% if (D(1) >= 0)
%     if (mod(Phi, 2*pi) >= D(1) + P) && (mod(Phi, 2*pi) <= D(2) - P)
%         F = 0;
%     else if (abs(mod(Phi, 2*pi) - D(1)) <= P) || (abs(mod(Phi, 2*pi) - D(2)) <= P)
%         F = 0.5;
%         else
%             F = 1;
%         end
%     end
% else if (D(1) >= -pi/2)
%         if (mod(Phi, 2*pi) >= 0) && (mod(Phi, 2*pi) <= D(2)) || (mod(Phi, 2*pi) >= 2*pi + D(1)) && (mod(Phi, 2*pi) <= 2*pi)
%             F = 0;
%         else
%             F = 1;
%         end
%     else
%        F = 0; 
%     end
% end


    
end

