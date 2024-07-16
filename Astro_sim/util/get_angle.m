function  [angle] = get_angle(fdA,fdB,constrain)
%whether constrain the angle to [-pi/2,pi/2]
      angle= 1/2*(atan(real(fdA)./real(fdB))+atan(imag(fdA)./imag(fdB)));
      if constrain
          angle(angle>pi/2)= angle(angle>pi/2)-pi;
          angle(angle<-pi/2)= angle(angle<-pi/2)+pi;
      end
end