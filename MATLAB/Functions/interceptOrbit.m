function [ deltaVStart, deltaVEnd, vIntersectOrbit ] = interceptOrbit( posStart, vStart, posEnd, vEnd, deltaTime, orbitTypeDebris, mu, tolerance, nMax )

   [ vRequiredStart, vIntersectOrbit ] = Lambert( posStart, posEnd, deltaTime, orbitTypeDebris, mu, tolerance, nMax );
   
   % Change in velocity required
   deltaVStart = vRequiredStart - vStart;
   deltaVEnd = vEnd - vIntersectOrbit;

end