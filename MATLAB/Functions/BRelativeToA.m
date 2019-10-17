
function [rLVLH_Rel, vLVLH_Rel, aLVLH_Rel] = BRelativeToA( rECI_A, vECI_A, rECI_B, vECI_B )

   hECI_A = cross( rECI_A, vECI_A );
   hECINorm_A = norm( hECI_A );

   i_unitVectorMovingFrame = ( rECI_A / norm( rECI_A ));
   k_unitVectorMovingFrame = ( hECI_A / hECINorm_A );
   j_unitVectorMovingFrame = cross( k_unitVectorMovingFrame, i_unitVectorMovingFrame );


   QmatECItoLVLH = [i_unitVectorMovingFrame';
                    j_unitVectorMovingFrame';
                    k_unitVectorMovingFrame'];

   angularVelocityLVLH = hECI_A / (norm( rECI_A ))^2;
   angularAccelerationLVLH = -2 * (dot( vECI_A, rECI_A ) / (norm( rECI_A ))^2) * angularVelocityLVLH;

   aECI_A = -(muEarth / norm(rECI_A)^3) * rECI_A;
   aECI_B = -(muEarth / norm(rECI_B)^3) * rECI_B;


   rECI_Rel = rECI_B - rECI_A;
   vECI_Rel = vECI_B - vECI_A - cross( angularVelocityLVLH, rECI_Rel );
   aECI_Rel = aECI_B - aECI_A - cross( angularAccelerationLVLH, rECI_Rel ) - cross( angularVelocityLVLH, cross( angularVelocityLVLH, rECI_Rel )) - 2 * cross( angularVelocityLVLH, vECI_Rel );


   rLVLH_Rel = QmatECItoLVLH * rECI_Rel;
   vLVLH_Rel = QmatECItoLVLH * vECI_Rel;
   aLVLH_Rel = QmatECItoLVLH * aECI_Rel;

end