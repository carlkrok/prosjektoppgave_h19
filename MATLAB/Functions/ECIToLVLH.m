function [ QmatECItoLVLH ] = ECIToLVLH( rECI, vECI )

   hECI = cross( rECI, vECI );
   hECINorm = norm( hECI );

   i_unitVectorMovingFrame = ( rECI / norm( rECI ));
   k_unitVectorMovingFrame = ( hECI / hECINorm );
   j_unitVectorMovingFrame = cross( k_unitVectorMovingFrame, i_unitVectorMovingFrame );


   QmatECItoLVLH = [i_unitVectorMovingFrame';
                    j_unitVectorMovingFrame';
                    k_unitVectorMovingFrame'];


end