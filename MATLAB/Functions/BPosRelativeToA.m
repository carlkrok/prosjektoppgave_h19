function rLVLH_Rel = BPosRelativeToA( rECI_A, vECI_A, rECI_B )

   hECI_A = cross( rECI_A, vECI_A );
   hECINorm_A = norm( hECI_A );

   i_unitVectorMovingFrame = ( rECI_A / norm( rECI_A ));
   k_unitVectorMovingFrame = ( hECI_A / hECINorm_A );
   j_unitVectorMovingFrame = cross( k_unitVectorMovingFrame, i_unitVectorMovingFrame );


   QmatECItoLVLH = [i_unitVectorMovingFrame';
                    j_unitVectorMovingFrame';
                    k_unitVectorMovingFrame'];

   rECI_Rel = rECI_B - rECI_A;

   rLVLH_Rel = QmatECItoLVLH * rECI_Rel;

end
