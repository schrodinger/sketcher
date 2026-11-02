// CRDGEN-428 uses the replacement certified on 2026-10-20, without Orn.
export const bicyclicPeptide =
  'CHEM1{[CCC(N1CN(C(CC)=O)CN(C(CC*)=O)C1)=O |$_R1;;;;;;;;;;_R2;;;;;;;_R3;;;;$|]}|PEPTIDE1{A.R.D.C.P.L.V.N.P.L.C.L.H.P.G.W.T.C.[dV].C}$PEPTIDE1,CHEM1,4:R3-1:R1|PEPTIDE1,CHEM1,11:R3-1:R2|CHEM1,PEPTIDE1,1:R3-18:R3$$$V2.0';

// Derived from CRDGEN-434's three certification examples, preserving the chain
// lengths, backbone/side-chain connections, and capped ends. Long ChEMBL
// symbols replace entries outside the existing database fixture; rendering
// tests also require the resulting beads to be wider than standard ones.
export const largeLabelPeptides = [
  'PEPTIDE1{[EtCO_dPhg].[Me_Cha].[Me_Lys(Boc)].[dD]}|PEPTIDE2{[EtCO_dPhg].T.[dS].[Me_Phg].C.[dPhe(2-NO2)].[dPhe(3-NO2)].[meR].[dPhe(2-Cl)].[Aib_OMe]}|PEPTIDE3{[dPhe(2-Me)].[Me_Hph].[dLys(Me)].[Phe(4-Me)].G.[Asp(OtBu)]}$PEPTIDE1,PEPTIDE3,4:R3-1:R1|PEPTIDE2,PEPTIDE1,5:R3-4:R2$$$V2.0',
  'PEPTIDE1{[Me_Nva].[Lys(Boc)].C.[Me_Spg].[Tyr(tBu)].T.[Glu(OtBu)].[Me_Phg].[meM].[Wil].[Me_1-Nal].[meR].[dPhe(2-Cl)].[dPhe(3-Cl)]}|PEPTIDE2{[dLys(Me)].[dTyr(Bzl)].[meL]}$PEPTIDE1,PEPTIDE2,3:R3-1:R1$$$V2.0',
  'PEPTIDE1{A.[Phe(4-Cl)].[Me_Abu].[Tyr(tBu)].[dThr(Bzl)].E.[Me_Lys(Me)2].[Lys(Cbz)].[dSer(tBu)].P.[Me_lalloI].[Me_3-Pal].[dTyr(Bzl)].D.[Phe(2-I)].[Me_Spg].[dLys(2-MP)].[Me_dalloT].[Lys(Me)2].[dF]}|PEPTIDE2{A.[meA].[dLys(Me)].[Cya].[dPhe(4-Me)].[Phe(4-NH2)]}|PEPTIDE3{[Me_lalloT].[dH].[Lys(BrAc)2]}$PEPTIDE1,PEPTIDE2,14:R3-1:R1|PEPTIDE1,PEPTIDE3,6:R3-1:R1$$$V2.0',
];
