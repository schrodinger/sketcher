#pragma once

#include <string>
#include <vector>

static const std::vector<std::string> INVALID_EXAMPLES{
    /* un numbered connection points 2 */
    R"(RNA1{R(C)[OP([*])([*])=O]}$$$$V2.0)",
    /* repeating chem monomer not valid */ R"(CHEM1{[Az]'3'}$$$$V2.0)",
    /* dot before and after the branch (invalid notation) */
    R"(RNA1{P.R.(U).P.R(T)}$$$$)",
    /* rna to rna strand both ends incorrectly fully specified */
    R"(RNA1{*}|RNA2{R(A)P.R(C)P}$RNA1,RNA2,1:R1-1:R1$$$V2.0)",
    /* dot before and after the branch invalid notation */
    R"(RNA1{P.R.(U).P.R(T)}$$$$)",
    /* unknown peptide peptide unknown peptide side incorrectly specified */
    R"(PEPTIDE1{*}|PEPTIDE2{A.C}$PEPTIDE1,PEPTIDE2,1:R1-?:?$$$V2.0)",
    /* rna wildcard to rna strand rna wildcard end fully specified */
    R"(RNA1{*}|RNA2{R(A)P.R(C)P}$RNA1,RNA2,1:R1-?:?$$$V2.0)",
    /* dot before the branch */
    R"(RNA1{R.([NC1=NC(N([*:1])C=C1)=O])P}$$$$)",
    /* repeating chem monomer   not valid */ R"(CHEM1{[Az]'3'}$$$$V2.0)",
    /* NOTE: UNSUPPORTED: repeating partial rna monomer not valid */ // R"(RNA1{R'3'}$$$$V2.0)",
    /* blob peptide blob specified peptide specified */
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,1:R1-1:R1$$$V2.0)",
    /* dot before and after the inline backbone before branch */
    R"(RNA1{P.[O[C@H]1[C@H]([*:1])O[C@H](CO[*:2])[C@H]1O[*:3]].(U)P.R(T)}$$$$V2.0)",
    /* use of wildcard multiple chems */
    R"(CHEM1{*}|CHEM2{[SMCC]}$CHEM2,CHEM1,1:R1-1:R1$$$V2.0)",
    /* blob peptide   blob specified and peptide unspecified */
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,1:R1-?:?$$$V2.0)",
    /* unknown peptide peptide peptide side fully specified */
    R"(PEPTIDE1{*}|PEPTIDE2{A.C}$PEPTIDE1,PEPTIDE2,1:?-?:R1$$$V2.0)",
    /* no dots */
    R"(RNA1{P[O[C@H]1[C@H]([*:1])O[C@H](CO[*:2])[C@H]1O[*:3]](U)P.R(T)}$$$$V2.0)",
    /* dot before and after the inline backbone and before branch */
    R"(RNA1{P.[O[C@H]1[C@H]([*:1])O[C@H](CO[*:2])[C@H]1O[*:3]].(U)P.R(T)}$$$$V2.0)",
    /* dot after the in line backbone and before the branch */
    R"(RNA1{P[O[C@H]1[C@H]([*:1])O[C@H](CO[*:2])[C@H]1O[*:3]].(U)P.R(T)}$$$$V2.0)",
    /* blob peptide   blob r group specified and peptide unspecified */
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,?:R1-?:?$$$V2.0)",
    /* blob peptide blob unspecified peptide r group specified */
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,?:?-?:R1$$$V2.0)",
    /* un numbered connection points */
    R"(RNA1{[O[C@H]1[C@H]([*])O[C@H](CO[*])[C@H]1O[*]](C)P}$$$$V2.0)",
    /* un numbered connection points */
    R"(RNA1{R(C)[OP([*])([*])=O]}$$$$V2.0)",
    /* dot before the branch invalid notation */
    R"(RNA1{P.R.(U)P.R(T)}$$$$)",
    /* dot after the in line backbone and before the branch */
    R"(RNA1{P[O[C@H]1[C@H]([*:1])O[C@H](CO[*:2])[C@H]1O[*:3]].(U)P.R(T)}$$$$V2.0)",
    /* blob peptide   blob unspecified and peptide r group specified */
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,?:?-?:R1$$$V2.0)",
    /* dot before the branch (invalid notation) */
    R"(RNA1{P.R.(U)P.R(T)}$$$$)",
    /* blob peptide   blob location specified and peptide unspecified */
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,1:?-?:?$$$V2.0)",
    /* unnumbered connection point no brackets 5 */
    R"(PEPTIDE1{A.[*C(=O)[C@H](C)N(*)C].A}$$$$)",
    /* dot before and after the branch */
    R"(RNA1{R.([NC1=NC(N([*:1])C=C1)=O]).P}$$$$)",
    /* blob peptide   blob specified and peptide specified */
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,1:R1-1:R1$$$V2.0)",
    /* in line helm with incorrectly numbered r groups */
    R"(PEPTIDE1{A.[[*:4]C(=O)[C@H](C)N([*:1])C].A}$$$$V2.0)",
    /* dot before the branch */
    R"(RNA1{R.([NC1=NC(N([*:1])C=C1)=O])P}$$$$)",
    /* rna wildcard to rna strand both ends (incorrectly) fully specified */
    R"(RNA1{*}|RNA2{R(A)P.R(C)P}$RNA1,RNA2,1:R1-1:R1$$$V2.0)",
    /* chems connected to each other incorrectly joined */
    R"(CHEM1{[Az].[MCC]}$$$$)",
    /* unnumbered connection point brackets 3 */
    R"(PEPTIDE1{A.[C(=O)[C@H](C)N([*])C]}$$$$)",
    /* unnumbered connection point no brackets */
    R"(PEPTIDE1{[C(=O)[C@H](C)N(*)C].P}$$$$)",
    /* unnumbered connection point brackets 2 */
    R"(PEPTIDE1{[[*]C(=O)[C@H](C)N([*])C].P}$$$$)",
    /* blob peptide blob location specified peptide unspecified */
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,1:?-?:?$$$V2.0)",
    /* chems connected to each other   incorrectly joined */
    R"(CHEM1{[Az].[MCC]}$$$$)",
    /* unspecified peptide connection specified chem */
    R"(CHEM1{*}|PEPTIDE1{A.C}$CHEM1,PEPTIDE1,1:R1-?:?$$$V2.0)",
    /* unnumbered connection point no brackets 4 */
    R"(PEPTIDE1{A.[*C(=O)[C@H](C)N(*)C]}$$$$)",
    /* unspecified peptide connection and specified chem */
    R"(CHEM1{*}|PEPTIDE1{A.C}$CHEM1,PEPTIDE1,1:R1-?:?$$$V2.0)",
    /* unnumbered connection point brackets 5 */
    R"(PEPTIDE1{A.[[*]C(=O)[C@H](C)N([*])C].A}$$$$)",
    /* dot before and after the branch */
    R"(RNA1{R.([NC1=NC(N([*:1])C=C1)=O]).P}$$$$)",
    /* blob peptide blob specified peptide unspecified */
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,1:R1-?:?$$$V2.0)",
    /* unnumbered connection point no brackets 3 */
    R"(PEPTIDE1{A.[C(=O)[C@H](C)N(*)C]}$$$$)",
    /* blob peptide blob r group specified peptide unspecified */
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,?:R1-?:?$$$V2.0)",
    /* specified peptide connection specified chem */
    R"(CHEM1{*}|PEPTIDE1{A.C}$CHEM1,PEPTIDE1,1:R1-1:R1$$$V2.0)",
    /* unknown peptide peptide and peptide side fully specified */
    R"(PEPTIDE1{*}|PEPTIDE2{A.C}$PEPTIDE1,PEPTIDE2,1:?-?:R1$$$V2.0)",
    /* unnumbered connection point brackets 4 */
    R"(PEPTIDE1{A.[[*]C(=O)[C@H](C)N([*])C]}$$$$)",
    /* use of multiple chems */
    R"(CHEM1{*}|CHEM2{[SMCC]}$CHEM2,CHEM1,1:R1-1:R1$$$V2.0)",
    /* Unsupported version */
    "CHEM1{*}$$$$V3.0",
    /* Connections section with | prefix */
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$|BLOB1,PEPTIDE1,?:?-1:R1$$$V2.0)",
    /* Smiles monomer in mixture */
    R"(RNA1{R.([NC1=NC(N([*:1])C=C1)=O],C,G).P}$$$$)",
    R"(RNA1{R.([NC1=NC(N([*:1])C=C1)=O]+C+G).P}$$$$)",
    R"(RNA1{R.([NC1=NC(N([*:1])C=C1)=O]+C,G).P}$$$$)",
    /* union and exclusive monomer lists in mixture */
    "PEPTIDE1{(A+G,L).L}$$$$V2.0",
    "PEPTIDE1{(A+G+P,L).L}$$$$V2.0",
    "PEPTIDE1{(A:0.4+G+P,L).L}$$$$V2.0",
    "PEPTIDE1{(A:0.4+G+P,L).L}$$$$V2.0",
    "PEPTIDE1{(A:?+G+P,L).L}$$$$V2.0",
};

static const std::vector<std::string>
    VALID_EXAMPLES{
        /* Monomer lists */
        R"(PEPTIDE1{A.(A+G+C).G.C}$$$$V2.0)",
        R"(PEPTIDE1{A.(A,G,L).G.C}$$$$V2.0)",
        R"(PEPTIDE1{A.(A:1+G:2+C:3).G.C}$$$$V2.0)",
        R"(PEPTIDE1{A.(A:1,G:9).G.C}$$$$V2.0)",
        R"(PEPTIDE1{A.A.(X,*).A.A.C.D.D.E.E}$$$$V2.0)",
        R"(PEPTIDE1{A.A.A.A.(A+G:?).A.A.A.A.(A+G+L).A.A.A.A.A.A.A.A.A.A.A.C.D.D.D.D.D.D.D.D.D.D.D.D.D}$$$$V2.0)",
        R"(PEPTIDE1{A.A.A.A.(A:1.1,G:69.5,W:25.5,[Nal]:3.9).A.A.C.D.D.E.E}$$$$V2.0)",
        R"(PEPTIDE1{[fmoc].D.(V,[dV]).Y.A}$$$$V2.0)",
        /* Extended annotations */
        R"(BLOB1{Bead}$$${"Name":"Gold particle conjugated with peptides","Load":26}$V2.0)",
        R"(PEPTIDE1{C.C}|BLOB1{Gold Particle}"Au10, Diameter:10nm"$PEPTIDE1,BLOB1,C:R3-?:?$G1(PEPTIDE1:20-34+BLOB1)${"Name":"Gold particle conjugated with peptides","Load":26}$V2.0)",
        R"(PEPTIDE1{L.V.A}$$${"PEPTIDE1":{"ChainType":"hc"}}$V2.0)",
        R"(RNA1{R(A)P.R(C)P.R(G)P}|RNA2{P.R(C)P.R(G)P.R(T)}$RNA1,RNA2,2:pair-9:pair|RNA1,RNA2,5:pair-6:pair|RNA1,RNA2,8:pair-3:pair$${"RNA1":{"strandtype":"ss"}}$V2.0)",
        R"(RNA1{R(A)P.R(C)P.R(G)}$$${"my chain":"my annotation"}$V2.0)",
        /* Monomer repetitions */
        R"(PEPTIDE1{A'11'}$$$$V2.0)",
        R"(PEPTIDE1{A'2'}$$$$V2.0)",
        R"(PEPTIDE1{A.(G.A.C)'3'.A}$$$$V2.0)",
        R"(PEPTIDE1{A.G'70-100'"repeatingMonomer"}$$$$V2.0)",
        R"(PEPTIDE1{A.G.A.C'5'.A}$$$$V2.0)",
        R"(PEPTIDE1{A.G.A.C.A'5-30'}$$$$V2.0)",
        R"(PEPTIDE1{[Bal]'11'}$$$$V2.0)",
        R"(RNA1{R(A)P.(R(N)P)'4'.(R(G)P)'3-7'}$$$$V2.0)",
        R"(RNA1{R(A)P.R(C)P.(R(T)P)'2'.R(G)}$$$$V2.0)",
        R"(RNA1{R(C)P.(R(U)P)'3'.R(A)P}$$$$V2.0)",
        /* Inline smiles */
        "CHEM1{[c1ccccc1]}$$$$V2.0",
        R"(CHEM1{[[*:1]OCCCCC=C]}$$$$)",
        R"(CHEM1{[[*:1]OCCCCC=C]}|CHEM2{[[*:2]C(=O)[C@H](C)N([*:1])CC]}$CHEM1,CHEM2,1:R1-1:R1$$$V2.0)",
        R"(CHEM1{[[*:2]C(=O)[C@H](C)N([*:1])CC]}$$$$)",
        "PEPTIDE1{C.C([O[*:1]])C}$$$$V2.0",
        R"(PEPTIDE1{A.[C(=O)[C@H](C)N([*:1])C]}$$$$)",
        R"(PEPTIDE1{A.[[*:2]C(=O)[C@H](C)N([*:1])C].A}$$$$)",
        R"(PEPTIDE1{A.[[*:2]C(=O)[C@H](C)N([*:1])C]}$$$$)",
        R"(PEPTIDE1{[C(=O)[C@H](C)N([*:1])C]}$$$$V2.0)",
        R"(PEPTIDE1{[[*:2]C(=O)[C@H](C)N([*:1])C]}$$$$)",
        R"(RNA1{P.[O[C@H]1[C@H]([*:1])O[C@H](CO[*:2])[C@H]1O[*:3]](U)P.R(T)}$$$$V2.0)",
        R"(RNA1{R(C)[OP([*:1])([*:2])=O]}$$$$V2.0)",
        R"(RNA1{[O[C@H]1[C@H]([*:3])O[C@H](CO[*:1])[C@H]1O[*:2]](C)P}$$$$V2.0)",
        /* Polymer groups */
        R"(PEPTIDE1{A.C.D}"Polymer 1"|PEPTIDE2{D.A.C}"Polymer 2"$$G3(PEPTIDE1,PEPTIDE2)$$V2.0)",
        R"(PEPTIDE1{A.C.D}"Polymer 1"|PEPTIDE2{D.A.C}"Polymer 2"$$G3(PEPTIDE1:2,PEPTIDE2:5)$$V2.0)",
        R"(PEPTIDE1{A.C.D}|PEPTIDE2{D.A.C}$$G3(PEPTIDE1+PEPTIDE2)$$V2.0)",
        R"(PEPTIDE1{A.C.D}|PEPTIDE2{D.A.C}$$G3(PEPTIDE1,PEPTIDE2)$$V2.0)",
        R"(PEPTIDE1{A.C.D}|PEPTIDE2{D.A.C}$$G3(PEPTIDE1:2+PEPTIDE2:1)$$V2.0)",
        R"(RNA1{R(A)P.R(C)P.R(G)}|RNA2{R(A)P.R(C)P.R(G)}$$G3(RNA1+RNA2)$$V2.0)",
        R"(RNA1{R(A)P.R(C)P.R(G)}|RNA2{R(A)P.R(C)P.R(G)}$$G3(RNA1,RNA2)$$V2.0)",
        /* Simple polymers */
        R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$PEPTIDE1,BLOB1,1:R1-?:?$$$V2.0)",
        R"(CHEM1{*}$$$$V2.0)",
        R"(CHEM1{[SMCC]}|PEPTIDE1{L.M}|RNA1{R(C)P.R(A)P}$RNA1,PEPTIDE1,6:R2-1:R1|PEPTIDE1,CHEM1,2:R2-1:R1$$$)",
        R"(CHEM1{[SMPEG2]}|PEPTIDE1{A}|RNA1{R(A)P}$$$$)",
        R"(CHEM1{[sDBl]}$$$$)",
        R"(PEPTIDE1{A.*.G.C}$$$$V2.0)",
        R"(PEPTIDE1{A.L.C}$$$$)",
        R"(PEPTIDE1{A.X.G.C}$$$$V2.0)",
        R"(PEPTIDE1{A.[*C(=O)[C@H](C)N(*)C |$_R2;;;;;;_R1;;;$].A}$$$$)",
        R"(PEPTIDE1{A.[*C(=O)[C@H](C)N(*)C |$_R2;;;;;;_R1;;;$]}$$$$)",
        R"(PEPTIDE1{A.[C(=O)[C@H](C)N([*])C |$;;;;;;;_R1;;;$]}$$$$)",
        R"(PEPTIDE1{A.[C[C@@H](C=O)N(C)(*) |$;;;;;;_R1$|]}$$$$)",
        R"(PEPTIDE1{A.[[*]C(=O)[C@H](C)N([*])C |$_R2;;;;;;_R1;;;$].A}$$$$)",
        R"(PEPTIDE1{A.[[*]C(=O)[C@H](C)N([*])C |$_R2;;;;;;_R1;;;$]}$$$$)",
        R"(PEPTIDE1{A.[meA].C}$$$$)",
        R"(PEPTIDE1{L.V.A}|PEPTIDE2{L.V.A}$$$$)",
        R"(PEPTIDE1{[*C(=O)[C@H](C)N(*)C |$_R2;;;;;;_R1;;;$]}$$$$)",
        R"(PEPTIDE1{[C(=O)[C@H](C)N(*)C |$_;;;;;;;_R1;;;$]}$$$$)",
        R"(PEPTIDE1{[C(=O)[C@H](C)N([*])C |$_;;;;;;;_R1;;;$]}$$$$)",
        R"(PEPTIDE1{[C[C@@H](N*)C(*)=O |$;;;_R1;;_R2;;$|].[O=C(*)[C@@H](C(C)C)N* |$;;_R2;;;;;;_R1;$|].[O=C(*)[C@@H](CC(C)C)N* |$;;_R2;;;;;;;_R1;$|]}$$$$V2.0)",
        R"(PEPTIDE1{[[*]C(=O)[C@H](C)N([*])C |$_R2;;;;;;_R1;;;$]}$$$$)",
        R"(RNA1{*.R(A)P.R(C)P}$$$$V2.0)",
        R"(RNA1{P.R(U)P.R(T)P.P.P}$$$$)",
        R"(RNA1{P.R(U)P.R(T)}$$$$)",
        R"(RNA1{R(A"mutation")P.R(U)P.R(G)P}$$$$V2.0)",
        R"(RNA1{R(A)P.R(N)P.R(C)P.R(C)P.R(C)}$$$$V2.0)",
        R"(RNA1{R(U)P.R(T)P.R(G)P.R(C)P.R(A)}$$$$)",
        R"(RNA1{R(U)P.R(T)P.R(G)P.R(C)P.R([daA])}$$$$)",
        R"(RNA1{[O[C@H]1[C@H](*)O[C@H](CO*)[C@H]1O* |$;;;;;;R3;;;;;R1;;;R2$](C)P}$$$$V2.0)",
        R"(RNA1{[O[C@H]1[C@H]([*])O[C@H](CO[*])[C@H]1O[*] |$;;;;;;R3;;;;;R1;;;R2$](C)P}$$$$V2.0)",
        R"(PEPTIDE1{D.I.Q.M.T.Q.S.P.S.S.L.S.A.S.V.G.D.R.V.T.I.T.C.R.A.S.Q.D.V.N.T.A.V.A.W.Y.Q.Q.K.P.G.K.A.P.K.L.L.I.Y.S.A.S.F.L.Y.S.G.V.P.S.R.F.S.G.S.R.S.G.T.D.F.T.L.T.I.S.S.L.Q.P.E.D.F.A.T.Y.Y.C.Q.Q.H.Y.T.T.P.P.T.F.G.Q.G.T.K.V.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C.E.V.Q.L.V.E.S.G.G.G.L.V.Q.P.G.G.S.L.R.L.S.C.A.A.S.G.F.N.I.K.D.T.Y.I.H.W.V.R.Q.A.P.G.K.G.L.E.W.V.A.R.I.Y.P.T.N.G.Y.T.R.Y.A.D.S.V.K.G.R.F.T.I.S.A.D.T.S.K.N.T.A.Y.L.Q.M.N.S.L.R.A.E.D.T.A.V.Y.Y.C.S.R.W.G.G.D.G.F.Y.A.M.D.Y.W.G.Q.G.T.L.V.T.V.S.S.A.S.T.K.G.P.S.V.F.P.L.A.P.S.S.K.S.T.S.G.G.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.Q.T.Y.I.C.N.V.N.H.K.P.S.N.T.K.V.D.K.K.V.E.P.P.K.S.C.D.K.T.H.T.C.P.P.C.P.A.P.E.L.L.G.G.P.S.V.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.V.T.C.V.V.V.D.V.S.H.E.D.P.E.V.K.F.N.W.Y.V.D.G.V.E.V.H.N.A.K.T.K.P.R.E.E.Q.Y.N.S.T.Y.R.V.V.S.V.L.T.V.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.V.S.N.K.A.L.P.A.P.I.E.K.T.I.S.K.A.K.G.Q.P.R.E.P.Q.V.Y.T.L.P.P.S.R.D.E.L.T.K.N.Q.V.S.L.T.C.L.V.K.G.F.Y.P.S.D.I.A.V.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.V.L.D.S.D.G.S.F.F.L.Y.S.K.L.T.V.D.K.S.R.W.Q.Q.G.N.V.F.S.C.S.V.M.H.E.A.L.H.N.H.Y.T.Q.K.S.L.S.L.S.P.G.K.E.V.Q.L.V.E.S.G.G.G.L.V.Q.P.G.G.S.L.R.L.S.C.A.A.S.G.F.N.I.K.D.T.Y.I.H.W.V.R.Q.A.P.G.K.G.L.E.W.V.A.R.I.Y.P.T.N.G.Y.T.R.Y.A.D.S.V.K.G.R.F.T.I.S.A.D.T.S.K.N.T.A.Y.L.Q.M.N.S.L.R.A.E.D.T.A.V.Y.Y.C.S.R.W.G.G.D.G.F.Y.A.M.D.Y.W.G.Q.G.T.L.V.T.V.S.S.A.S.T.K.G.P.S.V.F.P.L.A.P.S.S.K.S.T.S.G.G.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.Q.T.Y.I.C.N.V.N.H.K.P.S.N.T.K.V.D.K.K.V.E.P.P.K.S.C.D.K.T.H.T.C.P.P.C.P.A.P.E.L.L.G.G.P.S.V.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.V.T.C.V.V.V.D.V.S.H.E.D.P.E.V.K.F.N.W.Y.V.D.G.V.E.V.H.N.A.K.T.K.P.R.E.E.Q.Y.N.S.T.Y.R.V.V.S.V.L.T.V.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.V.S.N.K.A.L.P.A.P.I.E.K.T.I.S.K.A.K.G.Q.P.R.E.P.Q.V.Y.T.L.P.P.S.R.D.E.L.T.K.N.Q.V.S.L.T.C.L.V.K.G.F.Y.P.S.D.I.A.V.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.V.L.D.S.D.G.S.F.F.L.Y.S.K.L.T.V.D.K.S.R.W.Q.Q.G.N.V.F.S.C.S.V.M.H.E.A.L.H.N.H.Y.T.Q.K.S.L.S.L.S.P.G.K.D.I.Q.M.T.Q.S.P.S.S.L.S.A.S.V.G.D.R.V.T.I.T.C.R.A.S.Q.D.V.N.T.A.V.A.W.Y.Q.Q.K.P.G.K.A.P.K.L.L.I.Y.S.A.S.F.L.Y.S.G.V.P.S.R.F.S.G.S.R.S.G.T.D.F.T.L.T.I.S.S.L.Q.P.E.D.F.A.T.Y.Y.C.Q.Q.H.Y.T.T.P.P.T.F.G.Q.G.T.K.V.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C}$$$$V2.0)",
        /* Monomer annotations */
        R"(PEPTIDE1{*"IL6"}$$$$V2.0)",
        R"(PEPTIDE1{A.A.C"mutation".D.E.E}$$$$V2.0)",
        R"(RNA1{R(A)P.R"mutation"(U)P.R(G)P}$$$$V2.0)",
        R"(RNA1{R(A)P.R(G)P"mutation".R(U)P}$$$$V2.0)",
        R"(RNA1{R(A"mutation")P.R(U)P.R(G)P}$$$$V2.0)",
        /* Polymer annotations */
        R"(BLOB1{BEAD}"Animated Polystyrene"$$$$V2.0)",
        R"(BLOB1{Bead}"Aminated Polystyrene"|PEPTIDE1{A.G.T}$$$$V2.0)",
        R"(CHEM1{[A6OH]}"Test annotation"$$$$V2.0)",
        R"(PEPTIDE1{A.C.D.D.E}"HC"|PEPTIDE2{G.C.S.S.S.P.K.K.V.K}"LC"$$$$V2.0)",
        R"(RNA1{R(A)P.R(C)P.R(G)P}"Test annotation"$$$$V2.0)",
        /* Connections */
        R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,?:?-1:R1$$$V2.0)",
        R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$PEPTIDE1,BLOB1,1:R1-?:?$$$V2.0)",
        R"(CHEM1{[PEG2]}|PEPTIDE1{W.N.D.[Pen].G.[Orn].D.A.D.G.S.G.[Cap]}$PEPTIDE1,PEPTIDE1,13:R2-4:R3|PEPTIDE1,CHEM1,1:R1-1:R1$$$)",
        R"(CHEM1{[SMCC]}|CHEM2{[EG]}$CHEM1,CHEM2,1:R1-1:R1$$$)",
        R"(CHEM1{[SMCC]}|PEPTIDE1{A.A.A.C.D.D.E.E}$PEPTIDE1,CHEM1,C:R3-1:R1"Specific Conjugation"$$$V2.0)",
        R"(CHEM1{[SMCC]}|PEPTIDE1{L.M}|RNA1{R(C)P.R(A)P}$RNA1,PEPTIDE1,6:R2-1:R1|PEPTIDE1,CHEM1,2:R2-1:R1$$$)",
        R"(CHEM1{[SMPEG2]}|CHEM2{[SMCC]}|CHEM3{[PEG2]}|CHEM4{[MCC]}$CHEM4,CHEM3,1:R1-1:R1|CHEM2,CHEM1,1:R1-1:R1$$$)",
        R"(CHEM1{[[R1]OCCCCC=C]}|CHEM2{[[R2]C(=O)[C@H](C)N([R1])CC]}$CHEM1,CHEM2,1:R1-1:R1$$$V2.0)",
        R"(CHEM1{[hxy]"Annotation"}|PEPTIDE1{A.C.V.L.L}$CHEM1,PEPTIDE1,1:R1-1:R1$$$V2.0)",
        R"(PEPTIDE1{A.A.A.A.A.C.D.D.D.D.D.D.D.E.E.E}|CHEM1{[PEG2]}$PEPTIDE1,CHEM1,C:R3-1:R1$G1(PEPTIDE1+CHEM1:2.5)$$V2.0)",
        R"(PEPTIDE1{A.A.A.A.A.C.D.D.D.D.D.E.E.E.E}|PEPTIDE2{G.G.G.C.S.S.S.S.S.P.P.K.K.K.K}|CHEM1{[PEG2]}|CHEM2{[SMCC]}$PEPTIDE1,CHEM1,C:R3-1:R1|PEPTIDE2,CHEM2,C:R3-1:R1$G1(PEPTIDE1+CHEM1:2.5)|G2(PEPTIDE2+CHEM2:1.5)$$V2.0)",
        R"(PEPTIDE1{A.A.A.C.C}|PEPTIDE2{D.L.L.L.V.V.V}|PEPTIDE3{C.D.D}"Test annotation"$PEPTIDE1,PEPTIDE3,5:R2-1:R1|PEPTIDE3,PEPTIDE2,3:R2-1:R1$$$V2.0)",
        R"(PEPTIDE1{A.R.C.A.A.K.T.C.D.A}$PEPTIDE1,PEPTIDE1,8:R3-3:R3$$$)",
        R"(PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0)",
        R"(PEPTIDE1{D.I.Q.M.T.Q.S.P.S.S.L.S.A.S.V.G.D.R.V.T.I.T.C.R.A.S.Q.D.V.N.T.A.V.A.W.Y.Q.Q.K.P.G.K.A.P.K.L.L.I.Y.S.A.S.F.L.Y.S.G.V.P.S.R.F.S.G.S.R.S.G.T.D.F.T.L.T.I.S.S.L.Q.P.E.D.F.A.T.Y.Y.C.Q.Q.H.Y.T.T.P.P.T.F.G.Q.G.T.K.V.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C.E.V.Q.L.V.E.S.G.G.G.L.V.Q.P.G.G.S.L.R.L.S.C.A.A.S.G.F.N.I.K.D.T.Y.I.H.W.V.R.Q.A.P.G.K.G.L.E.W.V.A.R.I.Y.P.T.N.G.Y.T.R.Y.A.D.S.V.K.G.R.F.T.I.S.A.D.T.S.K.N.T.A.Y.L.Q.M.N.S.L.R.A.E.D.T.A.V.Y.Y.C.S.R.W.G.G.D.G.F.Y.A.M.D.Y.W.G.Q.G.T.L.V.T.V.S.S.A.S.T.K.G.P.S.V.F.P.L.A.P.S.S.K.S.T.S.G.G.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.Q.T.Y.I.C.N.V.N.H.K.P.S.N.T.K.V.D.K.K.V.E.P.P.K.S.C.D.K.T.H.T.C.P.P.C.P.A.P.E.L.L.G.G.P.S.V.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.V.T.C.V.V.V.D.V.S.H.E.D.P.E.V.K.F.N.W.Y.V.D.G.V.E.V.H.N.A.K.T.K.P.R.E.E.Q.Y.N.S.T.Y.R.V.V.S.V.L.T.V.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.V.S.N.K.A.L.P.A.P.I.E.K.T.I.S.K.A.K.G.Q.P.R.E.P.Q.V.Y.T.L.P.P.S.R.D.E.L.T.K.N.Q.V.S.L.T.C.L.V.K.G.F.Y.P.S.D.I.A.V.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.V.L.D.S.D.G.S.F.F.L.Y.S.K.L.T.V.D.K.S.R.W.Q.Q.G.N.V.F.S.C.S.V.M.H.E.A.L.H.N.H.Y.T.Q.K.S.L.S.L.S.P.G.K.E.V.Q.L.V.E.S.G.G.G.L.V.Q.P.G.G.S.L.R.L.S.C.A.A.S.G.F.N.I.K.D.T.Y.I.H.W.V.R.Q.A.P.G.K.G.L.E.W.V.A.R.I.Y.P.T.N.G.Y.T.R.Y.A.D.S.V.K.G.R.F.T.I.S.A.D.T.S.K.N.T.A.Y.L.Q.M.N.S.L.R.A.E.D.T.A.V.Y.Y.C.S.R.W.G.G.D.G.F.Y.A.M.D.Y.W.G.Q.G.T.L.V.T.V.S.S.A.S.T.K.G.P.S.V.F.P.L.A.P.S.S.K.S.T.S.G.G.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.Q.T.Y.I.C.N.V.N.H.K.P.S.N.T.K.V.D.K.K.V.E.P.P.K.S.C.D.K.T.H.T.C.P.P.C.P.A.P.E.L.L.G.G.P.S.V.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.V.T.C.V.V.V.D.V.S.H.E.D.P.E.V.K.F.N.W.Y.V.D.G.V.E.V.H.N.A.K.T.K.P.R.E.E.Q.Y.N.S.T.Y.R.V.V.S.V.L.T.V.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.V.S.N.K.A.L.P.A.P.I.E.K.T.I.S.K.A.K.G.Q.P.R.E.P.Q.V.Y.T.L.P.P.S.R.D.E.L.T.K.N.Q.V.S.L.T.C.L.V.K.G.F.Y.P.S.D.I.A.V.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.V.L.D.S.D.G.S.F.F.L.Y.S.K.L.T.V.D.K.S.R.W.Q.Q.G.N.V.F.S.C.S.V.M.H.E.A.L.H.N.H.Y.T.Q.K.S.L.S.L.S.P.G.K.D.I.Q.M.T.Q.S.P.S.S.L.S.A.S.V.G.D.R.V.T.I.T.C.R.A.S.Q.D.V.N.T.A.V.A.W.Y.Q.Q.K.P.G.K.A.P.K.L.L.I.Y.S.A.S.F.L.Y.S.G.V.P.S.R.F.S.G.S.R.S.G.T.D.F.T.L.T.I.S.S.L.Q.P.E.D.F.A.T.Y.Y.C.Q.Q.H.Y.T.T.P.P.T.F.G.Q.G.T.K.V.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C}$$$$V2.0)",
        R"(PEPTIDE1{D.I.V.M.T.Q.S.P.L.S.L.P.V.T.P.G.E.P.A.S.I.S.C.R.S.S.Q.S.L.L.Y.S.I.G.Y.N.Y.L.D.W.Y.L.Q.K.S.G.Q.S.P.Q.L.L.I.Y.L.G.S.N.R.A.S.G.V.P.D.R.F.S.G.S.G.S.G.T.D.F.T.L.K.I.S.R.V.E.A.E.D.V.G.F.Y.Y.C.M.Q.A.L.Q.T.P.Y.T.F.G.Q.G.T.K.L.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C}|PEPTIDE2{E.V.Q.L.V.E.S.G.G.G.L.E.Q.P.G.G.S.L.R.L.S.C.A.G.S.G.F.T.F.R.D.Y.A.M.T.W.V.R.Q.A.P.G.K.G.L.E.W.V.S.S.I.S.G.S.G.G.N.T.Y.Y.A.D.S.V.K.G.R.F.T.I.S.R.D.N.S.K.N.T.L.Y.L.Q.M.N.S.L.R.A.E.D.T.A.V.Y.Y.C.A.K.D.R.L.S.I.T.I.R.P.R.Y.Y.G.L.D.V.W.G.Q.G.T.T.V.T.V.S.S.A.S.T.K.G.P.S.V.F.P.L.A.P.C.S.R.S.T.S.E.S.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.K.T.Y.T.C.N.V.D.H.K.P.S.N.T.K.V.D.K.R.V.E.S.K.Y.G.P.P.C.P.P.C.P.A.P.E.F.L.G.G.P.S.V.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.V.T.C.V.V.V.D.V.S.Q.E.D.P.E.V.Q.F.N.W.Y.V.D.G.V.E.V.H.N.A.K.T.K.P.R.E.E.Q.F.N.S.T.Y.R.V.V.S.V.L.T.V.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.V.S.N.K.G.L.P.S.S.I.E.K.T.I.S.K.A.K.G.Q.P.R.E.P.Q.V.Y.T.L.P.P.S.Q.E.E.M.T.K.N.Q.V.S.L.T.C.L.V.K.G.F.Y.P.S.D.I.A.V.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.V.L.D.S.D.G.S.F.F.L.Y.S.R.L.T.V.D.K.S.R.W.Q.E.G.N.V.F.S.C.S.V.M.H.E.A.L.H.N.H.Y.T.Q.K.S.L.S.L.S.L.G}|PEPTIDE3{E.V.Q.L.V.E.S.G.G.G.L.E.Q.P.G.G.S.L.R.L.S.C.A.G.S.G.F.T.F.R.D.Y.A.M.T.W.V.R.Q.A.P.G.K.G.L.E.W.V.S.S.I.S.G.S.G.G.N.T.Y.Y.A.D.S.V.K.G.R.F.T.I.S.R.D.N.S.K.N.T.L.Y.L.Q.M.N.S.L.R.A.E.D.T.A.V.Y.Y.C.A.K.D.R.L.S.I.T.I.R.P.R.Y.Y.G.L.D.V.W.G.Q.G.T.T.V.T.V.S.S.A.S.T.K.G.P.S.V.F.P.L.A.P.C.S.R.S.T.S.E.S.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.K.T.Y.T.C.N.V.D.H.K.P.S.N.T.K.V.D.K.R.V.E.S.K.Y.G.P.P.C.P.P.C.P.A.P.E.F.L.G.G.P.S.V.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.V.T.C.V.V.V.D.V.S.Q.E.D.P.E.V.Q.F.N.W.Y.V.D.G.V.E.V.H.N.A.K.T.K.P.R.E.E.Q.F.N.S.T.Y.R.V.V.S.V.L.T.V.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.V.S.N.K.G.L.P.S.S.I.E.K.T.I.S.K.A.K.G.Q.P.R.E.P.Q.V.Y.T.L.P.P.S.Q.E.E.M.T.K.N.Q.V.S.L.T.C.L.V.K.G.F.Y.P.S.D.I.A.V.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.V.L.D.S.D.G.S.F.F.L.Y.S.R.L.T.V.D.K.S.R.W.Q.E.G.N.V.F.S.C.S.V.M.H.E.A.L.H.N.H.Y.T.Q.K.S.L.S.L.S.L.G}|PEPTIDE4{D.I.V.M.T.Q.S.P.L.S.L.P.V.T.P.G.E.P.A.S.I.S.C.R.S.S.Q.S.L.L.Y.S.I.G.Y.N.Y.L.D.W.Y.L.Q.K.S.G.Q.S.P.Q.L.L.I.Y.L.G.S.N.R.A.S.G.V.P.D.R.F.S.G.S.G.S.G.T.D.F.T.L.K.I.S.R.V.E.A.E.D.V.G.F.Y.Y.C.M.Q.A.L.Q.T.P.Y.T.F.G.Q.G.T.K.L.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C}$PEPTIDE1,PEPTIDE1,23:R3-93:R3|PEPTIDE2,PEPTIDE2,372:R3-430:R3|PEPTIDE2,PEPTIDE2,266:R3-326:R3|PEPTIDE3,PEPTIDE3,266:R3-326:R3|PEPTIDE3,PEPTIDE3,372:R3-430:R3|PEPTIDE4,PEPTIDE4,23:R3-93:R3|PEPTIDE3,PEPTIDE3,152:R3-208:R3|PEPTIDE2,PEPTIDE3,231:R3-231:R3|PEPTIDE2,PEPTIDE3,234:R3-234:R3|PEPTIDE3,PEPTIDE4,139:R3-219:R3|PEPTIDE2,PEPTIDE1,139:R3-219:R3|PEPTIDE2,PEPTIDE2,152:R3-208:R3|PEPTIDE2,PEPTIDE2,22:R3-96:R3|PEPTIDE3,PEPTIDE3,22:R3-96:R3|PEPTIDE4,PEPTIDE4,139:R3-199:R3|PEPTIDE1,PEPTIDE1,139:R3-199:R3$$$)",
        R"(PEPTIDE1{E.I.V.L.T.Q.S.P.A.T.L.S.L.S.P.G.E.R.A.T.L.S.C.R.A.S.K.G.V.S.T.S.G.Y.S.Y.L.H.W.Y.Q.Q.K.P.G.Q.A.P.R.L.L.I.Y.L.A.S.Y.L.E.S.G.V.P.A.R.F.S.G.S.G.S.G.T.D.F.T.L.T.I.S.S.L.E.P.E.D.F.A.V.Y.Y.C.Q.H.S.R.D.L.P.L.T.F.G.G.G.T.K.V.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C}|PEPTIDE2{Q.V.Q.L.V.Q.S.G.V.E.V.K.K.P.G.A.S.V.K.V.S.C.K.A.S.G.Y.T.F.T.N.Y.Y.M.Y.W.V.R.Q.A.P.G.Q.G.L.E.W.M.G.G.I.N.P.S.N.G.G.T.N.F.N.E.K.F.K.N.R.V.T.L.T.T.D.S.S.T.T.T.A.Y.M.E.L.K.S.L.Q.F.D.D.T.A.V.Y.Y.C.A.R.R.D.Y.R.F.D.M.G.F.D.Y.W.G.Q.G.T.T.V.T.V.S.S.A.S.T.K.G.P.S.V.F.P.L.A.P.C.S.R.S.T.S.E.S.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.K.T.Y.T.C.N.V.D.H.K.P.S.N.T.K.V.D.K.R.V.E.S.K.Y.G.P.P.C.P.P.C.P.A.P.E.F.L.G.G.P.S.V.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.V.T.C.V.V.V.D.V.S.Q.E.D.P.E.V.Q.F.N.W.Y.V.D.G.V.E.V.H.N.A.K.T.K.P.R.E.E.Q.F.N.S.T.Y.R.V.V.S.V.L.T.V.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.V.S.N.K.G.L.P.S.S.I.E.K.T.I.S.K.A.K.G.Q.P.R.E.P.Q.V.Y.T.L.P.P.S.Q.E.E.M.T.K.N.Q.V.S.L.T.C.L.V.K.G.F.Y.P.S.D.I.A.V.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.V.L.D.S.D.G.S.F.F.L.Y.S.R.L.T.V.D.K.S.R.W.Q.Q.G.N.V.F.S.C.S.V.M.H.E.A.L.H.N.H.Y.T.Q.K.S.L.S.L.S}|PEPTIDE3{E.V.Q.L.V.E.S.G.G.G.L.E.Q.P.G.G.S.L.R.L.S.C.A.G.S.G.F.T.F.R.D.Y.A.M.T.W.V.R.Q.A.P.G.K.G.L.E.W.V.S.S.I.S.G.S.G.G.N.T.Y.Y.A.D.S.V.K.G.R.F.T.I.S.R.D.N.S.K.N.T.L.Y.L.Q.M.N.S.L.R.A.E.D.T.A.V.Y.Y.C.A.K.D.R.L.S.I.T.I.R.P.R.Y.Y.G.L.D.V.W.G.Q.G.T.T.V.T.V.S.S.A.S.T.K.G.P.S.V.F.P.L.A.P.C.S.R.S.T.S.E.S.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.K.T.Y.T.C.N.V.D.H.K.P.S.N.T.K.V.D.K.R.V.E.S.K.Y.G.P.P.C.P.P.C.P.A.P.E.F.L.G.G.P.S.V.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.V.T.C.V.V.V.D.V.S.Q.E.D.P.E.V.Q.F.N.W.Y.V.D.G.V.E.V.H.N.A.K.T.K.P.R.E.E.Q.F.N.S.T.Y.R.V.V.S.V.L.T.V.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.V.S.N.K.G.L.P.S.S.I.E.K.T.I.S.K.A.K.G.Q.P.R.E.P.Q.V.Y.T.L.P.P.S.Q.E.E.M.T.K.N.Q.V.S.L.T.C.L.V.K.G.F.Y.P.S.D.I.A.V.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.V.L.D.S.D.G.S.F.F.L.Y.S.R.L.T.V.D.K.S.R.W.Q.E.G.N.V.F.S.C.S.V.M.H.E.A.L.H.N.H.Y.T.Q.K.S.L.S.L.S.L.G}|PEPTIDE4{D.I.V.M.T.Q.S.P.L.S.L.P.V.T.P.G.E.P.A.S.I.S.C.R.S.S.Q.S.L.L.Y.S.I.G.Y.N.Y.L.D.W.Y.L.Q.K.S.G.Q.S.P.Q.L.L.I.Y.L.G.S.N.R.A.S.G.V.P.D.R.F.S.G.S.G.S.G.T.D.F.T.L.K.I.S.R.V.E.A.E.D.V.G.F.Y.Y.C.M.Q.A.L.Q.T.P.Y.T.F.G.Q.G.T.K.L.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C}$PEPTIDE2,PEPTIDE2,22:R3-96:R3|PEPTIDE3,PEPTIDE3,372:R3-430:R3|PEPTIDE3,PEPTIDE3,266:R3-326:R3|PEPTIDE4,PEPTIDE4,139:R3-199:R3|PEPTIDE1,PEPTIDE1,23:R3-92:R3|PEPTIDE3,PEPTIDE4,139:R3-219:R3|PEPTIDE3,PEPTIDE3,152:R3-208:R3|PEPTIDE2,PEPTIDE2,367:R3-425:R3|PEPTIDE3,PEPTIDE3,22:R3-96:R3|PEPTIDE2,PEPTIDE2,261:R3-321:R3|PEPTIDE2,PEPTIDE3,226:R3-231:R3|PEPTIDE2,PEPTIDE3,229:R3-234:R3|PEPTIDE1,PEPTIDE1,138:R3-198:R3|PEPTIDE2,PEPTIDE1,134:R3-218:R3|PEPTIDE4,PEPTIDE4,23:R3-93:R3|PEPTIDE2,PEPTIDE2,147:R3-203:R3$$$)",
        R"(PEPTIDE1{L.M.P.Q.R.S.T}$PEPTIDE1,PEPTIDE1,7:R2-1:R1$$$)",
        R"(PEPTIDE1{N.P.F.V.L.P.[dV]}$PEPTIDE1,PEPTIDE1,7:R2-1:R1$$$)",
        R"(PEPTIDE1{Q.V.Q.L.Q.Q.P.G.A.D.L.V.M.P.G.A.P.V.K.L.S.C.L.A.S.G.Y.I.F.T.S.S.W.I.N.W.V.K.Q.R.P.G.R.G.L.E.W.I.G.R.I.D.P.S.D.G.E.V.H.Y.N.Q.D.F.K.D.K.A.T.L.T.V.D.K.S.S.S.T.A.Y.I.Q.L.N.S.L.T.S.E.D.S.A.V.Y.Y.C.A.R.G.F.L.P.W.F.A.D.W.G.Q.G.T.L.V.T.V.S.A.A.S.T.K.G.P.S.V.F.P.L.A.P.S.S.K.S.T.S.G.G.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.Q.T.Y.I.C.N.V.N.H.K.P.S.N.T.K.V.D.K.K.V.E.P.K.S.C}|PEPTIDE2{N.I.V.M.T.Q.S.P.K.S.M.Y.V.S.I.G.E.R.V.T.L.S.C.K.A.S.E.N.V.D.T.Y.V.S.W.Y.Q.Q.K.P.E.Q.S.P.K.L.L.I.Y.G.A.S.N.R.Y.T.G.V.P.D.R.F.T.G.S.G.S.A.T.D.F.T.L.T.I.S.S.V.Q.A.E.D.L.A.D.Y.H.C.G.Q.S.Y.N.Y.P.F.T.F.G.S.G.T.K.L.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C}$PEPTIDE1,PEPTIDE1,22:R3-96:R3|PEPTIDE2,PEPTIDE1,214:R3-220:R3|PEPTIDE1,PEPTIDE1,144:R3-200:R3|PEPTIDE2,PEPTIDE2,23:R3-88:R3|PEPTIDE2,PEPTIDE2,134:R3-194:R3$$$)",
        R"(PEPTIDE1{Q.V.Q.L.Q.Q.S.G.G.E.L.A.K.P.G.A.S.V.K.V.S.C.K.A.S.G.Y.T.F.S.S.F.W.M.H.W.V.R.Q.A.P.G.Q.G.L.E.W.I.G.Y.I.N.P.R.S.G.Y.T.E.Y.N.E.I.F.R.D.K.A.T.M.T.T.D.T.S.T.S.T.A.Y.M.E.L.S.S.L.R.S.E.D.T.A.V.Y.Y.C.A.S.F.L.G.R.G.A.M.D.Y.W.G.Q.G.T.T.V.T.V.S.S}|PEPTIDE2{D.I.Q.M.T.Q.S.P.S.S.L.S.A.S.V.G.D.R.V.T.I.T.C.R.A.S.Q.D.I.S.N.Y.L.A.W.Y.Q.Q.K.P.G.K.A.P.K.L.L.I.Y.Y.T.S.K.I.H.S.G.V.P.S.R.F.S.G.S.G.S.G.T.D.Y.T.F.T.I.S.S.L.Q.P.E.D.I.A.T.Y.Y.C.Q.Q.G.N.T.F.P.Y.T.F.G.Q.G.T.K.V.E.I.K}$PEPTIDE2,PEPTIDE2,23:R3-88:R3|PEPTIDE1,PEPTIDE1,22:R3-96:R3$$$)",
        R"(PEPTIDE1{Q.V.Q.L.V.Q.S.G.V.E.V.K.K.P.G.A.S.V.K.V.S.C.K.A.S.G.Y.T.F.T.N.Y.Y.M.Y.W.V.R.Q.A.P.G.Q.G.L.E.W.M.G.G.I.N.P.S.N.G.G.T.N.F.N.E.K.F.K.N.R.V.T.L.T.T.D.S.S.T.T.T.A.Y.M.E.L.K.S.L.Q.F.D.D.T.A.V.Y.Y.C.A.R.R.D.Y.R.F.D.M.G.F.D.Y.W.G.Q.G.T.T.V.T.V.S.S}|PEPTIDE2{E.I.V.L.T.Q.S.P.A.T.L.S.L.S.P.G.E.R.A.T.L.S.C.R.A.S.K.G.V.S.T.S.G.Y.S.Y.L.H.W.Y.Q.Q.K.P.G.Q.A.P.R.L.L.I.Y.L.A.S.Y.L.E.S.G.V.P.A.R.F.S.G.S.G.S.G.T.D.F.T.L.T.I.S.S.L.E.P.E.D.F.A.V.Y.Y.C.Q.H.S.R.D.L.P.L.T.F.G.G.G.T.K.V.E.I.K}$PEPTIDE2,PEPTIDE2,23:R3-92:R3|PEPTIDE1,PEPTIDE1,22:R3-96:R3$$$)",
        R"(RNA1{R(A)P.R(C)P.R(G)P}|RNA2{P.R(C)P.R(G)P.R(T)}$RNA1,RNA2,2:pair-9:pair|RNA1,RNA2,5:pair-6:pair|RNA1,RNA2,8:pair-3:pair$$$)",
        R"(RNA1{R(U)P.R(T)P}|RNA2{P.R(A)P.R(A)}$RNA1,RNA2,2:pair-6:pair|RNA1,RNA2,5:pair-3:pair$$$)",
        /* Nonstandard monomer names */
        "PEPTIDE1{[Phe_3Cl]}$$$$V2.0",
        "PEPTIDE1{[Phe-3Cl]}$$$$V2.0",
        "PEPTIDE1{[D-1Nal]}$$$$V2.0",
        "PEPTIDE1{[D-Phe_4F]}$$$$V2.0",
        "RNA1{R.[Phe_3Cl].P}$$$$V2.0",
        "RNA1{R.[Phe-3Cl].P}$$$$V2.0",
        "RNA1{R.[D-1Nal].P}$$$$V2.0",
        "RNA1{R.[D-Phe_4F].P}$$$$V2.0",
        "CHEM1{[Phe_3Cl]}$$$$V2.0",
        "CHEM1{[Phe-3Cl]}$$$$V2.0",
        "CHEM1{[D-1Nal]}$$$$V2.0",
        "CHEM1{[D-Phe_4F]}$$$$V2.0",
    };

static const std::vector<std::string> UNSUPPORTED_EXAMPLES{
    /* ambiguous non-blob connections */
    R"(RNA1{*}|RNA2{R(A)P.R(C)P}$RNA1,RNA2,?:?-1:R1$$$V2.0)",
    R"(PEPTIDE1{A.C.D.E.E}|PEPTIDE2{G.C.S.S.P.P.K}|CHEM1{[SMCC]}$G1,CHEM1,?:?-?:?$G1(PEPTIDE1+PEPTIDE2)$$V2.0)",
    R"(PEPTIDE1{A.C.D.E}|PEPTIDE2{G.C.S.P.P.K}|CHEM1{*}$G1,CHEM1,C:R3-1:?$G1(PEPTIDE1,PEPTIDE2)$$V2.0)",
    R"(CHEM1{*}|PEPTIDE1{A.C}$CHEM1,PEPTIDE1,?:?-1:R1$$$V2.0)",
    R"(PEPTIDE1{A}|PEPTIDE2{A.C}$PEPTIDE1,PEPTIDE2,1:R1-?:?$$$V2.0)",
    R"(CHEM1{[SMCC]}|PEPTIDE1{A.A.A.C.A.C.D.D.D.E.E.E.E.E}$PEPTIDE1,CHEM1,C:?-?:?$$$V2.0)",
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,?:?-1:?$$$V2.0)",
    R"(PEPTIDE1{A.C}|RNA1{R}$RNA1,PEPTIDE1,1:R2-?:?$$$V2.0)",
    R"(CHEM1{[Az]}|PEPTIDE1{A.A.A.A.A.A.A.A.A}|PEPTIDE2{A.C.D.D.D.D.D.D.D.D.D.D.D.D.D.D.D.D.D}|PEPTIDE3{A.A.A.A.A.A.A.A.A.A.A}"Test Annotation"$PEPTIDE1,PEPTIDE3,9:R2-1:R1|PEPTIDE3,PEPTIDE2,11:R2-1:R1|PEPTIDE3,CHEM1,?:?-1:R1$$$V2.0)",
    R"(CHEM1{[SMCC]}|PEPTIDE1{A.A.A.C.A.C.D.D.D.E.E.E.E.E}$PEPTIDE1,CHEM1,C:?-?:?$$$V2.0)",
    R"(PEPTIDE1{A.C}|RNA1{*}$RNA1,PEPTIDE1,?:?-1:R2$$$V2.0)",
    R"(PEPTIDE1{A.C}|RNA1{P}$RNA1,PEPTIDE1,1:R2-?:?$$$V2.0)",
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,?:?-?:?$$$V2.0)",
    R"(PEPTIDE1{A.A.A.A.A.C.D.D.D.E.E.E.E.E}|PEPTIDE2{G.G.G.C.S.S.S.P.P.P.P.K.K.K.K.K}|CHEM1{[SMCC]}$PEPTIDE1,CHEM1,?:?-?:?$$$V2.0)",
    R"(CHEM1{[Az]}|PEPTIDE1{A.A.A.A.A.A.A.A.A}|PEPTIDE2{A.C.D.D.D.D.D.D.D.D.D.D.D.D.D.D.D.D.D}|PEPTIDE3{A.A.A.A.A.A.A.A.A.A.A}"Test Annotation"$PEPTIDE1,PEPTIDE3,9:R2-1:R1|PEPTIDE3,PEPTIDE2,11:R2-1:R1|PEPTIDE3,CHEM1,?:?-1:R1$$$V2.0)",
    R"(CHEM1{*}|RNA1{R(C)P.R(C)P}|RNA2{R(A)P.R(G)P.R(C)P}$G3,CHEM1,?:?-1:?$G3(RNA2+RNA1)$$V2.0)",
    R"(CHEM1{*}|PEPTIDE1{A.C}$CHEM1,PEPTIDE1,?:?-1:R1$$$V2.0)",
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,?:?-?:?$$$V2.0)",
    R"(CHEM1{[SMCC]}|RNA1{R(A)P.R(C)P.R(G)}$CHEM1,RNA1,1:R1-?:?"my annotation"$$$V2.0)",
    R"(CHEM1{*}|PEPTIDE1{A.C}$CHEM1,PEPTIDE1,?:?-?:?$$$V2.0)",
    R"(PEPTIDE1{A.C}|RNA1{R}$RNA1,PEPTIDE1,1:R2-?:?$$$V2.0)",
    R"(PEPTIDE1{A.A.A.A.A.C.D.D.D.E.E.E.E.E}|PEPTIDE2{G.G.G.C.S.S.S.P.P.P.P.K.K.K.K.K}|CHEM1{[SMCC]}$PEPTIDE1,CHEM1,?:?-?:?$$$V2.0)",
    R"(PEPTIDE1{A.C.D.E.E}|PEPTIDE2{G.C.S.S.P.P.K}|CHEM1{[SMCC]}$G1,CHEM1,?:?-?:?$G1(PEPTIDE1+PEPTIDE2)$$V2.0)",
    R"(PEPTIDE1{A.C}|RNA1{*}$RNA1,PEPTIDE1,?:?-1:R2$$$V2.0)",
    R"(RNA1{*}|RNA2{R(A)P.R(C)P}$RNA1,RNA2,?:?-1:R1$$$V2.0)",
    R"(PEPTIDE1{A.C}|RNA1{P}$RNA1,PEPTIDE1,1:R2-?:?$$$V2.0)",
    R"(CHEM1{[SMCC]}|RNA1{R(A)P.R(C)P.R(G)}$CHEM1,RNA1,1:R1-?:?"my annotation"$$$V2.0)",
    R"(CHEM1{*}|PEPTIDE1{A.C}$CHEM1,PEPTIDE1,?:?-?:?$$$V2.0)",
    R"(PEPTIDE1{A}|PEPTIDE2{A.C}$PEPTIDE1,PEPTIDE2,1:R1-?:?$$$V2.0)",
    R"(PEPTIDE1{A.C.D.E}|PEPTIDE2{G.C.S.P.P.K}|CHEM1{*}$G1,CHEM1,C:R3-1:?$G1(PEPTIDE1,PEPTIDE2)$$V2.0)",
    R"(CHEM1{*}|RNA1{R(C)P.R(C)P}|RNA2{R(A)P.R(G)P.R(C)P}$G3,CHEM1,?:?-1:?$G3(RNA2+RNA1)$$V2.0)",
    R"(BLOB1{Bead}|PEPTIDE1{A.G.T}$BLOB1,PEPTIDE1,?:?-1:?$$$V2.0)",
    R"(PEPTIDE1{D.I.Q.M.T.Q.S.P.S.S.L.S.A.S.V.G.D.R.V.T.I.T.C.R.A.S.Q.D.V.N.T.A.V.A.W.Y.Q.Q.K.P.G.K.A.P.K.L.L.I.Y.S.A.S.F.L.Y.S.G.V.P.S.R.F.S.G.S.R.S.G.T.D.F.T.L.T.I.S.S.L.Q.P.E.D.F.A.T.Y.Y.C.Q.Q.H.Y.T.T.P.P.T.F.G.Q.G.T.K.V.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C}|PEPTIDE2{E.V.Q.L.V.E.S.G.G.G.L.V.Q.P.G.G.S.L.R.L.S.C.A.A.S.G.F.N.I.K.D.T.Y.I.H.W.V.R.Q.A.P.G.K.G.L.E.W.V.A.R.I.Y.P.T.N.G.Y.T.R.Y.A.D.S.V.K.G.R.F.T.I.S.A.D.T.S.K.N.T.A.Y.L.Q.M.N.S.L.R.A.E.D.T.A.V.Y.Y.C.S.R.W.G.G.D.G.F.Y.A.M.D.Y.W.G.Q.G.T.L.V.T.V.S.S.A.S.T.K.G.P.S.V.F.P.L.A.P.S.S.K.S.T.S.G.G.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.Q.T.Y.I.C.N.V.N.H.K.P.S.N.T.K.V.D.K.K.V.E.P.P.K.S.C.D.K.T.H.T.C.P.P.C.P.A.P.E.L.L.G.G.P.S.V.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.V.T.C.V.V.V.D.V.S.H.E.D.P.E.V.K.F.N.W.Y.V.D.G.V.E.V.H.N.A.K.T.K.P.R.E.E.Q.Y.N.S.T.Y.R.V.V.S.V.L.T.V.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.V.S.N.K.A.L.P.A.P.I.E.K.T.I.S.K.A.K.G.Q.P.R.E.P.Q.V.Y.T.L.P.P.S.R.D.E.L.T.K.N.Q.V.S.L.T.C.L.V.K.G.F.Y.P.S.D.I.A.V.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.V.L.D.S.D.G.S.F.F.L.Y.S.K.L.T.V.D.K.S.R.W.Q.Q.G.N.V.F.S.C.S.V.M.H.E.A.L.H.N.H.Y.T.Q.K.S.L.S.L.S.P.G.K}|PEPTIDE3{E.V.Q.L.V.E.S.G.G.G.L.V.Q.P.G.G.S.L.R.L.S.C.A.A.S.G.F.N.I.K.D.T.Y.I.H.W.V.R.Q.A.P.G.K.G.L.E.W.V.A.R.I.Y.P.T.N.G.Y.T.R.Y.A.D.S.V.K.G.R.F.T.I.S.A.D.T.S.K.N.T.A.Y.L.Q.M.N.S.L.R.A.E.D.T.A.V.Y.Y.C.S.R.W.G.G.D.G.F.Y.A.M.D.Y.W.G.Q.G.T.L.V.T.V.S.S.A.S.T.K.G.P.S.V.F.P.L.A.P.S.S.K.S.T.S.G.G.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.Q.T.Y.I.C.N.V.N.H.K.P.S.N.T.K.V.D.K.K.V.E.P.P.K.S.C.D.K.T.H.T.C.P.P.C.P.A.P.E.L.L.G.G.P.S.V.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.V.T.C.V.V.V.D.V.S.H.E.D.P.E.V.K.F.N.W.Y.V.D.G.V.E.V.H.N.A.K.T.K.P.R.E.E.Q.Y.N.S.T.Y.R.V.V.S.V.L.T.V.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.V.S.N.K.A.L.P.A.P.I.E.K.T.I.S.K.A.K.G.Q.P.R.E.P.Q.V.Y.T.L.P.P.S.R.D.E.L.T.K.N.Q.V.S.L.T.C.L.V.K.G.F.Y.P.S.D.I.A.V.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.V.L.D.S.D.G.S.F.F.L.Y.S.K.L.T.V.D.K.S.R.W.Q.Q.G.N.V.F.S.C.S.V.M.H.E.A.L.H.N.H.Y.T.Q.K.S.L.S.L.S.P.G.K}|PEPTIDE4{D.I.Q.M.T.Q.S.P.S.S.L.S.A.S.V.G.D.R.V.T.I.T.C.R.A.S.Q.D.V.N.T.A.V.A.W.Y.Q.Q.K.P.G.K.A.P.K.L.L.I.Y.S.A.S.F.L.Y.S.G.V.P.S.R.F.S.G.S.R.S.G.T.D.F.T.L.T.I.S.S.L.Q.P.E.D.F.A.T.Y.Y.C.Q.Q.H.Y.T.T.P.P.T.F.G.Q.G.T.K.V.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C}|CHEM1{[SMCC]}|CHEM2{[PEG2]}$CHEM1,CHEM2,1:R2-1:R1|PEPTIDE1,PEPTIDE1,23:R3-88:R3|PEPTIDE2,PEPTIDE2,371:R3-429:R3|PEPTIDE4,PEPTIDE4,134:R3-194:R3|PEPTIDE1,PEPTIDE1,134:R3-194:R3|PEPTIDE3,PEPTIDE3,265:R3-325:R3|PEPTIDE3,PEPTIDE4,224:R3-214:R3|PEPTIDE2,PEPTIDE3,233:R3-233:R3|PEPTIDE2,PEPTIDE2,265:R3-325:R3|PEPTIDE2,PEPTIDE2,147:R3-203:R3|PEPTIDE2,PEPTIDE3,230:R3-230:R3|PEPTIDE2,PEPTIDE1,224:R3-214:R3|PEPTIDE3,PEPTIDE3,147:R3-203:R3|PEPTIDE2,PEPTIDE2,22:R3-96:R3|PEPTIDE3,PEPTIDE3,22:R3-96:R3|PEPTIDE3,PEPTIDE3,371:R3-429:R3|PEPTIDE4,PEPTIDE4,23:R3-88:R3|G1,CHEM1,K:R3-1:R1$G1(PEPTIDE1+PEPTIDE2+PEPTIDE3+PEPTIDE4)$$V2.0)",
    R"(PEPTIDE1{D.I.Q.M.T.Q.S.P.S.S.L.S.A.S.V.G.D.R.V.T.I.T.C.R.A.S.Q.D.V.N.T.A.V.A.W.Y.Q.Q.K.P.G.K.A.P.K.L.L.I.Y.S.A.S.F.L.Y.S.G.V.P.S.R.F.S.G.S.R.S.G.T.D.F.T.L.T.I.S.S.L.Q.P.E.D.F.A.T.Y.Y.C.Q.Q.H.Y.T.T.P.P.T.F.G.Q.G.T.K.V.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C}|PEPTIDE2{E.V.Q.L.V.E.S.G.G.G.L.V.Q.P.G.G.S.L.R.L.S.C.A.A.S.G.F.N.I.K.D.T.Y.I.H.W.V.R.Q.A.P.G.K.G.L.E.W.V.A.R.I.Y.P.T.N.G.Y.T.R.Y.A.D.S.V.K.G.R.F.T.I.S.A.D.T.S.K.N.T.A.Y.L.Q.M.N.S.L.R.A.E.D.T.A.V.Y.Y.C.S.R.W.G.G.D.G.F.Y.A.M.D.Y.W.G.Q.G.T.L.V.T.V.S.S.A.S.T.K.G.P.S.V.F.P.L.A.P.S.S.K.S.T.S.G.G.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.Q.T.Y.I.C.N.V.N.H.K.P.S.N.T.K.V.D.K.K.V.E.P.P.K.S.C.D.K.T.H.T.C.P.P.C.P.A.P.E.L.L.G.G.P.S.V.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.V.T.C.V.V.V.D.V.S.H.E.D.P.E.V.K.F.N.W.Y.V.D.G.V.E.V.H.N.A.K.T.K.P.R.E.E.Q.Y.N.S.T.Y.R.V.V.S.V.L.T.V.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.V.S.N.K.A.L.P.A.P.I.E.K.T.I.S.K.A.K.G.Q.P.R.E.P.Q.V.Y.T.L.P.P.S.R.D.E.L.T.K.N.Q.V.S.L.T.C.L.V.K.G.F.Y.P.S.D.I.A.V.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.V.L.D.S.D.G.S.F.F.L.Y.S.K.L.T.V.D.K.S.R.W.Q.Q.G.N.V.F.S.C.S.V.M.H.E.A.L.H.N.H.Y.T.Q.K.S.L.S.L.S.P.G.K}|PEPTIDE3{E.V.Q.L.V.E.S.G.G.G.L.V.Q.P.G.G.S.L.R.L.S.C.A.A.S.G.F.N.I.K.D.T.Y.I.H.W.V.R.Q.A.P.G.K.G.L.E.W.V.A.R.I.Y.P.T.N.G.Y.T.R.Y.A.D.S.V.K.G.R.F.T.I.S.A.D.T.S.K.N.T.A.Y.L.Q.M.N.S.L.R.A.E.D.T.A.V.Y.Y.C.S.R.W.G.G.D.G.F.Y.A.M.D.Y.W.G.Q.G.T.L.V.T.V.S.S.A.S.T.K.G.P.S.V.F.P.L.A.P.S.S.K.S.T.S.G.G.T.A.A.L.G.C.L.V.K.D.Y.F.P.E.P.V.T.V.S.W.N.S.G.A.L.T.S.G.V.H.T.F.P.A.V.L.Q.S.S.G.L.Y.S.L.S.S.V.V.T.V.P.S.S.S.L.G.T.Q.T.Y.I.C.N.V.N.H.K.P.S.N.T.K.V.D.K.K.V.E.P.P.K.S.C.D.K.T.H.T.C.P.P.C.P.A.P.E.L.L.G.G.P.S.V.F.L.F.P.P.K.P.K.D.T.L.M.I.S.R.T.P.E.V.T.C.V.V.V.D.V.S.H.E.D.P.E.V.K.F.N.W.Y.V.D.G.V.E.V.H.N.A.K.T.K.P.R.E.E.Q.Y.N.S.T.Y.R.V.V.S.V.L.T.V.L.H.Q.D.W.L.N.G.K.E.Y.K.C.K.V.S.N.K.A.L.P.A.P.I.E.K.T.I.S.K.A.K.G.Q.P.R.E.P.Q.V.Y.T.L.P.P.S.R.D.E.L.T.K.N.Q.V.S.L.T.C.L.V.K.G.F.Y.P.S.D.I.A.V.E.W.E.S.N.G.Q.P.E.N.N.Y.K.T.T.P.P.V.L.D.S.D.G.S.F.F.L.Y.S.K.L.T.V.D.K.S.R.W.Q.Q.G.N.V.F.S.C.S.V.M.H.E.A.L.H.N.H.Y.T.Q.K.S.L.S.L.S.P.G.K}|PEPTIDE4{D.I.Q.M.T.Q.S.P.S.S.L.S.A.S.V.G.D.R.V.T.I.T.C.R.A.S.Q.D.V.N.T.A.V.A.W.Y.Q.Q.K.P.G.K.A.P.K.L.L.I.Y.S.A.S.F.L.Y.S.G.V.P.S.R.F.S.G.S.R.S.G.T.D.F.T.L.T.I.S.S.L.Q.P.E.D.F.A.T.Y.Y.C.Q.Q.H.Y.T.T.P.P.T.F.G.Q.G.T.K.V.E.I.K.R.T.V.A.A.P.S.V.F.I.F.P.P.S.D.E.Q.L.K.S.G.T.A.S.V.V.C.L.L.N.N.F.Y.P.R.E.A.K.V.Q.W.K.V.D.N.A.L.Q.S.G.N.S.Q.E.S.V.T.E.Q.D.S.K.D.S.T.Y.S.L.S.S.T.L.T.L.S.K.A.D.Y.E.K.H.K.V.Y.A.C.E.V.T.H.Q.G.L.S.S.P.V.T.K.S.F.N.R.G.E.C}|CHEM1{[SMCC]}|CHEM2{[PEG2]}$CHEM1,CHEM2,1:R2-1:R1|PEPTIDE1,PEPTIDE1,23:R3-88:R3|PEPTIDE2,PEPTIDE2,371:R3-429:R3|PEPTIDE4,PEPTIDE4,134:R3-194:R3|PEPTIDE1,PEPTIDE1,134:R3-194:R3|PEPTIDE3,PEPTIDE3,265:R3-325:R3|PEPTIDE3,PEPTIDE4,224:R3-214:R3|PEPTIDE2,PEPTIDE3,233:R3-233:R3|PEPTIDE2,PEPTIDE2,265:R3-325:R3|PEPTIDE2,PEPTIDE2,147:R3-203:R3|PEPTIDE2,PEPTIDE3,230:R3-230:R3|PEPTIDE2,PEPTIDE1,224:R3-214:R3|PEPTIDE3,PEPTIDE3,147:R3-203:R3|PEPTIDE2,PEPTIDE2,22:R3-96:R3|PEPTIDE3,PEPTIDE3,22:R3-96:R3|PEPTIDE3,PEPTIDE3,371:R3-429:R3|PEPTIDE4,PEPTIDE4,23:R3-88:R3|G1,CHEM1,K:R3-1:R1$G1(PEPTIDE1+PEPTIDE2+PEPTIDE3+PEPTIDE4)|G2(CHEM1:3.5+G1:1)$$V2.0)",
    /* missing monomers */
    R"(PEPTIDE1{A.(A,_).G.C}$$$$V2.0)",
    /* HELMV1 Annotations */
    R"(BLOB1{Bead}$$$BLOB1{Something:Here}|BLOB2{Something:There}$)",
};
