#define BOOST_TEST_MODULE rdkit_extensions_monomer_coordgen

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <rdkit/GraphMol/RWMol.h>
#include "schrodinger/rdkit_extensions/helm/monomer_coordgen.h"

#include <vector>
#include <string>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"

using helm::helm_to_rdkit;
using schrodinger::rdkit_extensions::compute_monomer_mol_coords;

namespace bdata = boost::unit_test::data;

typedef std::vector<std::vector<double>> coords_t;
BOOST_TEST_DONT_PRINT_LOG_VALUE(coords_t)
typedef std::pair<std::string, coords_t> mol_data_t;
BOOST_TEST_DONT_PRINT_LOG_VALUE(mol_data_t)

BOOST_DATA_TEST_CASE(
    TestComputeMonomerCoords,
    bdata::make(std::vector<mol_data_t>{
        {"CHEM1{[[*:1]N]}|PEPTIDE1{A.C.R.K}$PEPTIDE1,CHEM1,4:R2-1:R1|PEPTIDE1,"
         "PEPTIDE1,2:R3-4:R3$$$V2.0",
         {{0, 0}, {-3, -3}, {-1.5, -3}, {0, -3}, {0, -1.5}}},

        // Linear polymers
        {"RNA1{R(U)P.R(T)P.R(G)P.R(C)P.R(A)}$$$$",
         {{0, 0},
          {0, -1.5},
          {1.5, 0},
          {3, 0},
          {3, -1.5},
          {4.5, 0},
          {6, 0},
          {6, -1.5},
          {7.5, 0},
          {9, 0},
          {9, -1.5},
          {10.5, 0},
          {12, 0},
          {12, -1.5}}},
        {"PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$", {{0, 0}, {1.5, 0}}},

        {"PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0",
         {{0, 0}, {1.5, 0}, {3, 0}, {4.5, 0}, {1.5, -1.5}, {3, -1.5}}},

        {"RNA1{R(A)P.R(G)P.R(C)P.R(U)P.R(C)P.R(C)P.R(C)}|RNA2{R(U)P.R(G)P.R("
         "G)P.R(G)P.R(G)P.R(A)P.R(G)}$RNA1,RNA2,17:pair-11:pair|RNA1,RNA2,"
         "20:pair-8:pair|RNA1,RNA2,14:pair-14:pair|RNA1,RNA2,11:pair-17:"
         "pair|RNA1,RNA2,8:pair-20:pair$$$",
         {{0, 0},       {0, -1.5},    {1.5, 0},     {3, 0},       {3, -1.5},
          {4.5, 0},     {6, 0},       {6, -1.5},    {7.5, 0},     {9, 0},
          {9, -1.5},    {10.5, 0},    {12, 0},      {12, -1.5},   {13.5, 0},
          {15, 0},      {15, -1.5},   {16.5, 0},    {18, 0},      {18, -1.5},
          {24, -4.5},   {24, -3},     {22.5, -4.5}, {21, -4.5},   {21, -3},
          {19.5, -4.5}, {18, -4.5},   {18, -3},     {16.5, -4.5}, {15, -4.5},
          {15, -3},     {13.5, -4.5}, {12, -4.5},   {12, -3},     {10.5, -4.5},
          {9, -4.5},    {9, -3},      {7.5, -4.5},  {6, -4.5},    {6, -3}}},

        {"PEPTIDE1{[Sar].[dR].[dC].[dY].[dC].[dH].[dP].[dF]}|PEPTIDE2{[Sar].["
         "dR].[dC].[dY].[dC].[dH].[dP].[dF]}$PEPTIDE1,PEPTIDE2,5:R3-3:R3|"
         "PEPTIDE1,PEPTIDE2,3:R3-5:R3$$$",
         {{0, 0},
          {1.5, 0},
          {3, 0},
          {4.5, 0},
          {6, 0},
          {7.5, 0},
          {9, 0},
          {10.5, 0},
          {9, -1.5},
          {7.5, -1.5},
          {6, -1.5},
          {4.5, -1.5},
          {3, -1.5},
          {1.5, -1.5},
          {0, -1.5},
          {-1.5, -1.5}}},

        // Snaking linear polymers
        {"PEPTIDE1{K.W.L.N.A.L.L.H.H.G.L}$$$$V2.0",
         {{0, 0},
          {1.5, 0},
          {3, 0},
          {4.5, 0},
          {6, 0},
          {6.75, -1.299038},
          {6, -2.598076},
          {4.5, -2.598076},
          {3, -2.598076},
          {1.5, -2.598076},
          {0, -2.598076}}},

        {"PEPTIDE1{K.W.L.N.A.L.L.H.H.G.L.N.C.A.K.G.V.L.A.W.L.N.A.L.L.H}$$$$V2."
         "0",
         {{0, 0},
          {1.5, 0},
          {3, 0},
          {4.5, 0},
          {6, 0},
          {7.5, 0},
          {9, 0},
          {10.5, 0},
          {11.25, -1.299038},
          {10.5, -2.598076},
          {9, -2.598076},
          {7.5, -2.598076},
          {6, -2.598076},
          {4.5, -2.598076},
          {3, -2.598076},
          {1.5, -2.598076},
          {0, -2.598076},
          {-0.75, -3.897114},
          {0, -5.196152},
          {1.5, -5.196152},
          {3, -5.196152},
          {4.5, -5.196152},
          {6, -5.196152},
          {7.5, -5.196152},
          {9, -5.196152},
          {10.5, -5.196152}}},

        // Double stranded nucleic acids
        {"RNA1{R(C).P.R(A).P.R(T).P}|RNA2{R(C).P.R(A).P.R(T).P}$RNA1,RNA2,1:R1-"
         "5:R1|RNA1,RNA2,5:R1-1:R1$$$V2.0",
         {{0, 0},
          {-1.5, 0},
          {1.5, 0},
          {3, 0},
          {3, -1.5},
          {4.5, 0},
          {6, 0},
          {6, -1.5},
          {7.5, 0},
          {3, -4.5},
          {4.5, -4.5},
          {1.5, -4.5},
          {0, -4.5},
          {0, -3},
          {-1.5, -4.5},
          {-3, -4.5},
          {-3, -3},
          {-4.5, -4.5}}},

        {"PEPTIDE1{R.A.Y}|RNA1{r(G).p.r(A).p.r(U).p.r(C).p.[mR](U).p.[mR](U).p."
         "r(A).p}|RNA2{p.r(U).p.r(A).p.r(A).p.r(G).p.r(A).p.r(U).p.r(C).p.r(A)."
         "p.r(A)}$RNA1,RNA2,2:pair-21:pair|RNA1,RNA2,5:pair-18:pair|RNA1,RNA2,"
         "8:pair-15:pair|RNA1,RNA2,11:pair-12:pair|RNA1,RNA2,14:pair-9:pair|"
         "RNA1,RNA2,17:pair-6:pair|RNA1,RNA2,20:pair-3:pair|RNA1,PEPTIDE1,21:"
         "R2-1:R1$$$V2.0",
         {{0, 0},        {1.5, 0},      {3, 0},        {-19.5, -1.5},
          {-19.5, -3},   {-18, -1.5},   {-16.5, -1.5}, {-16.5, -3},
          {-15, -1.5},   {-13.5, -1.5}, {-13.5, -3},   {-12, -1.5},
          {-10.5, -1.5}, {-10.5, -3},   {-9, -1.5},    {-7.5, -1.5},
          {-7.5, -3},    {-6, -1.5},    {-4.5, -1.5},  {-4.5, -3},
          {-3, -1.5},    {-1.5, -1.5},  {-1.5, -3},    {0, -1.5},
          {0, -6},       {-1.5, -6},    {-1.5, -4.5},  {-3, -6},
          {-4.5, -6},    {-4.5, -4.5},  {-6, -6},      {-7.5, -6},
          {-7.5, -4.5},  {-9, -6},      {-10.5, -6},   {-10.5, -4.5},
          {-12, -6},     {-13.5, -6},   {-13.5, -4.5}, {-15, -6},
          {-16.5, -6},   {-16.5, -4.5}, {-18, -6},     {-19.5, -6},
          {-19.5, -4.5}, {-21, -6},     {-22.5, -6},   {-22.5, -4.5},
          {-24, -6},     {-25.5, -6},   {-25.5, -4.5}}},

        // branching polymer connections
        {"PEPTIDE1{A.C.G.C.G}|PEPTIDE2{K.L}|PEPTIDE3{S.T}$PEPTIDE1,PEPTIDE2,2:"
         "R2-1:R1|PEPTIDE1,PEPTIDE3,5:R1-1:R3$$$V2.0",
         {{0, 0},
          {1.5, 0},
          {3, 0},
          {4.5, 0},
          {6, 0},
          {1.5, -1.5},
          {3, -1.5},
          {6, -1.5},
          {7.5, -1.5}}},

        // Cyclic polymers
        {"PEPTIDE1{A.Y.V.[Orn].L.[dF].P.F.[dF].N}$PEPTIDE1,PEPTIDE1,10:R2-1:R1$"
         "$$",
         {{2.427051, 0},
          {1.963525, -1.426585},
          {0.75, -2.308263},
          {-0.75, -2.308263},
          {-1.963525, -1.426585},
          {-2.427051, 0},
          {-1.963525, 1.426585},
          {-0.75, 2.308263},
          {0.75, 2.308263},
          {1.963525, 1.426585}}},

        {"PEPTIDE1{G.C.C.S.L.P.R.C.A.L.N.C}$PEPTIDE1,PEPTIDE1,2:R3-8:R3|"
         "PEPTIDE1,PEPTIDE1,3:R3-12:R3$$$",
         {{0, 0},
          {1.5, 0},
          {1.933013, -0.75},
          {1.5, -1.5},
          {0, -1.5},
          {-1.06066, -0.43934},
          {-1.06066, 1.06066},
          {0, 2.12132},
          {1.5, 2.12132},
          {2.799038, 1.37132},
          {3.549038, 0.072282},
          {3.549038, -1.427718}}},

        {"PEPTIDE1{C.K.G.K.G.A.K.C.S.R.L.M.Y.D.C.C.T.G.S.C.R.S.G.K.C}$PEPTIDE1,"
         "PEPTIDE1,1:R3-16:R3|PEPTIDE1,PEPTIDE1,8:R3-20:R3|PEPTIDE1,PEPTIDE1,"
         "15:R3-25:R3$$$",
         {{0, 0},
          {1.569218, -0.571149},
          {2.848456, -1.644557},
          {3.68342, -3.090756},
          {3.9734, -4.735313},
          {3.68342, -6.379871},
          {2.848456, -7.82607},
          {1.569218, -8.899478},
          {0, -9.470627},
          {-1.385819, -8.896602},
          {-2.446479, -7.835942},
          {-3.020505, -6.450122},
          {-3.020505, -4.950122},
          {-2.446479, -3.564303},
          {-1.385819, -2.503643},
          {0, -1.929618},
          {1.15119, -2.484002},
          {1.947838, -3.482967},
          {2.232158, -4.728656},
          {1.947838, -5.974345},
          {1.15119, -6.97331},
          {0, -7.527694},
          {-1.06066, -6.467034},
          {-1.06066, -4.967034},
          {0, -3.906374}}},

        {"PEPTIDE1{C.F.Y.S.W.G.N.Y.W.S.Y.Y.G.W.C.[am]}|PEPTIDE2{C.F.Y.S.W.G.N."
         "Y.W.S.Y.Y.G.W.C.[am]}|PEPTIDE3{C.W.N.Y.Y.W.G.Y.S.F.S.G.W.Y.C.G.[am]}$"
         "PEPTIDE1,PEPTIDE1,1:R3-15:R3|PEPTIDE1,PEPTIDE1,1:R1-16:R2|PEPTIDE1,"
         "PEPTIDE2,3:R1-8:R2|PEPTIDE2,PEPTIDE2,1:R3-15:R3|PEPTIDE2,PEPTIDE2,1:"
         "R1-16:R2|PEPTIDE2,PEPTIDE3,11:R3-14:R3|PEPTIDE3,PEPTIDE3,1:R3-15:R3|"
         "PEPTIDE3,PEPTIDE3,1:R1-17:R2$$$",
         {{3.295433, 1.467221},    {3.607301, 0},
          {3.295433, -1.467221},   {2.413755, -2.680747},
          {1.114717, -3.430747},   {-0.377066, -3.58754},
          {-1.80365, -3.124014},   {-2.918368, -2.120318},
          {-3.528473, -0.75},      {-3.528473, 0.75},
          {-2.918368, 2.120318},   {-1.80365, 3.124014},
          {-0.377066, 3.58754},    {1.114717, 3.430747},
          {2.413755, 2.680747},    {3.905538, 2.83754},
          {9.509234, -10.142301},  {9.821102, -8.675079},
          {9.509234, -7.207858},   {8.627556, -5.994332},
          {7.328518, -5.244332},   {5.836735, -5.08754},
          {4.41015, -5.551065},    {3.295433, -6.554761},
          {2.685328, -7.925079},   {2.685328, -9.425079},
          {3.295433, -10.795397},  {4.41015, -11.799093},
          {5.836735, -12.262619},  {7.328518, -12.105826},
          {8.627556, -11.355826},  {10.119339, -11.512619},
          {5.476149, -15.882937},  {5.788017, -17.350158},
          {5.476149, -18.81738},   {4.594471, -20.030905},
          {3.295433, -20.780905},  {1.80365, -20.937698},
          {0.377066, -20.474172},  {-0.737652, -19.470477},
          {-1.347757, -18.100158}, {-1.347757, -16.600158},
          {-0.737652, -15.22984},  {0.377066, -14.226144},
          {1.80365, -13.762619},   {3.295433, -13.919411},
          {4.594471, -14.669411},  {5.807997, -13.787734},
          {6.689675, -15.001259}}},

        // CHEM polymers
        {"CHEM1{[SMCC]}|PEPTIDE1{L.M}|RNA1{R(C)P.R(A)P}$RNA1,PEPTIDE1,6:R2-1:"
         "R1|PEPTIDE1,CHEM1,2:R2-1:R1$$$",
         {{0, 0},
          {-1.5, -1.5},
          {0, -1.5},
          {-6, -3},
          {-6, -4.5},
          {-4.5, -3},
          {-3, -3},
          {-3, -4.5},
          {-1.5, -3}}},

        {"RNA1{R(A)P.R(A)P.R(G)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P}|RNA2{R(A)P.R(A)"
         "P.R(G)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P}|CHEM1{[sDBl]}$RNA2,CHEM1,24:"
         "R2-1:R3|RNA1,CHEM1,24:R2-1:R2$$$",
         {{0, 0},       {0, -1.5},    {1.5, 0},     {3, 0},       {3, -1.5},
          {4.5, 0},     {6, 0},       {6, -1.5},    {7.5, 0},     {9, 0},
          {9, -1.5},    {10.5, 0},    {12, 0},      {12, -1.5},   {13.5, 0},
          {15, 0},      {15, -1.5},   {16.5, 0},    {18, 0},      {18, -1.5},
          {19.5, 0},    {21, 0},      {21, -1.5},   {22.5, 0},    {0, -4.5},
          {0, -6},      {1.5, -4.5},  {3, -4.5},    {3, -6},      {4.5, -4.5},
          {6, -4.5},    {6, -6},      {7.5, -4.5},  {9, -4.5},    {9, -6},
          {10.5, -4.5}, {12, -4.5},   {12, -6},     {13.5, -4.5}, {15, -4.5},
          {15, -6},     {16.5, -4.5}, {18, -4.5},   {18, -6},     {19.5, -4.5},
          {21, -4.5},   {21, -6},     {22.5, -4.5}, {24, -2.25}}},

        {"RNA1{R(A)P.R(A)P.R(G)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P}|RNA2{R(A)P.R(A)"
         "P.R(G)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P}|CHEM1{[sDBl]}$RNA2,CHEM1,1:R2-"
         "1:R3|RNA1,CHEM1,1:R2-1:R2$$$",
         {{0, 0},       {-1.5, 0},    {1.5, 0},     {3, 0},       {3, -1.5},
          {4.5, 0},     {6, 0},       {6, -1.5},    {7.5, 0},     {9, 0},
          {9, -1.5},    {10.5, 0},    {12, 0},      {12, -1.5},   {13.5, 0},
          {15, 0},      {15, -1.5},   {16.5, 0},    {18, 0},      {18, -1.5},
          {19.5, 0},    {21, 0},      {21, -1.5},   {22.5, 0},    {0, -4.5},
          {0, -6},      {1.5, -4.5},  {3, -4.5},    {3, -6},      {4.5, -4.5},
          {6, -4.5},    {6, -6},      {7.5, -4.5},  {9, -4.5},    {9, -6},
          {10.5, -4.5}, {12, -4.5},   {12, -6},     {13.5, -4.5}, {15, -4.5},
          {15, -6},     {16.5, -4.5}, {18, -4.5},   {18, -6},     {19.5, -4.5},
          {21, -4.5},   {21, -6},     {22.5, -4.5}, {-1.5, -2.25}}},

        {"RNA1{R(A)P.R(A)P.R(G)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P}|RNA2{R(A)P.R(A)"
         "P.R(G)P.R(G)P.R(C)P.R(U)P.R(A)P.R(A)P}|CHEM1{[sDBl]}$RNA2,CHEM1,12:"
         "R2-1:R3|RNA1,CHEM1,12:R2-1:R2$$$",
         {{0, 0},       {0, -1.5},    {1.5, 0},     {3, 0},       {3, -1.5},
          {4.5, 0},     {6, 0},       {6, -1.5},    {7.5, 0},     {9, 0},
          {9, -1.5},    {10.5, 0},    {12, 0},      {12, -1.5},   {13.5, 0},
          {15, 0},      {15, -1.5},   {16.5, 0},    {18, 0},      {18, -1.5},
          {19.5, 0},    {21, 0},      {21, -1.5},   {22.5, 0},    {0, -4.5},
          {0, -6},      {1.5, -4.5},  {3, -4.5},    {3, -6},      {4.5, -4.5},
          {6, -4.5},    {6, -6},      {7.5, -4.5},  {9, -4.5},    {9, -6},
          {10.5, -4.5}, {12, -4.5},   {12, -6},     {13.5, -4.5}, {15, -4.5},
          {15, -6},     {16.5, -4.5}, {18, -4.5},   {18, -6},     {19.5, -4.5},
          {21, -4.5},   {21, -6},     {22.5, -4.5}, {10.5, -2.25}}},

        // Hairpin polymers
        {"RNA1{R(G)P.R(G)P.R(C)P.R(A)P.R(C)P.R(U)P.R(U)P.R(C)P.R(G)P.R(G)P.R(U)"
         "P.R(G)P.R(C)P.R(C)}$RNA1,RNA1,11:pair-32:pair|RNA1,RNA1,5:pair-38:"
         "pair|RNA1,RNA1,14:pair-29:pair|RNA1,RNA1,8:pair-35:pair|RNA1,RNA1,2:"
         "pair-41:pair$$$",
         {{-14.3458, 2.07818},  {-14.3458, 0.578184},  {-12.8458, 2.07818},
          {-11.3458, 2.07818},  {-11.3458, 0.578184},  {-9.84578, 2.07818},
          {-8.34578, 2.07818},  {-8.34578, 0.578184},  {-6.84578, 2.07818},
          {-5.34578, 2.07818},  {-5.34578, 0.578184},  {-3.84578, 2.07818},
          {-2.34578, 2.07818},  {-2.34578, 0.578184},  {-1.11131, 2.93028},
          {0.377754, 3.11109},  {0.558559, 4.60015},   {1.78028, 2.57918},
          {2.77496, 1.45641},   {4.10315, 2.1535},     {3.13394, 3.83796e-16},
          {2.77496, -1.45641},  {4.10315, -2.1535},    {1.78028, -2.57918},
          {0.377754, -3.11109}, {0.558559, -4.60015},  {-1.11131, -2.93028},
          {-2.34578, -2.07818}, {-2.34578, -0.578184}, {-3.84578, -2.07818},
          {-5.34578, -2.07818}, {-5.34578, -0.578184}, {-6.84578, -2.07818},
          {-8.34578, -2.07818}, {-8.34578, -0.578184}, {-9.84578, -2.07818},
          {-11.3458, -2.07818}, {-11.3458, -0.578184}, {-12.8458, -2.07818},
          {-14.3458, -2.07818}, {-14.3458, -0.578184}}},

        {"RNA1{R(G)P.R(G)P.R(C)P.R(A)P.R(C)P.R(U)P.R(G)P.R(G)P.R(U)P.R(G)P.R(C)"
         "P.R(C)}$RNA1,RNA1,11:pair-26:pair|RNA1,RNA1,5:pair-32:pair|RNA1,RNA1,"
         "14:pair-23:pair|RNA1,RNA1,8:pair-29:pair|RNA1,RNA1,2:pair-35:pair$$$",
         {{-12.3789, 2.635},   {-12.3789, 1.135},   {-10.8789, 2.635},
          {-9.37886, 2.635},   {-9.37886, 1.135},   {-7.87886, 2.635},
          {-6.37886, 2.635},   {-6.37886, 1.135},   {-4.87886, 2.635},
          {-3.37886, 2.635},   {-3.37886, 1.135},   {-1.87886, 2.635},
          {-0.378856, 2.635},  {-0.378856, 1.135},  {1.10588, 2.42153},
          {2.2395, 1.43924},   {3.50138, 2.2502},   {2.6621, 3.26013e-16},
          {2.2395, -1.43924},  {3.50138, -2.2502},  {1.10588, -2.42153},
          {-0.378856, -2.635}, {-0.378856, -1.135}, {-1.87886, -2.635},
          {-3.37886, -2.635},  {-3.37886, -1.135},  {-4.87886, -2.635},
          {-6.37886, -2.635},  {-6.37886, -1.135},  {-7.87886, -2.635},
          {-9.37886, -2.635},  {-9.37886, -1.135},  {-10.8789, -2.635},
          {-12.3789, -2.635},  {-12.3789, -1.135}}},
        // Cyclic structures
        {"CHEM1{[ClAc]}|PEPTIDE1{F.C.C.C.C.C.C.C.C.C.C.G}|CHEM2{[NH2]}$CHEM1,"
         "PEPTIDE1,1:R1-9:R3|CHEM1,PEPTIDE1,1:R2-1:R1|PEPTIDE1,CHEM2,10:R2-1:"
         "R1$$$V2.0",
         {{1.963526, 1.426585},
          {0.75, 2.308263},
          {-0.75, 2.308263},
          {-1.963525, 1.426585},
          {-2.427051, 0},
          {-1.963526, -1.426585},
          {-0.75, -2.308263},
          {0.75, -2.308263},
          {1.963525, -1.426585},
          {2.427051, 0},
          {3.927051, 0},
          {5.427051, 0},
          {6.927051, 0},
          {3.927051, -3.808263}}},
        // unconnected chains
        {R"(PEPTIDE1{A.C.D.D.E}"HC"|PEPTIDE2{G.C.S.S.S.P.K.K.V.K}"LC"$$$$V2.0)",
         {{0, 0},
          {1.5, 0},
          {3, 0},
          {4.5, 0},
          {6, 0},
          {0, -6.5},
          {1.5, -6.5},
          {3, -6.5},
          {4.5, -6.5},
          {6, -6.5},
          {7.5, -6.5},
          {9, -6.5},
          {10.5, -6.5},
          {12, -6.5},
          {13.5, -6.5}}}}),

    test_data)
{
    std::string helm = test_data.first;
    coords_t monomer_coords = test_data.second;
    auto mol = helm::helm_to_rdkit(helm);
    schrodinger::rdkit_extensions::compute_monomer_mol_coords(*mol);

    auto& conformer = mol->getConformer();
    bool failed = false;
    std::string expected_str = "{";
    std::string actual_str = "{";

    // Helper to convert double to string without unnecessary trailing
    // zeros
    auto fmt = [](double x) -> std::string {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6) << x;
        std::string s = oss.str();
        // Remove trailing zeros
        s.erase(s.find_last_not_of('0') + 1, std::string::npos);
        if (!s.empty() && s.back() == '.')
            s.pop_back(); // remove trailing dot
        return s;
    };

    for (auto monomer : mol->atoms()) {
        auto monomer_idx = monomer->getIdx();
        if (monomer_idx >= monomer_coords.size()) {
            failed = true;
            BOOST_FAIL("Monomer has "
                       << mol->getNumAtoms() << " atoms, but only "
                       << monomer_coords.size() << " sets of expected coords");
            break;
        }
        auto expected_coords = monomer_coords[monomer_idx];
        auto actual_coords = conformer.getAtomPos(monomer_idx);

        expected_str += "{" + fmt(expected_coords[0]) + ", " +
                        fmt(expected_coords[1]) + "}, ";
        actual_str +=
            "{" + fmt(actual_coords[0]) + ", " + fmt(actual_coords[1]) + "}, ";

        if (std::abs(expected_coords[0] - actual_coords[0]) >= 0.01 ||
            std::abs(expected_coords[1] - actual_coords[1]) >= 0.01) {
            failed = true;
        }
    }

    // Remove trailing comma and space
    if (expected_str.size() > 2)
        expected_str.pop_back(), expected_str.pop_back();
    if (actual_str.size() > 2)
        actual_str.pop_back(), actual_str.pop_back();
    expected_str += "}";
    actual_str += "}";

    BOOST_CHECK_MESSAGE(!failed, "Mismatch for HELM: "
                                     << helm
                                     << "\nExpected coords: " << expected_str
                                     << "\nComputed coords: " << actual_str);
}

BOOST_AUTO_TEST_SUITE(TestMonomerCoordgenCheckCoords)

static RDKit::RWMol
make_molecule_with_coords(const std::vector<std::pair<double, double>>& xy,
                          const std::vector<std::pair<int, int>>& bonds = {})
{
    auto mol = RDKit::RWMol();

    // Add atoms
    for (size_t i = 0; i < xy.size(); ++i)
        mol.addAtom(new RDKit::Atom(), true, true);

    // Add bonds
    for (auto [a, b] : bonds)
        mol.addBond(a, b);

    auto conf = new RDKit::Conformer(mol.getNumAtoms());
    for (size_t i = 0; i < xy.size(); ++i)
        conf->setAtomPos(i, RDGeom::Point3D(xy[i].first, xy[i].second, 0.0));

    mol.addConformer(conf);

    return mol;
}

BOOST_AUTO_TEST_CASE(DetectsClashingAtoms)
{

    // Two atoms very close together
    RDKit::ROMol mol = make_molecule_with_coords({{0.0, 0.0}, {0.01, 0.01}});

    BOOST_CHECK(!schrodinger::rdkit_extensions::has_no_clashes(
        mol)); // should detect clash
}

BOOST_AUTO_TEST_CASE(NoClashForSeparatedAtoms)
{
    RDKit::ROMol mol = make_molecule_with_coords({{0.0, 0.0}, {10.0, 10.0}});

    BOOST_CHECK(schrodinger::rdkit_extensions::has_no_clashes(mol)); // no clash
}

BOOST_AUTO_TEST_CASE(IntersectingSegments)
{
    RDGeom::Point3D p1(0, 0, 0), p2(2, 2, 0);
    RDGeom::Point3D q1(0, 2, 0), q2(2, 0, 0);
    BOOST_CHECK(
        schrodinger::rdkit_extensions::segments_intersect(p1, p2, q1, q2));
}

BOOST_AUTO_TEST_CASE(NonIntersectingSegments)
{
    RDGeom::Point3D p1(0, 0, 0), p2(1, 0, 0);
    RDGeom::Point3D q1(0, 1, 0), q2(1, 1, 0);
    BOOST_CHECK(
        !schrodinger::rdkit_extensions::segments_intersect(p1, p2, q1, q2));
}

BOOST_AUTO_TEST_CASE(CollinearTouchingSegments)
{
    RDGeom::Point3D p1(0, 0, 0), p2(2, 0, 0);
    RDGeom::Point3D q1(2, 0, 0), q2(3, 0, 0);
    BOOST_CHECK(schrodinger::rdkit_extensions::segments_intersect(
        p1, p2, q1, q2)); // endpoint touch
}

BOOST_AUTO_TEST_CASE(DistanceParallelSegments)
{
    RDGeom::Point3D p1(0, 0, 0), p2(1, 0, 0);
    RDGeom::Point3D q1(0, 1, 0), q2(1, 1, 0);
    double d = schrodinger::rdkit_extensions::compute_distance_between_segments(
        p1, p2, q1, q2);
    BOOST_CHECK_CLOSE(d, 1.0, 1e-4);
}

BOOST_AUTO_TEST_CASE(DistanceCrossingSegments)
{
    RDGeom::Point3D p1(0, 0, 0), p2(1, 1, 0);
    RDGeom::Point3D q1(0, 1, 0), q2(1, 0, 0);
    double d = schrodinger::rdkit_extensions::compute_distance_between_segments(
        p1, p2, q1, q2);
    BOOST_CHECK_SMALL(d, 1e-6);
}

BOOST_AUTO_TEST_CASE(BondCrossingDetected)
{
    // Make an "X"-shaped molecule with two crossing bonds
    RDKit::ROMol mol = make_molecule_with_coords(
        {
            {0.0, 0.0}, // 0
            {1.0, 1.0}, // 1
            {0.0, 1.0}, // 2
            {1.0, 0.0}  // 3
        },
        {
            {0, 1}, // bond 0–1
            {2, 3}  // bond 2–3 (crosses the first)
        });

    BOOST_CHECK(!schrodinger::rdkit_extensions::has_no_bond_crossings(
        mol)); // bonds intersect → should detect crossing
}

BOOST_AUTO_TEST_CASE(NoBondCrossing)
{
    // Make two parallel bonds that don't intersect
    RDKit::ROMol mol = make_molecule_with_coords(
        {
            {0.0, 0.0}, // 0
            {1.0, 0.0}, // 1
            {0.0, 1.0}, // 2
            {1.0, 1.0}  // 3
        },
        {
            {0, 1}, // horizontal lower
            {2, 3}  // horizontal upper, no intersection
        });

    BOOST_CHECK(schrodinger::rdkit_extensions::has_no_bond_crossings(
        mol)); // no crossing → should pass
}

BOOST_AUTO_TEST_CASE(TouchingBondsNoCross)
{
    // Two bonds that are nearly touching at one endpoint.
    // Will still be considered crossing because their distance is less
    // than the threshold
    RDKit::ROMol mol = make_molecule_with_coords(
        {
            {0.0, 0.0}, // 0
            {1.0, 0.0}, // 1
            {1.1, 0.0}, // 2
            {2.0, 1.0}  // 3
        },
        {
            {0, 1}, // first bond
            {2, 3}  // second bond
        });

    BOOST_CHECK(!schrodinger::rdkit_extensions::has_no_bond_crossings(mol));
}

BOOST_AUTO_TEST_CASE(IsGeometricallyRegularRing2D_RegularPolygon)
{
    // Create a regular hexagon (6-sided polygon), should be detected as
    // geometrically regular
    std::vector<std::pair<double, double>> hex_coords = {
        {1.0, 0.0},  {0.5, std::sqrt(3) / 2},   {-0.5, std::sqrt(3) / 2},
        {-1.0, 0.0}, {-0.5, -std::sqrt(3) / 2}, {0.5, -std::sqrt(3) / 2}};
    RDKit::ROMol hex_mol = make_molecule_with_coords(
        hex_coords, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}});
    std::vector<int> ring_atoms = {0, 1, 2, 3, 4, 5};
    BOOST_CHECK(schrodinger::rdkit_extensions::is_geometrically_regular_ring_2d(
        hex_mol, ring_atoms));

    // distort one vertex. The ring should no longer be considered regular
    hex_coords[0] = {1.5, 0.0};
    RDKit::ROMol distorted_hex_mol = make_molecule_with_coords(
        hex_coords, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}});
    BOOST_CHECK(
        !schrodinger::rdkit_extensions::is_geometrically_regular_ring_2d(
            distorted_hex_mol, ring_atoms));
    // Equilateral triangle, should be regular
    std::vector<std::pair<double, double>> tri_coords = {
        {0.0, 1.0}, {-std::sqrt(3) / 2, -0.5}, {std::sqrt(3) / 2, -0.5}};
    RDKit::ROMol tri_mol =
        make_molecule_with_coords(tri_coords, {{0, 1}, {1, 2}, {2, 0}});
    ring_atoms = {0, 1, 2};
    BOOST_CHECK(schrodinger::rdkit_extensions::is_geometrically_regular_ring_2d(
        tri_mol, ring_atoms));
}

BOOST_AUTO_TEST_CASE(FindAllConnectedMonomersOutsideRing)
{
    // check that we can count the number of monomers attached to each ring
    // monomer of a structure with two large rings.
    std::string helm =
        "PEPTIDE1{F}|PEPTIDE2{[dI].[Trp_Ome].[Asp_OMe].[Cys_Bn].[meG].[Phe_3Cl]"
        ".[dD].T.[dI].T.[dK].[aG].[3Pal].[xiIle].[meD].[Ala_tBu]}|PEPTIDE3{L.["
        "Pro_4Me3OH].S.[NMe2Abz].Q.[3Pal].[xiIle].[D-Hyp].[Ala_tBu].[dI].[Trp_"
        "Ome].[Asp_OMe].N.[meG].[Phe_34diCl].[Phe_34diCl]}$PEPTIDE2,PEPTIDE2,"
        "16:R2-1:R1|PEPTIDE3,PEPTIDE3,16:R2-1:R1|PEPTIDE3,PEPTIDE2,10:R3-1:R3|"
        "PEPTIDE1,PEPTIDE2,1:R2-9:R3$$$V2.0";
    auto mol = helm::helm_to_rdkit(helm);

    schrodinger::rdkit_extensions::compute_full_ring_info(*mol);
    auto rings = mol->getRingInfo()->atomRings();
    BOOST_REQUIRE(rings.size() == 2);

    // number of monomers attached to each ring monomer. This counts all atoms
    // that are encountered when walking the tree starting from each ring atom
    // and avoiding other atoms in the same ring.
    std::vector<std::vector<unsigned int>> ring_attachment_sizes = {
        {0u, 0u, 0u, 0u, 0u, 0u, 0u, 1u, 0u, 0u, 0u, 0u, 0u, 0u, 0u,
         16u}, // first ring
        {0u, 0u, 0u, 0u, 0u, 0u, 0u, 17u, 0u, 0u, 0u, 0u, 0u, 0u, 0u,
         0u}}; // second ring

    for (size_t ring_idx = 0; ring_idx < rings.size(); ++ring_idx) {
        const std::vector<int>& ring = rings[ring_idx];
        BOOST_REQUIRE(ring.size() == ring_attachment_sizes[ring_idx].size());
        for (size_t pos = 0; pos < ring.size(); ++pos) {
            int ring_atom_idx = ring[pos];
            // Find all connected monomers outside the ring.
            auto attached_monomers = schrodinger::rdkit_extensions::
                find_all_connected_monomers_outside_ring(*mol, ring_atom_idx,
                                                         ring);
            BOOST_CHECK_EQUAL(attached_monomers.size(),
                              ring_attachment_sizes[ring_idx][pos]);
        }
    }
}

BOOST_AUTO_TEST_CASE(TestPositionInCoilsScore)
{
    // Single-chain HELM, no branches/rings
    std::string helm = "PEPTIDE1{A.C.D.E.F.G.H.I.A.B.C.D.E.F.G}$$$$";
    auto mol = helm::helm_to_rdkit(helm);

    // Create explicit turns to simulate coiling layout
    std::vector<schrodinger::rdkit_extensions::TurnInfo> turns = {
        {1, 1, true, 0}, {3, 3, false, 0}, {7, 6, true, 0}};

    std::vector<double> scores = {0.5,     1.5,     2.5,     3.25,    3.5,
                                  3.75,    4.5,     5.14286, 5.28571, 5.42857,
                                  5.57143, 5.71429, 5.85714, 6,       6.5};

    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
        const double score =
            schrodinger::rdkit_extensions::get_position_in_coils_double(
                static_cast<int>(i), *mol, turns);
        BOOST_CHECK_CLOSE(score, scores[i], 1e-4);
    }
}

BOOST_AUTO_TEST_SUITE_END()
