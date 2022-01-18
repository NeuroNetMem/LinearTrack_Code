{\rtf1\ansi\ansicpg1252\cocoartf2636
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 LINEAR TRACK WORKFLOW: \
\
3 Directories are assumed: \
- Raw data\
- GLM Location \
- Processed Results (Dropbox, GitHub or whatever)\
\
\
0. COMPUTE GAMMA COUPLINGS 
In code_gammas 
- Basins and layer-resolved theta-gamma coupling AND layer-resolved gamma power 
- Gamma coupling to animal speed and locations on the track 


1. COMPUTE GLM\
In code_glm: \
- Compute_PhasePosition_GLM_1D\
\
2 COMPUTE CELL PROPERTIES\
In code_properties:\
- Move_Files\
- Compute_ThetaScores\
- Compute_PlaceFields\
- Compute_CellProperties \
- Define_ThetaLimits \
\
3 SPIKE TIMING MEASURES \
In code_spikes: \
- BurstAnalysis %Computes phase-spiking properties \
- CrossAnalysis % Computes pairwise-spike coordination \
- SingleTrial_PP % Computes measure for single trial phase precession \
- SpikeGamma_Analysis % Probability of spikes conditioned on dominating gamma \
- GammaPhase_Analysis % Probability of gamma dominating conditioned on theta phase \
\
4 PROCESS RESULTS AND PLOTS \
in code_plotting (should depend only on data available in Processed Results folder)\
\
- At the moment follows Guardamagna 2022 structure \
\
\
\
\
}