

A quantitative model of the switch cycle of an archaeal flagellar motor and
its sensory control, Nutsch et al, Biophys. J. 2005 (
[16192281](http://www.ncbi.nlm.nih.gov/pubmed/16192281) ) and del Rosario et
al, IET Syst. Biol. 2007 (
[17708428](http://www.ncbi.nlm.nih.gov/pubmed/17708428) ). This is the non-
cyclic model for spontaneous simulations used in creating Figure 5C of del
Rosario 2007. The value plotted in the figure is ks*A_43(t)/max(ks*A_43(t)),
where ks is a model parameter. In the figure, the asymmetric model with 10%
and 50% increase in parameter R_cw are compared with data for spontaneous,
repellent and attractant stimuli.

  

There are 5 SBML models provided:

spontaneous simulations (all light parameters Iuv, Ibl and Ior are zero)

repellent dark (sensor via SRII but with Ibl = 0)

repellent light (sensor via SRII)

attractant dark (sensor via SRI but with Iuv=Ior=0)

attractant light (sensor via SRI)

The 5 SBML files are "symmetric" models since the parameters in the clockwise
and counter-clockwise directions are equal. For the asymmetric simulations in
Figure 5c, parameter R_cw must be increased 10% and 50%.

  

We provide the following Matlab code to plot figure 5c using the Systems
Biology Toolbox:

\---- begining of Matlab code

myodeoptions = odeset('AbsTol', 1e-10, 'RelTol', 1e-8);

  

dark1 = linspace(0, 2, 1000);

lighton = linspace(2, 2+0.02, 100);

dark2 = linspace(2+0.02, 60, 1000);

%for spontaneous, no need to separate since all dark

tspan = [dark1,lighton(2:end),dark2(2:end)];

  

sbmodspont = SBmodel('Nutsch2005_phototaxis_noncyc_spont.xml');

sbmodrepdark = SBmodel('Nutsch2005_phototaxis_noncyc_rep_dark.xml');

sbmodreplight = SBmodel('Nutsch2005_phototaxis_noncyc_rep_light.xml');

sbmodattdark = SBmodel('Nutsch2005_phototaxis_noncyc_att_dark.xml');

sbmodattlight = SBmodel('Nutsch2005_phototaxis_noncyc_att_light.xml');

  

R_cw_nominal = SBparameters(sbmodspont, 'R_cw');

R_cw_10inc = R_cw_nominal + 0.1*R_cw_nominal;

R_cw_50inc = R_cw_nominal + 0.5*R_cw_nominal;

  

sbmodspont10 = SBparameters(sbmodspont, 'R_cw', R_cw_10inc);

sbmodspont50 = SBparameters(sbmodspont, 'R_cw', R_cw_50inc);

clear sbmodspont

  

sbmodrepdark10 = SBparameters(sbmodrepdark, 'R_cw', R_cw_10inc);

sbmodreplight10 = SBparameters(sbmodreplight, 'R_cw', R_cw_10inc);

sbmodrepdark50 = SBparameters(sbmodrepdark, 'R_cw', R_cw_50inc);

sbmodreplight50 = SBparameters(sbmodreplight, 'R_cw', R_cw_50inc);

clear sbmodrepdark sbmodreplight

  

sbmodattdark10 = SBparameters(sbmodattdark, 'R_cw', R_cw_10inc);

sbmodattlight10 = SBparameters(sbmodattlight, 'R_cw', R_cw_10inc);

sbmodattdark50 = SBparameters(sbmodattdark, 'R_cw', R_cw_50inc);

sbmodattlight50 = SBparameters(sbmodattlight, 'R_cw', R_cw_50inc);

clear sbmodattdark sbmodattlight

  

%Asymmetric Spontaneous Simulations, 10% increase in parameter R_cw

sboutput_spont_Rcw10 = SBsimulate(sbmodspont10, 'ode15s', tspan, [],
myodeoptions);

%Asymmetric Spontaneous Simulations, 50% increase in parameter R_cw

sboutput_spont_Rcw50 = SBsimulate(sbmodspont50, 'ode15s', tspan, [],
myodeoptions);

  

%Asymmetric Repellent Simulations, 10% increase in parameter R_cw

sboutput_rep_dark1_Rcw10 = SBsimulate(sbmodrepdark10, 'ode15s', dark1, [],
myodeoptions);

initcondafterdark1 = sboutput_rep_dark1_Rcw10.statevalues(end,:);

  

sboutput_rep_light_Rcw10 = SBsimulate(sbmodreplight10, 'ode15s', lighton,
initcondafterdark1, myodeoptions);

initcondafterlight = sboutput_rep_light_Rcw10.statevalues(end,:);

  

sboutput_rep_dark2_Rcw10 = SBsimulate(sbmodrepdark10, 'ode15s', dark2,
initcondafterlight, myodeoptions);

  

%Asymmetric Repellent Simulations, 50% increase in parmaeter R_cw

sboutput_rep_dark1_Rcw50 = SBsimulate(sbmodrepdark50, 'ode15s', dark1, [],
myodeoptions);

initcondafterdark1 = sboutput_rep_dark1_Rcw50.statevalues(end,:);

  

sboutput_rep_light_Rcw50 = SBsimulate(sbmodreplight50, 'ode15s', lighton,
initcondafterdark1, myodeoptions);

initcondafterlight = sboutput_rep_light_Rcw50.statevalues(end,:);

  

sboutput_rep_dark2_Rcw50 = SBsimulate(sbmodrepdark50, 'ode15s', dark2,
initcondafterlight, myodeoptions);

  

%Asymmetric Attractant Simulations, 10% increase in parameter R_cw

sboutput_att_dark1_Rcw10 = SBsimulate(sbmodattdark10, 'ode15s', dark1, [],
myodeoptions);

initcondafterdark1 = sboutput_att_dark1_Rcw10.statevalues(end,:);

  

sboutput_att_light_Rcw10 = SBsimulate(sbmodattlight10, 'ode15s', lighton,
initcondafterdark1, myodeoptions);

initcondafterlight = sboutput_att_light_Rcw10.statevalues(end,:);

  

sboutput_att_dark2_Rcw10 = SBsimulate(sbmodattdark10, 'ode15s', dark2,
initcondafterlight, myodeoptions);

  

%Asymmetric Attractant Simulations, 50% increase in parameter R_cw

sboutput_att_dark1_Rcw50 = SBsimulate(sbmodattdark50, 'ode15s', dark1, [],
myodeoptions);

initcondafterdark1 = sboutput_att_dark1_Rcw50.statevalues(end,:);

  

sboutput_att_light_Rcw50 = SBsimulate(sbmodattlight50, 'ode15s', lighton,
initcondafterdark1, myodeoptions);

initcondafterlight = sboutput_att_light_Rcw50.statevalues(end,:);

  

sboutput_att_dark2_Rcw50 = SBsimulate(sbmodattdark50, 'ode15s', dark2,
initcondafterlight, myodeoptions);

  

A44cwindex = stateindexSB(sbmodspont10, 'A_cw43');

A44ccwindex = stateindexSB(sbmodspont10, 'A_ccw43');

ks_cw = SBparameters(sbmodspont10, 'ks_cw');

ks_cc = SBparameters(sbmodspont10, 'ks_cc');

  

yfig5cspontRcw10 = (sboutput_spont_Rcw10.statevalues(:, A44cwindex)*ks_cw +
sboutput_spont_Rcw10.statevalues(:, A44ccwindex)*ks_cc) / ...

max(sboutput_spont_Rcw10.statevalues(:, A44cwindex)*ks_cw +
sboutput_spont_Rcw10.statevalues(:, A44ccwindex)*ks_cc);

yfig5cspontRcw50 = (sboutput_spont_Rcw50.statevalues(:, A44cwindex)*ks_cw +
sboutput_spont_Rcw50.statevalues(:, A44ccwindex)*ks_cc) / ...

max(sboutput_spont_Rcw50.statevalues(:, A44cwindex)*ks_cw +
sboutput_spont_Rcw50.statevalues(:, A44ccwindex)*ks_cc);

  

A44cwindex = stateindexSB(sbmodrepdark10, 'A_cw43');

A44ccwindex = stateindexSB(sbmodrepdark10, 'A_ccw43');

ks_cw = SBparameters(sbmodrepdark10, 'ks_cw');

ks_cc = SBparameters(sbmodrepdark10, 'ks_cc');

  

tfig5crepRcw10 = [sboutput_rep_dark1_Rcw10.time(:)',
sboutput_rep_light_Rcw10.time(:)', sboutput_rep_dark2_Rcw10.time(:)'];

A44cwrepRcw10 = [sboutput_rep_dark1_Rcw10.statevalues(:, A44cwindex);

sboutput_rep_light_Rcw10.statevalues(:, A44cwindex);

sboutput_rep_dark2_Rcw10.statevalues(:, A44cwindex)];

A44ccwrepRcw10 = [sboutput_rep_dark1_Rcw10.statevalues(:, A44ccwindex);

sboutput_rep_light_Rcw10.statevalues(:, A44ccwindex);

sboutput_rep_dark2_Rcw10.statevalues(:, A44ccwindex)];

yfig5crepRcw10 = (A44cwrepRcw10*ks_cw + A44ccwrepRcw10*ks_cc) / ...

max(A44cwrepRcw10*ks_cw + A44ccwrepRcw10*ks_cc);

  

tfig5crepRcw50 = [sboutput_rep_dark1_Rcw50.time(:)',
sboutput_rep_light_Rcw50.time(:)', sboutput_rep_dark2_Rcw50.time(:)'];

A44cwrepRcw50 = [sboutput_rep_dark1_Rcw50.statevalues(:, A44cwindex);

sboutput_rep_light_Rcw50.statevalues(:, A44cwindex);

sboutput_rep_dark2_Rcw50.statevalues(:, A44cwindex)];

A44ccwrepRcw50 = [sboutput_rep_dark1_Rcw50.statevalues(:, A44ccwindex);

sboutput_rep_light_Rcw50.statevalues(:, A44ccwindex);

sboutput_rep_dark2_Rcw50.statevalues(:, A44ccwindex)];

yfig5crepRcw50 = (A44cwrepRcw50*ks_cw + A44ccwrepRcw50*ks_cc) / ...

max(A44cwrepRcw50*ks_cw + A44ccwrepRcw50*ks_cc);

  

A44cwindex = stateindexSB(sbmodattdark10, 'A_cw43');

A44ccwindex = stateindexSB(sbmodattdark10, 'A_ccw43');

ks_cw = SBparameters(sbmodattdark10, 'ks_cw');

ks_cc = SBparameters(sbmodattdark10, 'ks_cc');

  

tfig5cattRcw10 = [sboutput_att_dark1_Rcw10.time(:)',
sboutput_att_light_Rcw10.time(:)', sboutput_att_dark2_Rcw10.time(:)'];

A44cwattRcw10 = [sboutput_att_dark1_Rcw10.statevalues(:, A44cwindex);

sboutput_att_light_Rcw10.statevalues(:, A44cwindex);

sboutput_att_dark2_Rcw10.statevalues(:, A44cwindex)];

A44ccwattRcw10 = [sboutput_att_dark1_Rcw10.statevalues(:, A44ccwindex);

sboutput_att_light_Rcw10.statevalues(:, A44ccwindex);

sboutput_att_dark2_Rcw10.statevalues(:, A44ccwindex)];

yfig5cattRcw10 = (A44cwattRcw10*ks_cw + A44ccwattRcw10*ks_cc) / ...

max(A44cwattRcw10*ks_cw + A44ccwattRcw10*ks_cc);

  

tfig5cattRcw50 = [sboutput_att_dark1_Rcw50.time(:)',
sboutput_att_light_Rcw50.time(:)', sboutput_att_dark2_Rcw50.time(:)'];

A44cwattRcw50 = [sboutput_att_dark1_Rcw50.statevalues(:, A44cwindex);

sboutput_att_light_Rcw50.statevalues(:, A44cwindex);

sboutput_att_dark2_Rcw50.statevalues(:, A44cwindex)];

A44ccwattRcw50 = [sboutput_att_dark1_Rcw50.statevalues(:, A44ccwindex);

sboutput_att_light_Rcw50.statevalues(:, A44ccwindex);

sboutput_att_dark2_Rcw50.statevalues(:, A44ccwindex)];

yfig5cattRcw50 = (A44cwattRcw50*ks_cw + A44ccwattRcw50*ks_cc) / ...

max(A44cwattRcw50*ks_cw + A44ccwattRcw50*ks_cc);

  

figure

plot(tfig5crepRcw10, yfig5crepRcw10, 'y', 'linewidth', 2)

hold on

plot(tfig5crepRcw50, yfig5crepRcw50, 'k')

legend('R_{cw} increased 10%', 'R_{cw} increased 50%')

plot(sboutput_spont_Rcw10.time, yfig5cspontRcw10, 'y', 'linewidth', 2)

plot(sboutput_spont_Rcw50.time, yfig5cspontRcw50, 'k')

plot(tfig5cattRcw10, yfig5cattRcw10, 'y', 'linewidth', 2)

plot(tfig5cattRcw50, yfig5cattRcw50, 'k')

grid on

myaxis = axis; axis([0 60 0 1.2])

text(1, 1.1, 'repellent'); text(10, 1.1, 'spontaneous'); text(25, 1.1,
'attractant')

xlabel('time, s'); ylabel('reversals per time interval (1/s)')

  

\---- end of Matlab code

This model originates from BioModels Database: A Database of Annotated
Published Models. It is copyright (c) 2005-2011 The BioModels.net Team.  
To the extent possible under law, all copyright and related or neighbouring
rights to this encoded model have been dedicated to the public domain
worldwide. Please refer to [CC0 Public Domain
Dedication](http://creativecommons.org/publicdomain/zero/1.0/) for more
information.

In summary, you are entitled to use this encoded model in absolutely any
manner you deem suitable, verbatim, or with modification, alone or embedded it
in a larger context, redistribute it, commercially or not, in a restricted way
or not..  
  
To cite BioModels Database, please use: [Li C, Donizelli M, Rodriguez N,
Dharuri H, Endler L, Chelliah V, Li L, He E, Henry A, Stefan MI, Snoep JL,
Hucka M, Le Novère N, Laibe C (2010) BioModels Database: An enhanced, curated
and annotated resource for published quantitative kinetic models. BMC Syst
Biol., 4:92.](http://www.ncbi.nlm.nih.gov/pubmed/20587024)

