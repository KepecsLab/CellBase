function [all_spikes, all_waves]  = TrodeLoadingEngine_NLX_1eMinus4(fname)
% if nlx ntt file is loaded in units of 10^-4 s (see David Redish's C++
% file)
[all_spikes, all_waves] = LoadTT_NeuralynxNT(fname);
all_spikes= all_spikes*1e-4;