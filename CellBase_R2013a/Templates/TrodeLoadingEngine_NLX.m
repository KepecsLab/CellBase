function [all_spikes, all_waves]  = TrodeLoadingEngine_NLX(fname)

[all_spikes, all_waves] = LoadTT_NeuralynxNT(fname);
all_spikes= all_spikes*10^-4;