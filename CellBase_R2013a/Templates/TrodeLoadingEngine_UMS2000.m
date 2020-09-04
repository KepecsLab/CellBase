function [all_spikes, all_waves]  = TrodeLoadingEngine_UMS2000(fname)

[all_spikes, all_waves] = LoadTT_UMS2000(fname);
all_spikes= all_spikes*1;