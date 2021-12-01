# Beckr-Hickl_SPC_Depth_px
Takes SDT files from a B&amp;H SPC card, X-correlates or fits and produces depth images.

Requires the use of sdts where the histograms are in pixel order - i.e using the 'FiFo image' mode or 'scan sync in' mode.

To use:<br />
- Select sdt using top button.<br />
- If gating required, use the two boxes to enter your values. Leave the end gate at '0' for no gating. Recommended to use no gating intially to see histogram.<br />
- If no external IRF sdt file is available select the 'Automatically determine IRF' box to take the brightest pixel (highest peak) and use that as the IRF for the x-correlation.<br />
- 'Subtract back reflections' - NOT TESTED - takes an external sdt file where the histogram is with no target returns. Will then subtract this from the histograms prior to the processing. May deal with back reflections/noise in the system.<br />
- 'Fit instead of X-corr' - Will take each histogram and fit an exponentially modified gaussian to the data. SLOW! Don't use for big images.<br />
- 'Correct aspect in graph?' - Takes the scene size scanned and corrects the aspect of the image. <br />
- 'Largest value cutoff?' - If there are odd large values in the plot can cut them off <br />
- 'Remove empty pixels?' - If the peak returns in a pixel are just noise (i.e ~1 for SNSPD) then the fitting/X-Corr won't work so the pixels are left blank if selected.<br />
- 'Bin the data?' - When using small bin sizes (i.e 200/400fs) we can gain SNR by binning the data at the expense of resolution. Note that if you enter a number that the bins cannot be equally divided into it won't work. Use 2/4/8/etc. <br />
- 'Plot it' - Plots the data. Will return a histogram of the brightest pixel to allow tuning of the gating, a 2d image of the peak bins, a 3d image of the peak bins and a pixel counts map.<br />

After selection of parameters and 'Plot it' is selected the brightest histogram will be plotted with a fitted EMG so you can review the gating etc before committing to the full processing.<br />

Also in file:<br />
- A cmd_line version