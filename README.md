# Beckr-Hickl_SPC_Depth_px
Takes SDT files from a B&amp;H SPC card, X-correlates and produces depth images.

Requires the use sdts where the histograms are in pixel order - i.e using the 'FiFo image' mode or 'scan sync in' mode.

To use:<br />
- Select sdt using top button.<br />
- If gating required, use the two boxes to enter your values. Leave the end gate at '0' for no gating. Recommended to use no gating intially to see histogram.<br />
- If no external IRF sdt file is available select the 'Automatically determine IRF' box to take the brightest pixel (highest peak) and use that as the IRF for the x-correlation.<br />
- 'Subtract back reflections' - NOT TESTED - takes an external sdt file where the histogram is with no target returns. Will then subtract this from the histograms prior to the processing. May deal with back reflections/noise in the system.<br />
- 'Plot it' - Plots the data. Will return a histogram of the brightest pixel to allow tuning of the gating, a 2d image of the peak bins, a 3d image of the peak bins and a pixel counts map.<br />

Also in file:<br />
- A cmd_line version