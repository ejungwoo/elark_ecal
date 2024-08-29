RUN=113 root -q -b -l make_histograms.C
RUN=113 root -q -b -l run_energy_correction.C
RUN=113 root -q -b -l run_slope_correction.C
RUN=113 root -q -b -l make_c1_histograms.C
RUN=113 root -q -b -l run_ballistic_correction.C
RUN=113 root -q -b -l make_c2_histograms.C

RUN=303 root -q -b -l make_histograms.C
RUN=303 root -q -b -l run_energy_correction.C
RUN=303 root -q -b -l run_slope_correction.C
RUN=303 root -q -b -l make_c1_histograms.C
RUN=303 root -q -b -l run_ballistic_correction.C
RUN=303 root -q -b -l make_c2_histograms.C

#root summary.C
