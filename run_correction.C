#include "run_slope_correction.h"
#include "run_ballistic_correction.h"
#include "run_energy_calibration.h"

void run_correction(int iCorrection, bool runViewer=true)
{
    if (iCorrection==0) run_energy_calibration(runViewer);
    if (iCorrection==1) run_slope_correction(runViewer);
    if (iCorrection==2) run_ballistic_correction(runViewer);
}
