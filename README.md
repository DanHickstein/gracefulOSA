# gracefulOSA

Software for collecting spectra using a Yokogowa/Ando spectrum analyzer (OSA) while rotating a waveplate using a Thorlabs K10CR1. So far, this has been tested with the following spectrum analyzers:

ANDO AQ-6315A
Yokogawa AQ6375 (set to AQ6317 compatibility mode)

This should work with any of the modern Yokogawa OSAs that have the AQ6317 legacy compatiblity mode. To activate this mode press System >> More >> GPIB Settings >> Command Format >> AQ6317.

This software assumes GPIB communication. It is necessary to manually set the GPIB address in the script to match the optical spectrum analyzer. 
