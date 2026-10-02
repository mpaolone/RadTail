# RadTail
calculation of the radiative tail

1) build raveshift with "make"
2) edit the runRaveshift.pl to include the beam energy and angle you want to calculate the radiate tail.
3) run the script with "perl runRaveshift.pl" to create the needed inputFF2.dat file
4) edit the calcRadTail.C script.  The top portion of the script lists all the values you should change to match your desired calculation
4) run the script with root