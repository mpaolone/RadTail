#This script will run raveshift and convert the output for use in the calcRadTail.C script.  Here you must set the beam-energy and angle equal to the same values used in the calcRadTail.C file.  See that file for more details.

my $beames = 400.0;
my $ang = 90.0;

my $model_opt = 0; #0 for Sum of Gaussians, 1 for Forrier-Bessel version-A, or 2 for Forrier-Bessel version-B.

my $minAng = 5.0;
my $maxAng = $ang + 5.0;

open(RAVIN,">ravin.txt");
       if($model_opt == 0){
               print RAVIN "12C  SOG Model
               6,12
               300,23
               1.20
               0.0, .016690
               0.4, .05032
               1.0, .128621
               1.3, .180515
               1.7, .219097
               2.3, .278416
               2.7, .058779
               3.5  .057817
               4.3, .007739
               5.4, .002001
               6.7, .000007
               1,1
               $beames,$minAng,$maxAng,0.05,0,0,0,
               ,,,,,,,,
               ~    ";
       }elsif($model_opt == 1){
               print RAVIN "12C  FBS  Model
               6,12
               100,18
               8.0
               1.5737e-02
               3.8897e-02
               3.7085e-02
               1.4795e-02
               -4.4831e-03
               -1.0057e-02
               -6.8696e-03
               2.8813e-03
               -7.7229e-04
               6.6908e-05
               1.0636e-04
               -3.6864e-05
               -5.0135e-06
               9.4550e-06
               -4.7687e-06
               0.0000e+00
               0.0000e+00
               1,1
               $beames,$minAng,$maxAng,0.05,0,0,0,
               ,,,,,,,,
               ~    ";
       }elsif($model_opt == 2){
               print RAVIN "12C  FBS  Model
               6,12
               100,18
               8.0
               1.5721e-02
               3.8732e-02
               3.6808-02
               1.4671-02
               -4.43277-03
               -9.7752-03
               -6.8908e-03
               2.7631e-03
               -6.3568e-04
               7.1809e-05
               1.8441e-04
               7.5066-05
               5.1069e-06
               1.4308e-05
               2.3170e-06
               6.8465e-07
               0.0000e+00
               1,1
               $beames,$minAng,$maxAng,0.05,0,0,0,
               ,,,,,,,,
               ~    ";
       }

       close(RAVIN);
	   
 system("./raveshift ravin.txt >& trash");
 #cleanup
 #system("rm ravin.txt");
 #system("rm trash");
 
 open(OUTF,">inputFF2.dat");
 open(RAVOUT,"output.dat");

 my $startRead = 0;

 while(<RAVOUT>){
	 chomp();
	 if($_ eq "       ----------------------------------------------------------------"){
		 if($startRead){
			 $startRead = 0;
		 }else{
			 $startRead = 1;
			 $n++;
		 }
	 }
	 if($startRead){
		 if($_ ne "       ----------------------------------------------------------------"){
			 my @a = split(/\s+/);
			 print OUTF sprintf($a[2]."  ".$a[6]."\n");
		 }
	 }
 }
 close(RAVOUT);
 system("rm output.dat");
 close(OUTF);
 

