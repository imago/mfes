#print "Hello world from vol2stl!\n";

$num_args = $#ARGV + 1;
if ($num_args != 1) {
  print "\nUsage: perl vol2stl.pl <NETGEN-volFile>\n";
  exit;
}

$filePath=$ARGV[0];
  
my @tri_idx;
my @points;
my $nrSurfaces = 0;
my $nrPoints = 0;
open (MYFILE, $filePath);
 while ($line = <MYFILE>) {
  #   print substr($line, 0, 7)."\n";
    if (substr($line, 0, 7) eq "surface"){
	$nrSurfaces = <MYFILE>;
#	print "#surfaces: $nrSurfaces\n";
	for ($i = 0; $i < $nrSurfaces; $i++){
	    $nextline = <MYFILE>;
	    @tokens = split(" ", $nextline);
	    push(@tri_idx, $tokens[5]);
	    push(@tri_idx, $tokens[6]);
	    push(@tri_idx, $tokens[7]);
#	    print "$tokens[5] $tokens[6] $tokens[7]\n";
	}
    }

    if (substr($line, 0, 6) eq "points"){
	$nrPoints = <MYFILE>;
#	print "#points: $nrPoints\n";
	for ($i = 0; $i < $nrPoints; $i++){
	    $nextline = <MYFILE>;
	    @tokens = split(" ", $nextline);
	    push(@points, $tokens[0]);
	    push(@points, $tokens[1]);
	    push(@points, $tokens[2]);
#	    print "$tokens[0] $tokens[1] $tokens[2]\n";
	}
    }

 }

print "solid object\n";
for ($i = 0; $i < $nrSurfaces; $i++){
    $idx = $tri_idx[($i*3)]-1;
    $idy = $tri_idx[($i*3)+1]-1;
    $idz = $tri_idx[($i*3)+2]-1;

    $v1x = $points[($idx*3)];
    $v1y = $points[($idx*3)+1];
    $v1z = $points[($idx*3)+2];

    $v2x = $points[($idy*3)];
    $v2y = $points[($idy*3)+1];
    $v2z = $points[($idy*3)+2];

    $v3x = $points[($idz*3)];
    $v3y = $points[($idz*3)+1];
    $v3z = $points[($idz*3)+2];

    print "facet normal 0 0 0\n";
    print "outer loop\n";
    print "    vertex $v1x $v1y $v1z\n";
    print "    vertex $v2x $v2y $v2z\n";
    print "    vertex $v3x $v3y $v3z\n";
    print "endloop\n";
    print "endfacet\n";
}

print "endsolid object\n"
