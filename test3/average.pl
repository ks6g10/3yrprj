#!/usr/bin/perl
my %val = ();
my %count = ();
my %max = ();
my %min = ();
my $total = 0;
my @aa = ();
my %aaa = (); 
open(FILE,"output");
$min{"1"} = 999999999;
$min{"2"} = 999999999;
$min{"3"} = 999999999;
$min{"4"} = 999999999;
$min{"5"} = 999999999;
$min{"6"} = 999999999;
$min{"7"} = 999999999;
$min{"8"} = 999999999;
$min{"9"} = 999999999;
$min{"10"} = 999999999;
while(<FILE>) {
    if($_ =~ /l([0-9]*)..([0-9]*)/g) {
	if($1) {
	    push( @{ $aaa{ $1 } }, int($2) );
	    #$val{$1} += $2;	    
	    #$count{$1} += 1;
	    if($2 > $max{$1}) {
		$max{$1} =$2; 
	    }
	    if($2 < $min{$1}) {
		$min{$1} =$2; 
	    }
	}
    }
}
close(FILE);
my @sorted;# = sort { $a->[1] <=> $b->[1] } @aa;
print "x\t min\t lq\t median\t uq\t max\n";
for(sort keys %aaa) {
    my $len = @{$aaa{$_}};
#    print $_."\t";
#    $aaa{$_} = sort(@{$aaa{$_}});
    my @array = sort {$a <=> $b} @{$aaa{$_}};
    $aaa{$_} = @array;
    my $means = meany(@array);
    my $medians = median(@array);
    my $max = max(@array);
    my $min = min(@array);
    my $lq = lowerqt(@array);
    my $uq = upperqt(@array);
    print "$_\t $min\t $lq\t $medians\t $uq\t $max\n";
#    print "mean\t $means \t median\t $medians\t max\t $max\t min\t $min\t uq\t $uq\t lq\t $lq\n";
#    push(@sorted,@tmpsort);
}


for (keys %val) {
    $val{$_} = $val{$_}/$count{$_};
    $total += $val{$_};
}

foreach $value (sort keys %val)
{
     print "$value\t$val{$value} \tpercentage ". 100*($val{$value}/$total)." \tmin\t$min{$value}\tmax\t$max{$value}\n";
}

sub lowerqt {
    my @vals = @_;
    my $len = @vals;
    my $l = $len*(25/100);
    if($l != int($i)) {
	return $vals[int($l)];
    } else {
	return int(($vals[int($l)]+$vals[int($l)+1])/2);

    }

}


sub upperqt {
    my @vals = @_;
    my $len = @vals;
    my $l = int($len*(75/100));
    return $vals[$l];
}

sub meany {
    my @vals = @_;
    my $count  = @vals;
    my $average = 0;
    for (@vals) {
	$average += $_;
    }

    return int($average/$count);
}

sub max {
    my @vals = @_; 
    my $len = @vals;
    return ($vals[$len-1]);
}

sub min {
    my @vals = @_; 
    return ($vals[0]);
}

sub median
{
    my @vals = @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return int(($vals[int($len/2)-1] + $vals[int($len/2)])/2);
    }
}
