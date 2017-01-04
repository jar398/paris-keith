use strict;
package CellCycle;

# use for CellCycle lmy $dataFileHome="./data.";
$dataFileHome=$ENV{'KRMLNM'/data} if (defined $ENV{'KRMLNM'});


my %llToCycle=();
my %cycleToLl=('_any'=>{});

my $init=0;
sub _init
{
    return if ($init);
    $init=1;
    open(IN,"$dataFileHome/StanfordCellCycle.cat");
    $_=<IN>;
    while ($_ = <IN>)
    {
	chomp;
	my @lls=split(/\t/,$_);
	my $label=shift(@lls);
	my ($phase)=($label=~/: ([^ ]+)/);
	$cycleToLl{$phase}={};
	foreach my $ll(@lls)
	{
	    $cycleToLl{$phase}->{$ll}=1;
	    $cycleToLl{'_any'}->{$ll}=1;
	    $llToCycle{$ll}=$phase;
	}
    }
}
sub phases
{
    &_init();
    return grep(/^[A-Z]/,keys %cycleToLl);
}
sub llForPhase
{
    my ($phase)=@_;
    &_init();
    return keys %{$cycleToLl{$phase}};
}
sub cycleForLl
{
    my ($ll)=@_;
    &_init();
    return $llToCycle{$ll};
}
sub getCycling
{
    my @lls=@_;
    &_init();
    my %positives=();
    foreach my $ll(@lls)
    {
	next unless (defined $llToCycle{$ll});
	$positives{$ll}=$llToCycle{$ll};
    }
    return \%positives;
}
1;



