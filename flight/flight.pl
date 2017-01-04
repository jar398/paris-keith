#!/mpi/tools/bin/perl
use strict;

use KrHomologene;

my %genes=();
my %symb=();
my %syn=();
open(IN,"aliases.fly");
while ($_ = <IN>)
{
    chomp;
    my ($syn,$id,$symb,$desc)=split(/\t/);
    $genes{$id}={'syn'=>{},'id'=>$id,'symb'=>$symb,'desc'=>$desc} unless (ref $genes{$id});
    
    $symb{$symb}=$genes{$id};
    $genes{$id}->{'syn'}->{$syn}=1;
    $syn{$syn}=$genes{$id};
}
open(OUT,">flight.sigs");
print OUT join("\t","FlightHomologs","LOCUSLINK"),"\n";
foreach my $arg(@ARGV)
{
    open(IN,$arg);
    my ($prefix)=split(/\./,$arg);
    my ($pmid,$title)=();
    my %sig=();
    my %humSig=();
    while ($_ = <IN>)
    {
	chomp;
	$pmid=$1 if (/^PMID: ?([0-9]+)/);
	$title=$1 if (/^SCREEN: ?(.*)/);
	if (/(FBgn[0-9]+|HDC[0-9]+)/)
	{
	    my $code=$1; my $origCode=$code;
	    $code=~tr/A-Z/a-z/;
	    my $methCode=undef;
	    my ($row,$cellLine,$id,$fbgn,$symb)=split(/\t/);
	    $symb=~s/ *$//;

	    my $id=undef;
	    if (ref $syn{$code})
	    {
	       $id=$syn{$code}->{'id'}; $methCode="C";
	    }
	    elsif (ref $symb{$symb})
	    {
		$id=$symb{$symb}->{'id'}; $methCode="S";
	    }
	    my @humanLl=();
	    unless (defined $id)
	    {
		$id="NoID"; 
	    }
	    else
	    {
		@humanLl=KrHomologene::getHomologs($id,"human");
		foreach my $ll(@humanLl)
		{
		    $humSig{$ll}++;
		}
	    }
	    my $humCount=scalar(@humanLl);
	    print join("\t",$prefix,$pmid,$origCode,$symb,$id,$methCode,$humCount,join(',',sort {$a<=>$b} @humanLl),$title),"\n";
	    
	}
    }
    if (scalar(keys %humSig)>5)
    {
	print OUT join("\t","Flight-$prefix","Homologs of Drosphila RNAi screen hits; $title (PMID:$pmid)",sort {$a<=>$b} keys %humSig),"\n";
    }
}

