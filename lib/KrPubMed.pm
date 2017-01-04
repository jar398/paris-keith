#!/mpi/tools/bin/perl

use strict;
use LWP::Simple;
package KrPubMedEntry;

my $defaultPm=new KrPubMed("KrPubMed");

sub _getWin
{
    my ($string,$pos)=@_;
    my $winSt=$pos-20;
    
    my $win=substr($string,$winSt,41);
    return $win;
}

sub new
{
    my ($type,$pmid,$text)=@_;
    my $this={'pmid'=>$pmid, 'AuthorList'=>[] };
    bless $this,$type;
    my @lines=();
    if (ref $text)
    {
	@lines=@{$text};
    }
    else
    {
	@lines=split(/\n/,$text);
    }
    if (scalar(@lines)==0)
    {
	my $hash=$defaultPm->cachedFetch($pmid);
	$this->{'lines'}=$hash->{$pmid};
    }
    else
    {
	$this->{'lines'}=\@lines;
    }
    $this->_parse();
    return $this;
}

sub guessNonhumanSpecies
{
    my ($this)=@_;
    my %speciesCand=();
    foreach my $key('Title','Abstract')
    {
	my $text=$this->{$key};
	$text=~s/yeast (one|two|three).hybrid//g;
	foreach my $cand ($text=~/ (human|aplysia|dictyostelium|rice|mouse|mice|murine|arabidopsis|candida|drosophila|zebrafish|[a-z]+bacillus|strepto[a-z]+|physarum|xenopus|[a-z]?ovine|escherichia|[a-z]+coccus|chloroplast|[a-z]+bacter|soybean|tomato|neurospora|maize|aspergillus|caenorhabditis|nematode|virus|hiv|c. ?elegans|rats?|porcine|pigs|limulus|subtilis|bacillus)[\.\, ]/ig)
	{
	    $cand="mouse" if ($cand=~/mouse|murine|mice/i);
	    $speciesCand{$cand}=1;
	}
	if ($text=~/yeast|saccharomyces|pombe|schizosaccharomyces/i)
	{	    
	    my @yeasts=();
	    if ($text=~/fission yeast|pombe|schizosaccharomyces/i)
	    {
		push(@yeasts,"pombe");
	    }
	    if ($text=~/[^a-z]saccharomyces|baker.s yeast|budding yeast|cerevisiae/i)
	    {
		push(@yeasts,"cerevisiae");
	    }
	    push(@yeasts,"yeast") unless (scalar(@yeasts)>0);
	    foreach my $yeast(@yeasts)
	    {
		$speciesCand{$yeast}=1;
	    }
	}
	last if (scalar(keys %speciesCand)>0);
    }
    return \%speciesCand;
}

sub candPhosphoSites
{
    my ($this)=@_;
    my %where=();
    my %aaList=();    
    my %posList=();
    foreach my $key('Abstract','Title')
    {
	my $field=$this->{$key};
	$field=~s/(Ser|Thr|Tyr)-([0-9]+) and -/$1-$2 and $1-/;
	$field=~s/T[47][\- ](PNK|poly|[RD]NA)//g;
	$field=~s/T7-[a-z]+//;
	$field=~s/T84 .?(colon|cell l[a-z]+|human|monolayer|epithelial)//;
	$field=~s/(human|polarized|lines?),? T84//;
	$field=~s/[pP]hage T[47]//g;
	$field=~s/T47-D//;
	$field=~s/T[0-9]+ phage//g;
	$field=~s/S[0-9]+ cell//g;
	$field=~s/S100[^A-CE-Z]//;
	$field=~s/3T3//g;
	$field=~s/.ibosomal S6 //g;
	$field=~s/.ibosomal protein S6 //g;
	$field=~s/.yclin T[12]//;
	$field=~s/S6 kinase//g;
	$field=~s/S6[Kk]//g;
	$field=~s/RPS6//ig;
	$field=~s/LY[0-9]+//g;
	$field=~s/[STY][0-9]+ cells//g;
	$field=~s/Y27632//g;
	while ( $field=~/(ser|serine|tyr|tyrosine|thr|threonine)s?([- \(]|residues? ?)?([0-9][ ,\)\/0-9andor]+)/ig
		##  ||$this->{$key}=~/[STY][0-9]+/
	      )
	{
	    my ($aaLong,$posList)=($1,$3);
	    my $stringPos=pos $field;
	    
	    
	    foreach my $pos(split(/[^0-9]/,$posList))
	    {
		next if ($pos<1);
		$aaLong=~s/^tyr/yr/i;
		my $aa=substr($aaLong,0,1); $aa=~tr/a-z/A-Z/;
		$where{$key}++;
		$aaList{$aa}={} unless (ref $aaList{$aa});
		$aaList{$aa}->{$pos}=&_getWin($field,$stringPos);
	    }
	}
	while ($field=~/[^A-Za-z0-9]([STY])([0-9]+)[^A-Za-z0-9]/g)
	{
	    my ($aa,$pos)=($1,$2);
	    next if ($pos<1);
	    my $stringPos=pos $field;
	    pos $field=$stringPos-length($aa.$pos);
	    $where{$key}++;
	    $aaList{$aa}={} unless (ref $aaList{$aa});
	    $aaList{$aa}->{$pos}=&_getWin($field,$stringPos);
	}
    }
    foreach my $aa(keys %aaList)
    {
	foreach my $pos(keys %{$aaList{$aa}})
	{
	    $posList{$pos}={} unless (ref $posList{$pos});
	    $posList{$pos}->{$aa}=$aaList{$aa}->{$pos};
	}
    }
    return (\%aaList,\%where,\%posList);
}

my $matchTerms=0;


sub _initializeMatchTerms
{
    return if (ref $matchTerms);
    $matchTerms={};
    foreach my $row(
		    [1000,"immunoprecip","modif","subunit","ligand","associated","intermolecular",
		     "activat","repress","regulat"],
		    [2000,"bind","phosphorylat","complex","ligand",
		     "affinity","associate[ds]","effector","recruit",
		     "modif","cleavage","target","hetero[a-z]{1,5}imer",
		     "adaptor","interact","partner","inhibit","scaffold","tether","co.?(activat|repress)"
		    ],
		    [3000,"ternary","(di|tri|tetra|penta|hexa)mer"],
		    [5000,
		     "two.hybrid","co.?immunoprecip","ternary",
		     "hetero[a-z]{1,5}imer","multi-?(protein|subunit)"		     
		    ])
    {
	my @terms=@{$row};
	my $score=shift(@terms);
	foreach my $term(@terms)
	{
	    $matchTerms->{$term}=$score;
	}
    }
}
sub getYear
{
    my ($this)=@_;
    return  $this->{'Year'};
}
sub getJournal
{
    my ($this)=@_;
    return  $this->{'Source'};
}

sub getPmid
{
    my ($this)=@_;
    return $this->{'pmid'};
}
sub getTitle
{
    my ($this)=@_;
    return $this->{'Title'};
}
sub getAbstract
{
    my ($this)=@_;
    return $this->{'Abstract'};
}
sub interactScore
{
    my ($this)=@_;    
    unless (defined $this->{'scored'})
    {
	&_initializeMatchTerms();
	$this->{'scores'}={} unless (ref $this->{'scores'});
	my $scores=$this->{'scores'};
	
	$scores->{'misc'}-=1000 unless ($this->{'Language'}=~/eng/i);
	$scores->{'misc'}-=5000 unless (length($this->{'Abstract'})>10);
	$scores->{'misc'}-=10000 if ($this->{'Title'}=~/Prediction of the coding sequences of unidentified human genes/);
	$scores->{'misc'}-=5000 if ($this->{'Title'}=~/pseudogene|haplotype|mutation|auto.?(antigen|antibod)/i);
	

	# term-based positive scores
	# high-scoring
	foreach my $term(keys %{$matchTerms})
	{
	    
	    foreach my $key("Title","Abstract")
	    {
		if ($this->{$key}=~/$term([a-z]*) ([^ ]+)/i)
		{
		    my $follow=$1;
		    my $inc=$matchTerms->{$term};
		    $inc*=2 if ($follow=~/^(with|to|by)/);
		    if ($key eq 'Title')
		    {
			$inc*=5;
			$inc*=3 if ($scores->{'queryTitle'}>0);
		    }
		    $scores->{'interact'}+=$inc;
		}
	    }
	}
	# lower scoring
	foreach my $term("signal","ubiquitin","proteasom",)
	{
	    foreach my $key("Title","Abstract")
	    {
		if ($this->{$key}=~/$term/)
		{
		    $scores->{'misc'}+=300;
		}
	    }
	}

	# organism scores
	
	foreach my $key("Title","Abstract")
	{
	    
	    my ($n,$p)=();
	    my $oScore=1000; $oScore*=4 if ($key eq 'Title');
	    if ($this->{$key}=~/yeast two.hybrid/) {}
	    elsif ($this->{$key}=~/bacteri|plasmodium|trypanos|candida|zebrafish|drosophila|xenopus|laevis|melanogaster|maize|thaliana|elegans|yeast|escherichia|cerevisae|yeast|saccharomyces|coli|pombe|aspergillus|arabidopsis|flower|maize|floral|chicken|bacteriophage|chloroplast|hepatitis [a-z]+ virus/i)
	    {
		$n=1;
	    }
	    
	    if ($this->{$key}=~/human|mouse|mammal|bovine|[^a-z]rat[^a-z]/i)
	    {
		$p=1;
	    }
	    if ($p)
	    {
		$scores->{'organism'}+=$oScore;
	}
	    else
	    {
		$scores->{'organism'}-=$oScore;
	    }
	}
	if ($this->{'Title'}=~/prediction of the coding/i)
	{
	    $scores->{'interact'}-=3000;
	}
	if (
	    ($this->{'Title'}=~/gene|chromosom/i && $this->{'Title'}=~/mapping|location/i) ||
	    ($this->{'Title'}=~/histocompatibility|genomic|translocation|polymorphism|deletion|insertin|chromosom [0-9XY]/i))
	{
	    $scores->{'interact'}-=5000;
	}
	my $score=0;
	foreach my $value(values %{$scores})
	{
	    $score+=$value;
	}
	$this->{'scored'}=$score;
    }
    return $this->{'scored'};
}

sub batchNew
{
    my @pmids=@_;
    shift(@pmids) if (ref $pmids[0]);
    my %entries=();
    my $textHash=$defaultPm->cachedFetch(@pmids);
    foreach my $pmid(keys %{$textHash})
    {
	$entries{$pmid}=new KrPubMedEntry($pmid,$textHash->{$pmid});;
    }
    return \%entries;
}

sub _parse
{
    my ($this)=@_;
    my $mode="";
    my $lastTag="";
    unless (ref $this->{'lines'})
    {
	die "Bad parse for $this->{'pmid'}\n";
    }
    foreach my $line(@{$this->{'lines'}})
    {
	if ($line=~/\<([A-Z0-9]+)/i)
	{
   	    $lastTag=$1;
	}
	if ($mode eq "Journal")
	{
	    if ($line=~/<(Year|Volume)>([0-9]+)/)
	    {
		$this->{$1}=$2 if ($2>0);
	    }
	    if ($line=~/<\/Journal>/)
	    {
		$mode="";
	    }
	    if ($line=~/<MedlineDate>([0-9]{4})/)
	    {
		$this->{'Year'}=$1 unless ($this->{'Year'}>0);
	    }
	}
	elsif ($mode eq 'Author')
	{
	    if ($line=~/<LastName>(.*)<\/LastName>/)
	    {
		push(@{$this->{'AuthorList'}},$1);
	    }
	    $mode="" if ($line=~/<\/AuthorList>/);	    
	}
	else
	{
	    if ($line=~/<Journal>/)
	    {
		$mode="Journal";
	    }
	    elsif ($line=~/<Language>([a-z]+)/)
	    {
		$this->{'Language'}=$1;
	    }
	    elsif ($line=~/<AuthorList/)
	    {
		$this->{'AuthorList'}=[];
		$mode="Author";
	    }
	    elsif ($line=~/<Med/)
	    {
		if ($line=~/<MedlineTA>(.*)<\/MedlineTA>/)
		{
		    $this->{'Source'}=$1;
		    if ($this->{'Source'}=~/Rev|Curr Op|Bioessays/)
		    {
			$this->{'Type'}='Review?' unless (defined $this->{'Type'});
		    }
		}
		if ($line=~/<MedlinePgn>(.*)<\/MedlinePgn>/)
		{
		    $this->{'Pages'}=$1;
		}
            }
	    elsif ($line=~/<AbstractText>(.*)<\/AbstractText/)
	    {
		$this->{'Abstract'}=$1;
		if ($this->{'Abstract'}=~/(i|current|we|this|short|will|present) review/i)
		{
		    $this->{'Type'}='Review?';
		}
	    }
	    elsif ($line=~/<PublicationType/ && /Review/)
	    {
		$this->{'Type'}="Review";
	    }
	    elsif ($line=~/<ArticleTitle>(.*)<\/ArticleTitle/)
	    {
		$this->{'Title'}=$1;
		$this->{'Type'}="Review?" if ($this->{'Title'}=~/\(Review\)/);
	    }
	}
    }
    $this->{'Author1'}=$this->{'AuthorList'}->[0];
}

package KrPubMed;

my $spliceSize=120;
my $maxSplice=int($spliceSize*1.2);

my $eUtilBase="http://eutils.ncbi.nlm.nih.gov/entrez/eutils";


sub batchNew
{
my @pmids=@_;
  return KrPubMedEntry::batchNew(@pmids);
}

sub new
{
    # it is recommended that each tool have a distinct name!
    my $type=shift;
    my $this={};
    bless $this,$type;
    $this->{'tool'}=shift;
    $this->{'tool'}="KrPubMed" unless (defined $this->{'tool'});
    $this->{'user'}=$ENV{'USER'}; $this->{'user'}="millennium" unless (defined $this->{'user'});
    $this->{'idtag'}="\&email=$this->{'user'}\@mpi.com\&tool=$this->{'tool'}";

    $defaultPm=$this if (!ref($defaultPm) || $defaultPm->{'idtag'} eq "KrPubMed");
    return $this;
}

sub fetchByPmid
{
    my ($this,$pmid)=@_;
    return $this->_queryNcbi("efetch","&id=$pmid");
}

# mindate= maxdate=  must be paired. Dates can be years
# reldate=90  search last 90 days

sub getEntry
{
    my ($this,$pmid)=@_;
    my $text=$this->fetchByPmid($pmid);
    return new KrPubMedEntry($pmid,$text);
}
sub pmidsForQuery
{
    # requires properly formatted query
    my ($this,$query,$pmidExclude)=@_;
    $pmidExclude={} unless (ref $pmidExclude);
    my $result=$this->_queryNcbi("esearch",$query);
    open(PFQ,">krpubmed.pmidsforquery.out");
    print PFQ $query,"\n\n";
    print PFQ $result;
    close PFQ;
    my $inIds=();
    my @lines=split(/\n/,$result);
    my @pmids=();
    foreach my $line(@lines)
    {
	if ($line=~/\<IdList\>/)
	{
	    $inIds=1;
	}
	elsif ($line=~/\<.IdList/)
	{
	    $inIds=0;
	}
	if ($inIds && 
	    $line=~/\<Id.([0-9]+)/)
	{
	    my $pmid=$1;
	    push(@pmids,$pmid) unless (defined $pmidExclude->{$pmid});
	}
    }
    return @pmids;
}

sub _queryNcbi
{
    my ($this,$action,$query)=@_;

    my $time=time();
    if ($this->{'lastTime'}>0)
    {
	# code to obey NCBI guidelines
	my $delta=$time-$this->{'lastTime'};
	if ($delta<3)
	{
	    my $s=3-$delta;
	    print STDERR "# sleeping $s s before querying $action with $query\n";
	    sleep($s);
	}
    }
    $this->{'lastTime'}=$time;
    $query=~tr/ /+/;
    
    my $queryString=join('',$eUtilBase,'/',$action,'.fcgi?term=',$query,$this->{'idtag'},'&db=pubmed','&retmode=xml');
    $queryString=~s/term=term=/term=/;
 #print STDERR $queryString,"\n";
    $this->{'lastQuery'}=$queryString;
    return LWP::Simple::get($queryString);
}

my $pmXmlCacheDir="$ENV{'HOME'}/medxml" if ( -d "$ENV{'HOME'}/medxml");
sub _xmlFilePath
{
    my ($pmid)=@_;
    my ($last)=($pmid=~/([0-9])$/);
    return "$pmXmlCacheDir/$last/$pmid.pm.xml";
}
sub _getCached
{
    my @parms=@_;
    shift(@parms) if (ref $parms[0]);
    my $pmid=shift(@parms);
    chomp($pmid);
    die join(" ","Bad PMID! $pmid",@parms) unless ($pmid=~/^[0-9]+$/);
    my $medXmlFile=&_xmlFilePath($pmid);
    my @lines=();
    if ( -f $medXmlFile )
    {
	open(MEDXML,$medXmlFile);
	while (my $line = <MEDXML>)
	{
	    push(@lines,$line);
	}
    }
    else
    {
	print STDERR join(" ","KrPubMed::_getCached",@parms,"attempted on non-existant $medXmlFile\n");
    }
    return @lines;
}

sub pmidsForProtGis
{
    # currently protein only
    my @gis=@_;
    my $molType="protein";
    my $this;
    if (ref $gis[0])
    {
	$this=shift(@gis);
    }
    else
    {
	$this=$defaultPm;
    }
    my %pmids=();
    while (scalar(@gis)>0)
    {
	my @queries=();
	if (scalar(@gis)<$maxSplice)
	{
	    @queries=@gis; @gis=();
	}
	else
	{
	    @queries=splice(@gis,0,$spliceSize);
	}
	my $text=$this->_queryNcbi("elink","dbfrom=protein\&db=pubmed\&id=".join(',',@queries));
	my $on=0;
	
	foreach my $line(split(/\n/,$text))
	{
	    $on=1 if ($line=~/DbTo/);
	    next unless ($on);
	    if ($line=~/<Id>([0-9]+)/)
	    {
		$pmids{$1}++;
	    }
	}
    }
    return keys %pmids;
}

sub cachedFetch
{
    my @pmidList=@_;
# print "cF: ",join(' ',@pmidList),"\n";
    my $this;
    if (ref $pmidList[0])
    {
	$this=shift(@pmidList);
    }
    else
    {
	$this=$defaultPm;
    }
    my %pmids=();
    my @retrievals=();
    my $i=0;
    foreach my $pmid(@pmidList)
    {
	next if (ref $pmids{$pmid});  # already loaded
	my $medXmlFile=&_xmlFilePath($pmid);
	
	my @lines=();
	if ( -f $medXmlFile )
	{
	    my $p=$pmid;
	    @lines=_getCached($p,"KR",$i);
	    $pmids{$pmid}=\@lines;
	    # don't try to consolidate with identical below; this fixes bad retrieves!
	}
#	print "# pmid=$pmid\n";
	if (scalar(@lines)<2)  # if bad retrieve, @lines too short!
	{
	    push(@retrievals,$pmid);
#	print "# pmid2=$pmid\n";
	}
	$i++;
    }
#    print join(" ","retrievals:",@retrievals),"\n";
    while (scalar(@retrievals)>0)
    {
	my @queries;
	if (scalar(@retrievals)>$maxSplice)
	{
	    @queries=splice(@retrievals,0,$spliceSize);
	}
	else
	{
	    @queries=@retrievals;
	    @retrievals=();
	}
	my @lines=split(/\n/,$this->fetchByPmid(join(',',@queries)));
	foreach my $line(@lines)
	{
	    if ($line=~/\<PMID\>([0-9]+)/)
	    {
		my $pmid=$1;
		$pmid=~s/^0+//; # leading zero kill
		my $filePath=&_xmlFilePath($pmid);
		open(MEDXML,">$filePath");
#		print "#\twriting\t$filePath\n";
		print MEDXML "</PubmedArticle>\n";
	    }
	    chomp($line);
	    print MEDXML $line,"\n";
	    close MEDXML if ($line=~/<\/PubmedArticle>/);
	}
	close MEDXML;
	foreach my $pmid(@queries)
	{
	    $pmids{$pmid}=[ $this->_getCached($pmid,"LR") ];
	}
    }
    close MEDXML;
    return \%pmids;
}


1;






