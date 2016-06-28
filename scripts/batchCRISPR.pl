#!/usr/bin/env_perl
use lib "/home/rqge/perl5lib/lib";

use Bio::SearchIO; 
use Getopt::Std;
#use Bio::Tools::Run::StandAloneBlast;

#-----------------------------------------------------

$StarTime=time();

$tmpdir="datas";

if ( -d $tmpdir )
{;} 
else
{mkdir $tmpdir;}

##$fnadir="fnas";
##if ( -d $fnadir )
##{;} 
##else
##{mkdir $fnadir;}

$lookfilelike="fna";

&find_fileindir("./",$lookfilelike);

@typelists=();
$filelike="gff3|fna.PILERCROUT|fna.CRTOUT|gff";
#$filelike="sameRLen.gff3";
@typelists=split(/\|/,$filelike);
$typenumb=@typelists;
$resultflag=0;
for (my $typecount=0;$typecount<$typenumb;$typecount++)
{
  $resultflag=0;
  $resultF=$typelists[$typecount]."\."."result";
  #write into result file:total result file
  @allresult=();
  @allresult=find_result_fileindir("./",$typelists[$typecount]);
  if ($allresult[0] ne 0)
  {
    open(ResultF,">$resultF")||die "$!\n";
    print ResultF @allresult;
    close(ResultF);
  }
}

###old version :no dr infor
for (my $typecount=0;$typecount<$typenumb;$typecount++)
{
  $resultflag=1;
  $resultF=$typelists[$typecount]."\."."oldresult";
  #write into result file:total result file
  @allresult=();
  @allresult=find_result_fileindir("./",$typelists[$typecount]);
  if ($allresult[0] ne 0)
  {
    open(ResultF,">$resultF")||die "$!\n";
    print ResultF @allresult;
    close(ResultF);
  }
}


$EndTime=time();
$AllTime=$EndTime-$StarTime;
print "&&&&&&&&&CRISPRdigger-pilercr-CRT-results&&&&&&&&&&&&The  total result program use the all time minutes:$AllTime s &&&&&&&&&&&&&&&\n ";
exit(0);


sub find_fileindir
{
  local($dir,$filelike) = @_;
  opendir(DIR,"$dir"|| die "can't open this $dir");
  local @files =readdir(DIR);
  closedir(DIR);
  for $file (@files){
    next if($file=~m/\.$/ || $file =~m/\.\.$/);
    if ($file =~/\.($filelike)$/i)
    {         
      system("perl CRISPRdigger.pl -i $dir\/$file \n");   
	  #system("pilercr -in $dir\/$file -out $dir\/$file\.PILERCROUT  -noinfo -minrepeat 23 -maxrepeat 55 -minspacer 10 -maxspacer 120 \n");  
      #system("java -cp /home/rqge/bin/CRT1.2-CLI.jar crt -minRL 23 -maxRL 55 -minSL 10 -maxSL 120 $dir\/$file $dir\/$file\.CRTOUT \n");
      
#     system ("cp $dir\/$file .\/ \n");
#      system(" perl CRISPRFinder.pl $file $file\.gffout \n ");
      ##system("rm -r .\/*.fna \n");
      ##system("rm -r .\/$file \n");
    }
    elsif(-d "$dir/$file"){
            find_fileindir("$dir/$file",$filelike );
    }      
  }
}

sub find_result_fileindir
{
  local($dir,$filelike) = @_;
  opendir(DIR,"$dir"|| die "can't open this $dir");
  local @files =readdir(DIR);
  closedir(DIR);
  for $file (@files){
    next if($file=~m/\.$/ || $file =~m/\.\.$/);
    if ($file =~/\.($filelike)$/i)
    {
      #print "$dir\/$file \n";      
      # system("cp $dir\/$file $putindir\/$file ");
      $curfiledir="$dir/$file";
      # system("perl resultcount.pl -i $dir\/$file \n");
       #every result file statistics
       @oneresult=();
       if ($resultflag eq 0)
       {
	@oneresult=readresult($curfiledir);
       }
       elsif($resultflag eq 1)
       {
	@oneresult=old_readresult($curfiledir);
       } 
       if ($oneresult[0] ne 0)
       {
	push @allresult,@oneresult;
       }
       ##total result file:detail file        
    }
    elsif(-d "$dir/$file"){
            find_result_fileindir("$dir/$file",$filelike );
    }    
  }
  return @allresult; 
}

sub readresult
{	
	my ($filedirname) =@_;
	##@filepath=split('\/',$filedirname);
	##$filename=$filepath[$#filepath];
	@resultStrs=();
	if ($filedirname=~/\.(gff3)$/)
	{
	  #count CRISPR details: all numbers,beg-end in every CRISPR
	  # NC_***  ProgrammeName CRISPRnumbs beg1 end1 SP1num beg2 end2 SP2num ...
	  
	  @resultStrs=mygff3out($filedirname);
	  return @resultStrs;
	}
	if($filedirname=~/\.(PILERCROUT)$/)
	{ 
	  @resultStrs=pilercrout($filedirname);
	  return @resultStrs;
	}
	
	if($filedirname=~/\.(CRTOUT)$/)
	{ 
	  @resultStrs=crtout($filedirname);
	  return @resultStrs;
	}
	elsif($filedirname=~/\.(gff)$/i)
	{
	  
	  @resultStrs=finderout($filedirname);
	  return @resultStrs;
	}
	else
	{
	  return 0;
	}			
}

sub old_readresult
{	
	my ($filedirname) =@_;
	##@filepath=split('\/',$filedirname);
	##$filename=$filepath[$#filepath];
	@resultStrs=();
	if ($filedirname=~/\.(gff3)$/)
	{
	  #count CRISPR details: all numbers,beg-end in every CRISPR
	  # NC_***  ProgrammeName CRISPRnumbs beg1 end1 SP1num beg2 end2 SP2num ...
	  
	  @resultStrs=old_mygff3out($filedirname);
	  return @resultStrs;
	}
	if($filedirname=~/\.(PILERCROUT)$/)
	{ 
	  @resultStrs=old_pilercrout($filedirname);
	  return @resultStrs;
	}
	
	if($filedirname=~/\.(CRTOUT)$/)
	{ 
	  @resultStrs=old_crtout($filedirname);
	  return @resultStrs;
	}
	elsif($filedirname=~/\.(gff)$/i)
	{
	  
	  @resultStrs=old_finderout($filedirname);
	  return @resultStrs;
	}
	else
	{
	  return 0;
	}			
}


#count CRISPR details: all numbers,beg-end in every CRISPR ;;;add all DR positions information
# Genome_name	NC_***  ProgrammeName CRISPRnumbs beg1 end1 SP1num   dr1b dr1e;dr2b dr2e;       beg2 end2 SP2num ...
#Thermococcus barophilus MP chromosome	NC_014804.1 mygff3out  4   398320 398761 6 579477 580049 8 1334169  1334873 10  1530306 1532069 25
## the result file like::
##gi|212223144|ref|NC_011529.1|	LibRepeatMasker	CRISPR	294121	295760	196	+	.	ID=294121;Name=0;CRISPRnum=0;StrLen=1640;aSPLen=37;SPnum=24;aDRLen=29;isCRISPR=1
## gi|212223144|ref|NC_011529.1| LibRepeatMasker Repeats 294121	294145	196	+	.	ID=21;Parent=294121;Name=0|CRISPRnum=0
sub mygff3out
{
  @resultgff3=();
  my ($filename) =@_;
  $crisprnumb=0;
  #@locations=();
  $programname="CRISPRdigger";
  @getname=();
  $crisprloc=();
  %crisprId=();
  %repeatId=();
  %sortCRISPR=();
  
  open(GF3, "$filename")||die "$!\n";
  ##every crispr like:: |beg1 end1 SP1num;dr1b,dr1e;dr2b,dr2e |beg2 end2 SP2num;dr21b,dr22e;dr31b,dr32e;  
  while(<GF3>)
  {
     chomp;
     next if /^#/;  # Allow embedded comments.
     $body =$_;
     $body=~s/^\s+//g;
     my @gffstr=();
     @gffstr =split(/\s+/,$body);
     $IdCRISPR=();
     $isCRISPR=0;
     ($isCRISPR)=($gffstr[$#gffstr]=~/isCRISPR\=(\d+)/);
       
     if (($gffstr[2] eq "CRISPR")&&($isCRISPR eq 1))
     {
       @getname=split('\|',$gffstr[0]);
       ($IdCRISPR)=($gffstr[$#gffstr]=~/ID\=(\d+)/);
       $crisprnumb++;
       ($strspnum)=($gffstr[$#gffstr]=~/SPnum\=(\d+)/);
       $begposition=$gffstr[3];
       $endposition= $gffstr[4];
       $crisprloc="\|".$begposition."\t".$endposition."\t".$strspnum."\;";
       $crisprId{$IdCRISPR}=$crisprloc;  
     }
     elsif($gffstr[2] eq "Repeats")
     {
        my ($parentId)=($gffstr[$#gffstr]=~/Parent\=(\d+)\;/);
        my $tempdr=$gffstr[3]."\t".$gffstr[4];
        if (exists $repeatId{$parentId})
        { 
	  $repeatId{$parentId}=$repeatId{$parentId}."\;".$tempdr;
        }
	else
	{
	  $repeatId{$parentId}=$tempdr;
	}  
      }
  }
  close(GF3);
  if ($crisprnumb>0)
  {
    $truename=();  
    $gename=$getname[3];
    $truename=fdnamefromacc($gename); 
    $allCRISPRs=();
    
  #  foreach $CRId(sort keys %crisprId)
    foreach $CRId(sort keys %crisprId)
    {
      foreach $DrId(keys %repeatId)
      {
	if ($CRId eq $DrId)
	{
	  $oneCRISPR=$crisprId{$CRId}.$repeatId{$DrId};
	  my ($seqloc)=($crisprId{$CRId}=~/\|(\d+)\t/);
	  $sortCRISPR{$seqloc}= $oneCRISPR ;
	  last;
	}	
      }
##      $allCRISPRs.=$oneCRISPR;
    }
    
#    foreach $DRId(sort keys %sortCRISPR)
     foreach $lockey(sort{ $a <=> $b} keys%sortCRISPR)
    {
      $allCRISPRs.=$sortCRISPR{$lockey}; 
    }
    
    push @resultgff3,$truename."\t".$gename."\t".$programname."\t".$crisprnumb.$allCRISPRs."\n";
    ##push @resultgff3,$truename."\t".$gename."\t".$programname."\t".$crisprnumb.$crisprloc."\n";
    return  @resultgff3 ;
  }
  else
  {
    return 0;
  }
}

#count CRISPR details: all numbers,beg-end in every CRISPR
# Genome_name	NC_***  ProgrammeName CRISPRnumbs beg1 end1 SP1num  beg2 end2 SP2num ...
#Thermococcus barophilus MP chromosome	NC_014804.1 mygff3out  4   398320 398761 6 579477 580049 8 1334169  1334873 10  1530306 1532069 25
sub old_mygff3out
{
  @resultgff3=();
 my ($filename) =@_;
  $crisprnumb=0;
  #@locations=();
  $programname="CRISPRdigger";
  @getname=();
  $crisprloc=();
  %sortCRISPR=();
  
  open(GF3, "$filename")||die "$!\n";
  
  while(<GF3>)
  {
     chomp;
     next if /^#/;  # Allow embedded comments.
     $body =$_;
     $body=~s/^\s+//g;
     @gffstr=();
     @gffstr =split(/\s+/,$body);
     ($isCRISPR)=($gffstr[$#gffstr]=~/isCRISPR\=(\d+)/);
       
     if (($gffstr[2] eq "CRISPR")&&($isCRISPR eq 1))
     {
       $crisprnumb++;
       ($strspnum)=($gffstr[$#gffstr]=~/SPnum\=(\d+)/);
       $begposition=$gffstr[3];
       $endposition= $gffstr[4];
       $oneCR="\|".$begposition."\t".$endposition."\t".$strspnum;
       $sortCRISPR{$begposition}= $oneCR;
     # $crisprloc.=$oneCR;
      #push @locations, "\t".$begposition."\t".$endposition."\t".$strspnum;
     }
     else
     {
       next;
     } 
  }
  close(GF3);
  if ($crisprnumb>0)
  {
    $allCRISPR=();
    $truename=();
    @getname=split('\|',$gffstr[0]);
    $gename=$getname[3];
    $truename=fdnamefromacc($gename);
    
    # foreach $DRId(sort keys %sortCRISPR)
     foreach $keyloc(sort{$a <=> $b} keys%sortCRISPR)
    {
      $allCRISPR.=$sortCRISPR{$keyloc}; 
    }
    
    push @resultgff3,$truename."\t".$gename."\t".$programname."\t".$crisprnumb.$allCRISPR."\n";
    ##push @resultgff3,$truename."\t".$gename."\t".$programname."\t".$crisprnumb.$crisprloc."\n";
    return  @resultgff3 ;
  }
  else
  {
    return 0;
  }
}

##piler-CR CRISPR output contents
##count CRISPR details: all numbers,beg-end in every CRISPR ;;;add all DR positions information
## Genome_name	NC_***  ProgrammeName CRISPRnumbs beg1 end1 SP1num   dr1b dr1e;dr2b dr2e;       beg2 end2 SP2num ...
###  |beg1 end1 SP1num;dr1b,dr1e;dr2b,dr2e |beg2 end2 SP2num;dr21b,dr22e;dr31b,dr32e;
###Thermococcus sp	NC_015865.1	mygff3out	5|131167	133209	3;131167,131195;131235,131263;131301,131329|1907425	1908326	3;1907425,1907454;1907492,1907521;1907560,1907589
##Thermococcus barophilus MP chromosome	NC_014804.1 mygff3out  4   398320 398761 6 579477 580049 8 1334169  1334873 10  1530306 1532069 25
sub pilercrout
{
  @resultgff3=();
  my ($filename) =@_;
  my $crisprnumb=0;
  #@locations=();
  $programname="pilerCR";
  my @getname=();
  my $onecrispr=(); 
  open(CR, "$filename")||die "$!\n";
  my $geneflag=0;
  while(<CR>)
  {
     chomp;
     if ($_=~/(\d+) putative CRISPR arrays found/)
     {
      $crisprnumb=$1;
     } 
     if ($_=~/SUMMARY BY SIMILARITY/)
     {
      last;
     }   
     next if /^#/;  # Allow embedded comments.      
     
     $body =$_;
     $body=~s/^\s+//g;   
     if ($body=~/Array\s+\d+/)
     {
	my ($parentId)=($body=~/Array\s+(\d+)/);
	my $begflag=0;
	my $spnumb=-1;
	my $drs=();
	my $lineflag=0;
	while (<CR>)
	{
	  chomp;
	  $body =$_;
	  $body=~s/^\s+//g;
	   if ($body=~/==========/)
	  {$lineflag++;}
	  if ($lineflag eq 2)
	  {
	    $endpos=$gffstr[0]+$gffstr[1]-1;
	    $onecrispr.="\|".$begpos."\t".$endpos."\t".$spnumb.$drs;
	    last;
	  }	  
	  @gffstr=();
	  @gffstr =split(/\s+/,$body); 
	  if ((/^>/)&&($geneflag eq 0))
	  {
	   $geacc=();
	   $truename=();
	   $geacc=$_;
	   
	   #$geacc=~/ref\|(NC_\d+\.\d+)\|\s+(.*?)(\,\s+|\s+)complete genome/i;
	   #$gename=$1;
	   #$truename=$2;
           
           ($gename)=($geacc=~/ref\|(NC_\d+\.\d+)\|/i);
          #???($truename)=($geacc=~/NC_\d+\.\d+\|\s+((\w+\s+)*\w+)\,/i);
          $truename=fdnamefromacc($gename);
           
	   $geneflag=1;
	   next;
	  }  	  
	  if ($gffstr[0]=~/\d+$/)
	  {
	    if ($begflag eq 0)
	    {
	      $begpos=$gffstr[0];
	      $begflag=1;
	    }
	    $begdr=$gffstr[0];
	    $endr=$gffstr[0]+$gffstr[1]-1;
	    $spnumb++;
	    $drs.="\;".$begdr."\t".$endr; 	  
	  }	  
	}
      } 
  }
  close(CR);
  if ($crisprnumb>0)
  {
    ##push @resultgff3,$gename."\t".$programname."\t".$crisprnumb.@locations."\n";
    push @resultgff3,$truename."\t".$gename."\t".$programname."\t".$crisprnumb.$onecrispr."\n";
    return  @resultgff3 ;
  }
  else
  {
    return 0;
  } 
}
#piler-CR CRISPR output contents
#count CRISPR details: all numbers,beg-end in every CRISPR
# Genome_name	NC_***  ProgrammeName CRISPRnumbs beg1 end1 SP1num beg2 end2 SP2num ...
#Thermococcus barophilus MP chromosome	NC_014804.1 mygff3out  4   398320 398761 6 579477 580049 8 1334169  1334873 10  1530306 1532069 25
sub old_pilercrout
{
  @resultgff3=();
  my ($filename) =@_;
  my $crisprnumb=0;
  #@locations=();
  $programname="pilerCR";
  my @getname=();
  my $crisprloc=();
  
  open(CR, "$filename")||die "$!\n";
  $starflag=0;
  while(<CR>)
  {
     chomp;
     if ($_=~/SUMMARY BY POSITION/)
     {
      $starflag=1;
      next;
     }
     
     next if /^#/||($starflag eq 0);  # Allow embedded comments. 
     ##  From "SUMMARY BY POSITION" begin countnumb
     
     if (/^>/)
     {
      $geacc=();
      $geacc=$_;
      $truename=(); 
      ($gename)=($geacc=~/ref\|(NC_\d+\.\d+)\|/i);
      #???($truename)=($geacc=~/NC_\d+\.\d+\|\s+((\w+\s+)*\w+)\,/i);
      $truename=fdnamefromacc($gename);
      next;
     }
     if (length($_)<5)
     {next;}
     my $body =$_;
     $body=~s/^\s+//g;
     my @gffstr=();
     @gffstr =split(/\s+/,$body);
     if ($gffstr[0]=~/\d+$/)
     {
      $crisprnumb++;
      $begposition=$gffstr[2];
      $endposition=$begposition+$gffstr[3];
      $strspnum=$gffstr[4]-1;
      $crisprloc.="\|".$begposition."\t".$endposition."\t".$strspnum;
     }   
  }
  close(CR);
  if ($crisprnumb>0)
  {
    ##push @resultgff3,$gename."\t".$programname."\t".$crisprnumb.@locations."\n";
    push @resultgff3,$truename."\t".$gename."\t".$programname."\t".$crisprnumb.$crisprloc."\n";
    return  @resultgff3 ;
  }
  else
  {
    return 0;
  }
  
}

#CRT CRISPR output contents
sub crtout
{
  @resultgff3=();
  my ($filename) =@_;
  my $crisprnumb=0;
  #@locations=();
  $programname="CRT";
  my $onecrispr=();
  my $gename=();
  my $truename=();
  my $CRISPRId=0;
  open(CR, "$filename")||die "$!\n";
  my $geneflag=0;
  while(<CR>)
  {
     chomp;
     next if /^#/;  # Allow embedded comments.
     $body =$_;
     $body=~s/^\s+//g;   
     if ($body=~/ORGANISM/)
     {
#       #($gename)=($body=~/ref\|(NC_\d+\.\d+)\|/i);
#       ##$truename=fdnamefromacc($gename);
#       #($truename)=($body=~/NC_\d+\.\d+\|\s+(.*?)[\,\d+|\d+]complete genome/i);      
#	$body=~/ref\|(NC_\d+\.\d+)\|\s+(.*?)(\,\s+|\s+)complete genome/i;
#	$gename=$1;
#	$truename=$2;
        
          my @lines=();
        @lines=split(/\|/,$body);
        $gename=$lines[3];
        
	#??? $truename=$2;
	$truename=fdnamefromacc($gename);
        
     }  	    
     if ($body=~/CRISPR\s+\d+/)
     {
	($CRISPRId)=($body=~/CRISPR\s+(\d+)/);
	my $begflag=0;
	my $spnumb=-1;
	my $drs=();
	my $lineflag=0;
	my ($begpos)=($body=~/Range:\s+(\d+)\s+\-/);
	my ($endpos)=($body=~/Range:\s+\d+\s+\-\s+(\d+)/);	  
	while (<CR>)
	{
	  chomp;
	  $body =$_;
	  $body=~s/^\s+//g;
	   if ($body=~/--------/)
	  {$lineflag++;}
	  if ($lineflag eq 2)
	  {
	    $onecrispr.="\|".$begpos."\t".$endpos."\t".$spnumb.$drs;
	    last;
	  }	  
	  if ($body=~/^\d+\s+/)
	  {
	    $spnumb++;
	    ($begdr)=($body=~/^(\d+)\s+/);
	    ($drlen)=($body=~/\[\s+(\d+)\,/);
	    $endr=$begdr+$drlen-1;
	    $drs.="\;".$begdr."\t".$endr; 
	  }	  
	}
      } 
  }
  close(CR);
  
  $crisprnumb=$CRISPRId;
  if ($crisprnumb>0)
  {
    ##push @resultgff3,$gename."\t".$programname."\t".$crisprnumb.@locations."\n";
    push @resultgff3,$truename."\t".$gename."\t".$programname."\t".$crisprnumb.$onecrispr."\n";
    return  @resultgff3 ;
  }
  else
  {
    return 0;
  }   
}

##CRT CRISPR output contents
#CRT CRISPR output contents
#count CRISPR details: all numbers,beg-end in every CRISPR
# Genome_name	NC_***  ProgrammeName CRISPRnumbs beg1 end1 SP1num beg2 end2 SP2num ...
#Thermococcus barophilus MP chromosome	NC_014804.1 mygff3out  4   398320 398761 6 579477 580049 8 1334169  1334873 10  1530306 1532069 25
sub old_crtout
{
  @resultgff3=();
  my ($filename) =@_;
  my $crisprnumb=0;
  #@locations=();
  $programname="CRT";
  my $onecrispr=();
  my $gename=();
  my $truename=();
  my $CRISPRId=0;
  open(CR, "$filename")||die "$!\n";
  my $geneflag=0;
  while(<CR>)
  {
     chomp;
     next if /^#/;  # Allow embedded comments.
     $body =$_;
     $body=~s/^\s+//g;   
     if ($body=~/ORGANISM/)
     {
       ##($gename)=($body=~/ref\|(NC_\d+\.\d+)\|/i);
       ###$truename=fdnamefromacc($gename);
       ##($truename)=($body=~/NC_\d+\.\d+\|\s+(.*?)[\,\d+|\d+]complete genome/i);      
	#$body=~/ref\|(NC_\d+\.\d+)\|\s+(.*?)(\,\s+|\s+)complete genome/i;
	#$gename=$1;
        
        my @lines=();
        @lines=split(/\|/,$body);
        $gename=$lines[3];
        
	#??? $truename=$2;
	$truename=fdnamefromacc($gename);
     }  	    
     if ($body=~/CRISPR\s+\d+/)
     {
	($CRISPRId)=($body=~/CRISPR\s+(\d+)/);
	my $begflag=0;
	my $spnumb=-1;
	my $drs=();
	my $lineflag=0;
	my ($begpos)=($body=~/Range:\s+(\d+)\s+\-/);
	my ($endpos)=($body=~/Range:\s+\d+\s+\-\s+(\d+)/);	  
	while (<CR>)
	{
	  chomp;
	  $body =$_;
	  $body=~s/^\s+//g;
	   if ($body=~/--------/)
	  {$lineflag++;}
	  if ($lineflag eq 2)
	  {
	    $onecrispr.="\|".$begpos."\t".$endpos."\t".$spnumb;
	    last;
	  }	  
	  if ($body=~/^\d+\s+/)
	  {
	    $spnumb++;
	  }	  
	}
      } 
  }
  close(CR);
  
  $crisprnumb=$CRISPRId;
  if ($crisprnumb>0)
  {
    ##push @resultgff3,$gename."\t".$programname."\t".$crisprnumb.@locations."\n";
    push @resultgff3,$truename."\t".$gename."\t".$programname."\t".$crisprnumb.$onecrispr."\n";
    return  @resultgff3 ;
  }
  else
  {
    return 0;
  }   
}

##fine bacteria name from NC_number
sub old_fdnamefromacc
{
  my ($accid)=@_;
  ##truename and Accession reference table
  $findfile="bacteria.summary.txt";
  @nameacc=();
  $truename=();
  open(FN, "$findfile")||die "$!\n";
  
  while(<FN>)
  {
     chomp;
     next if /^#/;  # Allow embedded comments.
     @allacc=();
     my $oneline =$_;
     $oneline=~s/^\s+//g;
     @nameacc =split(/\t/,$oneline);
     @allacc=split(/\./,$nameacc[0]);
     if (($accid=~/($allacc[0])/i)||($accid eq $nameacc[0]))    
     {
      $truename=$nameacc[5];
      return $truename;
     }    
  }
  close(FN);
  return 0;
}
##fine bacteria name from NC_number
sub fdnamefromacc
{
  my ($accid)=@_;
  ##truename and Accession reference table
  $findfile="GoldenCRISPRs";
  @nameacc=();
  $truename=();
  open(FN, "$findfile")||die "$!\n";
  
  while(<FN>)
  {
     chomp;
     next if /^#/;  # Allow embedded comments.
     @allacc=();
     my $oneline =$_;
     $oneline=~s/^\s+//g;
     @nameacc =split(/\t/,$oneline);
     $ncname=$nameacc[1];
     if (($accid=~/($ncname)/i)||($accid eq $ncname))    
     {
      $truename=$nameacc[0];
      return $truename;
     }    
  }
  close(FN);
  return 0;
}
#count CRISPR details: all numbers,beg-end in every CRISPR ;;;add all DR positions information
# Genome_name	NC_***  ProgrammeName CRISPRnumbs beg1 end1 SP1num   dr1b dr1e;dr2b dr2e;       beg2 end2 SP2num ...
#Thermococcus barophilus MP chromosome	NC_014804.1 mygff3out  4   398320 398761 6 579477 580049 8 1334169  1334873 10  1530306 1532069 25
## the result file like::
##gi|212223144|ref|NC_011529.1|	LibRepeatMasker	CRISPR	294121	295760	196	+	.	ID=294121;Name=0;CRISPRnum=0;StrLen=1640;aSPLen=37;SPnum=24;aDRLen=29;isCRISPR=1
## gi|212223144|ref|NC_011529.1| LibRepeatMasker Repeats 294121	294145	196	+	.	ID=21;Parent=294121;Name=0|CRISPRnum=0
sub finderout
{
  @resultgff3=();
  my ($filename) =@_;
  $crisprnumb=0;
  #@locations=();
  $programname="CRISPRFinder";
  @getname=();
  $crisprloc=();
  %crisprId=();
  %repeatId=();
  %sortCRISPR=();
  
  open(GF3, "$filename")||die "$!\n";
  ##every crispr like:: |beg1 end1 SP1num;dr1b,dr1e;dr2b,dr2e |beg2 end2 SP2num;dr21b,dr22e;dr31b,dr32e;  
  while(<GF3>)
  {
     chomp;
     next if /^#/;  # Allow embedded comments.
     $body =$_;
     $body=~s/^\s+//g;
     my @gffstr=();
     @gffstr =split(/\s+/,$body);
     $IdCRISPR=();  
       
     if ($gffstr[2] eq "CRISPR")
     {
       @getname=split('\|',$gffstr[0]);
       ($IdCRISPR)=($gffstr[$#gffstr]=~/ID\=(.*?)$/);
       $crisprnumb++;
       ($strspnum)=($gffstr[$#gffstr]=~/Number_of_spacers\=(\d+)\;/);
       $begposition=$gffstr[3];
       $endposition= $gffstr[4];
       $crisprloc="\|".$begposition."\t".$endposition."\t".$strspnum."\;";
       $crisprId{$IdCRISPR}=$crisprloc;  
     }
     elsif($gffstr[2] eq "CRISPRdr")
     {
        my ($parentId)=($gffstr[$#gffstr]=~/Parent\=(.*?)$/);
        my $tempdr=$gffstr[3]."\,".$gffstr[4];
        if (exists $repeatId{$parentId})
        { 
	  $repeatId{$parentId}=$repeatId{$parentId}."\;".$tempdr;
        }
	else
	{
	  $repeatId{$parentId}=$tempdr;
	}  
      }
  }
  close(GF3);
  if ($crisprnumb>0)
  {
    $truename=();  
    $gename=$getname[3];
    $truename=fdnamefromacc($gename); 
    $allCRISPRs=();
    
  #  foreach $CRId(sort keys %crisprId)
    foreach $CRId(sort keys %crisprId)
    {
      foreach $DrId(keys %repeatId)
      {
	if ($CRId eq $DrId)
	{
	  $oneCRISPR=$crisprId{$CRId}.$repeatId{$DrId};
	  my ($seqloc)=($crisprId{$CRId}=~/\|(\d+)\t/);
	  $sortCRISPR{$seqloc}= $oneCRISPR ;
	  last;
	}	
      }
##      $allCRISPRs.=$oneCRISPR;
    }
    
#    foreach $DRId(sort keys %sortCRISPR)
     foreach $lockey(sort{ $a <=> $b} keys%sortCRISPR)
    {
      $allCRISPRs.=$sortCRISPR{$lockey}; 
    }
    
    push @resultgff3,$truename."\t".$gename."\t".$programname."\t".$crisprnumb.$allCRISPRs."\n";
    ##push @resultgff3,$truename."\t".$gename."\t".$programname."\t".$crisprnumb.$crisprloc."\n";
    return  @resultgff3 ;
  }
  else
  {
    return 0;
  }
}

#count CRISPR details: all numbers,beg-end in every CRISPR
# Genome_name	NC_***  ProgrammeName CRISPRnumbs beg1 end1 SP1num  beg2 end2 SP2num ...
#Thermococcus barophilus MP chromosome	NC_014804.1 mygff3out  4   398320 398761 6 579477 580049 8 1334169  1334873 10  1530306 1532069 25
sub old_finderout
{
  @resultgff3=();
 my ($filename) =@_;
  $crisprnumb=0;
  #@locations=();
  $programname="CRISPRFinder";
  @getname=();
  $crisprloc=();
  %sortCRISPR=();
  
  open(GF3, "$filename")||die "$!\n";
  
  while(<GF3>)
  {
     chomp;
     next if /^#/;  # Allow embedded comments.
     $body =$_;
     $body=~s/^\s+//g;
     @gffstr=();
     @gffstr =split(/\s+/,$body);
    # ($isCRISPR)=($gffstr[$#gffstr]=~/isCRISPR\=(\d+)/);
       
     if ($gffstr[2] eq "CRISPR")
     {
       @getname=split('\|',$gffstr[0]);
       $crisprnumb++;
       ($strspnum)=($gffstr[$#gffstr]=~/Number_of_spacers\=(\d+)\;/);
       $begposition=$gffstr[3];
       $endposition= $gffstr[4];
       $oneCR="\|".$begposition."\t".$endposition."\t".$strspnum;
       $sortCRISPR{$begposition}= $oneCR;
     # $crisprloc.=$oneCR;
      #push @locations, "\t".$begposition."\t".$endposition."\t".$strspnum;
     }
     else
     {
       next;
     } 
  }
  close(GF3);
  if ($crisprnumb>0)
  {
    $allCRISPR=();
    $truename=();
  #  @getname=split('\|',$gffstr[0]);
    $gename=$getname[3];
    $truename=fdnamefromacc($gename);
    
    # foreach $DRId(sort keys %sortCRISPR)
     foreach $keyloc(sort{$a <=> $b} keys%sortCRISPR)
    {
      $allCRISPR.=$sortCRISPR{$keyloc}; 
    }
    
    push @resultgff3,$truename."\t".$gename."\t".$programname."\t".$crisprnumb.$allCRISPR."\n";
    ##push @resultgff3,$truename."\t".$gename."\t".$programname."\t".$crisprnumb.$crisprloc."\n";
    return  @resultgff3 ;
  }
  else
  {
    return 0;
  }
}

#delete null line in the file
##deletenulline($cmpresults);
sub old_deletenulline
{
    my ($infile)=@_;
    open(AF, "$infile")||die "$!\n";
    $tempfile="$infile.temp";
    open (TF,">$tempfile")||  die "Cannot create file"; 
    
    while(<AF>)
    {	
	my $len=length("\n");
	chomp;
	$tline=$_;
	$tline=~s/[\r\n]+//g;
	#print $_;
	#if ((!/^[\r\n]*$/)&&(length($tline) > 2))
	if ((length($tline) > 2))
	{print TF "$tline \n";}
    }
    close(AF);
    close(TF);
    #my $rm = "rm $dir/$inputrefseqfile"; system($rm);
    system("rm $infile ");
    rename("$tempfile", "$infile");
    
}

##deletenulline($cmpresults);
sub deletenulline
{
    my ($infile)=@_;
    
    open(AF, "$infile")||die "$!\n";
    my $tempfile="$infile.temp";
    open (TF,">$tempfile")||  die "Cannot create file:$tempfile"; 
    
    while(<AF>)
    {	
	#chomp;
	$tline=$_;
	$tline=~s/[\r\n]+//g;
	#print $_;
	#if ((!/^[\r\n]*$/)&&(length($tline) > 2))
	
	if ((length($tline) > 2))
	{
		print TF "$tline \n";
	}
    }
    
    close(TF);   
    close(AF);
    #my $rm = "rm $dir/$inputrefseqfile"; system($rm);
    system("rm -f $infile \n");
    rename("$tempfile", "$infile");
}

