#!/usr/bin/env_perl
use lib "/home/rqge/perl5lib/lib";
#
##use lib "/home/smile/perl5lib/lib";
#Author: grayspring
#use strict;
#use Bio::Seq;
use Bio::SearchIO; 
use Bio::SeqIO;
use Getopt::Std;
use Bio::SimpleAlign;
use  Bio::Tools::Run::Alignment::Clustalw;


#use Bio::Tools::Run::StandAloneBlast;

#-----------------------------------------------------
getopts("i:B:d:t:s:p:R:g");

$Inputfile   = defined $opt_i ? $opt_i : "";
#$output    = defined $opt_o ? $opt_o : "";

$SPMinDRratio=defined $opt_l ? $opt_l : 0.5;
$SPMaxDRratio=defined $opt_m ? $opt_m : 3.0;
$MinSP     = defined $opt_g ? $opt_g : 10;
$MaxSP     = defined $opt_p ? $opt_p : 120;
$DiffTempSim = defined $opt_t ? $opt_t : 0.80;#0.7-0.9
$SPsimilary = defined $opt_s ? $opt_s : 0.50;#0.4-0.8
$SameRratio = defined $opt_R ? $opt_R : 0.50;#0.4-0.8
$blastsameR = defined $opt_B ? $opt_B : 0.35;#0.2-0.7
$MaxDR     = defined $opt_d ? $opt_d: 47 ;
$MinDR     = defined $opt_r ? $opt_r : 10;

##read the fasta data ,sub the first line,then substr

usuage() if((!$Inputfile)||($Help));
# -----------------------------------------------------

$StarTime=time();

#open(IF, "$Input")||die "$!\n";
#考虑增加去掉文件路径代码
#增加去掉文件路径代码
@pathfile=split(/\//,$Inputfile);
$filefullname=$pathfile[$#pathfile];
$filepath=();
for ($pathI=0;$pathI<$#pathfile;$pathI++)
{	
    $filepath=$filepath.$pathfile[$pathI]."\/";		
}

@filename =split( /\./,$filefullname);

#@filename =split( /\./,$Input);
#$Infile = $filename[0];

my @DRstr='';
my @lastgff='';
my $subfilenum=0;
$outfile=$filepath.$filename[0];
$DRout="$outfile.dr";
open(DRoutF,">>$DRout")||die "$!\n";
push @DRstr,"#sequenceName","Start","End","oldID","ID","SPdist","CRISPRnum","Score","DRstring","SPstring\n";
print DRoutF join("\t", @DRstr);

$SPout= "$outfile.sp";
open(spoutF,">>$SPout")||die "$!\n";
#push @spgff,"#SeqenceID","Source","Type","Start","End","Score","Strand","Phase","Rtype","CRISPRnum","SPnum","spLen","isame","SPstring\n";
#print spoutF join("\t", @spgff);

$gffoutput= "$outfile.gff3";
open(crisproutF,">>$gffoutput")||die "$!\n";
push @lastgff,"#SeqID","Source","Type","Start","End","Score","Strand","Phase","Attributes\n";
print crisproutF join("\t", @lastgff);


my $bodystr;
my $egSeqObj = new Bio::SeqIO( -file   => "$Inputfile", -format => "fasta");
while ( my $egSeqO = $egSeqObj->next_seq )
{
	$tline = $egSeqO->seq; $bodystr = "\U$tline\E"; $Inputid  = $egSeqO->id;
	
	$subfilenum++;
	my @inputsets=split(/\s+/,$Inputid);
	my $subfile=$inputsets[0];
	my $Input="subfile"."_".$subfilenum;
	
	$Infile = $filepath.$filename[0].'_'.$Input;
	
	open (IF,">$Input")||die "Error in opening the file: $Input\n";
    print IF "\>$subfile\n"; 
	print IF "$bodystr\n";
    close(IF);

	$fnastrlen=length($bodystr);
	$default_l=int( 1 + log($fnastrlen) / log(4) );
	$extendeachsides= $MaxDR-$default_l;

	#if (stat("$Infile\_rep_filt_stg1.fasta")<0)
	if ((-e "$Infile\_rep_filt_stg1.fasta" )&&(-s "$Infile\_rep_filt_stg1.fasta" >40 ))
	{
		print  "file $Infile\_rep_filt_stg1.fasta exist\n";
		print " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n";
		
	}
	else
	{
		system ("build_lmer_table -tandem 60 -min 2 -sequence $Input -freq $Infile\.freq ");
		system ("RepeatScout -sequence $Input -output $Infile\_rep.fasta -freq $Infile\.freq -L $extendeachsides -goodlength 20 -minthresh 2 -tandemdist 60");
		system ("cat $Infile\_rep.fasta | /home/rqge/bin/RepeatScout/filter-stage-1.prl > $Infile\_rep_filt_stg1.fasta ");
		
	}
		
	my @stg1args=stat("$Infile\_rep_filt_stg1.fasta");
	##$stg1size = -s "$Infile\_rep_filt_stg1.fasta";
	##if ($stg1size >0)
	if ($stg1args[7] >0)
	{
		##delete len(DR)>50 string in ***_filt_stg1.fasta
		open(Stg1,"$Infile\_rep_filt_stg1.fasta")||die "$!\n";
		$stg1out="$Infile\.stg1.fasta";
		open (Stout,">$stg1out")||die "$!\n";
		
		while (<Stg1>)
		{
			chomp;
			next if /^#/;  # Allow embedded comments.
			if(/^>/) 
			{
				 $Name =$_;	
			}
			else 
			{
				my $tempRstr=$_;
				if ((length($tempRstr)<= $MaxDR)&&(length($tempRstr)>= $MinDR))
				{
					print Stout $Name."\n";
					print Stout $tempRstr."\n";				
				}		
			}		
		}	
		close(Stg1);
		close(Stout);
		
		##if ((-e "$Infile\_rep_filt_stg2_thresh2.fasta" )&&(-s "$Infile\_rep_filt_stg2_thresh2.fasta" >40 ))
		if ((-e "$Input.out" )&&(-s "$Input.out" >40 ))
		{
			print  "file $Input.out exist\n";
		}
		else
		{
			system ("RepeatMasker -cutoff 120 -s -lib $stg1out $Input ");		
			$repeatime=time();
			$onerepeatime=$repeatime-$StarTime;
			print "***********repeatmaker.out ********running first time RepeatMasker ,the program use the time (s):$onerepeatime s ************** \n ";		
				
		} 
		
		####next step : judge the spacer :0.6~2 times length(DR)
				
			$output= "$Input.out";
			$dataF="$output.gff3.data";
			open(DatF,">$dataF")||die "$!\n";
			open(outP,"$output")||die "$!\n";
			##131 severs use this command
			system ("grep \"Unspecified\" $output > $dataF");

		##system ("grep \"Unknown\" $output > $dataF");
		close(DatF);
		close(outP);	
		
		#mkdir datas for putting the tempt results
		$tmpdir="datas";
		if ( -d $tmpdir )
		{;} 
		else

		{mkdir $tmpdir;}
		
		$samelendata="$dataF.samelen";
		open(DatF,"$dataF")||die "$!\n";
		open(SameData,">$samelendata")||die "$!\n";
		my @line;
		my $linenum=0;
		%allsameRs=();
		@similR=();
	#	$templatefile="$Infile\_rep_filt_stg2_thresh2.fasta";
		$templatefile="$Infile\.stg1.fasta";
		##find the silmlarly R
		%allsameRs=sameRtemplate($templatefile);		
		
		#@similR=values %allsameRs;
		### find the different direction--Reverse complement string in template R
		#将$templatefile自身比对找出存在的反相互补序列
		%reversestr=();
		%reversestr=ReverseCmpTemp($templatefile,%allsameRs);
			
		###每次读两行进行分析，赋值最后一列，做为新的R类型		
		while (<DatF>)
		{   
		    chomp;
		    next if /^#/;  # Allow embedded comments.
		    if ($linenum==0)
		    {
			$linenum++;
			$fstnewRstr =$_;
			$fstnewRstr=~s/^\s+//;
			@fstRstr=split(/\s+/,$fstnewRstr);
			($fstRvalue)=($fstRstr[9]=~/^R\=(\d+)$/);
			$fstdirect=$fstRstr[8];
			$fstRkey=gethashkey($fstRvalue,%allsameRs);
			#$fstrvRkey=gethashkey($fstRvalue,%reversestr);
			$fstrvRkey=getpositiveR($fstdirect,$fstRvalue,%reversestr);
			
			if ($fstrvRkey>=0)
			{
			     $fstRstr[15]=$fstrvRkey;
			}
			elsif ($fstRkey>=0)
			{
			    $fstRstr[15]=$fstRkey;	
			}
			else
			{$fstRstr[15]= $fstRvalue;}
			
			push @line,[@fstRstr];     
		    }
		    else
		    {
			$sndnewRstr =$_;
			$sndnewRstr=~s/^\s+//;
			@sndRstr=split(/\s+/,$sndnewRstr);
			$linenum++;
			($sndRvalue)=($sndRstr[9]=~/^R\=(\d+)$/);
			$snddirect=$sndRstr[8];
			$sndRkey=gethashkey($sndRvalue,%allsameRs);
			#$srvRkey=gethashkey($sndRvalue,%reversestr);
			if ($sndRkey>=0)
			{$sndrvRkey=getpositiveR($snddirect,$sndRkey,%reversestr);}
			else
			{$sndrvRkey=getpositiveR($snddirect,$sndRvalue,%reversestr);}
			
			## 判断两行的R类型是否一致或类似或反向互补
			if ((($fstRvalue==$sndRvalue)&&( $fstdirect == $snddirect))||((gethashkey($fstRvalue,%allsameRs)>=0)&&(gethashkey($fstRvalue,%allsameRs)==gethashkey($sndRvalue,%allsameRs))&&( $fstdirect == $snddirect))||(($fstdirect != $snddirect)&&(($sndrvRkey>0)&&(($sndrvRkey==$fstrvRkey)))))
			{
				$sndRstr[15]=$fstRstr[15];
			}
			elsif ($sndrvRkey>=0)
			{
			     $sndRstr[15]=$sndrvRkey;
			}
			elsif ($sndRkey>=0)
			{
			    $sndRstr[15]=$sndRkey;	
			}
			else
			{$sndRstr[15]= $sndRvalue;}
			
			push @line,[@sndRstr];
		        @fstRstr=@sndRstr;
		        ($fstRvalue)=($fstRstr[9]=~/^R\=(\d+)$/);
		        $fstdirect=$fstRstr[8];
		        $fstRkey=gethashkey($fstRvalue,%allsameRs);
			#$fstrvRkey=gethashkey($fstRvalue,%reversestr);
			#$fstrvRkey=getpositiveR($fstdirect,$fstRvalue,%reversestr);		
			if ($fstRkey>=0)
			{$fstrvRkey=getpositiveR($snddirect,$fstRkey,%reversestr);}
			else
			{$fstrvRkey=getpositiveR($fstdirect,$fstRvalue,%reversestr);}
			
		    }
		   	    
		   # print @line;
		}
		close(DatF);											
		
		#foreach(sort{($a->[9]=~/^R\=(\d+)$/) <=> ($b->[9]=~/^R\=(\d+)$/) or $a->[5] <=> $b->[5]}@line)
		foreach(sort{$a->[15] cmp $b->[15] or $a->[5] <=> $b->[5]}@line)
		{  
		    print SameData join("\t",@$_)."\n";
		}
		close(SameData);
		
		#analyse the R gap distance
		$sameRLendata="$samelendata.sameRLen";
		#$sameseqstr="$output.seqstr";
		open(SameData,"$samelendata")||die "$!\n";
		open(SameLenData,">$sameRLendata")||die "$!\n";
			#open(SameSeqData,">$sameseqstr")||die "$!\n";
			my @firstfields='';
			my @fields='';
			my @annotation='';
			my @allgff='';
			
			my @sameRstring='';
			my @result='';
		#@Rdist;
		#$nameR;
		my $samelennum=0;
		
		push @annotation,"#score","perc div","del","ins.","sequence","Start","End","left","Strand","repeat","class/family","begin","end","left","ID","other","dist","identify\n";
		print SameLenData join("\t", @annotation);
		
		$alloutput= "$samelendata.gff";
		open(alloutF,">$alloutput")||die "$!\n";
		push @allgff,"#SeqID","Source","Type","Start","End","Score","Strand","Phase","Attributes\n";
		print alloutF join("\t", @allgff);
		
			# $DRout="$Infile.sameRLen.dr";
			# open(DRoutF,">$DRout")||die "$!\n";
			# push @DRstr,"sequenceName","Start","End","oldID","ID","SPdist","CRISPRnum","Score","DRstring","SPstring\n";
			# print DRoutF join("\t", @DRstr);
			
			# $SPout= "$Infile.sameRLen.sp";
			# open(spoutF,">$SPout")||die "$!\n";
			# #push @spgff,"#SeqenceID","Source","Type","Start","End","Score","Strand","Phase","Rtype","CRISPRnum","SPnum","spLen","isame","SPstring\n";
			# #print spoutF join("\t", @spgff);
			
			# $gffoutput= "$Infile.sameRLen.gff3";
			# open(crisproutF,">$gffoutput")||die "$!\n";
			# push @lastgff,"#SeqID","Source","Type","Start","End","Score","Strand","Phase","Attributes\n";
			# print crisproutF join("\t", @lastgff);
			
		while (<SameData>)
		{
		    $dist=-1;
		    chomp;
		    next if /^#/;  # Allow embedded comments.
		    $body =$_;
		    #$body=~s/^\n//;
		    
		    #out.gff3 files
		    $body=~s/^\s+//g;
		    my @bodystring =split(/\s+/,$body);
		    if ($bodystring[0]=='')
		    {
			shift(@fields);
		    }		
		    $gffall[0]= $bodystring[4];
		    $gffall[1]= "LibRepeatMasker";
		    $gffall[2]= "Repeats";
		    $gffall[3]= $bodystring[5]; 
		    $gffall[4]= $bodystring[6];
		    $gffall[5]= $bodystring[0];
		    $gffall[6]= $bodystring[8];
		    if ($bodystring[8] eq "C")
		      { $gffall[6]="-";}    
		    $gffall[7]= "."; 
		    my ($nameR)=($bodystring[9]=~/^R\=(\d+)$/);
		    $length=$bodystring[6]-$bodystring[5];
		    my $idstr="ID=".$bodystring[14].";"."Name=".$nameR.";"."Len=".$length;
		    $gffall[8]= $idstr;    
		    print alloutF join("\t", @gffall)."\n";      
		    
		    if ($samelennum <1)
		    {
			#@sameRstring='';
			$samelennum++;
			@firstfields =split(/\s+/,$body);
			$firstfields[16]=$dist;
			$oneRstring= join("\t", @firstfields);
			
			push @sameRstring,$oneRstring;
		    }
		    else
		    {
			@fields =split(/\s+/,$body);
			##the same R
					if (($firstfields[15] eq $fields[15] )&&($firstfields[4] eq $fields[4]))
			{    	    
			    $dist=$fields[5]-$firstfields[6];
			    $fields[16]=$dist;
			    #push @fields, $dist;
			    ($nameR)=($fields[15]=~/^R\=(\d+)$/);	    
			    $samelennum++;
			    
			    if ( $fields[16]>=-1)
			    {
				$oneRstring= join("\t", @fields);
				push @sameRstring,$oneRstring;
				@firstfields=@fields;
			    }
			    ##else???存在两行记录位置有交叉区，如何取舍其中一个？与DR相似性高的一行留下
			    
			}						
			else
			{
			    #the same R can be dealed with in here :$Rdist[$nameR]
			    #first to sort,then count the distance
			    $Rnum=@sameRstring;
			    if ($Rnum>1)
			    {
				##add code to deal with the interval R value 
				@result=findnearnum(@sameRstring);
				print SameLenData @result;				
				#count out the spaces and CRISPR gff3 files
				my($DR_ref,$space_ref,$crispr_ref)=&findspaceCRISPR(@result);
				@DRoneR=@$DR_ref;
				@SPoneR=@$space_ref;
				@CRISPRoneR=@$crispr_ref;
				print DRoutF @DRoneR;
				print spoutF @SPoneR;
				print crisproutF @CRISPRoneR;
						
			    }
			    
			    @sameRstring=();
			    @result=();
			    $samelennum=1;  
			    $dist=-1;
			    @firstfields=@fields;
			    $firstfields[16]=$dist;
			    $oneRstring= join("\t", @firstfields);
			    push @sameRstring,$oneRstring;	    
			}			
		    }   	 
		}
			#the last R type strings 
			@result=findnearnum(@sameRstring);
			print SameLenData @result;
				
			#count out the spaces and CRISPR gff3 files
			my($DR_ref,$space_ref,$crispr_ref)=&findspaceCRISPR(@result);
			@DRoneR=@$DR_ref;
			@SPoneR=@$space_ref;
			@CRISPRoneR=@$crispr_ref;
			print DRoutF @DRoneR;
			print spoutF @SPoneR;
			print crisproutF @CRISPRoneR;		
			
		close(SameData);
		close(SameLenData);
		#close(SameSeqData);
			# close (DRoutF);
			# close(spoutF);
			# close(crisproutF);
		close(alloutF);
	#}
}

###更新gff文件内容将存在overlap的合并（初期R类型可能不一致，但实际比较类似，需要将不同类型R存在overlap的CRISPRs合并）
system("rm  -rf $tmpdir \n");
system("rm -f $Input\.* $Infile\_* $Infile\.stg* $Infile\.freq");
system("rm  -f formatdb.log \n");

#$output= "$Input.out";
#$dataF="$output.gff3.data";
#$samelendata="$dataF.samelen";
#$sameRLendata="$samelendata.sameRLen";

##$gffoutput= "$sameRLendata.gff3";


##deletefalsegff($gffoutput);
	system("rm -f $Input");
}
close (DRoutF);
close(spoutF);
close(crisproutF);


$EndTime=time();
$AllTime=$EndTime-$StarTime;
print "The program  one time use the all time minutes:$AllTime s \n ";



#--------------------------------------------------------------------------------------------------
sub usuage {
    print "\nHi, need some help?\n";
    print STDERR <<"    _EOT_";

    Usage :CRISPRdigger.pl <options> <specification> <default>
     
	\$Inputfile   = defined $opt_i ? $opt_i : "";
	\$SPMinDRratio=defined $opt_l ? $opt_l : 0.5;
	\$SPMaxDRratio=defined $opt_m ? $opt_m : 3.0;
	\$MinSP     = defined $opt_g ? $opt_g : 10;
	\$MaxSP     = defined $opt_p ? $opt_p : 120;
	\$DiffTempSim = defined $opt_t ? $opt_t : 0.80;#0.7-0.9
	\$SPsimilary = defined $opt_s ? $opt_s : 0.50;#0.4-0.8
	\$SameRratio = defined $opt_R ? $opt_R : 0.50;#0.4-0.8
	\$blastsameR = defined $opt_B ? $opt_B : 0.35;#0.2-0.7
	\$MaxDR     = defined $opt_d ? $opt_d: 47 ; 

    _EOT_
    exit(1)
}

#----------------------------------------------------------------------------------------------

##delete false gff lines
sub deletefalsegff
{
    my ($infile)=@_;
    my (@gfflines,@tailgff,$GDNCconts);
    open(AF, "$infile")||die "can't open this $infile\n";
    my $tempfile="$infile.temp";
    open (TF,">$tempfile")||  die "Cannot create file:$tempfile\n"; 
    my ($NCname)=($infile=~/(NC_\d*)./);
    
    while(<AF>)
    {	
	#chomp;
	my $tline=$_;
	$tline=~s/[\r\n]+//g;
	@gfflines=split(/\t/,$tline);
	@tailgff=split(/;/,$gfflines[8]);
	#($NCcont)=($gfflines[0]=~/(NC_*)./);
	
	if (($gfflines[2] eq "CRISPR")&&( $tailgff[$#tailgff] eq "isCRISPR=1"))
	{	
		
	   my $notinGD=0;
	   my $questionCRISPR=0;
	   $begloc=$gfflines[3];
	   $endloc=$gfflines[4];
	   
	  # $notinGD=notinGDdb($begloc,$endloc,$gfflines[8],$GDNCconts);
		##compare DRs,SPs similary to modify isCRISPR=1 or isCRISPR=2
		
		$drfilename="$NCname";
		$questionCRISPR=cmpDRSP($tline,$drfilename);
		
		if ($questionCRISPR eq 1)
		{		
		  $_=$tline;
		  s/isCRISPR=1/isCRISPR=3/;
		  $tline=$_;
		}
		##
		print TF "$tline \n";
	}
	else
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

###


##delete false gff lines referring the GoldenDB
sub old_deletefalsegff
{
    my ($infile)=@_;
    my (@gfflines,@tailgff,$GDNCconts);
    open(AF, "$infile")||die "can't open this $infile\n";
    my $tempfile="$infile.temp";
    open (TF,">$tempfile")||  die "Cannot create file:$tempfile\n"; 
    my ($NCname)=($infile=~/(NC_\d*)./);
    
    ##my $GDNCconts=NClinesGD($NCname);
     $GDdbname="GoldenCRISPRs";
     $dir="./../";
     $GDNCconts= findNCline($dir,$GDdbname,$NCname);
    
    while(<AF>)
    {	
	#chomp;
	my $tline=$_;
	$tline=~s/[\r\n]+//g;
	@gfflines=split(/\t/,$tline);
	@tailgff=split(/;/,$gfflines[8]);
	#($NCcont)=($gfflines[0]=~/(NC_*)./);
	
	if (($gfflines[2] eq "CRISPR")&&( $tailgff[$#tailgff] eq "isCRISPR=1"))
	{
	   my $notinGD=0;
	   my $questionCRISPR=0;
	   $begloc=$gfflines[3];
	   $endloc=$gfflines[4];
	   
	   $notinGD=notinGDdb($begloc,$endloc,$gfflines[8],$GDNCconts);
	   
	   if( $notinGD eq 1)
	   {
		##compare DRs,SPs similary to modify isCRISPR=1 or isCRISPR=2
		
		$drfilename="$NCname";	
		$questionCRISPR=cmpDRSP($tline,$drfilename);
		print "$$$$$$$$$$ questionCRISPR ******* $questionCRISPR*****\n";
		
		if ($questionCRISPR eq 1)
		{		
		  $_=$tline;
		  s/isCRISPR=1/isCRISPR=3/;
		  $tline=$_;
		}
		##
		print TF "$tline \n";
	   }
	   else
	   {
		print TF "$tline \n";
	   }
	}
	else
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

###
sub find_DRSPfile
{
  local($dir,$DRfilelike,$locbeg,$locend) = @_;
  opendir(DIR,"$dir"|| die "can't open this $dir");
  local @files =readdir(DIR);
  closedir(DIR);
  for $file (@files){
    next if($file=~m/\.$/ || $file =~m/\.\.$/);
    if ($file =~/($DRfilelike.*?\.dr)/i)
    {
     
      $curfiledir="$dir/$file";
     
      ($DRsumAT,$DRsumGC,$SPsumAT,$SPsumGC,$DRstrs_ref,$SPstrs_ref)=readDRSP($curfiledir,$locbeg,$locend);
     
    }
    elsif(-d "$dir/$file"){
            find_DRSPfile("$dir/$file",$DRfilelike,$locbeg,$locend );
    }    
  }
  return ($DRsumAT,$DRsumGC,$SPsumAT,$SPsumGC,$DRstrs_ref,$SPstrs_ref); 
}

sub readDRSP
{
     my ($DRpath,$locbeg,$locend)=@_;
     my (@DRstr,@SPstr);
     $DRsumAT=0;
     $DRsumGC=0;
     $SPsumAT=0;
     $SPsumGC=0;
     
     open(DrF, "$DRpath")||die "can't open this $DRpath\n";
     while (<DrF>)
     {
	my $DRline=$_;
	$DRline=~s/[\r\n]+//g;
	@DRSPlines=split(/\t/,$DRline);
	
	if (($DRSPlines[1] >= $locbeg)&&($DRSPlines[2] <= $locend))
	{
	  my $DRsnote=$DRSPlines[1]."_".$DRSPlines[2];
	  my $SPsnote=$DRSPlines[2]."_".$DRSPlines[1];
	  my $DRstr=$DRSPlines[8];
	  my $SPstr=$DRSPlines[9];
	  
	  ($DRATpercent,$DRGCpercent)= ATGCPercentage($DRstr);
	  ($SPATpercent,$SPGCpercent)= ATGCPercentage($SPstr);
	  $DRsumAT+=$DRATpercent;
	  $DRsumGC+=$DRGCpercent;
	  $SPsumAT+=$SPATpercent;
	  $SPsumGC+=$SPGCpercent;
	  
	  #my $seq_DRobj=Bio::Seq->new(-seq => "$DRstr", -display_id => "$DRSPlines[0]", -desc => "$DRsnote" );
	  #push @DRstr,$seq_DRobj."\n";
	  #
	  #my $seq_SPobj=Bio::Seq->new(-seq => "$SPstr", -display_id => "$DRSPlines[0]", -desc => "$SPsnote" );
	  #push @SPstr,$seq_SPobj."\n";
	  #
	  push @DRstr, ">".$DRSPlines[0].$DRsnote."\n";
	  push @DRstr, $DRstr."\n";
	  
	  push @SPstr, ">".$DRSPlines[0].$SPsnote."\n";
	  push @SPstr, $SPstr."\n";   
	}
     }
     pop @SPstr;
     pop @SPstr;
     close(DrF);
     return ($DRsumAT,$DRsumGC,$SPsumAT,$SPsumGC,\@DRstr,\@SPstr);	
}

###compare DRs,SPs similiary in *.dr files
##questionable CRISPR return 1,else return 0
sub cmpDRSP
{
  my ($tline,$DRfilelike)=@_;
  
  $dir="./";
  
  my @gfflines=split(/\t/,$tline);
  my $CRbeg=$gfflines[3];
  my $CRend=$gfflines[4];
  my (@DRstrs,@SPstrs,$DRsim,$SPsim);
  my ($SPnum)=($tline=~/SPnum=(\d+);/);
  my $DRnum=$SPnum+1;
  my $DRaveAT=0;
  my $DRaveGC=0;
  my $SPaveAT=0;
  my $SPaveGC=0;
  my $DRsumAT=0;
  my $DRsumGC=0;
  my $SPsumAT=0;
  my $SPsumGC=0;

   ($DRsumAT,$DRsumGC,$SPsumAT,$SPsumGC,$DRstrs_ref,$SPstrs_ref)=find_DRSPfile($dir,$DRfilelike,$CRbeg,$CRend);
  
     @DRstrs=@$DRstrs_ref;
     @SPstrs=@$SPstrs_ref;
     ###$DRsim (0-100)
     $DRfile=$gfflines[3]."_".$CRbeg."_".$CRend;
     $SPfile=$gfflines[3]."_".$CRend."_".$CRbeg;
     open DRF,"> $DRfile" or die "can't open this $DRfile\n";
     open SPF,"> $SPfile" or die "can't open this $SPfile\n";
     print DRF @DRstrs;
     print SPF @SPstrs;
     close DRF;
     close SPF;
     
     $DRsim=checkDRSPAlign($DRfile);
     $SPsim=checkDRSPAlign($SPfile);
      
     print " !!!!!!!!!!!!!!DRsim;;SPsim  _________ $DRsim;;$SPsim \n";
     
     $DRaveAT= sprintf("%0.3f" ,$DRsumAT/$DRnum);
     $DRaveGC= sprintf("%0.3f" ,$DRsumGC/$DRnum);
     
     $SPaveAT= sprintf("%0.3f" ,$SPsumAT/($DRnum-1));
     $SPaveGC= sprintf("%0.3f" ,$SPsumGC/($DRnum-1));
     
     print " #######DRaveAT;DRaveGC;;;SPaveAT;SPaveGC  @@@@@@@ $DRaveAT;;$DRaveGC ;;$SPaveAT;;$SPaveGC\n";

     `rm -r $DRfile `;
     `rm -r $SPfile `;
     ############
     ###########
     ########
     
     ##(!(($SPsim<$SPsimilary)&&($DRsim>$DiffTempSim)))
    if ((($DRsim < 40)&&($DRsim > 0)&&($DRnum < 10))||(($SPsim > 60)&&($SPsim <= 100)&&($DRnum < 10))||(($DRsim < (100*$blastsameR))&&($DRsim > 0))
	||((($SPsim>(100*$SPsimilary))||($DRsim<=(100*$SameRratio)))&&($DRaveAT > 0.7)&&($DRaveAT <= 1)&&($SPaveAT > 0.7)&&($SPaveAT <= 1))
	||((($SPsim>(100*$SPsimilary))||($DRsim<=(100*$SameRratio)))&&($DRaveGC > 0.7)&&($DRaveGC <= 1)&&($SPaveGC > 0.7)&&($SPaveGC <= 1))
	||(($DRaveGC > 0.75)&&($DRaveGC <= 1)&&($SPaveGC > 0.75)&&($SPaveGC <= 1))||(($DRaveAT > 0.75)&&($DRaveAT <= 1)&&($SPaveAT > 0.75)&&($SPaveAT <= 1)))
     {
	return 1;
     }
     else
     {
	return 0;
     }
  		
}

###Calulate the AT/GC percentage
sub ATGCPercentage
{
        my ($Seq) = @_;
        return 0 if not $Seq;
        $Seq =~ tr/atcg/ATCG/;
        my @aSplitSeq = split(//,$Seq);
        my %hCountBase = ( "G" => 0, "T" => 0,
                       "A" => 0, "C" => 0 );   
        # Count and for cal for base frequence
        $hCountBase{ $_ }++ foreach ( @aSplitSeq );
        $GCPercentage = ( $hCountBase{ G } + $hCountBase{ C } ) / length( $Seq  );
	$ATPercentage = ( $hCountBase{ A } + $hCountBase{ T } ) / length( $Seq  );
       
        return ($ATPercentage,$GCPercentage);
}

##short CRISPRs DR compare to define is or questionable CRISPR
###Fasta formate sequence
sub checkDRSPAlign
{
 my($file) = @_;
 #my(@seq_array) = @_;
 #$seq_array_ref = \@seq_array;
 my $ident=0;
use Bio::SimpleAlign;

 eval
 {
 	#$ENV{CLUSTALDIR} = '/home/username/clustalw1.8/';
 	#use Bio::Tools::Run::Alignment::Clustalw;
	#use Bio::AlignIO::clustalw;
	#open STDOUT, "|tee stdout >/dev/null 2>&1";
	my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM');
	
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
	
	#my $factory = Bio::AlignIO::clustalw->new(@params);
      #my $factory = clustalw2->new(@params);

     	#  Pass the factory a list of sequences to be aligned.
	
	# my $aln = Bio::SimpleAlign->new();
	
      my $aln = $factory->align($file); # $aln is a SimpleAlign object.
     
      ##my $aln = $factory->align($seq_array_ref); # $aln is a SimpleAlign object.
  	$ident = $aln->percentage_identity;

 };
	##DRs' similiary 
	##if($ident >= 50) {return 1;} else {return 0;}
	return $ident;
}


## get contents about some NCname in GDdb (GoldenCRISPRs)
sub NClinesGD
{
  my ($NCname)=@_;
  ##my (@NClinGD);
  $GDdbname="GoldenCRISPRs";
  my ($findGDline,$GDpath);
  
  $dir="./../../";
  $GDpath=findpath($dir,$GDdbname); 
  print "&&&&&&&&& GDpath&&&&&&&&&& $GDpath \n";
  
  if (length($GDpath)>0)
  {
     open(GdF, "$GDpath")||die "can't open this $GDpath\n";
     while (<GdF>)
     {
	my $GDline=$_;
	$GDline=~s/[\r\n]+//g;
	@GDCRlines=split(/\t/,$GDline);
	if ($GDCRlines[1] eq $NCname)
	{
	  ## push @NClinGD,$GDline."\n";
	  $findGDline=$GDline;	  
	  return $findGDline;
	}
     }
     close(GdF);	    
  }
  eles
  {
    return $findGDline;
    print "$GDdbname does not exist!";
  }		
}

###find the file path from dir:./../../
sub findNCline
{
  my ($dirpath,$filelike,$NCname)=@_;
  my($filepath,$NCline);
  ##$dir ="./../../";
  
   opendir(DIR,"$dirpath"|| die "can't open this $dirpath");
  local @files =readdir(DIR);
  closedir(DIR);
  for $file (@files)
  {
    next if($file=~m/\.$/ || $file =~m/\.\.$/);
    if ($file =~/($filelike)/i)
    {         
      # $filepath=`pwd`;
       #return $filepath;
       $filepath="$dirpath/$file";
       
       $GDNCconts=getGDNC($filepath,$NCname);
       return ($GDNCconts);
    }
    elsif(-d "$dir/$file"){
            findNCline("$dirpath/$file",$filelike,$NCname );
    }      
  }
   return ($GDNCconts);
  #return $dirpath;
}

###
sub getGDNC
{
  my ($filename,$NCname)=@_;
  my ($findGDline,$GDline);
  open(GdF, "$filename")||die "can't open this $filename\n";
  while (<GdF>)
 { 
   $GDline=$_;
   $GDline=~s/[\r\n]+//g;
   @GDCRlines=split(/\t/,$GDline);
   if ($GDCRlines[1] eq $NCname)
   {
     ## push @NClinGD,$GDline."\n";
     $findGDline=$GDline;
     
     return ($findGDline);
   }
 }
 close(GdF);	    
 
 #return ($findGDline);
 print "$GDdbname does not exist in GDdb!";	
	
}


###some gff3 line exist in or not in GDdatabase
###exist return 0; no exist return 1
##param:($NCcont,$begloc,$endloc)
###line in gff::  gi|407472532|ref|NC_018664.1| LibRepeatMasker CRISPR	704762	706990	222 - . ID=704762;Name=0;CRISPRnum=0;StrLen=2229;aSPLen=36;SPnum=33;aDRLen=30;isCRISPR=1
###line in GDstand::  Acidimicrobium ferrooxidans DSM 10331 NC_013124 CRISPRGolden	2 |996601 998276	27|1000880 1003348	40
sub notinGDdb
{
  my ($gffbeg,$gffend,$oneCRinfor,$NClinesGD)= @_;
  my $noinGD=1;
  my @GDeachCR=split(/\|/,$NClinesGD);
  my $CRnum=$#GDeachCR;
  ##look for GDdb 	
  my @CRinfor=split(/;/,$oneCRinfor);
  my $aSPlen=$CRinfor[4];
  my $SPnum =$CRinfor[5];
  my $aDRlen=$CRinfor[6];
  my $midgffloc=($gffbeg+$gffend)/2;
  my $leftgffloc = $gffbeg+2*$aDRlen+$aSPlen;
  my $rightgffloc= $gffend-2*$aDRlen-$aSPlen;
   
  for (my $CRloop=1;$CRloop<=$CRnum;$CRloop++)
  {
    my @eachCR=split(/\t/,$GDeachCR[$CRloop]);
    my $GDbeg=$eachCR[0];
    my $GDend=$eachCR[1];
    if ((($midgffloc>$GDbeg)&&($midgffloc<$GDend))||(($leftgffloc>$GDbeg)&&($leftgffloc<$GDend))||(($rightgffloc>$GDbeg)&&($rightgffloc<$GDend)))
    {
	$noinGD=0;
	return $noinGD;
    }	
  }
  return $noinGD;	
}


#find the nearest data sets in one group R (include all lines)(include add lost DR)
sub findnearnum
{
    my (@Rstring)=@_;
    $Rstrnum=@Rstring;
    $presentnum=1;
    $samenum=0;
    $countsame=1;   
    @resultstr=();
    while($presentnum<$Rstrnum-1)
    {
	#the first element[16]=-1,so  begin in $presentnum=1
	@beforestr=split /\t/,$Rstring[$presentnum-1];
	@onestr=split /\t/,$Rstring[$presentnum];
	$oneDRlen=$onestr[6]-$onestr[5];
	$DRSPlen=$oneDRlen+$onestr[16] ;
	@anotherstr=split /\t/,$Rstring[$presentnum+1];
	$beforestr[17]=0;
	$onestr[17]=0;
	$anotherstr[17]=0;
	$isflag=0;
	
	#while((((abs($onestr[16]/$anotherstr[16]-1)<0.2)||(abs($onestr[16]-$anotherstr[16]))<10))&&($Rstrnum-1>$presentnum++)&&($onestr[16]>0))
	#while((((abs($onestr[16]/$anotherstr[16]-1)<0.2)||(abs($onestr[16]-$anotherstr[16]))<10))&&($onestr[16]>0))
	while( ($Rstrnum-1>$presentnum++)&&($anotherstr[16]!=0)&&($onestr[16]!=0))	
	{
		if (($onestr[16]<$MaxSP)&&(($anotherstr[16])>0.5*$oneDRlen||($onestr[16])>0.5*$oneDRlen)
		    &&(($anotherstr[16])<2.5*$oneDRlen||($onestr[16])<2.5*$oneDRlen)&&
		((abs($onestr[16]/$anotherstr[16]-1)<0.2)||(abs($onestr[16]-$anotherstr[16])<15))&&($onestr[16]>$MinSP)&&($anotherstr[16]>$MinSP))
		{
		      
		  if ($anotherstr[16]<0)
		  {
		      @anotherstr=split /\t/,$Rstring[$presentnum+1];
		      next;
		  }	    
		  $isflag=1;	
		  if ($samenum==0)
		  {
		      $beforestr[17]=$countsame;
		      $onestr[17]=$countsame;
		      $anotherstr[17]=$countsame;
		      #update the DR SP length
		      $oneDRlen=$onestr[6]-$onestr[5];
		      $DRSPlen=$oneDRlen+$onestr[16] ;
		      #$countsame++;
		      $befstr= join ("\t",@beforestr)."\n";
		      push @resultstr,$befstr;
		  }
		  else
		  {
		      $onestr[17]=$countsame;
		      $anotherstr[17]=$countsame;
		      #update the DR SP length
		      $oneDRlen=$onestr[6]-$onestr[5];
		      $DRSPlen=$oneDRlen+$onestr[16] ;
		      #$countsame++;
		  }
		  
		  $samenum++;
		  $firstr= join ("\t",@onestr)."\n";
		  push @resultstr,$firstr;
		  @beforestr=split /\t/,$Rstring[$presentnum-1];
		  @onestr=@anotherstr;
		  @anotherstr=split /\t/,$Rstring[$presentnum+1];
		}
		##interval one DR condition
		elsif(($onestr[16]<$MaxSP)&&($anotherstr[16]>$DRSPlen)&&($anotherstr[16]<2*$DRSPlen))
		{
		    $isflag=1;
		    #中间存在间隔DR情况处理:@onestr__gap__@anotherstr:四舍五入
		      $beginlocDR=int(($anotherstr[5]-$onestr[5])/2+0.5)+$onestr[5];
		      $DRlen=int(($onestr[6]-$onestr[5]+$anotherstr[6]-$anotherstr[5])/2+0.5+1);
		      $anlaystr=substr($bodystr,$beginlocDR-1,$DRlen);
		      $DRstring=substr($bodystr,$onestr[5]-1,$onestr[6]-$onestr[5]+1);
		      
		      ($astr,$fstr,$samelen)=TirSimilar($DRstring,$anlaystr);
		      ($cmp_result,$same_similary)=Str_Compare($astr,$fstr);
		      if ($same_similary> $SameRratio )
		      {
			#add new DR into @resultstr
			@newDRcontents=();
			@newDRcontents=@onestr;
			$newDRcontents[5]=$beginlocDR;
			$newDRcontents[6]=$beginlocDR+$DRlen;
			$newDRcontents[16]=$newDRcontents[5]-$onestr[6];
			$newDRcontents[17]=$countsame;
			$beforestr[17]=$countsame;
			$onestr[17]=$countsame;
		        $anotherstr[17]=$countsame;
			$anotherstr[16]=$anotherstr[5]-$newDRcontents[6];
			if ($beforestr[16]== -1)
			{
			    $beforestr[17]=$countsame;
			    $befstr=join ("\t",@beforestr)."\n";
			    push @resultstr,$befstr;
			    $samenum++;
			}		
			$firstr= join ("\t",@onestr)."\n";
			$newstr= join ("\t",@newDRcontents)."\n";	
			push @resultstr,$firstr;
			push @resultstr,$newstr;
			$samenum++;
			$samenum++;
			@beforestr=@newDRcontents;
		        @onestr=@anotherstr;
		        @anotherstr=split /\t/,$Rstring[$presentnum+1];						
		      }
		      else
		      {		
			last;
		      }
		}
		##@beforestr is the first string, condition:@beforestr__gap_____@onestr 
		elsif(($anotherstr[16]<$MaxSP)&&($onestr[16]>$DRSPlen)&&($onestr[16]<2*$DRSPlen))
		{
		      $isflag=1;
		      #中间存在间隔DR情况处理:@beforestr__gap_____@onestr 
		      $beginlocDR=int(($onestr[5]-$beforestr[5])/2+0.5)+$beforestr[5];
		      $DRlen=int(($onestr[6]-$onestr[5]+$anotherstr[6]-$anotherstr[5])/2+0.5+1);
		      $anlaystr=substr($bodystr,$beginlocDR-1,$DRlen);
		      $DRstring=substr($bodystr,$onestr[5]-1,$onestr[6]-$onestr[5]+1);
		      
		      ($astr,$fstr,$samelen)=TirSimilar($DRstring,$anlaystr);
		      ($cmp_result,$same_similary)=Str_Compare($astr,$fstr);
		      if ($same_similary> $SameRratio )
		      {
			#add new DR into @resultstr
			@newDRcontents=();
			@newDRcontents=@onestr;
			$newDRcontents[5]=$beginlocDR;
			$newDRcontents[6]=$beginlocDR+$DRlen;
			$newDRcontents[16]=$newDRcontents[5]-$beforestr[6];
			$newDRcontents[17]=$countsame;
			$beforestr[17]=$countsame;
			$onestr[17]=$countsame;
		        $anotherstr[17]=$countsame;
			$onestr[16]=$onestr[5]-$newDRcontents[6];
			if ($beforestr[16]== -1)
			{
			    $beforestr[17]=$countsame;
			    $befstr=join ("\t",@beforestr)."\n";
			    push @resultstr,$befstr;
			    $samenum++;
			}
							
			$newstr= join ("\t",@newDRcontents)."\n";
			$firstr= join ("\t",@onestr)."\n";			
			push @resultstr,$newstr;
			push @resultstr,$firstr;
			$samenum++;
			$samenum++;						
			@beforestr=@onestr;
		        @onestr=@anotherstr;
		        @anotherstr=split /\t/,$Rstring[$presentnum+1];
		      }
		      else
		      {
			last;
		      }
		}	
		else
		{last;}
	}
    		
	##the other probable crispr
	if((($isflag eq 0)&&($presentnum==2))||(($isflag eq 0)&&($beforestr[16]>0)))
	#if($isflag eq 0)
	{
	    $befstr= join ("\t",@beforestr)."\n";
	    push @resultstr,$befstr;
	    #$presentnum++;	    
	}
	else
	{
	    if ($onestr[16]>0)
	    {
		    $countsame++;	
		    $firstr= join ("\t",@onestr)."\n";
		    push @resultstr,$firstr;
	    }
	    $presentnum++;
	}    	      	    	    
	$samenum=0;
	#$presentnum++;
    }
        #the last two lines ?????????????????
	if(($presentnum==$Rstrnum-1)||($presentnum==$Rstrnum))
	#if($presentnum==$Rstrnum)
	{
	    ##@temptstr=@onestr;
	    #@onestr=@anotherstr;
	    #@anotherstr=split /\t/,$Rstring[$presentnum];
	    #$onestr[17]=0;
	    @beforestr=split /\t/,$Rstring[$Rstrnum-3];
	    @onestr=split /\t/,$Rstring[$Rstrnum-2];
	    @anotherstr=split /\t/,$Rstring[$Rstrnum-1];
	    if (!($onestr[17]))
	    {
		    $onestr[17]=0;
	    }		
	    $anotherstr[17]=0; 		
	    if(($onestr[16]>2)&&($onestr[16]<$MaxSP)&&((abs(($onestr[16]/$anotherstr[16])-1)<0.2)||((abs($onestr[16]-$anotherstr[16])<15))))
	    {
		$onestr[17]=$countsame;
		$anotherstr[17]=$countsame;   		
	    }
		    
	    if(($isflag eq 0)||($presentnum==$Rstrnum-1))
	    {
		    $firstr= join ("\t",@onestr)."\n";
		    push @resultstr,$firstr;
	    }
	     $anostr= join ("\t",@anotherstr)."\n";
	     push @resultstr,$anostr;
	}	      	
    return  @resultstr;
}



#find the nearest data sets in one group R (include all lines)
#sub findnearnumold
#{
#    my (@Rstring)=@_;
#    $Rstrnum=@Rstring;
#    $presentnum=1;
#    $samenum=0;
#    $countsame=1;
#    
#    @resultstr=();
#    while($presentnum<$Rstrnum-1)
#    {
#	#the first element[16]=-1,so  begin in $presentnum=1
#	@beforestr=split /\t/,$Rstring[$presentnum-1];
#	@onestr=split /\t/,$Rstring[$presentnum];
#	@anotherstr=split /\t/,$Rstring[$presentnum+1];
#	$beforestr[17]=0;
#	$onestr[17]=0;
#	$anotherstr[17]=0;
#	$isflag=0;	
#	#while((((abs($onestr[16]/$anotherstr[16]-1)<0.2)||(abs($onestr[16]-$anotherstr[16]))<10))&&($Rstrnum-1>$presentnum++)&&($onestr[16]>0))
#	while( ($anotherstr[16]!=0)&&($onestr[16]!=0)&&($onestr[16]<$MaxSP*1.2)&&((abs($onestr[16]/$anotherstr[16]-1)<0.2)||(abs($onestr[16]-$anotherstr[16])<15)) &&($onestr[16]>10)&&($Rstrnum-1>$presentnum++))
#	#while((((abs($onestr[16]/$anotherstr[16]-1)<0.2)||(abs($onestr[16]-$anotherstr[16]))<10))&&($onestr[16]>0))
#	{
#	    if ($anotherstr[16]<0)
#	    {
#		@anotherstr=split /\t/,$Rstring[$presentnum+2];
#		next;
#	    }
#	    
#	    $isflag=1;	
#	    if ($samenum==0)
#	    {
#		$beforestr[17]=$countsame;
#		$onestr[17]=$countsame;
#		$anotherstr[17]=$countsame;
#		#$countsame++;
#		$befstr= join ("\t",@beforestr)."\n";
#		push @resultstr,$befstr;
#	    }
#	    else
#	    {
#		$onestr[17]=$countsame;
#		$anotherstr[17]=$countsame;
#		#$countsame++;
#	    }
#	    
#	    $samenum++;
#	    $firstr= join ("\t",@onestr)."\n";
#	    push @resultstr,$firstr;
#	    @onestr=@anotherstr;
#	    @anotherstr=split /\t/,$Rstring[$presentnum+1];
#	}
#	
#	#    if($samenum eq 0)
#	#    {
#	#	$befstr= join ("\t",@beforestr)."\n";
#	#	push @resultstr,$befstr;
#	#    }
#	#    else
#	#    {
#	#	$firstr= join ("\t",@onestr)."\n";
#	#        push @resultstr,$firstr;
#	#	#$presentnum++;
#	#    }
#	
#	if((($isflag eq 0)&&($presentnum==1))||(($isflag eq 0)&&($beforestr[16]>0)))
#	{
#	    $befstr= join ("\t",@beforestr)."\n";
#	    push @resultstr,$befstr;
#	}
#	 else
#	{
#	    if ($onestr[16]>0)
#	    {
#		    $countsame++;	
#		    $firstr= join ("\t",@onestr)."\n";
#		    push @resultstr,$firstr;
#	    }
#	    $presentnum++;
#	}
#				
#	$samenum=0;
#	$presentnum++;	    
#    } 
#     #the last two lines ?????????????????
#	    if(($presentnum==$Rstrnum-1)||($presentnum==$Rstrnum))
#	     #if($presentnum==$Rstrnum)
#	    {
#		##@temptstr=@onestr;
#		#@onestr=@anotherstr;
#		#@anotherstr=split /\t/,$Rstring[$presentnum];
#		#$onestr[17]=0;
#		@beforestr=split /\t/,$Rstring[$Rstrnum-3];
#		@onestr=split /\t/,$Rstring[$Rstrnum-2];
#		@anotherstr=split /\t/,$Rstring[$Rstrnum-1];
#		if (!($onestr[17]))
#		{
#			$onestr[17]=0;
#		}		
#		$anotherstr[17]=0; 		
#		if(($onestr[16]>2)&&($onestr[16]<$MaxSP*1.2)&&((abs(($onestr[16]/$anotherstr[16])-1)<0.2)||((abs($onestr[16]-$anotherstr[16])<15))))
#	        {
#		    $onestr[17]=$countsame;
#		    $anotherstr[17]=$countsame;   		
#		}		
#		if(($isflag eq 0)||($presentnum==$Rstrnum-1))
#		{
#			$firstr= join ("\t",@onestr)."\n";
#			push @resultstr,$firstr;
#		}
#		 $anostr= join ("\t",@anotherstr)."\n";
#		 push @resultstr,$anostr;
#	}	      	
#    return  @resultstr;
#}

#find CRISPRs in one kind R 
sub crisprsameR
{
   (@sameRstr)=@_;
   $Rnum=@sameRstr;
   $numstr=0;
   @oneRgff=();
   $CRISPRnum=0;
   
    while ($numstr<$Rnum-1)
    #if (($fields[17]==1)&&($fields[16]<$MaxSP*1.2))
    {       
        $exist=0;
	$spacelen=0;
	$spacenums=0;
	$aveDRlen=0;
	@outgff=();
	
	#($nameR)=($fields[15]=~/^R\=(\d+)$/);
	@fstrs=split (/\n/,$sameRstr[$numstr]);
	$fstr=$fstrs[0];
	#($sndstr)=($sameRstr[$numstr+1]=~/(.+)\n $/);	
	@sndstrs=split (/\n/,$sameRstr[$numstr+1]);
	$sndstr=$sndstrs[0];
	
	@Ffields =split(/\t/,$fstr);
	@Sfields=split(/\t/,$sndstr);
	if ($Ffields[17] ne 0)
	{
		$outgff[0]= $Ffields[4];
		$outgff[1]= "LibRepeatMasker";
		$outgff[2]= "Repeats";
		$outgff[3]= $Ffields[5]; 
		#$outgff[4]= $fields[6];#modify
		$outgff[4]= $Ffields[6];
		$aveDRlen=$Ffields[6]-$Ffields[5]+1;
		$outgff[5]= $Ffields[0];
		$outgff[6]= $Ffields[8];
		if ($Ffields[8] eq "C")
		  { $outgff[6]="-";}      
		$outgff[7]= ".";	
		 $nameR=$Ffields[15];
	       # my ($nameR)=($fields[9]=~/^R\=(\d+)$/);
		$ID=$Ffields[14];		
        }
	while (($Sfields[17] ne 0)&&($Ffields[17] eq $Sfields[17])&&($Sfields[16]<$MaxSP))
	{
		$exist++;
		$numstr++;
		$spacelen+=$Sfields[16]-1;
		$aveDRlen+=$Sfields[6]-$Sfields[5]+1;
		$outgff[4]= $Sfields[6];#modify		
		@sndstrs=split (/\n/,$sameRstr[$numstr+1]);
		$sndstr=$sndstrs[0];		
		@Sfields=split(/\t/,$sndstr);
	}	
	if ($exist ne 0)
	{		
		$StrLen=$outgff[4]-$outgff[3]+1;
		$avespace=int($spacelen/$exist);
		$aveDRlen=int($aveDRlen/($exist+1));
		$spacenums=$exist;
		my $idstr="ID=".$ID.";"."R=".$nameR.";"."CRISPRnum=".$CRISPRnum.";"."StrLen=".$StrLen.";"."aSPLen=".$avespace.";"."SPnum=".$spacenums.";"."aDRLen=".$aveDRlen.";"."isCRISPR=".1;
		#my $idstr="ID=".$ID.";".$nameR.";"."StrLen=".$StrLen;
		$outgff[8]= $idstr;
		$CRISPRnum++;
		push @oneRgff ,join("\t", @outgff)."\n"; 		
	}
	#else
	$numstr++;	
    }
    
	return @oneRgff;	
}

#find spaces in one kind R
sub findspaces
{
   (@sameRstr)=@_;
   $Rnum=@sameRstr;
   $numstr=0;
   @oneRspace=();
   $CRISPRnum=0;
   @outgff=();
   #@onespaces=();
   $simil=0;
   
    while ($numstr<$Rnum-1)
    #if (($fields[17]==1)&&($fields[16]<$MaxSP))
    {       
        $exist=0;
	@oneCRISPR=();
	@onespaces=();
	@spaces=();
	@oddspaces=();
	@evenspaces=();
	
	#($nameR)=($fields[15]=~/^R\=(\d+)$/);
	@fstrs=split (/\n/,$sameRstr[$numstr]);
	$fstr=$fstrs[0];
	#($sndstr)=($sameRstr[$numstr+1]=~/(.+)\n $/);	
	@sndstrs=split (/\n/,$sameRstr[$numstr+1]);
	$sndstr=$sndstrs[0];
	
	@Ffields =split(/\t/,$fstr);
	@Sfields=split(/\t/,$sndstr);
	if ($Ffields[17] ne 0)
	{
		$SPnum=0;
		$outgff[0]= $Ffields[4];
		$outgff[1]= "LibRepeatMasker";
		$outgff[2]= "Repeats";
		$outgff[3]= $Ffields[6]+1; 
		#$outgff[4]= $fields[6];#modify
		#$outgff[4]= $Ffields[6];	
		$outgff[5]= $Ffields[0];
		$outgff[6]= $Ffields[8];
		if ($Ffields[8] eq "C")
		  { $outgff[6]="-";}      
		$outgff[7]= ".";	
		 $nameR=$Ffields[15];
	       # my ($nameR)=($fields[9]=~/^R\=(\d+)$/);       			
        }
	while (($Sfields[17] ne 0)&&($Ffields[17] eq $Sfields[17])&&($Sfields[16]<$MaxSP))
	{	
		$numstr++;
		$outgff[3]= $Ffields[6]+1;	
		$outgff[4]= $Sfields[5]-1;#modify
		#repeat string is adjacent in some conditions
		if (($outgff[4] < $outgff[3])||($outgff[4] == $outgff[3]))
		{			
			@Ffields=@Sfields;
			@sndstrs=split (/\n/,$sameRstr[$numstr+1]);
			$sndstr=$sndstrs[0];		
			@Sfields=split(/\t/,$sndstr);
			if($exist>0)
			{$exist--;}
			next;
		 }
		
		$outgff[5]= $Ffields[0];
		$outgff[6]= $Ffields[8];
		if ($Ffields[8] eq "C")
		  { $outgff[6]="-";}      
		$outgff[7]= ".";	
		$outgff[8]=$nameR;
		$outgff[9]=$CRISPRnum;
		$outgff[10]=$SPnum++;
		$outgff[11]=$outgff[4]-$outgff[3]+1; #space length		
		$outgff[12]=0;                  ##is or not same spaces
		$outgff[13]=getsubstr($outgff[3],$outgff[4],$bodystr);
		
		my @seqnamestr=split(/\|/,$outgff[0]);
		
		$spaces[0]=">".$seqnamestr[3]."_".$outgff[8]."_".$outgff[9]."_".$outgff[10]."_".$outgff[11]."_".$outgff[12];
		$spaces[1]=$outgff[13];
		push @onespaces,join("\n",@spaces)."\n";	
		push @oneCRISPR ,join("\t", @outgff)."\n";
		
		if ($exist%2==0)
		{
		    $spaces[0]=">even"."_".$seqnamestr[3]."_".$outgff[8]."_".$outgff[9]."_".$outgff[10]."_".$outgff[11]."_".$outgff[12];
		    push @evenspaces,join("\n",@spaces)."\n";
		}
		else
		{
		    $spaces[0]=">odd"."_".$seqnamestr[3]."_".$outgff[8]."_".$outgff[9]."_".$outgff[10]."_".$outgff[11]."_".$outgff[12];
		    push @oddspaces,join("\n",@spaces)."\n";	
		}
		
		@Ffields=@Sfields;
		@sndstrs=split (/\n/,$sameRstr[$numstr+1]);
		$sndstr=$sndstrs[0];		
		@Sfields=split(/\t/,$sndstr);
		$exist++;
	}
	
	if ($exist ne 0)
	{
		#my $idstr="ID=".$ID.";"."Name=".$nameR.";"."StrLen=".$StrLen;
		##my $idstr="ID=".$ID.";".$nameR.";"."StrLen=".$StrLen;
		#$outgff[8]= $idstr;
		$simil=spacescomp(@onespaces);	
		 my $numline=@oneCRISPR;
		if ($simil > $SPsimilary)
		{ 
		   for (my $num=0;$num<$numline;$num++)
		   {
			@eachline=split (/\t/,$oneCRISPR[$num]);
			#if ($eachline[12]== 0)
			{$eachline[12]=1;}
			push @oneRspace ,join("\t", @eachline);
		   }   			
		}
		else
		{
			$oddflag=0;
			$evenflag=0;
			if ($exist >3)
			{
			    $similodd=spacescomp(@oddspaces);
			    if ($similodd > $SPsimilary)
			    {
				$oddflag=1;
			       for (my $num=1;$num<$numline;$num++,$num++)
			       {
				    @eachline=split (/\t/,$oneCRISPR[$num]);
				     if ($eachline[12]== 0)
				     {$eachline[12]=2;}
				     push @oneRspace ,join("\t", @eachline);
			       }   			
			    }
			}
			if ($exist >2)
			{
			    $simileven=spacescomp(@evenspaces);
			    if ($simileven > $SPsimilary)
			    {
				$evenflag=1;
				for (my $num=0;$num<$numline;$num++,$num++)
				{
				     @eachline=split (/\t/,$oneCRISPR[$num]);
				     if ($eachline[12]== 0)
				     {$eachline[12]=2;}
				     push @oneRspace ,join("\t", @eachline);
				}   			
			     }
			}
			
		     if (($evenflag==0)&&($oddflag == 0))
		     {
			push @oneRspace ,@oneCRISPR;
		     }
		     elsif ($evenflag==0)
		     {
			for (my $num=0;$num<$numline;$num++,$num++)
			{
				@eachline=split (/\t/,$oneCRISPR[$num]);				    
				$eachline[12]=3;
				push @oneRspace ,join("\t", @eachline);
			}   		
		     }
		     elsif ($oddflag==0)
		     {
			for (my $num=1;$num<$numline;$num++,$num++)
			{
				@eachline=split (/\t/,$oneCRISPR[$num]);
				$eachline[12]=3;
				push @oneRspace ,join("\t", @eachline);
			}  
		     }
		}
		$CRISPRnum++;			
	}
	$numstr++;	
    }  
	return @oneRspace;			
}

#find spaces and CRISPRs in one kind R
sub findspaceCRISPR
{
   (@sameRstr)=@_;
   $Rnum=@sameRstr;
   $numstr=0;
   @oneRspace=();
   @oneRDR=();
   @oneRgff=();
   $CRISPRnum=0;
   @spgff=();
   #@onespaces=();
   $simil=0;
   $isextendDR=0;
   $existCRISPR=0;
   
    while ($numstr<$Rnum-1)
    #if (($fields[17]==1)&&($fields[16]<$MaxSP))
    {       
        $exist=0;
	@oneCRISPR=();
	@onespaces=();
	@spaces=();
	@oddspaces=();
	@evenspaces=();	
	$spacelen=0;
	$spacenums=0;
	$isCRISPR=0;
	$aveDRlen=0;
	@outgff=();
	@oneDR=();
	
	$sameRstr[$numstr]=~s/[\r\n]+$//g;
	$sameRstr[$numstr+1]=~s/[\r\n]+$//g;
	chomp($sameRstr[$numstr]);
	chomp($sameRstr[$numstr+1]);
	@Ffields =split(/\t/,$sameRstr[$numstr]);
	@Sfields=split(/\t/,$sameRstr[$numstr+1]);
	if ($Ffields[17] ne 0)
	{
		#space length :all,odd,even space len		
		#spaces file contents
		$SPnum=0;
		$spgff[0]= $Ffields[4];
		$spgff[1]= "LibRepeatMasker";
		$spgff[2]= "Repeats";
		$spgff[3]= $Ffields[6]+1; 
		#$spgff[4]= $fields[6];#modify
		#$spgff[4]= $Ffields[6];	
		$spgff[5]= $Ffields[0];
		$spgff[6]= $Ffields[8];
		if ($Ffields[8] eq "C")
		  { $spgff[6]="-";}      
		$spgff[7]= ".";	
		 $nameR=$Ffields[15];
	        #avoid ID= child ID ;its value= first local
		#$ID=$Ffields[14];
		$ID=$Ffields[5];
	       
	       #CRISPR GFF3 file contents
	        $outgff[0]= $Ffields[4];
		$outgff[1]= "LibRepeatMasker";
		#$outgff[2]= "Repeats";
		$outgff[2]= "CRISPR";
		$outgff[3]= $Ffields[5]; 
		#$outgff[4]= $fields[6];#modify
		$outgff[4]= $Ffields[6];
		$aveDRlen=$Ffields[6]-$Ffields[5]+1;
		$outgff[5]= $Ffields[0];
		$outgff[6]= $Ffields[8];
		if ($Ffields[8] eq "C")
		  { $outgff[6]="-";}      
		$outgff[7]= ".";
		$tmpDRstr=substr($bodystr,$Ffields[5]-1,$Ffields[6]-$Ffields[5]+1);
		push @oneDR, $Ffields[4]."\t".$Ffields[5]."\t".$Ffields[6]."\t".$Ffields[14]."\t".$Ffields[15]."\t".$Ffields[16]."\t".$CRISPRnum."\t".$Ffields[0]."\t".$tmpDRstr."\n";	
        }
	while (($Sfields[17] ne 0)&&($Ffields[17] eq $Sfields[17])&&($Sfields[16]<$MaxSP)&&($Sfields[16]>$MinSP))
	{	
		$numstr++;
		$spgff[3]= $Ffields[6]+1;	
		$spgff[4]= $Sfields[5]-1;#modify
		#repeat string is adjacent in some conditions
		if (($spgff[4] < $spgff[3])||($spgff[4] == $spgff[3]))
		{			
			#pop @oneDR;
			@Ffields=@Sfields;
			##$sameRstr[$numstr+1]=~s/[\r\n]+$//g;
			@sndstrs=split (/\n/,$sameRstr[$numstr+1]);
			$sndstr=$sndstrs[0];
			$sndstr=~s/[\r\n]+$//g;
			@Sfields=split(/\t/,$sndstr);
			if($exist>0)
			{$exist--;}
			next;	
		 }
		
		$spgff[5]= $Ffields[0];
		$spgff[6]= $Ffields[8];
		if ($Ffields[8] eq "C")
		  { $spgff[6]="-";}      
		$spgff[7]= ".";	
		$spgff[8]=$nameR;
		$spgff[9]=$CRISPRnum;
		$spgff[10]=$SPnum++;
		$spgff[11]=$spgff[4]-$spgff[3]+1; #space length		
		$spgff[12]=0;                  ##is or not same spaces
		$spgff[13]=substr($bodystr,$spgff[3]-1,$spgff[11]);
		@seqnamestr=split(/\|/,$spgff[0]);	
		
		if (length($spgff[13])>2)
		{
			$spaces[0]=">".$seqnamestr[3]."_".$spgff[8]."_".$spgff[9]."_".$spgff[10]."_".$spgff[11]."_".$spgff[12];
			$spaces[1]=$spgff[13];
		}
		
		my $tempaveDRlen=sprintf("%.1f",$aveDRlen/($exist+1));
		if (($spgff[11]<$SPMaxDRratio*$tempaveDRlen)&&($spgff[11]>$SPMinDRratio*$tempaveDRlen)&& length($spaces[1])>5)
		{
			push @onespaces,join("\n",@spaces)."\n";	
			push @oneCRISPR ,join("\t", @spgff)."\n";
			###############################
			$tmpDRstr=substr($bodystr,$Sfields[5]-1,$Sfields[6]-$Sfields[5]+1);
			push @oneDR, $Sfields[4]."\t".$Sfields[5]."\t".$Sfields[6]."\t".$Sfields[14]."\t".$Sfields[15]."\t".$Sfields[16]."\t".$CRISPRnum."\t".$Sfields[0]."\t".$tmpDRstr."\n";		
			################	
			if ($exist%2==0)
			{
			    $spaces[0]=">even"."_".$seqnamestr[3]."_".$spgff[8]."_".$spgff[9]."_".$spgff[10]."_".$spgff[11]."_".$spgff[12];
			    push @evenspaces,join("\n",@spaces)."\n";
			}
			else
			{
			    $spaces[0]=">odd"."_".$seqnamestr[3]."_".$spgff[8]."_".$spgff[9]."_".$spgff[10]."_".$spgff[11]."_".$spgff[12];
			    push @oddspaces,join("\n",@spaces)."\n";	
			}
			#CRISPR gff3 files
			$spacelen+=$Sfields[16]-1;
			$aveDRlen+=$Sfields[6]-$Sfields[5]+1;
			$outgff[4]= $Sfields[6];#modify		
			#next line 
			@Ffields=@Sfields;
			@sndstrs=split (/\n/,$sameRstr[$numstr+1]);
			$sndstr=$sndstrs[0];
			$sndstr=~s/[\r\n]+$//g;
			@Sfields=split(/\t/,$sndstr);
			$exist++;
			
		}
		else
		{			
			#next line 
			#@Ffields=@Sfields;
			@sndstrs=split (/\n/,$sameRstr[$numstr+1]);
			$sndstr=$sndstrs[0];
			$sndstr=~s/[\r\n]+$//g;
			@Sfields=split(/\t/,$sndstr);			
		}			
		
	}	
	if ($exist ne 0)
	{
		$isCRISPR=1; #is or not a CRISPR
		if ($exist<2)
		{
		    $isCRISPR=2; #space num<2 is questional CRISPR
		}
		$simil=spacescomp(@onespaces);	
		 my $numline=@oneCRISPR;
		if (($simil > $SPsimilary)&&($simil<=1))
		{
		   $isCRISPR=0; #no CRISPR
		   ###the same space string do not display in sp file
		   for (my $num=0;$num<$numline;$num++)
		   {
			@eachline=split (/\t/,$oneCRISPR[$num]);
			#if ($eachline[12]== 0)
			{$eachline[12]=1;}
			push @oneRspace ,join("\t", @eachline);
		   }   			
		}
		else
		{
			$oddflag=0;
			$evenflag=0;
			if ($exist >3)
			{
			    $similodd=spacescomp(@oddspaces);
			    if (($similodd > $SPsimilary)&&($similodd<=1))
			    {
				$isCRISPR=2;#questional CRISPR
				$oddflag=1;
			       for (my $num=1;$num<$numline;$num++,$num++)
			       {
				    @eachline=split (/\t/,$oneCRISPR[$num]);
				     if ($eachline[12]== 0)
				     {$eachline[12]=2;}
				     push @oneRspace ,join("\t", @eachline);
			       }   			
			    }
			}
			if ($exist >2)
			{
			    $simileven=spacescomp(@evenspaces);
			    if (($simileven > $SPsimilary)&&($simileven<=1))
			    {
				$isCRISPR=2;#questional CRISPR
				$evenflag=1;
				for (my $num=0;$num<$numline;$num++,$num++)
				{
				     @eachline=split (/\t/,$oneCRISPR[$num]);
				     if ($eachline[12]== 0)
				     {$eachline[12]=2;}
				     push @oneRspace ,join("\t", @eachline);
				}   			
			     }
			}			
		     if (($evenflag==0)&&($oddflag == 0))
		     {
			push @oneRspace ,@oneCRISPR;
		     }
		     elsif ($evenflag==0)
		     {
			for (my $num=0;$num<$numline;$num++,$num++)
			{
				@eachline=split (/\t/,$oneCRISPR[$num]);				    
				$eachline[12]=3;
				push @oneRspace ,join("\t", @eachline);
			}   		
		     }
		     elsif ($oddflag==0)
		     {
			for (my $num=1;$num<$numline;$num++,$num++)
			{
				@eachline=split (/\t/,$oneCRISPR[$num]);
				$eachline[12]=3;
				push @oneRspace ,join("\t", @eachline);
			}  
		     }
		}		
		#CRISPR gff3 file contents
		if ($isCRISPR ne 0)
		{
			##the last DR's end position:$outgff[4],the first DR's begin position: $outgff[3]
			$StrLen=$outgff[4]-$outgff[3]+1;
			$avespace=int($spacelen/$exist);
			$aveDRlen=int($aveDRlen/($exist+1));
			$spacenums=$exist;			
			#look for both sides of the present CRISPR :is or not exist the similarly DR strings	
			$oneDRSPlen=$avespace+$aveDRlen;
			$oneDRstr=substr($bodystr,$outgff[3]-1,$aveDRlen);			
			#(seqname,$Rtype,benginlocal,cmpstrlen,cmpstr)
			#return one or two DR(beginlocal)
			@leftaddDRlocs=();
			@rightaddDRlocs=();
			##if use dynamic programming ,need the next 3 lines 
			#$multiR=2;
			#@leftaddDRlocs=blastDRstr($seqnamestr[3],$nameR,$outgff[3]-$multiR*$oneDRSPlen,$multiR,$oneDRSPlen,$oneDRstr);
			#@rightaddDRlocs=blastDRstr($seqnamestr[3],$nameR,$outgff[4]+int(0.8*$avespace),$multiR,$oneDRSPlen,$oneDRstr);
			###if use blast compare,need the next two command lines
			## length of fna file : $fnastrlen
			if ($outgff[3]-2*$oneDRSPlen-1>0)
			{
			    @leftaddDRlocs=blastDRstrold($seqnamestr[3],$nameR,$outgff[3]-2*$oneDRSPlen-1,2*$oneDRSPlen,$oneDRstr);
			}
			elsif ($outgff[3]-$oneDRSPlen-1>0)
			{
			    @leftaddDRlocs=blastDRstrold($seqnamestr[3],$nameR,$outgff[3]-$oneDRSPlen-1,$oneDRSPlen,$oneDRstr);
			}
			if ($outgff[4]+2*$oneDRSPlen<$fnastrlen)
			{
			    @rightaddDRlocs=blastDRstrold($seqnamestr[3],$nameR,$outgff[4]+int(0.8*$avespace)-1,2*$oneDRSPlen-int(0.5*$avespace),$oneDRstr);
			   ## @rightaddDRlocs=blastDRstrold($seqnamestr[3],$nameR,$outgff[4]+int($SPMinDRratio*$aveDRlen)-1,2*$oneDRSPlen-int(0.5*$avespace),$oneDRstr);
			}
			elsif($outgff[4]+$oneDRSPlen<$fnastrlen)
			{
			    @rightaddDRlocs=blastDRstrold($seqnamestr[3],$nameR,$outgff[4]+int(0.8*$avespace)-1,$oneDRSPlen-int(0.5*$avespace),$oneDRstr);
			   ## @rightaddDRlocs=blastDRstrold($seqnamestr[3],$nameR,$outgff[4]+int($SPMinDRratio*$aveDRlen)-1,$oneDRSPlen-int(0.5*$avespace),$oneDRstr);
			}
			
			#add new DR into DRrows
			#In old DRrows left side add DR
			#modify the add numbers lest that existing the overlap DRSP string

			$leftnum=0;
			$rightnum=0;
			$leftnum=@leftaddDRlocs;
			$rightnum=@rightaddDRlocs;
			if ($leftnum>1)
			{
				if (($leftaddDRlocs[0]-$leftaddDRlocs[1]>$aveDRlen)&&($leftaddDRlocs[0]-$leftaddDRlocs[1]>$avespace)&&($leftaddDRlocs[0]-$leftaddDRlocs[1]<$oneDRSPlen+$aveDRlen)||($leftaddDRlocs[1]-$leftaddDRlocs[0]>$aveDRlen)&&($leftaddDRlocs[1]-$leftaddDRlocs[0]>$avespace)&&($leftaddDRlocs[1]-$leftaddDRlocs[0]<$oneDRSPlen+$aveDRlen))
				{$leftnum=2;}
				else
				{ $leftnum=1;}
			}
			if($rightnum>1)
			{
				if (($rightaddDRlocs[0]-$rightaddDRlocs[1]>$aveDRlen)&&($rightaddDRlocs[0]-$rightaddDRlocs[1]>$avespace)&&($rightaddDRlocs[0]-$rightaddDRlocs[1]<$oneDRSPlen+$aveDRlen)||($rightaddDRlocs[1]-$rightaddDRlocs[0]>$aveDRlen)&&($rightaddDRlocs[1]-$rightaddDRlocs[0]>$avespace)&&($rightaddDRlocs[1]-$rightaddDRlocs[0]<$oneDRSPlen+$aveDRlen))
				{$rightnum=2;}
				else
				{$rightnum=1;}
			}
			$addleftDRnum=0;			
			$addrightDRnum=0;		
			if ($leftnum>1)
			{
			#blast result's sequence is unarranged
			    if ($leftaddDRlocs[0]>$leftaddDRlocs[1])
			    {
				$leftleftDRbeg=$leftaddDRlocs[1];
				$leftDRbeg=  $leftaddDRlocs[0];				
			    }
			    else
			    {
				$leftleftDRbeg=$leftaddDRlocs[0];
				$leftDRbeg=  $leftaddDRlocs[1];
			    }
			    $addleftDRnum=2;
			    
			     ##add  two sp contents
			    @llspgff=();
			    @llspaces=();
			    @lspgff=();
			    @lspaces=();
			    $llspgff[0]=$nameR;
			    $llspgff[1]=$CRISPRnum;
			    $llspgff[2]=$SPnum++;
			    $llspgff[3]=$leftDRbeg-$leftleftDRbeg-$aveDRlen; #space length		
			    $llspgff[4]=0;                  ##is or not same spaces
			    $llspgff[5]=substr($bodystr,$leftleftDRbeg+$aveDRlen-1,$llspgff[3]);	
			    $llspaces[0]=">".$seqnamestr[3]."_".$llspgff[0]."_".$llspgff[1]."_".$llspgff[2]."_".$llspgff[3]."_".$llspgff[4];
			    $llspaces[1]=$llspgff[5];
			    if (($llspgff[3]>$MinSP)&&($llspgff[3]<$SPMaxDRratio*$aveDRlen)&&($llspgff[3]>$SPMinDRratio*$aveDRlen))
			    {
				push @onespaces,join("\n",@llspaces)."\n";
				#add DR contents according to the *.dr format
				$tmpDRstr=substr($bodystr,$leftleftDRbeg,$aveDRlen+1);
				push @oneDR, $Ffields[4]."\t".($leftleftDRbeg)."\t".($leftleftDRbeg+$aveDRlen-1)."\t".(0)."\t".$Ffields[15]."\t".(0)."\t".$CRISPRnum."\t".$Ffields[0]."\t".$tmpDRstr."\n";			    
			    }
			    ##another spacer
			    $lspgff[0]=$nameR;
			    $lspgff[1]=$CRISPRnum;
			    $lspgff[2]=$SPnum++;
			    $lspgff[3]=$outgff[3]-$leftDRbeg-$aveDRlen; #space length		
			    $lspgff[4]=0;                  ##is or not same spaces
			    $lspgff[5]=substr($bodystr,$leftDRbeg+$aveDRlen-1,$lspgff[3]);	
			    $lspaces[0]=">".$seqnamestr[3]."_".$lspgff[0]."_".$lspgff[1]."_".$lspgff[2]."_".$lspgff[3]."_".$lspgff[4];
			    $lspaces[1]=$lspgff[5];
			    if (($lspgff[3]>$MinSP)&&($lspgff[3]<$SPMaxDRratio*$aveDRlen)&&($lspgff[3]>$SPMinDRratio*$aveDRlen))
			    {
				push @onespaces,join("\n",@lspaces)."\n";
				#add DR contents according to the *.dr format
				$tmpDRstr=substr($bodystr,$leftDRbeg,$aveDRlen+1);
				 push @oneDR, $Ffields[4]."\t".($leftDRbeg)."\t".($leftDRbeg+$aveDRlen-1)."\t".(0)."\t".$Ffields[15]."\t".(0)."\t".$CRISPRnum."\t".$Ffields[0]."\t".$tmpDRstr."\n";			    
			    }			    
			    #modify the CRISPR begin/end localtion
			    $outgff[3]=$leftleftDRbeg ;
			    if ($llspgff[3]<=0)
			    {
				$outgff[3]=$leftDRbeg ;
			    }
			    if(($llspgff[3]<=0)||($lspgff[3]<=0))
			    {$addleftDRnum=1;}
			}
			elsif($leftnum>0)
			{
			   
			   if(($leftaddDRlocs[0]+$aveDRlen)<($outgff[3]-$oneDRSPlen))
			   {
				$leftleftDRbeg=	$leftaddDRlocs[0];
				$leftDRbeg=  $leftaddDRlocs[0]+$oneDRSPlen;
				##next if ---modify value =2 
				$addleftDRnum=0;
				
				 ##add  two sp contents
				@llspgff=();
				@llspaces=();
				@lspgff=();
				@lspaces=();
				$llspgff[0]=$nameR;
				$llspgff[1]=$CRISPRnum;
				$llspgff[2]=$SPnum++;
				$llspgff[3]=$leftDRbeg-$leftleftDRbeg-$aveDRlen; #space length		
				$llspgff[4]=0;                  ##is or not same spaces
				$llspgff[5]=substr($bodystr,$leftleftDRbeg+$aveDRlen-1,$llspgff[3]);	
				$llspaces[0]=">".$seqnamestr[3]."_".$llspgff[0]."_".$llspgff[1]."_".$llspgff[2]."_".$llspgff[3]."_".$llspgff[4];
				$llspaces[1]=$llspgff[5];
				if ($llspgff[3]>$MinSP)
				{
					$addleftDRnum++;
					push @onespaces,join("\n",@llspaces)."\n";
					##add DR contents
					$tmpDRstr=substr($bodystr,$leftleftDRbeg,$aveDRlen+1);
					push @oneDR, $Ffields[4]."\t".($leftleftDRbeg)."\t".($leftleftDRbeg+$aveDRlen-1)."\t".(0)."\t".$Ffields[15]."\t".(0)."\t".$CRISPRnum."\t".$Ffields[0]."\t".$tmpDRstr."\n";
				}
				##another spacer
				$lspgff[0]=$nameR;
				$lspgff[1]=$CRISPRnum;
				$lspgff[2]=$SPnum++;
				$lspgff[3]=$outgff[3]-$leftDRbeg-$aveDRlen; #space length		
				$lspgff[4]=0;                  ##is or not same spaces
				$lspgff[5]=substr($bodystr,$leftDRbeg+$aveDRlen-1,$lspgff[3]);	
				$lspaces[0]=">".$seqnamestr[3]."_".$lspgff[0]."_".$lspgff[1]."_".$lspgff[2]."_".$lspgff[3]."_".$lspgff[4];
				$lspaces[1]=$lspgff[5];
				
				if (($lspgff[3]>$MinSP)&&($lspgff[3]<$SPMaxDRratio*$aveDRlen)&&($lspgff[3]>$SPMinDRratio*$aveDRlen))
				{
					$addleftDRnum++;
					push @onespaces,join("\n",@lspaces)."\n";
					##add DR contents
					$tmpDRstr=substr($bodystr,$leftDRbeg,$aveDRlen+1);
					push @oneDR, $Ffields[4]."\t".($leftDRbeg)."\t".($leftDRbeg+$aveDRlen-1)."\t".(0)."\t".$Ffields[15]."\t".(0)."\t".$CRISPRnum."\t".$Ffields[0]."\t".$tmpDRstr."\n";				
				}
				$outgff[3]=$leftaddDRlocs[0] ;
			   }
			   else
			   {
				$leftDRbeg=  $leftaddDRlocs[0];
				$addleftDRnum=0;
				
				##add  one spacer content
				@lspgff=();
				@lspaces=();
				$lspgff[0]=$nameR;
				$lspgff[1]=$CRISPRnum;
				$lspgff[2]=$SPnum++;
				$lspgff[3]=$outgff[3]-$leftDRbeg-$aveDRlen; #space length		
				$lspgff[4]=0;                  ##is or not same spaces
				$lspgff[5]=substr($bodystr,$leftDRbeg+$aveDRlen-1,$lspgff[3]);	
				$lspaces[0]=">".$seqnamestr[3]."_".$lspgff[0]."_".$lspgff[1]."_".$lspgff[2]."_".$lspgff[3]."_".$lspgff[4];
				$lspaces[1]=$lspgff[5];
				
				if (($lspgff[3]>$MinSP)&&($lspgff[3]<$SPMaxDRratio*$aveDRlen)&&($lspgff[3]>$SPMinDRratio*$aveDRlen))
				{
					$addleftDRnum++;
					push @onespaces,join("\n",@lspaces)."\n";
					
					#add DR content
					$tmpDRstr=substr($bodystr,$leftDRbeg,$aveDRlen+1);
					push @oneDR, $Ffields[4]."\t".($leftDRbeg)."\t".($leftDRbeg+$aveDRlen-1)."\t".(0)."\t".$Ffields[15]."\t".(0)."\t".$CRISPRnum."\t".$Ffields[0]."\t".$tmpDRstr."\n";
				}
				$outgff[3]=$leftaddDRlocs[0] ;
			   }
			}
			#old DRrows right side add DR	
			if ($rightnum>1)
			{
			#blast result's sequence is unarranged
			    if ($rightaddDRlocs[0]>$rightaddDRlocs[1])
			    {
				$rightDRbeg=$rightaddDRlocs[1];
				$rightrightDRbeg= $rightaddDRlocs[0];
			    }
			    else
			    {
				$rightDRbeg=$rightaddDRlocs[0];
				$rightrightDRbeg= $rightaddDRlocs[1];
			    }
			    $addrightDRnum=0;
				
				##add  two sp contents
				@rrspgff=();
				@rrspaces=();
				@rspgff=();
				@rspaces=();
				$rspgff[0]=$nameR;
				$rspgff[1]=$CRISPRnum;
				$rspgff[2]=$SPnum++;
				$rspgff[3]=$rightDRbeg-$outgff[4]+1; #space length		
				$rspgff[4]=0;                  ##is or not same spaces
				$rspgff[5]=substr($bodystr,$outgff[4],$rspgff[3]);	
				$rspaces[0]=">".$seqnamestr[3]."_".$rspgff[0]."_".$rspgff[1]."_".$rspgff[2]."_".$rspgff[3]."_".$rspgff[4];
				$rspaces[1]=$rspgff[5];
				if (($rspgff[3]>$MinSP)&&($rspgff[3]<$SPMaxDRratio*$aveDRlen)&&($rspgff[3]>$SPMinDRratio*$aveDRlen))
				{
				    $addrightDRnum++;
				    push @onespaces,join("\n",@rspaces)."\n";
				    ##add DR into *.dr  :@oneDR
				    $tmpDRstr=substr($bodystr,$rightDRbeg,$aveDRlen+1);
				    push @oneDR, $Ffields[4]."\t".($rightDRbeg)."\t".($rightDRbeg+$aveDRlen-1)."\t".(0)."\t".$Ffields[15]."\t".(0)."\t".$CRISPRnum."\t".$Ffields[0]."\t".$tmpDRstr."\n";			    
				}
				##another spacer
				$rrspgff[0]=$nameR;
				$rrspgff[1]=$CRISPRnum;
				$rrspgff[2]=$SPnum++;
				$rrspgff[3]=$rightrightDRbeg-$rightDRbeg-$aveDRlen; #space length		
				$rrspgff[4]=0;                  ##is or not same spaces
				$rrspgff[5]=substr($bodystr,$rightDRbeg+$aveDRlen-1,$rrspgff[3]);	
				$rrspaces[0]=">".$seqnamestr[3]."_".$llspgff[0]."_".$rrspgff[1]."_".$rrspgff[2]."_".$rrspgff[3]."_".$rrspgff[4];
				$rrspaces[1]=$rrspgff[5];
				if (($rrspgff[3]>$MinSP)&&($rrspgff[3]<$SPMaxDRratio*$aveDRlen)&&($rrspgff[3]>$SPMinDRratio*$aveDRlen))
				{
					$addrightDRnum++;
					push @onespaces,join("\n",@rrspaces)."\n";					    			    
					$tmpDRstr=substr($bodystr,$rightrightDRbeg,$aveDRlen+1);
					push @oneDR, $Ffields[4]."\t".($rightrightDRbeg)."\t".($rightrightDRbeg+$aveDRlen-1)."\t".(0)."\t".$Ffields[15]."\t".(0)."\t".$CRISPRnum."\t".$Ffields[0]."\t".$tmpDRstr."\n";
				}  
			    #modify the CRISPR begin/end localtion
			    $outgff[4]=$rightrightDRbeg+$aveDRlen ;
			    if ($rrspgff[3]<=0)
			    {
				$outgff[4]=$rightDRbeg+$aveDRlen  ;
			    }
			    if(($rrspgff[3]<=0)||($rspgff[3]<=0))
			    {$addrightDRnum=1;}	
			}
			elsif($rightnum>0)
			{   
			   if($rightaddDRlocs[0]>($outgff[4]+$oneDRSPlen))
			   {
				$rightrightDRbeg=$rightaddDRlocs[0];
				$rightDRbeg=  $rightaddDRlocs[0]-$oneDRSPlen;
				##next if modify value=2
				$addrightDRnum=0;
				
				##add  two sp contents
				@rrspgff=();
				@rrspaces=();
				@rspgff=();
				@rspaces=();
				$rspgff[0]=$nameR;
				$rspgff[1]=$CRISPRnum;
				$rspgff[2]=$SPnum++;
				$rspgff[3]=$rightDRbeg-$outgff[4]+1; #space length		
				$rspgff[4]=0;                  ##is or not same spaces
				$rspgff[5]=substr($bodystr,$outgff[4],$rspgff[3]);	
				$rspaces[0]=">".$seqnamestr[3]."_".$rspgff[0]."_".$rspgff[1]."_".$rspgff[2]."_".$rspgff[3]."_".$rspgff[4];
				$rspaces[1]=$rspgff[5];
				if (($rspgff[3]>$MinSP)&&($rspgff[3]<$SPMaxDRratio*$aveDRlen)&&($rspgff[3]>$SPMinDRratio*$aveDRlen))
				{
					$addrightDRnum++;
					push @onespaces,join("\n",@rspaces)."\n";
					#add DR contents
					$tmpDRstr=substr($bodystr,$rightDRbeg,$aveDRlen+1);
					push @oneDR, $Ffields[4]."\t".($rightDRbeg)."\t".($rightDRbeg+$aveDRlen-1)."\t".(0)."\t".$Ffields[15]."\t".(0)."\t".$CRISPRnum."\t".$Ffields[0]."\t".$tmpDRstr."\n";	
				}
				##another spacer
				$rrspgff[0]=$nameR;
				$rrspgff[1]=$CRISPRnum;
				$rrspgff[2]=$SPnum++;
				$rrspgff[3]=$rightrightDRbeg-$rightDRbeg-$aveDRlen; #space length		
				$rrspgff[4]=0;                  ##is or not same spaces
				$rrspgff[5]=substr($bodystr,$rightDRbeg+$aveDRlen-1,$rrspgff[3]);	
				$rrspaces[0]=">".$seqnamestr[3]."_".$llspgff[0]."_".$rrspgff[1]."_".$rrspgff[2]."_".$rrspgff[3]."_".$rrspgff[4];
				$rrspaces[1]=$rrspgff[5];				
				if (($rrspgff[3]>$MinSP)&&($rrspgff[3]<$SPMaxDRratio*$aveDRlen)&&($rrspgff[3]>$SPMinDRratio*$aveDRlen))
				{
					$addrightDRnum++;
					push @onespaces,join("\n",@rrspaces)."\n";	
					#add DR contents
					$tmpDRstr=substr($bodystr,$rightrightDRbeg,$aveDRlen+1);
					push @oneDR, $Ffields[4]."\t".($rightrightDRbeg)."\t".($rightrightDRbeg+$aveDRlen-1)."\t".(0)."\t".$Ffields[15]."\t".(0)."\t".$CRISPRnum."\t".$Ffields[0]."\t".$tmpDRstr."\n";
				}				
				$outgff[4]=$rightaddDRlocs[0]+$aveDRlen ;
			   }
			   else
			   {
			    $rightDRbeg=  $rightaddDRlocs[0];
			    $addrightDRnum=0;
			    
			    ##add  two sp contents
				@rspgff=();
				@rspaces=();
				$rspgff[0]=$nameR;
				$rspgff[1]=$CRISPRnum;
				$rspgff[2]=$SPnum++;
				$rspgff[3]=$rightDRbeg-$outgff[4]+1; #space length		
				$rspgff[4]=0;                  ##is or not same spaces
				$rspgff[5]=substr($bodystr,$outgff[4],$rspgff[3]);	
				$rspaces[0]=">".$seqnamestr[3]."_".$rspgff[0]."_".$rspgff[1]."_".$rspgff[2]."_".$rspgff[3]."_".$rspgff[4];
				$rspaces[1]=$rspgff[5];
				if (($rspgff[3]>$MinSP)&&($rspgff[3]<$SPMaxDRratio*$aveDRlen)&&($rspgff[3]>$SPMinDRratio*$aveDRlen))
				{
				    $addrightDRnum++;
				    push @onespaces,join("\n",@rspaces)."\n";	
				    #add one DR
				    $tmpDRstr=substr($bodystr,$rightDRbeg,$aveDRlen+1);
				    push @oneDR, $Ffields[4]."\t".($rightDRbeg)."\t".($rightDRbeg+$aveDRlen-1)."\t".(0)."\t".$Ffields[15]."\t".(0)."\t".$CRISPRnum."\t".$Ffields[0]."\t".$tmpDRstr."\n";
				    $outgff[4]=$rightaddDRlocs[0]+$aveDRlen ;
				}
			   }
			}						
			#add new DR infor into CRISPR file			
			$spacenums=$spacenums+$addleftDRnum+$addrightDRnum;
			if (($spacenums>2)&&($evenflag==0)&&($oddflag==0))
			{$isCRISPR=1;}	
			
			#add DR ;check SP is or not same?
			$simil=spacescomp(@onespaces);	
		        if (($simil > $SPsimilary)&&($simil<=1))
		        {
			  #$CRISPRnum--;
			  $isCRISPR=0; #no CRISPR
			  next;
		        }
			
			my $idstr="ID=".$ID.";"."Name=".$nameR.";"."CRISPRnum=".$CRISPRnum.";"."StrLen=".$StrLen.";"."aSPLen=".$avespace.";"."SPnum=".$spacenums.";"."aDRLen=".$aveDRlen.";"."isCRISPR=".$isCRISPR;			
			$outgff[8]= $idstr;
			push @oneRgff ,join("\t", @outgff)."\n";							
			#if addrightDRnum addleftDRnum >0, re-sort the same R			
			if (($leftnum>0)||($rightnum>0))
			{
				$isextendDR=1;
				@rowsDR=();
				foreach(@oneDR)
				{
				    @onedrline=();
				    @onedrline=split("\t",$_);
				    push @rowsDR,[@onedrline] ;
				}
				foreach(sort{$a->[4] cmp $b->[4] or $a->[1] <=> $b->[1]}@rowsDR)
				{			    
				    push @oneRDR,join("\t",@$_);
				}				
			}		
			#write the DR strings into new files			
			else
			{
				foreach(@oneDR)
				{
				    push @oneRDR, $_;
				}
			}
			
			if ($isCRISPR ne 0)
			{
			    $existCRISPR=1;
			}			
		}
		$CRISPRnum++;			
	}	
	$numstr++;
    }        
    #exist extend DRs
    my $oneRgffnumb=0;
    $oneRgffnumb=@oneRgff;
    if (($isextendDR==1)&&($existCRISPR ne 0))
    {
	#process the overlap DRs
	#modify the CRISPR contents if exist overlap DRs
	@modRgff=();
	@modRgff=modRgff3(@oneRgff);
      #according to the @modRgff to modify the DR contents @oneRDR
	@modRDR=();
	@modRDR=modoneRDR(@oneRDR);
	
	#add questional CRISPR content according to the CRISPRs (@modRgff)
	# add gff line to @modRgff AND add DRs to @modRDR
	my($RDR_ref,$Rgff_ref)=questionCRISPR(\@sameRstr,\@modRDR,\@modRgff);
	@addquesRDR=@$RDR_ref;
	@addquesRgff=@$Rgff_ref;		
	
	#update gff file contents (add DRs)
	@upgff=();
	@upgff=updategff(\@addquesRDR,\@addquesRgff);	
	
	#update @modRDR:add the information of SPs
	@modDRSP=();
	@modDRSP=getDRSPinfor(@addquesRDR);
	@modSP=();
	@modSP=getSPinfor(@addquesRDR);
	
	##return ( \@modRDR,\@oneRspace,\@modRgff );
	#return ( \@modDRSP,\@oneRspace,\@upgff);
	return ( \@modDRSP,\@modSP,\@upgff);
    }
    elsif (($oneRgffnumb ne 0)&&($existCRISPR ne 0))
    {
	#add questional CRISPR content according to the CRISPRs (@oneRgff)
	# add gff line to @oneRgff AND add DRs to @oneRDR
	my($RDR_ref,$Rgff_ref)=questionCRISPR(\@sameRstr,\@oneRDR,\@oneRgff);
	@addquesRDR=@$RDR_ref;
	@addquesRgff=@$Rgff_ref;	
	
	#update gff file contents
	 @upgff=();
	 @upgff=updategff(\@addquesRDR,\@addquesRgff);
	 #@upgff=updategff(\@oneRDR,\@oneRgff);
	
	#update @oneRDR:add the information of SPs
	 @DRSP=();
	 @DRSP=getDRSPinfor(@addquesRDR);
	 @SP=();
	 @SP=getSPinfor(@addquesRDR);
	 ##return ( \@oneRDR,\@oneRspace,\@oneRgff);
	#return ( \@DRSP,\@oneRspace,\@upgff );
	return ( \@DRSP,\@SP,\@upgff );
    }
    else
	{return 0;}
		
}

#add questional CRISPR
sub questionCRISPR
{
    my @tmpRDR=();
    my @tmpRgff=();
    my @searchstrs=();
    my ($sameR_ref,$DR_ref,$gff_ref)=@_;
    @searchstrs=@$sameR_ref;
    @tmpRDR=@$DR_ref;
    @tmpRgff=@$gff_ref;
    $strnums=@searchstrs;
    
    for (my $num=0;$num<$strnums;$num++)
    {
	 @fstDR=split("\t",$tmpRDR[$num]);
	 @sndDR=split("\t",$tmpRDR[$num+1]);
	 my $DRlen=$fstDR[6]-$fstDR[5]+1;
	 if (($fstDR[16]>0.5*$DRlen)&&(($fstDR[16]< 2.5*$DRlen)||($fstDR[16]<$MaxSP))&&($fstDR[17] eq 0))
	 {
	    $CRISPRnum=$fstDR[15];
	    $tmpDRstr=substr($bodystr,$fstDR[5],$fstDR[6]-$fstDR[5]+1);
	    push @tmpRDR, $fstDR[4]."\t".$fstDR[5]."\t".$fstDR[6]."\t".$fstDR[14]."\t".$fstDR[15]."\t".$fstDR[16]."\t".$CRISPRnum."\t".$fstDR[0]."\t".$tmpDRstr."\n";
	    $tmpDRstr=substr($bodystr,$sndDR[5],$sndDR[6]-$sndDR[5]+1);
	    push @tmpRDR, $sndDR[4]."\t".$sndDR[5]."\t".$sndDR[6]."\t".$sndDR[14]."\t".$sndDR[15]."\t".$sndDR[16]."\t".$CRISPRnum."\t".$sndDR[0]."\t".$tmpDRstr."\n";
	    my @outgff=();
	    $outgff[0]= $fstDR[4];
	    $outgff[1]= "LibRepeatMasker";
	    $outgff[2]= "CRISPR";
	    $outgff[3]= $fstDR[5]; 
	    #$outgff[4]= $fields[6];#modify
	    $outgff[4]= $sndDR[6];
	    $aveDRlength=int(($fstDR[6]-$fstDR[5]+1+$sndDR[6]-$sndDR[5]+1)/2);
	    if ($aveDRlength> $MaxDR)
	    {next;}
	    $aveSPlength=$sndDR[5]-$fstDR[6]-1;
	    $StrLen=$sndDR[6]-$fstDR[5]+1;
	    $outgff[5]= $fstDR[0];
	    $outgff[6]= $fstDR[8];
	    if ($fstDR[8] eq "C")
	    { $outgff[6]="-";}      
	    $outgff[7]= ".";
	    $spacenums=1;
	    $isCRISPR=2;
	    $nameR=$fstDR[15];
	    $ID=$fstDR[5];
	    my $idstr="ID=".$ID.";"."Name=".$nameR.";"."CRISPRnum=".$CRISPRnum.";"."StrLen=".$StrLen.";"."aSPLen=".$aveSPlength.";"."SPnum=".$spacenums.";"."aDRLen=".$aveDRlength.";"."isCRISPR=".$isCRISPR;			
	    $outgff[8]= $idstr;
	    push @tmpRgff ,join("\t", @outgff)."\n";
	    $num++;
	 }
	 else
	 {next;}	
    }	
    return  \@tmpRDR,\@tmpRgff;	
}

#the order of the (@tmpRDR,@tmpRgff) must be identical
sub updategff
{
    my @tmpRDR=();
    my @tmpRgff=();
    my ($DR_ref,$gff_ref)=@_;
    @tmpRDR=@$DR_ref;
    @tmpRgff=@$gff_ref;
    $gfflines=@tmpRgff;
    $DRlines=@tmpRDR;
    my $DRlinenum=0;
    
    for (my $gfflinenum=0;$gfflinenum<$gfflines;$gfflinenum++)
    {
	@analygff=split("\t",$tmpRgff[$gfflinenum]);
	@lastelemDR=split(";",$analygff[8]);
	for (;$DRlinenum<$DRlines;$DRlinenum++)
	{
	    @analyDR=split("\t",$tmpRDR[$DRlinenum]);
	   
	    if (($analyDR[1] >= $analygff[3])&&($analyDR[2] <= $analygff[4]))
	    {
		my @oneDRgff=();
		$oneDRgff[0]=$analyDR[0];
		$oneDRgff[1]=$analygff[1];
		$oneDRgff[2]="Repeats";
		$oneDRgff[3]=$analyDR[1];
		$oneDRgff[4]=$analyDR[2];		
		$oneDRgff[5]=$analyDR[7];		
		$oneDRgff[6]=$analygff[6];
		$oneDRgff[7]=$analygff[7];
		my $DRID= $analyDR[3];
		if ($DRID eq 0)
		{
		    $DRID=$analyDR[2];
		}
		my $DRname=$analyDR[4];
		my $CRISPRnum=$analyDR[6];
		my ($parent)=($lastelemDR[0]=~/ID\=(\d+)/);
		my $lastring="ID=".$DRID.";"."Parent=".$parent.";"."Name=".$DRname."|"."CRISPRnum=".$CRISPRnum."\n";		
		$oneDRgff[8]= $lastring ;				
		push @tmpRgff, join ("\t",@oneDRgff)   ;
		next;
	    }
	    else
	   { last;}	
	}		
    }	
    return @tmpRgff;	
}

##update @modRDR:add the information of SPs
##update @oneRDR:add the information of SPs
 #sequenceName	Start	End	oldID	ID	SPdist	CRISPRnum	DRstring
    #gi|386003090|ref|NC_017528.1|	3114518	3114555	1817	99	37	0	GGTTTCCGTCCCCTCTCGGGGTTTTGGGTCTGACGACG
     #gi|386003090|ref|NC_017528.1|	3115958	3115995	0	99	0	2	CATGCCGGGGCGGTTCAGGGTTTTGGGTCTGACGACT
sub getDRSPinfor
{
    my @DRinfor=@_;
    my $linenums=@DRinfor;
    @DRSPcontents=();
    #read two line to add SP contents
    
    for (my $line=0;$line<$linenums;$line++)
    {
	chomp($DRinfor[$line]);
	chomp($DRinfor[$line+1]);
	@fstDRstr=split("\t",$DRinfor[$line]);
	@sndDRstr=split("\t",$DRinfor[$line+1]);
	if (($fstDRstr[4]eq $sndDRstr[4])&&($fstDRstr[6]eq $sndDRstr[6]))
	{
	    my $SPstrlen=$sndDRstr[1]-$fstDRstr[2]-1;
	    my $SPstr=substr($bodystr,$fstDRstr[2],$SPstrlen);
	    push @fstDRstr,$SPstr;
	    push @DRSPcontents,join("\t",@fstDRstr)."\n";	
	}
	else
	{
	    push @DRSPcontents,join("\t",@fstDRstr)."\n";
	}
    }  
    return @DRSPcontents;	
}

##update @modRDR:display the information of SPs
##update @oneRDR:display the information of SPs
 #sequenceName	Start	End	oldID	ID	SPdist	CRISPRnum	DRstring
    #gi|386003090|ref|NC_017528.1|	3114518	3114555	1817	99	37	0	GGTTTCCGTCCCCTCTCGGGGTTTTGGGTCTGACGACG
     #gi|386003090|ref|NC_017528.1|	3115958	3115995	0	99	0	2	CATGCCGGGGCGGTTCAGGGTTTTGGGTCTGACGACT
sub getSPinfor
{
    my @DRinfor=@_;
    my $linenums=@DRinfor;
    @SPcontents=();
    @getname=();
    $spnum=0;
    #read two line to add SP contents
    
    for (my $line=0;$line<$linenums;$line++)
    {
	chomp($DRinfor[$line]);
	chomp($DRinfor[$line+1]);
	@fstDRstr=split("\t",$DRinfor[$line]);
	@sndDRstr=split("\t",$DRinfor[$line+1]);
	
	@getname=split('\|',$fstDRstr[0]);
	#$SPnotestr=">".$Input."|".$getname[3]."|".$fstDRstr[4]."|".$fstDRstr[6]."|".$spnum;
	$SPnotestr=">".$getname[3]."|".$fstDRstr[4]."|".$fstDRstr[6]."|".$spnum;
	if (($fstDRstr[4]eq $sndDRstr[4])&&($fstDRstr[6] eq $sndDRstr[6]))
	{
	    my $SPstrlen=$sndDRstr[1]-$fstDRstr[2]-1;
	    my $SPstr=substr($bodystr,$fstDRstr[2],$SPstrlen);
	    push @SPcontents,$SPnotestr."\n";
	    push @SPcontents,$SPstr."\n";
	    $spnum++;
	}
	elsif($fstDRstr[6] ne $sndDRstr[6])
	{
	    $spnum=0;	
	}
	
    }  
    return @SPcontents;	
}

#according to the @modRgff to modify the DR contents @oneRDR
sub modoneRDR
{
	my @oneRDR=@_;
	my %mergeimg=();
	my $existmerge=0;
	@modRDR=();
	$Rnums=@oneRDR;
	my $curRnum=1;
	@fstRDR=split("\t",$oneRDR[$curRnum-1]);
	@sndRDR=split("\t",$oneRDR[$curRnum]);
	$SPlen=$sndRDR[1]-$fstRDR[2];
	$DRlen=(($fstRDR[2]-$fstRDR[1])+($sndRDR[2]-$sndRDR[1]))/2;
	$oneDRSPlen=$SPlen+$DRlen;
	while ($curRnum < $Rnums)
	{
	    #distance between < one space length
	    if (((($sndRDR[1]-$fstRDR[2])< $SPlen*1.5)&&(($sndRDR[1]-$fstRDR[2])>$SPlen*0.6))||($fstRDR[6]==$sndRDR[6]))
	    {
		if ( ($fstRDR[6]ne $sndRDR[6]) &&(($sndRDR[1]-$fstRDR[2])< $SPlen*1.5)&&(($sndRDR[1]-$fstRDR[2])>$SPlen*0.6))
		{
		    $existmerge=1;
		    $mergeimg{$sndRDR[6]}=$mergeimg{$fstRDR[6]};
		}
		if ($existmerge eq 0)
		{
		    $mergeimg{$fstRDR[6]}=$fstRDR[6];
		}
		else
		{
		    $fstRDR[6]=$mergeimg{$fstRDR[6]};	
		}	
		push @modRDR,join("\t",@fstRDR);
		$curRnum++;
		@fstRDR=split("\t",$oneRDR[$curRnum-1]);
		#$fstRDR[6]=$mergeimg{$fstRDR[6]};
		@sndRDR=split("\t",$oneRDR[$curRnum]);
		#$sndRDR[6]=$mergeimg{$sndRDR[6]};
		next;
	    }
	    elsif(($fstRDR[6] ne $sndRDR[6])&&(($fstRDR[2]-$sndRDR[1])>$DRlen*0.6)&&(($fstRDR[2]-$sndRDR[1])<$oneDRSPlen*0.8))
	    {
	    #exist one overlap DR	    	   
	    #delete one overlap (sndDR delete),add two length DR
		$existmerge=1;
		$mergeimg{$sndRDR[6]}=$mergeimg{$fstRDR[6]};
		if ($sndRDR[3] ne 0)
		{
		    $sndRDR[6]=$mergeimg{$sndRDR[6]};
		    push @modRDR,join("\t",@sndRDR);
		}
		else
		{
		    $fstRDR[6]=$mergeimg{$fstRDR[6]};
		    push @modRDR,join("\t",@fstRDR);
		}	    
		$curRnum++;$curRnum++;
		@fstRDR=split("\t",$oneRDR[$curRnum-1]);
		@sndRDR=split("\t",$oneRDR[$curRnum]);
		next;		
	    }
	    elsif (($fstRDR[6] ne $sndRDR[6])&&(($fstRDR[1]-$sndRDR[1])>$oneDRSPlen*0.6)&&(($fstRDR[1]-$sndRDR[1])<$oneDRSPlen*1.5))
	    {
	    #exist two overlap DRs
	        $existmerge=1;
	        $mergeimg{$sndRDR[6]}=$mergeimg{$fstRDR[6]};
		my @thdRDR=();
		@thdRDR=split("\t",$oneRDR[$curRnum+1]);
		if ($sndRDR[3] ne 0)
		{
		    pop @modRDR;
		    $sndRDR[6]=$mergeimg{$sndRDR[6]};
		    push @modRDR,join("\t",@sndRDR);
		    $curRnum++;$curRnum++;
		}
		elsif ($thdRDR[3] ne 0)
		{
		    $thdRDR[6]=$mergeimg{$thdRDR[6]};
		    push @modRDR,join("\t",@thdRDR);
		    $curRnum++;$curRnum++;$curRnum++;
		}
		else
		{
		    $fstRDR[6]=$mergeimg{$fstRDR[6]};
		    push @modRDR,join("\t",@fstRDR);
		    $curRnum++;$curRnum++;$curRnum++;
		}	
		@fstRDR=split("\t",$oneRDR[$curRnum-1]);
		@sndRDR=split("\t",$oneRDR[$curRnum]);
		next;	
	    }
	    elsif (($fstRDR[6] ne $sndRDR[6])&&(($fstRDR[1]-$sndRDR[1])>$oneDRSPlen*1.5)&&(($fstRDR[1]-$sndRDR[1])<$oneDRSPlen*2.5))
	    {
	    #exist three overlap DRs
	        $existmerge=1;
	        $mergeimg{$sndRDR[6]}=$mergeimg{$fstRDR[6]};
	        my @thdRDR=();
		@thdRDR=split("\t",$oneRDR[$curRnum+1]);
		if ($sndRDR[3] ne 0)
		{
		    pop @modRDR;
		    pop @modRDR;
		    $sndRDR[6]=$mergeimg{$sndRDR[6]};
		    push @modRDR,join("\t",@sndRDR);
		    $curRnum++;$curRnum++;
		}
		elsif ($fstRDR[3] ne 0)
		{
		    $fstRDR[6]=$mergeimg{$fstRDR[6]};
		    push @modRDR,join("\t",@fstRDR);
		    $curRnum++;$curRnum++;$curRnum++;$curRnum++;
		}
		elsif ($thdRDR[3] ne 0)
		{
		    pop @modRDR;
		    $thdRDR[6]=$mergeimg{$thdRDR[6]};
		    push @modRDR,join("\t",@thdRDR);
		    $curRnum++;$curRnum++;$curRnum++;			
		}
		else
		{
		     $curRnum++;$curRnum++;$curRnum++;
		}	
		@fstRDR=split("\t",$oneRDR[$curRnum-1]);
		@sndRDR=split("\t",$oneRDR[$curRnum]);
		next;	
	    }
	    elsif (($fstRDR[6] ne $sndRDR[6])&&(($fstRDR[1]-$sndRDR[1])>$oneDRSPlen*2.5)&&(($fstRDR[1]-$sndRDR[1])<$oneDRSPlen*3.5))
	    {
	    #exist four overlap DRs
	        $existmerge=1;
	        $mergeimg{$sndRDR[6]}=$mergeimg{$fstRDR[6]};
	        my @befRDR=();
		@befRDR=split("\t",$oneRDR[$curRnum-2]);
	        my @thdRDR=();
		@thdRDR=split("\t",$oneRDR[$curRnum+1]);
		if ($befRDR[3] ne 0)
		{
		    $curRnum++;$curRnum++;$curRnum++;$curRnum++;					
		}
		elsif ($fstRDR[3] ne 0)
		{
		    $curRnum++;$curRnum++;$curRnum++;$curRnum++;$curRnum++;	
		}
		elsif ($sndRDR[3] ne 0)
		{
		    pop @modRDR; pop @modRDR; pop @modRDR;
		    $sndRDR[6]=$mergeimg{$sndRDR[6]};
		    push @modRDR,join("\t",@sndRDR);
		    $curRnum++;$curRnum++;	
		}
		else
		{
		    pop @modRDR;
		    $curRnum++;$curRnum++;$curRnum++;
		}
		@fstRDR=split("\t",$oneRDR[$curRnum-1]);
		@sndRDR=split("\t",$oneRDR[$curRnum]);
		next;
	    }
	    else
	    {
	    #the boundary of the two crisprnums
	        $existmerge=0;
	        $fstRDR[6]=$mergeimg{$fstRDR[6]};
	        push @modRDR,join("\t",@fstRDR);
		#$mergeimg{$fstRDR[6]}=$fstRDR[6];
		$curRnum++;
		@fstRDR=split("\t",$oneRDR[$curRnum-1]);
		@sndRDR=split("\t",$oneRDR[$curRnum]);
		next;   	
	    }
	}
	#add the last one into @modRDR
	$fstRDR[6]=$mergeimg{$fstRDR[6]};
	 push @modRDR,join("\t",@fstRDR);
	 
	return @modRDR;
}

#modify the CRISPR contents if exist overlap DRs
sub modRgff3 
{
	my @oneRgff=@_;
	@modRgff=();
	$Rgffnum=@oneRgff;
	my $oRnum=1;
	$isoverlap=0;	
	@fstRgff=split("\t",$oneRgff[$oRnum-1]);
	@sndRgff=split("\t",$oneRgff[$oRnum]);
	@sndexp=split(";",$sndRgff[$#sndRgff]);
	($sDRlen)=($sndexp[6]=~/\=(\d+)/);
	($sSPlen)=($sndexp[4]=~/\=(\d+)/);
	$sDRSPlen=$sDRlen+$sSPlen;
	### modified  the reference like :2.5,0.8 and so on .................
	while ($oRnum<$Rgffnum)
	{
	    @currentstr=();
	    $isoverlap=1;	  	    
	    if ((($sndRgff[3]-$fstRgff[4])< $sSPlen*1.5 )&& (($sndRgff[3]-$fstRgff[4])> 0 ))
	    #if ((($sndRgff[3]-$fstRgff[4])< $sSPlen*1.5 )&& (($sndRgff[3]-$fstRgff[4])> $sSPlen *0.6 ))
	    {
	    # combine one CRISPR:num:sp1+sp2+1,DR1 add DR2
		($fstspnum)=($fstRgff[$#fstRgff]=~/SPnum\=(\d+)/);
		($sndspnum)=($sndRgff[$#sndRgff]=~/SPnum\=(\d+)/);
		$tmpspnum=$fstspnum+$sndspnum+1;
		($fstCRISPRnums)=($fstRgff[$#fstRgff]=~/CRISPRnum\=(\d+(\,\d+)*)/);
		($sndCRISPRnums)=($sndRgff[$#sndRgff]=~/CRISPRnum\=(\d+)/);
		$tmpCRISPRnums=$fstCRISPRnums."\,".$sndCRISPRnums;		
		@currentstr=@fstRgff;
		$currentstr[4]=$sndRgff[4];
		$_=$currentstr[8];
		s/SPnum=.*;aDRLen/SPnum=$tmpspnum;aDRLen/;
		s/CRISPRnum=.*;StrLen/CRISPRnum=$tmpCRISPRnums;StrLen/;
		$currentstr[8]=$_;
	    }
	    elsif ((($fstRgff[4]-$sndRgff[3])< $sDRlen*1.5 )&& (($fstRgff[4]-$sndRgff[3])> 0 ))    
	    #elsif ((($fstRgff[4]-$sndRgff[3])< $sDRlen*1.5 )&& (($fstRgff[4]-$sndRgff[3])> $sDRlen*0.6 ))
	    {
	    #num:sp1+sp2 DR1+DR2-1
		($fstspnum)=($fstRgff[$#fstRgff]=~/SPnum\=(\d+)/);
		($sndspnum)=($sndRgff[$#sndRgff]=~/SPnum\=(\d+)/);
		$tmpspnum=$fstspnum+$sndspnum;		
		#$tmpCRISPRnums=($fstRgff[$#fstRgff]=~/CRISPRnum\=(\d+)(\,\d+)*\;/)."\,".($sndRgff[$#sndRgff]=~/CRISPRnum\=(\d+)/);		
		($fstCRISPRnums)=($fstRgff[$#fstRgff]=~/CRISPRnum\=(\d+(\,\d+)*)/);
		($sndCRISPRnums)=($sndRgff[$#sndRgff]=~/CRISPRnum\=(\d+)/);
		$tmpCRISPRnums=$fstCRISPRnums."\,".$sndCRISPRnums;
		@currentstr=@fstRgff;
		$currentstr[4]=$sndRgff[4];		
		$_=$currentstr[8];
		s/SPnum=.*;aDRLen/SPnum=$tmpspnum;aDRLen/;
		s/CRISPRnum=.*;StrLen/CRISPRnum=$tmpCRISPRnums;StrLen/;
		$currentstr[8]=$_;
	    }
	    elsif((($fstRgff[4]-$sndRgff[3])< 2*$sDRSPlen )&& (($fstRgff[4]-$sndRgff[3])> $sDRSPlen ))
	    #elsif((($fstRgff[3]-$sndRgff[3])< 1.5*$sDRSPlen )&& (($fstRgff[3]-$sndRgff[3])> 0.6* $sDRSPlen ))
	    {
	    #num:sp1+sp2-1, DR1 +DR2-2
		($fstspnum)=($fstRgff[$#fstRgff]=~/SPnum\=(\d+)/);
		($sndspnum)=($sndRgff[$#sndRgff]=~/SPnum\=(\d+)/);
		$tmpspnum=$fstspnum+$sndspnum-1;	
		#$tmpCRISPRnums=($fstRgff[$#fstRgff]=~/CRISPRnum\=(\d+)(\,\d+)*\;/)."\,".($sndRgff[$#sndRgff]=~/CRISPRnum\=(\d+)/);
		($fstCRISPRnums)=($fstRgff[$#fstRgff]=~/CRISPRnum\=(\d+(\,\d+)*)/);
		($sndCRISPRnums)=($sndRgff[$#sndRgff]=~/CRISPRnum\=(\d+)/);
		$tmpCRISPRnums=$fstCRISPRnums."\,".$sndCRISPRnums;
		@currentstr=@fstRgff;
		$currentstr[4]=$sndRgff[4];		
		$_=$currentstr[8];
		s/SPnum=.*;aDRLen/SPnum=$tmpspnum;aDRLen/;
		s/CRISPRnum=.*;StrLen/CRISPRnum=$tmpCRISPRnums;StrLen/;
		$currentstr[8]=$_;
	    }
	    elsif((($fstRgff[4]-$sndRgff[3])< 3*$sDRSPlen )&& (($fstRgff[4]-$sndRgff[3])> 2*$sDRSPlen ))
	    #elsif((($fstRgff[3]-$sndRgff[3])< 2.5*$sDRSPlen )&& (($fstRgff[3]-$sndRgff[3])> 1.6*$sDRSPlen ))
	    {
	    #num:sp1+sp2-2, DR1 +DR2-3
		($fstspnum)=($fstRgff[$#fstRgff]=~/SPnum\=(\d+)/);
		($sndspnum)=($sndRgff[$#sndRgff]=~/SPnum\=(\d+)/);
		$tmpspnum=$fstspnum+$sndspnum-2;
		($fstCRISPRnums)=($fstRgff[$#fstRgff]=~/CRISPRnum\=(\d+(\,\d+)*)/);
		($sndCRISPRnums)=($sndRgff[$#sndRgff]=~/CRISPRnum\=(\d+)/);
		$tmpCRISPRnums=$fstCRISPRnums."\,".$sndCRISPRnums;
		@currentstr=@fstRgff;
		$currentstr[4]=$sndRgff[4];		
		$_=$currentstr[8];
		s/SPnum=.*;aDRLen/SPnum=$tmpspnum;aDRLen/;
		s/CRISPRnum=.*;StrLen/CRISPRnum=$tmpCRISPRnums;StrLen/;
		$currentstr[8]=$_;
	    }
	     elsif((($fstRgff[4]-$sndRgff[3])< 4*$sDRSPlen ) && (($fstRgff[4]-$sndRgff[3])> 3*$sDRSPlen ))
	     #elsif((($fstRgff[3]-$sndRgff[3])< 3.5*$sDRSPlen )&& (($fstRgff[3]-$sndRgff[3])> 2.6*$sDRSPlen ))
	    {
	    #num:sp1+sp2-3, DR1 +DR2-4
		($fstspnum)=($fstRgff[$#fstRgff]=~/SPnum\=(\d+)/);
		($sndspnum)=($sndRgff[$#sndRgff]=~/SPnum\=(\d+)/);
		$tmpspnum=$fstspnum+$sndspnum-3;
		($fstCRISPRnums)=($fstRgff[$#fstRgff]=~/CRISPRnum\=(\d+(\,\d+)*)/);
		($sndCRISPRnums)=($sndRgff[$#sndRgff]=~/CRISPRnum\=(\d+)/);
		$tmpCRISPRnums=$fstCRISPRnums."\,".$sndCRISPRnums;
		@currentstr=@fstRgff;
		$currentstr[4]=$sndRgff[4];		
		$_=$currentstr[8];
		s/SPnum=.*;aDRLen/SPnum=$tmpspnum;aDRLen/;
		s/CRISPRnum=.*;StrLen/CRISPRnum=$tmpCRISPRnums;StrLen/;
		$currentstr[8]=$_;
	    }
	    else
	    {
		$isoverlap=0;
		push @modRgff,join("\t",@fstRgff);
		@fstRgff=@sndRgff;
		$oRnum++;          
		@sndRgff=split("\t",$oneRgff[$oRnum]);
		@sndexp=split(";",$sndRgff[$#sndRgff]);
		($sDRlen)=($sndexp[6]=~/\=(\d+)/);
		($sSPlen)=($sndexp[4]=~/\=(\d+)/);
		$sDRSPlen=$sDRlen+$sSPlen;
		next;
	    }	    
	    if ($isoverlap==1)
	    {
		@fstRgff=@currentstr;
		$oRnum++;
		@sndRgff=split("\t",$oneRgff[$oRnum]);
		@sndexp=split(";",$sndRgff[$#sndRgff]);
		($sDRlen)=($sndexp[6]=~/\=(\d+)/);
		($sSPlen)=($sndexp[4]=~/\=(\d+)/);
		$sDRSPlen=$sDRlen+$sSPlen;		
		next;
	    }						
	}
	#add the last one into modRgff
	push @modRgff,join("\t",@fstRgff);
	
	return @modRgff;	
}

#(seqname,Rvalue,benginlocal,cmpstrlen,cmpstr) blast
#return oen or two DR(begin,end)
sub blastDRstrold
{
    my ($seqname,$Rtype,$beginloc,$cmplen,$DRstr)=@_;
    $drlen=length($DRstr);
    my $cmpstr=substr($bodystr,$beginloc,$cmplen);
    $firstr=">".$seqname."_".$beginloc."_".$cmplen."\n";
    $tmpFstr=">".$seqname."_".$Rtype."_".$drlen."\n";
    
    $templateF=$seqname."_".$Rtype."_".$drlen."_".$beginloc;
    open (tmpF,">./datas/$templateF")||die "$!\n";
    #open (tmpF,">$templateF")||die "$!\n";
    print tmpF $tmpFstr;
    print tmpF $DRstr;
    close(tmpF);
    $cmpF=$seqname."_".$beginloc."_".$cmplen;
    open (IF,">./datas/$cmpF")||die "$!\n";
    #open (IF,">$cmpF")||die "$!\n";
    print IF $firstr;
    print IF $cmpstr;
    close(IF);
    
#    
	#    #system("/lustrefs/home/temp/rqge/bin/rmblast/bin/formatdb -i $templateF -p F");
	#    #system("/lustrefs/home/temp/rqge/bin/rmblast/bin/blastall -p blastn -i $cmpF -d $templateF -o $cmpF.out -m 8");
    
##system("/home/rqge/bin/ncbi-rmblastn-2.2.28/bin/formatdb -i ./datas/$templateF -p F");
##system("/bwfs/home/rqge/bin/ncbi-rmblastn-2.2.28/bin/blastall -p blastn -i ./datas/$cmpF -d ./datas/$templateF -W 7 -o ./datas/$cmpF.out -m 8");

system("formatdb -i ./datas/$templateF -p F");
system("blastall -p blastn -i ./datas/$cmpF -d ./datas/$templateF -W 7 -o ./datas/$cmpF.out -m 8");

##system("/bwfs/home/rqge/bin/rmblast-2.2.27/bin/blastall -p blastn -i ./datas/$cmpF -d ./datas/$templateF -o ./datas/$cmpF.out -m 8");
#    #
#   # system("/bwfs/home/rqge/bin/rmblast-2.2.27/bin/makeblastdb -in ./datas/$templateF -dbtype nucl");
#    #system("/bwfs/home/rqge/bin/rmblast-2.2.27/bin/blastn -query ./datas/$cmpF -db ./datas/$templateF -out ./datas/$cmpF.out -outfmt 6");
#    

	#    #system("/bwfs/home/rqge/bin/rmblast-2.2.27/bin/makeblastdb -in ./datas/$templateF -out ./datas/$templateF.out -dbtype nucl");
	#    #system("/bwfs/home/rqge/bin/rmblast-2.2.27/bin/blastn -query ./datas/$cmpF -db ./datas/$templateF.out -out ./datas/$cmpF.out -outfmt 6");
    
    open (outF,"./datas/$cmpF.out")||die "$!\n";
     #open (outF,"$cmpF.out")||die "$!\n";
    @addDRbegs=();
    while (<outF>)
    {
	my $body=$_;
	@oneline=split("\t",$body);
	$sameratio=$oneline[3]/$drlen;
	if ($sameratio> $blastsameR)
	{
	    $tempbeg=$oneline[8];
	    if ($oneline[8]>$oneline[9])
	    {
		$tempbeg=$oneline[9];
	    }
	    $DRbeg=$beginloc+$oneline[6]-$tempbeg;
	    push @addDRbegs,$DRbeg;
	}	
    }
   # #system("rm ./datas/$cmpF.*");
   # #system("rm ./datas/$templateF.*");
   ## unlink "./datas/$cmpF.*";
   # #unlink "./datas/$templateF.*";
    return @addDRbegs;	
}


#(seqname,Rvalue,benginlocal,cmpstrlen,cmpstr) blast
#return oen or two DR(begin,end)
#边界存在排序不一致情??故不能用blast
sub blastDRstr
{
    #my ($seqname,$Rtype,$beginloc,$cmplen,$DRstr)=@_;
    my ($seqname,$Rtype,$beginloc,$multiR,$oneDRSPlen,$DRstr)=@_;
    $drlen=length($DRstr);
    #my $cmpstr=substr($bodystr,$beginloc,$cmplen);
    @addDRbegs=();
    my $sameloc=0;
    for (my $loopnum=0;$loopnum<$multiR;$loopnum++ )
    {	
	$cmpstr=substr($bodystr,$beginloc+$sameloc+$loopnum*$oneDRSPlen,$drlen);
	($astr,$fstr,$samelen)=TirSimilar($DRstr,$cmpstr);
	($cmp_result,$same_similary)=Str_Compare($astr,$fstr);
	
	if ($same_similary>$SameRratio)
	{
	    $sameloc=0;
	    @cmpresult=split("",$cmp_result);
	    while($cmpresult[$sameloc] ne "*")
	    { $sameloc++;}
	    $DRbeg=$beginloc+$sameloc+$loopnum*$oneDRSPlen;
	    push @addDRbegs,$DRbeg;	
	}	
	#print $fstr;print "\n";
	#print $astr; print "\n";
	#print $cmp_result;print "\n";
	#print $same_similary;	
    }
    return @addDRbegs;
}


#clustalw the sameR spaces
sub spacescomp
{
    #change the format
    my @compstr=@_;
    ($tempfile)=($compstr[0]=~/^>(.*)\n/);
    open (IF,">./datas/$tempfile")||die "Error in opening the file: spacescomp $tempfile\n";
    #open (IF,">$tempfile")||die "$!\n";
    print IF @compstr; 
    close(IF);
    
    $dir="./datas";
    chdir($dir);
    deletenulline($tempfile); 
    chdir("..");
    
    system(" clustalw -INFILE=./datas/$tempfile -ALIGN -OUTFILE=./datas/$tempfile.multialign "); 
    open (clustoutF,"<./datas/$tempfile.multialign")||die "Error in opening the file: $tempfile.multialign!";
   
    #open (clustoutF,"$tempfile.multialign")||die "$!\n";
    @strbody=();
    $linenum=0;
    $samenum=0;
    $similarly=0;   
    #$strlen=1;
    $avestrlen=1;
    $cmpnum=0;
    while(<clustoutF>)
    {
	#$samenum++ if ($_=~/\*/g);
	next if /^#/;  # Allow embedded comments.
	$linenum++;
	$body =$_;
	$body=~s/^\s+//;
#	if (($linenum > 1)&&(length($body) ne 0)&&($body=~/A|T|C|G/i))
#	{
#	    @strbody=split(/\s+/,$body);
#            @strhead=split(/_/,$strbody[0]);
#	    if ($strlen < $strhead[$#strhead-1])
#            {
#                $strlen=$strhead[$#strhead-1];
#            }
#	}
	
	if (($linenum > 1)&&(length($body) ne 0)&&($body=~/A|T|C|G/i))
	{
	    @strbody=split(/\s+/,$body);
            @strhead=split(/_/,$strbody[0]);
	    $cmpnum++;
	    $avestrlen+=$strhead[$#strhead-1];    	    
	}    
        if ($body=~/\*/)
        {
            $coutnum=length($body);
            @samestr=split(//,$body);
            for (my $num=0;$num<$coutnum;$num++)
            {
                if ($samestr[$num] eq '*')
                {
                    $samenum++;
                }	
            }
        }       
    }
    if ($cmpnum>0)
    {
	$avestrlen=$avestrlen/$cmpnum;
    }
    $similarly=$samenum/$avestrlen;
    close(clustoutF); 
    ##delete the temp clustalw files
    #system("rm -f ./datas/$tempfile ./datas/$tempfile.multialign ./datas/*.dnd");
    #system("rm -f  ./datas/*.dnd");
    #unlink "./datas/$tempfile" ,"./datas/$tempfile.multialign" ,"./datas/*.dnd";
    #unlink  "./datas/*.dnd";
    
    return $similarly;	
}

#display the same gaps string contents in the same kind of R
sub samegapstring
{
   ($bodystr,@Rstr)=@_;
   $Rnum=@Rstr;
   $presentnum=0;
   @sequencestr=();
  
   while ($presentnum<$Rnum-2)
   {
	@onestr=split /\t/,$Rstr[$presentnum];
	if (($onestr[17]ne 0)&&($presentnum<$Rnum-2))
	{
	     push @sequencestr,"#".$Rstr[$presentnum];
	    
	     $oneresultstr=substr($bodystr,$onestr[5],$onestr[6]-$onestr[5]+1)."\n";
	     push @sequencestr,$oneresultstr;
	     
	     @beforstr=@onestr;	
	     @onestr=split /\t/,$Rstr[++$presentnum];
	     if($onestr[16]>$MaxSP)
	     {
		 pop @sequencestr;
		 pop @sequencestr;
	     } 
	    
	     while (($onestr[17]ne 0)&&($presentnum<$Rnum-1)&&($onestr[16]<$MaxSP)&&($onestr[17]eq $beforstr[17])&&($onestr[16]> -2))
	     {
		#@beforestr=split /\t/,$Rstr[$presentnum-1];
		$beginnum=$beforstr[6]+1;
		$longnum=$onestr[5]-$beforstr[6]-1;
		#if($longnum<0)
		#{
		#	$beginnum=$beforstr[5];
		#	$longnum=$onestr[5]-$beforstr[5];
		#	$oneresultstr="\n".substr($bodystr,$beginnum,$longnum)."\n";
		#	push @sequencestr,$oneresultstr;	
		#}
		if(($longnum<$MaxSP)&&($longnum>0))
		{
			$oneresultstr=substr($bodystr,$beginnum,$longnum)."\n";
			push @sequencestr,$oneresultstr;
		}
		
		$oneresultstr=substr($bodystr,$onestr[5],$onestr[6]-$onestr[5]+1)."\n";
		push @sequencestr,"##".$Rstr[$presentnum];
		push @sequencestr,$oneresultstr;
		
		@beforstr=@onestr;
		@onestr=split /\t/,$Rstr[++$presentnum];
		
	     }
	    # if (! @sequencestr)
	     {
		 push @sequencestr,"//"."\n";
	     }	     
	}
	else
	{
	    $presentnum++;
	 }	
   }
    
   return @sequencestr;	
}

#互补字符序列修改
sub DNA_reverser {
    my($Seq) = @_;
	$Seq = reverse $Seq;
	$Seq =~ tr/ACGTacgt/TGCAtgca/;
    return($Seq);
}
#比较两个相同长度字符串各个字符的异同
sub Str_Compare
{
    my ($Seq1,$Seq2) =@_;
    $prelen = length($Seq1);
    $afterlen = length($Seq2);
    @BP1 = split(//, $Seq1);
    @BP2 = split(//, $Seq2);
    $sum_same=0;
    $cmp_result='';
    if ($prelen != $afterlen){
    return "two strings are not same length";}
    for ($i=0;$i<$prelen;$i++)
    {
	if ($BP1[$i] eq $BP2[$i])
	{
	    $cmp_result.="*";
	    $sum_same++;
	}
	else
	{$cmp_result.="_";}	
    }
    
    $similary = ($sum_same)/($prelen);
    return ($cmp_result,$similary);
}

#compare two strings identify?
sub TirSimilar 
{
   my($Seq1, $Seq2) = @_;
   #my ($prestr,$afterstr)=@_;
   $prelen = length($Seq1);
   $afterlen = length($Seq2);
   @BP1 = split(//, $Seq1);
   @BP2 = split(//, $Seq2);
   $MATCH=1;
   $MISMATCH=-3;
   $GAP =-1;

   $equal_score=0;
   $left_score=0;
   $up_score=0;   
   my $prem=0;
   my $aftern=0;
   $same_len=0;      
   $score[0][0]=0;
   $point[0][0]=0; 

   for ( $i=1;$i<=$prelen;$i++)
   {
       $score[$i][0]=0;
       $point[$i][0]=3;
   }
   for ( $j=1;$j<=$afterlen;$j++)
   {
       $score[0][$j]=$j*$GAP;
       $point[0][$j]=2;
   }       
   for ( $i=1;$i<=$prelen;$i++)
   {
       $prech= $BP1[$i-1];
      for ( $j=1;$j<=$afterlen;$j++)
       {
          $afterch=$BP2[$j-1];
          if ($prech eq $afterch)
          {
             $equal_score=$score[$i-1][$j-1]+$MATCH;             
          }
          else
          {
              $equal_score=$score[$i-1][$j-1]+$MISMATCH;
          }          
          $left_score=$score[$i][$j-1]+$GAP;
          $up_score= $score[$i-1][$j]+$GAP;
          if ($equal_score>$left_score)
	  {
             if ($equal_score>$up_score)
             {
                $score[$i][$j]=$equal_score;
                $point[$i][$j]=1;
             }
             else
             {
                 $score[$i][$j]=$up_score;
                 $point[$i][$j]=3;    
             }
	  }
           else
             {
                if ($up_score >= $left_score) 
                {
                    $score[$i][$j]   = $up_score;
                    $point[$i][$j] = 3;
                }
                else 
                {
                    $score[$i][$j]   = $left_score;
                    $point[$i][$j] = 2;     
                }    
             }   
 #       //printf("score i,j:%d\t%d\t%d\n",i,j,point[i][j]);      
       }      
   }  
    $i--;
    $j--;
   while(($i > 0) || ($j > 0))
   {
       if (0==$point[$i][$j])
        {  last;}
       if ($point[$i][$j]==1)
       {
           $pre[$prem++]=$BP1[$i-1];
           $after[$aftern++]=$BP2[$j-1];
           $i--;
           $j--;
           if($pre[$prem] eq $after[$aftern])
           {
               $same_len++;
           }
       }
       if ($point[$i][$j]==2)
       {
           $pre[$prem++]='-';
           $after[$aftern++]=$BP2[$j-1];
           $j--;                     
       }
       if ($point[$i][$j]==3)
       {
           $pre[$prem++]=$BP1[$i-1];
           $after[$aftern++]='-';
           $i--;                
       }       
   } 
#string revert
   $prerev=reverse(@pre);
   $afterrev=reverse(@after);
   @pre = ();
   @after = ();

   #printf("compare result:\n   prestr:$prerev\n afterstr:$afterrev\n");  
  return ($prerev,$afterrev,$same_len);
}

#analyses the template R string--all hash templates
#$Infile\_rep_filt_stg2_thresh2.fasta
sub sameRtemplate
{
    my ($infile)=@_;
    %similR=();
    #%reversecmp=();
    @Rnames=();
    $Rnums=0;
    my $RName=();
    open(DiffR, "$infile")||die "$!\n";
    while(<DiffR>)
    {  
	chomp;
	if(/^>(\S+)/)
        {
	    my $allname=$1;
	    ($RName)=($allname=~/^R\=(\d+)/);
            $Rnames[$Rnums++]=$RName;
	}
        else
        {
	    $Seq{$RName}= $_;
            $Len{$RName}=length($_);
            $flag{$RName}=0;
	    #$reverseflag{$RName}=0;
	}    
    }
    close(DiffR);    
    for (my $loop=0;$loop<$Rnums-1;$loop++)
    {
        @someRrow=();
	my $fname=$Rnames[$loop];
	
        for (my $secloop=$loop+1;$secloop<$Rnums;$secloop++)
        {
           
            my $sname=$Rnames[$secloop];
            #if ((abs($Len{$fname}/$Len{$sname}-1)<0.3)&&($flag{$sname}==0))
	    if ($flag{$sname}==0)
            {
                #compare strings
                ($fstr,$astr,$samelen)=TirSimilar($Seq{$fname},$Seq{$sname});
		($cmp_result,$same_similary)=Str_Compare($fstr,$astr);
                #if similarly>0.85 into one group; $flag{$sname}=1; push($similR{$fname},$sname)
                if ($same_similary > $DiffTempSim)
                {
                   if ($flag{$fname}==0)
		   {
			push (@someRrow,$fname);
		        $flag{$fname}=1;
		    }
		   push (@someRrow,$sname);
                   $flag{$sname}=1;
                }		
            }
	    
#	   ##考虑在这里增加反相互补串的查找。。。 	   
#	   if ((abs($Len{$fname}/$Len{$sname}-1)<0.3)&&($reverseflag{$fname}==0)&&($reverseflag{$sname}==0))
#            {
#                #compare strings
#		my $reversefname=DNA_reverser($Seq{$fname});
#                ($fstr,$astr,$samelen)=TirSimilar($reversefname,$Seq{$sname});
#		($cmp_result,$same_similary)=Str_Compare($fstr,$astr);
#                #if similarly>0.85 into one group; $flag{$sname}=1; push($similR{$fname},$sname)
#                if ($same_similary > $DiffTempSim)
#                {
#                   $reversecmp{$fname}=$sname;
#		   $reversecmp{$sname}=$fname;
#                   $reverseflag{$sname}=1;
#		   $reverseflag{$fname}=1;
#                }		
#            }	    	    
        }
	my $rownum=@someRrow;	
	my $simRrows=join("\t",@someRrow);	
	if ($rownum>0)
	{$similR{$fname}=$simRrows;}
	#else {$similR{$fname}=$fname;}
    }
    return %similR;     
}

#analyses the template R string
#$Infile\_rep_filt_stg2_thresh2.fasta
sub sameRtemplatePart
{
    my ($infile)=@_;
    %similR=();
    @Rnames=();
    $Rnums=0;
    my $RName=();
    open(DiffR, "$infile")||die "$!\n";
    while(<DiffR>)
    {  
	chomp;
	if(/^>(\S+)/)
        {
	    my $allname=$1;
	    ($RName)=($allname=~/^R\=(\d+)/);
            $Rnames[$Rnums++]=$RName;
	}
        else
        {
	    $Seq{$RName}= $_;
            $Len{$RName}=length($_);
            $flag{$RName}=0;
	}    
    }
    close(DiffR);    
    for (my $loop=0;$loop<$Rnums-1;$loop++)
    {
        @someRrow=();
	my $fname=$Rnames[$loop];
        for (my $secloop=$loop+1;$secloop<$Rnums;$secloop++)
        {
           
            my $sname=$Rnames[$secloop];
            if ((abs($Len{$fname}/$Len{$sname}-1)<0.3)&&($flag{$fname}==0)&&($flag{$sname}==0))
            {
                #compare strings
                ($fstr,$astr,$samelen)=TirSimilar($Seq{$fname},$Seq{$sname});
		($cmp_result,$same_similary)=Str_Compare($fstr,$astr);
                #if similarly>0.85 into one group; $flag{$sname}=1; push($similR{$fname},$sname)
                if ($same_similary > $DiffTempSim)
                {
                   push (@someRrow,$sname);
                   $flag{$sname}=1;
                }		
            }
        }
	my $rownum=@someRrow;	
	my $simRrows=join("\t",@someRrow);	
	if ($rownum>0)
	{$similR{$fname}=$simRrows;}   
    }
    return %similR;     
}


##after R template update,check is or not exist the reverse complement strings
sub ReverseCmpTemp
{
   my ($infile,%newallR)=@_;
    %reverseR=();
    @Rnames=();
    $Rnums=0;
    my $RName=();
    open(DiffR, "$infile")||die "$!\n";
    while(<DiffR>)
    {  
	chomp;
	if(/^>(\S+)/)
        {
	    my $allname=$1;
	    ($RName)=($allname=~/^R\=(\d+)/);
            $Rnames[$Rnums++]=$RName;
	}
        else
        {
	    $Seq{$RName}= $_;
            $Len{$RName}=length($_);
            $flag{$RName}=0;
	}    
    }
    close(DiffR);    
    for (my $loop=0;$loop<$Rnums-1;$loop++)
    {
	my $fname=$Rnames[$loop];
	my $fkey=-1;
	   $fkey=gethashkey($fname,%newallR);
	if (exists($reverseR{$fkey}))
	{ next;  }
	$reversefstr=DNA_reverser($Seq{$fname});
	
        for (my $secloop=$loop+1;$secloop<$Rnums;$secloop++)
        {
           
            my $sname=$Rnames[$secloop];
	    my $skey=-1;
	       $skey=gethashkey($sname,%newallR);
	    if (exists($reverseR{$skey}))
		{ next;  }
            #if ((abs($Len{$fname}/$Len{$sname}-1)<0.3)&&($flag{$fname}==0)&&($flag{$sname}==0))
	     if (($flag{$fname}==0)&&($flag{$sname}==0))
            {
                #compare strings
                ($fstr,$astr,$samelen)=TirSimilar($reversefstr,$Seq{$sname});
		($cmp_result,$same_similary)=Str_Compare($fstr,$astr);
                #if similarly>0.85 into one group; $flag{$sname}=1; push($similR{$fname},$sname)
                if ($same_similary > $DiffTempSim)
                {
		   if (($fkey>=0)&&($skey>=0))
		   {
			$reverseR{$fkey}=$skey;
			#$reverseR{$skey}=$fkey;
			$flag{$sname}=1;
			$flag{$fname}=1;
		   }
		   elsif ($fkey>=0)
		   {
			$reverseR{$fkey}=$sname;
			#$reverseR{$sname}=$fkey;
			$flag{$sname}=1;
			$flag{$fname}=1;
		   }
		   elsif($skey>=0)
		   {
			#$reverseR{$fname}=$skey;
			$reverseR{$skey}=$fname;
			$flag{$sname}=1;
			$flag{$fname}=1;
		   }
		   else
		   {
			$reverseR{$fname}=$sname;
			#$reverseR{$sname}=$fname;
			$flag{$sname}=1;
			$flag{$fname}=1;
		   }
		   last;
                }		
            }
        }	  
    }
    return %reverseR;     	
}

##?ü祷袢ash key value##value中无重复
sub gethashkey
{
   my ($somevalue,%allR)=@_;  
   while (($key,$values)=each %allR)
  {
    @onevalues=split(/\t/,$values);
    for (my $num=0;$num<@onevalues;$num++)
    {
	if ($somevalue eq $onevalues[$num])
	{
	     return $key;
	}       
    }       
  }   
  return -1;	
}

###通过判断方向，及反向互补，获得赋值某R
sub getpositiveR
{
  my ($Rdir,$Rvalue,%reverseR)=@_;
  while (($key,$values)=each %reverseR)
  {
    if (($Rdir eq "+")&&(($Rvalue==$values)||($Rvalue==$key)))
    {
	return $Rvalue;
    }
    elsif ((($Rdir eq "-")||($Rdir eq "C"))&&($Rvalue==$values))
    {
	return $key;
    }
    elsif ((($Rdir eq "-")||($Rdir eq "C"))&&($Rvalue==$key))
    {
	return $values;
    }
  }   
   return -1;	
}

sub outogff
{
    my $analystr=@_;
    @fields =split(/\t/,$analystr);
    @outgff=();
	$outgff[0]= $fields[4];
	$outgff[1]= "LibRepeatMasker";
	$outgff[2]= "Repeats";
        $outgff[3]= $fields[5]; #begin
        $outgff[4]= $fields[6];#end
        $outgff[5]= $fields[0];
        $outgff[6]= $fields[8];
        if ($fields[8] eq "C")
          { $outgff[6]="-";}      
        $outgff[7]= ".";	
        my ($nameR)=($fields[9]=~/^R\=(\d+)$/);
        # my ($nameR)=($fields[15]=~/^R\=(\d+)$/);
	$ID=$fields[14];	
	$StrLen=$outgff[4]-$outgff[3]+1;
	my $idstr="ID=".$ID.";"."Name=".$nameR.";"."StrLen=".$StrLen;
	$outgff[8]= $idstr;	 
        $gffstr= join("\t", @outgff);      	
	return $gffstr;	
}

#delete null line in the file
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

